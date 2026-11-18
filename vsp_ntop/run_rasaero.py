"""RASAero II on an OpenVSP model of a rocket, driven from nTop.

    python -m vsp_ntop.run_rasaero SETTINGS_FILE

RASAero II (Rogers Aeroscience, free, Windows) has no command line. This driver

1. turns the OpenVSP model into a RASAero II design: the fuselage (``body``, or the model's only fuselage) sampled on
   its OpenVSP surface gives the nose (its length is where the radius first reaches the body's; its shape the
   closest of Conical, Tangent Ogive and Von Karman Ogive), the body tube and, when the radius falls at the tail, a
   boattail; the wings (``fins``, or every wing) sampled on their surfaces give one fin set: the count (every surface,
   mirrored copies included), and the exposed root chord, tip chord, span and sweep distance outboard of the body, the
   thickness and a hexagonal section with the leading- and trailing-edge bevel lengths. RASAero II takes trapezoidal
   fins of one kind only, so the fins must all be alike;
2. runs RASAero II (``RASAERO_PATH``: the folder holding ``RASAero II.exe``) through its window, with Microsoft Active
   Accessibility (``rasaero_ui``: the scripts of hypanel's RASAero comparison): opens the design and runs its Run
   Test at each total angle of attack the flight condition needs, once per altitude (the launch site altitude);
3. reads the Run Test files (normal force, axial force and centre of pressure against Mach number, every 0.01) and
   writes ``<case>.ntvsp.coefficients.csv`` at every condition of the flight condition: RASAero II is axisymmetric,
   so at (alpha, beta) the normal force is taken at the total angle acos(cos a cos b) along the crossflow, with the
   centre of pressure on the body axis, in the flight condition's reference (RASAero's coefficients are on the
   body's cross-section area). Also ``<case>.CDX1`` (the design) and ``<case>.ntvsp.rasaero_design.txt`` (what was
   fitted).

Settings: model, flight, body, fins, turbulent (1: all turbulent), surface (RASAero's Surface text), timeout_s,
rasaero_dir, mach_min, mach_max.
"""
from __future__ import annotations

import bisect
import math
import os
import re
import subprocess
import time
import xml.sax.saxutils as xs

from .common import DriverError, Log, Settings, main_wrapper, openvsp_dir, run_vspscript, vsp_path, write_csv
from .flight import Flight
from .run_body_lift import combine, ring_script, rings

INCH = 0.0254
UI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rasaero_ui")
HEADER_ALT = ("Mach", "Altitude", "AoA", "Beta", "CA", "CY", "CN", "Cl", "Cm", "Cn", "Sref", "Bref", "Cref", "Xref",
              "Yref", "Zref")


# --- geometry ------------------------------------------------------------------------------------------------------

def ogive(t, fineness):
    """Tangent ogive r / R at x / L = t for L / R = 2 fineness."""
    lr = 2.0 * fineness
    rho = (1.0 + lr * lr) / 2.0
    return math.sqrt(max(rho * rho - (lr * (1.0 - t)) ** 2, 0.0)) + 1.0 - rho


def von_karman(t):
    th = math.acos(1.0 - 2.0 * t)
    return math.sqrt((th - math.sin(2.0 * th) / 2.0) / math.pi)


def body_radius(radius):
    """The body tube's radius: the largest radius the body holds over at least a fifth of its stations (so a
    local bulge, such as a skinning overshoot near a base, does not count)."""
    ordered = sorted(radius)
    return ordered[int(0.8 * (len(ordered) - 1))]


def fit_nose(xs, radius):
    """(length, radius, shape name, RMS of r/R) of the nose: from the tip to where r first reaches 99% of the
    body tube's radius."""
    big = body_radius(radius)
    j = next((i for i, r in enumerate(radius) if r >= 0.999 * big), None)
    if j is None:
        j = next(i for i, r in enumerate(radius) if r >= 0.99 * big)
    x0, length = xs[0], xs[j] - xs[0]
    if length <= 0.0:
        raise DriverError("the body has no nose: its radius is largest at the front")
    fineness = length / (2.0 * big)
    pts = [((x - x0) / length, r / big) for x, r in zip(xs[:j + 1], radius[:j + 1])]
    shapes = {"Conical": lambda t: t, "Tangent Ogive": lambda t: ogive(t, fineness), "Von Karman Ogive": von_karman}
    rms = {n: math.sqrt(sum((f(t) - y) ** 2 for t, y in pts) / len(pts)) for n, f in shapes.items()}
    best = min(rms, key=rms.get)
    return length, big, best, rms


def wing_script(model, names, n_u=60, n_w=64):
    """Points on each wing's surfaces: per surface (mirrored copies included) n_u + 1 rings of n_w points."""
    pick = ([f'    {{ array<string> @ids = FindGeomsWithName( "{nm}" ); if ( ids.size() == 0 ) '
             f'{{ Print( "ERR no component named {nm}" ); return; }} chosen.insertLast( ids[0] ); }}' for nm in names]
            if names else
            ['    array<string> @all = FindGeoms();',
             '    for ( int i = 0; i < int( all.size() ); i++ ) { if ( GetGeomTypeName( all[i] ) == "Wing" ) '
             'chosen.insertLast( all[i] ); }'])
    return "\n".join([
        "void main()", "{", f'    ReadVSPFile( "{vsp_path(model)}" );', "    Update();", "    array<string> chosen;"]
        + pick + [
        '    if ( chosen.size() == 0 ) { Print( "ERR the model has no wing (fin)" ); return; }',
        "    for ( int g = 0; g < int( chosen.size() ); g++ ) {",
        "      int ns = GetTotalNumSurfs( chosen[g] );",
        "      for ( int sidx = 0; sidx < ns; sidx++ ) {",
        f"        for ( int a = 0; a <= {n_u}; a++ ) {{",
        '          string line = "NTVSP_WING " + g + "_" + sidx + " " + a;',
        f"          for ( int b = 0; b < {n_w}; b++ ) {{",
        f"            vec3d q = CompPnt01( chosen[g], sidx, a / {float(n_u)!r}, b / {float(n_w)!r} );",
        '            line += " " + formatFloat( q.x(), "", 0, 8 ) + " " + formatFloat( q.y(), "", 0, 8 ) + " " + '
        'formatFloat( q.z(), "", 0, 8 ); }',
        "          Print( line ); } } }",
        '    Print( "NTVSP_OK" );', "}"]) + "\n"


def wing_rings(text):
    out = {}
    for line in text.splitlines():
        p = line.split()
        if len(p) < 9 or p[0] != "NTVSP_WING":
            continue
        v = [float(t) for t in p[3:]]
        out.setdefault(p[1], []).append([v[i:i + 3] for i in range(0, len(v) - 2, 3)])
    return out


def fin_planform(ring_list, body_radius_at):
    """Exposed trapezoid of one fin surface: (root LE x, root chord, tip chord, span, sweep distance, thickness,
    LE bevel, TE bevel), outboard of the body."""
    def spread(ring):
        cy, cz = sum(q[1] for q in ring) / len(ring), sum(q[2] for q in ring) / len(ring)
        n = math.hypot(cy, cz) or 1.0
        h = [((q[1] - cy) * -cz + (q[2] - cz) * cy) / n for q in ring]
        return max(h) - min(h)

    # OpenVSP's wing surface starts and ends with its end caps: rings at the root (and tip) span whose section
    # shrinks to nothing. Keep one ring per span station, the fullest.
    by_rho = {}
    for ring in ring_list:
        rho = math.hypot(sum(q[1] for q in ring) / len(ring), sum(q[2] for q in ring) / len(ring))
        key = round(rho, 6)
        if key not in by_rho or spread(ring) > spread(by_rho[key]):
            by_rho[key] = ring
    stations = []
    for rho, ring in by_rho.items():
        le = min(ring, key=lambda q: q[0])
        te = max(ring, key=lambda q: q[0])
        stations.append((rho, le[0], te[0], ring))
    stations.sort(key=lambda s: s[0])
    rhos = [s[0] for s in stations]
    tip = stations[-1]
    # root: where the planform crosses the body surface (or the innermost station if the fin stops outside it)
    root = None
    for s in stations:
        if s[0] >= body_radius_at(s[1]) - 1e-9:
            root = s
            break
    if root is None:
        raise DriverError("a fin lies entirely inside the body")
    i = stations.index(root)
    rb = body_radius_at(root[1])
    if i > 0 and stations[i - 1][0] < rb < root[0]:
        a, b = stations[i - 1], root
        t = (rb - a[0]) / (b[0] - a[0])
        le_x, te_x = a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])
    else:
        rb, le_x, te_x = root[0], root[1], root[2]
    chord, tip_chord = te_x - le_x, tip[2] - tip[1]
    # thickness and bevels from the root ring: half-thickness along the chord
    ring = root[3]
    cy = sum(q[1] for q in ring) / len(ring)
    cz = sum(q[2] for q in ring) / len(ring)
    ny, nz = -cz, cy  # normal to the fin's plane in y-z: across the span direction
    nn = math.hypot(ny, nz) or 1.0
    ny, nz = ny / nn, nz / nn
    pts = [(q[0], (q[1] - cy) * ny + (q[2] - cz) * nz) for q in ring]
    n = 40
    x0, x1 = root[1], root[2]
    thick = []
    for k in range(n + 1):
        x = x0 + (x1 - x0) * k / n
        near = [h for xx, h in pts if abs(xx - x) <= (x1 - x0) / n]
        thick.append(max(near) - min(near) if len(near) >= 2 else 0.0)
    tmax = max(thick)
    cr = x1 - x0
    le_bevel = next((k for k, t in enumerate(thick) if t >= 0.95 * tmax), 0) * cr / n
    te_bevel = (n - max(k for k, t in enumerate(thick) if t >= 0.95 * tmax)) * cr / n
    return {"root_le_x": le_x, "root_chord": chord, "tip_chord": tip_chord, "span": tip[0] - rb,
            "sweep": tip[1] - le_x, "thickness": tmax, "le_bevel": le_bevel, "te_bevel": te_bevel, "root_radius": rb}


def design_xml(nose, body, fins, n_fins, altitude_ft, turbulent, surface):
    """A RASAero II design (.CDX1); lengths in inches."""
    n_len, radius, shape, _ = nose
    x_tip, x_end, boattail = body
    d = 2.0 * radius
    tube_len = (x_end - boattail[0] if boattail else x_end) - (x_tip + n_len)
    i = lambda v: f"{v / INCH:.6f}"  # noqa: E731
    fin_loc = x_end - fins["root_le_x"] if not boattail else (boattail[0] - fins["root_le_x"])
    parts = [f"""    <NoseCone>
      <PartType>NoseCone</PartType>
      <Length>{i(n_len)}</Length>
      <Diameter>{i(d)}</Diameter>
      <Shape>{shape}</Shape>
      <BluntRadius>0</BluntRadius>
      <Location>0</Location>
      <Color>Black</Color>
    </NoseCone>
    <BodyTube>
      <PartType>BodyTube</PartType>
      <Length>{i(tube_len)}</Length>
      <Diameter>{i(d)}</Diameter>
      <LaunchLugDiameter>0</LaunchLugDiameter>
      <LaunchLugLength>0</LaunchLugLength>
      <RailGuideDiameter>0</RailGuideDiameter>
      <RailGuideHeight>0</RailGuideHeight>
      <LaunchShoeArea>0</LaunchShoeArea>
      <Location>{i(n_len)}</Location>
      <Color>Black</Color>
      <BoattailLength>0</BoattailLength>
      <BoattailRearDiameter>0</BoattailRearDiameter>
      <BoattailOffset>0</BoattailOffset>
      <Overhang>0</Overhang>
      <Fin>
        <Count>{n_fins}</Count>
        <Chord>{i(fins['root_chord'])}</Chord>
        <Span>{i(fins['span'])}</Span>
        <SweepDistance>{i(fins['sweep'])}</SweepDistance>
        <TipChord>{i(fins['tip_chord'])}</TipChord>
        <Thickness>{i(fins['thickness'])}</Thickness>
        <LERadius>0</LERadius>
        <Location>{i(fin_loc)}</Location>
        <AirfoilSection>Hexagonal</AirfoilSection>
        <FX1>{i(fins['le_bevel'])}</FX1>
        <FX3>{i(fins['te_bevel'])}</FX3>
      </Fin>
    </BodyTube>"""]
    if boattail:
        parts.append(f"""    <BoatTail>
      <PartType>BoatTail</PartType>
      <Length>{i(boattail[1])}</Length>
      <Diameter>{i(d)}</Diameter>
      <RearDiameter>{i(boattail[2])}</RearDiameter>
      <Location>{i(n_len + tube_len)}</Location>
      <Color>Black</Color>
    </BoatTail>""")
    return f"""<RASAeroDocument>
  <FileVersion>2</FileVersion>
  <RocketDesign>
{chr(10).join(parts)}
    <Surface>{xs.escape(surface)}</Surface>
    <CP>0</CP>
    <ModifiedBarrowman>False</ModifiedBarrowman>
    <Turbulence>{turbulent}</Turbulence>
    <SustainerNozzle>0</SustainerNozzle>
    <Booster1Nozzle>0</Booster1Nozzle>
    <Booster2Nozzle>0</Booster2Nozzle>
    <UseBooster1>False</UseBooster1>
    <UseBooster2>False</UseBooster2>
    <Comments>written by nTop_vsp from an OpenVSP model</Comments>
  </RocketDesign>
  <LaunchSite>
    <Altitude>{altitude_ft:.3f}</Altitude>
    <Pressure>0</Pressure>
    <RodAngle>0</RodAngle>
    <RodLength>10</RodLength>
    <Temperature>59</Temperature>
    <WindSpeed>0</WindSpeed>
  </LaunchSite>
  <Recovery>
    <Altitude1>1000</Altitude1>
    <Altitude2>1000</Altitude2>
    <DeviceType1>None</DeviceType1>
    <DeviceType2>None</DeviceType2>
    <Event1>False</Event1>
    <Event2>False</Event2>
    <Size1>1</Size1>
    <Size2>1</Size2>
    <EventType1>None</EventType1>
    <EventType2>None</EventType2>
    <CD1>1.33</CD1>
    <CD2>1.33</CD2>
  </Recovery>
  <MachAlt />
</RASAeroDocument>
"""


# --- running RASAero II ----------------------------------------------------------------------------------------

def ps(script, *args, timeout=600):
    cmd = ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", os.path.join(UI, script)] + list(args)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if p.returncode:
        raise DriverError(f"{script}: {(p.stderr or p.stdout).strip()[-400:]}")
    return p.stdout


def rasaero_running():
    p = subprocess.run(["tasklist", "/FI", "IMAGENAME eq RASAero II.exe", "/NH"], capture_output=True, text=True)
    return "RASAero II.exe" in p.stdout


def start_rasaero(folder, log):
    exe = os.path.join(folder, "RASAero II.exe")
    if not os.path.isfile(exe):
        raise DriverError(f"RASAERO_PATH = {folder}: no 'RASAero II.exe' there (point it at the RASAero II folder)")
    if rasaero_running():
        return False
    log(f"starting {exe}")
    subprocess.Popen([exe], cwd=folder)
    for _ in range(60):
        time.sleep(1.0)
        if rasaero_running():
            time.sleep(4.0)
            return True
    raise DriverError("RASAero II did not start")


def read_run_test(path):
    """{Mach: (CN, CA, x_cp in metres from the nose tip, low-alpha x_cp)} from a Run Test file."""
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if re.match(r"\s*\d", line):
                rows.append([float(v) for v in line.split()])
    if not rows or len(rows) % 6:
        raise DriverError(f"{path}: not a RASAero II Run Test file ({len(rows)} rows)")
    table = {}
    for i in range(0, len(rows), 6):
        g = rows[i:i + 6]
        table[round(g[0][0], 2)] = (g[3][2], g[4][5], g[3][3] * INCH, g[1][2] * INCH)
    return table


def at_mach(table, mach):
    ms = sorted(table)
    if not ms[0] - 1e-9 <= mach <= ms[-1] + 1e-9:
        raise DriverError(f"Mach {mach:g} is outside RASAero II's range ({ms[0]:g} to {ms[-1]:g})")
    j = min(max(bisect.bisect_left(ms, mach), 1), len(ms) - 1)
    m0, m1 = ms[j - 1], ms[j]
    t = 0.0 if m1 == m0 else (mach - m0) / (m1 - m0)
    t = min(max(t, 0.0), 1.0)
    return tuple((1 - t) * a + t * b for a, b in zip(table[m0], table[m1]))


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = openvsp_dir(s)
    folder = s.text("rasaero_dir")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    model = s.text("model")
    if not model.lower().endswith(".vsp3"):
        model += ".vsp3"
    if not os.path.isfile(model):
        raise DriverError(f"OpenVSP model {model} not found")
    timeout = s.real("timeout_s", 3600)
    flight = Flight(s.text("flight"))
    lo, hi = s.real("mach_min", 0.0), s.real("mach_max", 1e9)
    conds = [c for c in flight.conditions if lo <= c["Mach"] < hi]
    if not conds:
        raise DriverError(f"the flight condition has no Mach number from {lo:g} up to {hi:g}")

    # 1. geometry from the OpenVSP model
    body_names = [n.strip() for n in s.text("body", "").split("|") if n.strip()]
    fin_names = [n.strip() for n in s.text("fins", "").split("|") if n.strip()]
    text = run_vspscript(vsp_dir, ring_script(model, body_names, 300, 64), work, case + ".body", log, timeout)
    parts = rings(text)
    if len(parts) != 1:
        raise DriverError(f"RASAero II needs one body: the model has {len(parts)} fuselages (name it in Body)")
    xs, areas, centres, _ = combine(parts, 300)
    radius = [math.sqrt(max(a, 0.0) / math.pi) for a in areas]
    nose = fit_nose(xs, radius)
    n_len, big, shape, rms = nose
    tail = [i for i in range(len(xs)) if xs[i] > xs[0] + n_len]
    drop = [i for i in tail if radius[i] < 0.95 * big]  # a real boattail, not a skinning wiggle at the base
    boattail = None
    if drop and drop[0] > tail[0]:
        start = xs[drop[0] - 1]
        boattail = (start, xs[-1] - start, 2.0 * radius[-1])
    body = (xs[0], xs[-1], boattail)
    axis_y = sum(c[0] for c in centres) / len(centres)
    axis_z = sum(c[1] for c in centres) / len(centres)

    def body_radius_at(x):
        if x <= xs[0] or x >= xs[-1]:
            return 0.0
        j = bisect.bisect_right(xs, x)
        t = (x - xs[j - 1]) / (xs[j] - xs[j - 1])
        return radius[j - 1] + t * (radius[j] - radius[j - 1])

    text = run_vspscript(vsp_dir, wing_script(model, fin_names), work, case + ".fins", log, timeout)
    surfaces = wing_rings(text)
    if not surfaces:
        raise DriverError("no fin surfaces found")
    plan = [fin_planform([[(q[0], q[1] - axis_y, q[2] - axis_z) for q in r] for r in rl], body_radius_at)
            for rl in surfaces.values()]
    ref = plan[0]
    for p in plan[1:]:  # measured where each fin meets the fitted round body: a noncircular body spreads them
        for k in ("root_chord", "tip_chord", "span", "sweep"):
            gap = abs(p[k] - ref[k]) / max(ref["root_chord"], 1e-9)
            if gap > 0.10:
                raise DriverError(f"the fins differ ({k}: {ref[k]:.4f} vs {p[k]:.4f} m): RASAero II takes one fin set")
            if gap > 0.01:
                log(f"warning: the fins differ ({k}: {ref[k]:.4f} vs {p[k]:.4f} m, {100 * gap:.1f}% of the root "
                    "chord); RASAero II takes one fin set: their mean")
    fins = {k: sum(p[k] for p in plan) / len(plan) for k in ref}
    n_fins = len(plan)
    summary = (f"nose {shape} (RMS of r/R: " + ", ".join(f"{k} {v:.3f}" for k, v in rms.items()) +
               f"), length {n_len:.4f} m, diameter {2 * big:.4f} m; body to x = {xs[-1]:.4f} m"
               + (f", boattail {boattail[1]:.4f} m to diameter {boattail[2]:.4f} m" if boattail else "")
               + f"; {n_fins} fins: root chord {fins['root_chord']:.4f}, tip chord {fins['tip_chord']:.4f}, span "
               f"{fins['span']:.4f}, sweep distance {fins['sweep']:.4f}, thickness {fins['thickness']:.4f}, bevels "
               f"{fins['le_bevel']:.4f} / {fins['te_bevel']:.4f} m, root leading edge at x = {fins['root_le_x']:.4f} m")
    log("RASAero II design: " + summary)
    open(os.path.join(work, case + ".ntvsp.rasaero_design.txt"), "w", encoding="utf-8").write(summary + "\n")

    # 2. run RASAero II: one design per altitude, Run Test at every total angle of attack needed
    def alpha_t(c):
        a, b = math.radians(c["AoA"]), math.radians(c["Beta"])
        return math.degrees(math.acos(max(-1.0, min(1.0, math.cos(a) * math.cos(b)))))

    # RASAero II's angle-of-attack box takes whole degrees only (it drops a decimal point: 8.94 becomes 894), so the
    # Run Test is made at every whole degree up to the largest total angle and interpolated between them.
    angles = list(range(0, int(math.ceil(max(alpha_t(c) for c in conds) - 1e-9)) + 1))
    if len(angles) < 2:
        angles = [0, 1]
    started = start_rasaero(folder, log)
    tables = {}
    turbulent = "True" if s.flag("turbulent", True) else "False"
    try:
        for alt in sorted({round(c["Altitude"], 3) for c in conds}):
            tag = f"_h{alt:g}"
            design = os.path.join(work, f"{case}{tag}.CDX1")
            open(design, "w", encoding="utf-8").write(
                design_xml(nose, body, fins, n_fins, alt / 0.3048, turbulent, s.text("surface", "Smooth (Zero Roughness)")))
            log(f"RASAero II: {design}, Run Test at alpha 0 to {angles[-1]} degrees")
            out = ps("open_design.ps1", "-File", design)
            if "All Turbulent" not in out:
                log("note: open_design did not report the turbulence check box: " + out.strip()[-200:])
            for line in out.strip().splitlines():
                log("  open: " + line.strip()[:200])
            out = ps("run_test.ps1", "-OutDir", work, "-Tag", f"{case}{tag}", "-Alphas",
                     ",".join(str(a) for a in angles), timeout=timeout)
            for line in out.strip().splitlines():
                log("  run test: " + line.strip()[:200])
            for a in angles:
                path = os.path.join(work, f"runtest{case}{tag}_alpha{a}.txt")
                if not os.path.isfile(path):
                    raise DriverError(f"RASAero II wrote no Run Test file for alpha {a:g} ({path})")
                tables[(alt, a)] = read_run_test(path)
    finally:
        if started:
            subprocess.run(["taskkill", "/IM", "RASAero II.exe", "/F"], capture_output=True)

    # 3. coefficients on the flight condition's reference
    sref, bref, cref, (xr, yr, zr) = flight.ref
    scale = math.pi * big * big / sref  # RASAero's coefficients are on the body's cross-section area
    rows = []
    for c in conds:
        at = alpha_t(c)
        alt = round(c["Altitude"], 3)
        j = min(int(at), angles[-1] - 1)
        t = at - j
        lo_, hi_ = at_mach(tables[(alt, j)], c["Mach"]), at_mach(tables[(alt, j + 1)], c["Mach"])
        cn = (1 - t) * lo_[0] + t * hi_[0]
        ca = (1 - t) * lo_[1] + t * hi_[1]
        # the moment (CN x_cp) is linear between the runs; at no normal force, the low-angle centre of pressure
        x_cp = ((1 - t) * lo_[0] * lo_[2] + t * hi_[0] * hi_[2]) / cn if abs(cn) > 1e-9 else lo_[3]
        cn, ca = cn * scale, ca * scale
        a, b = math.radians(c["AoA"]), math.radians(c["Beta"])
        vy, vz = -math.sin(b), math.sin(a) * math.cos(b)
        h = math.hypot(vy, vz)
        cy, cz = (cn * vy / h, cn * vz / h) if h > 1e-12 else (0.0, 0.0)
        rx, ry, rz = x_cp + xs[0] - xr, axis_y - yr, axis_z - zr  # x_cp is from the nose tip
        f = (ca, cy, cz)  # CA along +x (aft): drag
        mx = ry * f[2] - rz * f[1]
        my = rz * f[0] - rx * f[2]
        mz = rx * f[1] - ry * f[0]
        rows.append([c["Mach"], c["Altitude"], c["AoA"], c["Beta"], ca, cy, cz, mx / bref, my / cref, mz / bref,
                     sref, bref, cref, xr, yr, zr])
    write_csv(os.path.join(work, case + ".ntvsp.coefficients.csv"), HEADER_ALT, rows)
    log(f"{len(rows)} conditions written; RASAero II coefficients scaled by {scale:.4f} (body area over Sref)")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
