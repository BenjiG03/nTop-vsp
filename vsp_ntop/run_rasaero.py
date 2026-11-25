"""RASAero II on an OpenVSP model of a rocket, driven from nTop.

    python -m vsp_ntop.run_rasaero SETTINGS_FILE
"""
from __future__ import annotations

import os
import re
import subprocess
import time
import xml.sax.saxutils as xs

import numpy as np

from .common import DriverError, Log, Settings, main_wrapper, openvsp_dir, run_vspscript, vsp_path, write_csv
from .flight import Flight
from .run_body_lift import combine, ring_script, rings

INCH = 0.0254
UI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rasaero_ui")
HEADER_ALT = ("Mach", "Altitude", "AoA", "Beta", "CA", "CY", "CN", "Cl", "Cm", "Cn", "Sref", "Bref", "Cref", "Xref",
              "Yref", "Zref")


# --- geometry

def ogive(t, fineness):
    """Tangent ogive r / R at x / L = t for L / R = 2 fineness."""
    lr = 2.0 * fineness
    rho = (1.0 + lr * lr) / 2.0
    return np.sqrt(np.maximum(rho * rho - (lr * (1.0 - t)) ** 2, 0.0)) + 1.0 - rho


def von_karman(t):
    th = np.arccos(1.0 - 2.0 * t)
    return np.sqrt((th - np.sin(2.0 * th) / 2.0) / np.pi)


def body_radius(radius):
    """The body tube's radius: the largest radius the body holds over at least a fifth of its stations (so a
    local bulge, such as a skinning overshoot near a base, does not count)."""
    return np.sort(radius)[int(0.8 * (len(radius) - 1))]


def fit_nose(xs, radius):
    """(length, radius, shape name, RMS of r/R) of the nose: from the tip to where r first reaches 99.9% (else 99%)
    of the body tube's radius."""
    xs, radius = np.asarray(xs), np.asarray(radius)
    big = body_radius(radius)
    reach = np.flatnonzero(radius >= 0.999 * big)
    j = reach[0] if len(reach) else np.flatnonzero(radius >= 0.99 * big)[0]
    length = xs[j] - xs[0]
    if length <= 0.0:
        raise DriverError("the body has no nose: its radius is largest at the front")
    fineness = length / (2.0 * big)
    t, y = (xs[:j + 1] - xs[0]) / length, radius[:j + 1] / big
    shapes = {"Conical": t, "Tangent Ogive": ogive(t, fineness), "Von Karman Ogive": von_karman(t)}
    rms = {n: float(np.sqrt(np.mean((f - y) ** 2))) for n, f in shapes.items()}
    return length, big, min(rms, key=rms.get), rms


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
        out.setdefault(p[1], []).append(np.array(p[3:3 + 3 * ((len(p) - 3) // 3)], float).reshape(-1, 3))
    return out


def fin_planform(ring_list, body_radius_at):
    """Exposed trapezoid of one fin surface: (root LE x, root chord, tip chord, span, sweep distance, thickness,
    LE bevel, TE bevel), outboard of the body."""
    def plane(ring):
        """The ring's centre in (y, z) and the unit normal to the fin's plane there."""
        c = ring[:, 1:].mean(axis=0)
        return c, np.array([-c[1], c[0]]) / (np.hypot(*c) or 1.0)

    def spread(ring):
        c, n = plane(ring)
        h = (ring[:, 1:] - c) @ n
        return h.max() - h.min()

    # OpenVSP's wing surface starts and ends with its end caps
    by_rho = {}
    for ring in ring_list:
        key = round(float(np.hypot(*ring[:, 1:].mean(axis=0))), 6)
        if key not in by_rho or spread(ring) > spread(by_rho[key]):
            by_rho[key] = ring
    stations = sorted((rho, ring[:, 0].min(), ring[:, 0].max(), ring) for rho, ring in by_rho.items())
    tip = stations[-1]
    # root: where the planform crosses the body surface (or the innermost station if the fin stops outside it)
    i = next((k for k, st in enumerate(stations) if st[0] >= body_radius_at(st[1]) - 1e-9), None)
    if i is None:
        raise DriverError("a fin lies entirely inside the body")
    root = stations[i]
    rb = body_radius_at(root[1])
    if i > 0 and stations[i - 1][0] < rb < root[0]:
        a, b = stations[i - 1], root
        t = (rb - a[0]) / (b[0] - a[0])
        le_x, te_x = a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])
    else:
        rb, le_x, te_x = root[0], root[1], root[2]
    chord, tip_chord = te_x - le_x, tip[2] - tip[1]
    # thickness and bevels from the root ring
    c, nrm = plane(root[3])
    x, h = root[3][:, 0], (root[3][:, 1:] - c) @ nrm
    n = 40
    x0, x1 = root[1], root[2]
    thick = np.zeros(n + 1)
    for k, xk in enumerate(x0 + (x1 - x0) * np.arange(n + 1) / n):
        near = h[np.abs(x - xk) <= (x1 - x0) / n]
        thick[k] = np.ptp(near) if len(near) >= 2 else 0.0
    tmax = thick.max()
    full = np.flatnonzero(thick >= 0.95 * tmax)
    cr = x1 - x0
    le_bevel = full[0] * cr / n
    te_bevel = (n - full[-1]) * cr / n
    return {k: float(v) for k, v in (("root_le_x", le_x), ("root_chord", chord), ("tip_chord", tip_chord),
                                      ("span", tip[0] - rb), ("sweep", tip[1] - le_x), ("thickness", tmax),
                                      ("le_bevel", le_bevel), ("te_bevel", te_bevel), ("root_radius", rb))}


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


# --- running RASAero II

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
    values = np.array([table[m] for m in ms])
    return tuple(float(np.interp(mach, ms, values[:, k])) for k in range(values.shape[1]))


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
    radius = np.sqrt(np.maximum(areas, 0.0) / np.pi)
    nose = fit_nose(xs, radius)
    n_len, big, shape, rms = nose
    tail = np.flatnonzero(xs > xs[0] + n_len)
    drop = tail[radius[tail] < 0.95 * big]  # a real boattail
    boattail = None
    if len(drop) and drop[0] > tail[0]:
        start = xs[drop[0] - 1]
        boattail = (start, xs[-1] - start, 2.0 * radius[-1])
    body = (xs[0], xs[-1], boattail)
    axis_y, axis_z = centres.mean(axis=0)

    def body_radius_at(x):
        return float(np.interp(x, xs, radius)) if xs[0] < x < xs[-1] else 0.0

    text = run_vspscript(vsp_dir, wing_script(model, fin_names), work, case + ".fins", log, timeout)
    surfaces = wing_rings(text)
    if not surfaces:
        raise DriverError("no fin surfaces found")
    plan = [fin_planform([r - [0.0, axis_y, axis_z] for r in rl], body_radius_at) for rl in surfaces.values()]
    ref = plan[0]
    for p in plan[1:]:  # measured where each fin meets the fitted round body
        for k in ("root_chord", "tip_chord", "span", "sweep"):
            gap = abs(p[k] - ref[k]) / max(ref["root_chord"], 1e-9)
            if gap > 0.10:
                raise DriverError(f"the fins differ ({k}: {ref[k]:.4f} vs {p[k]:.4f} m): RASAero II takes one fin set")
            if gap > 0.01:
                log(f"warning: the fins differ ({k}: {ref[k]:.4f} vs {p[k]:.4f} m, {100 * gap:.1f}% of the root "
                    "chord); RASAero II takes one fin set: their mean")
    fins = {k: float(np.mean([p[k] for p in plan])) for k in ref}
    n_fins = len(plan)
    summary = (f"nose {shape} (RMS of r/R: " + ", ".join(f"{k} {v:.3f}" for k, v in rms.items()) +
               f"), length {n_len:.4f} m, diameter {2 * big:.4f} m; body to x = {xs[-1]:.4f} m"
               + (f", boattail {boattail[1]:.4f} m to diameter {boattail[2]:.4f} m" if boattail else "")
               + f"; {n_fins} fins: root chord {fins['root_chord']:.4f}, tip chord {fins['tip_chord']:.4f}, span "
               f"{fins['span']:.4f}, sweep distance {fins['sweep']:.4f}, thickness {fins['thickness']:.4f}, bevels "
               f"{fins['le_bevel']:.4f} / {fins['te_bevel']:.4f} m, root leading edge at x = {fins['root_le_x']:.4f} m")
    log("RASAero II design: " + summary)
    open(os.path.join(work, case + ".ntvsp.rasaero_design.txt"), "w", encoding="utf-8").write(summary + "\n")

    # 2. run RASAero II
    def alpha_t(c):
        a, b = np.radians(c["AoA"]), np.radians(c["Beta"])
        return float(np.degrees(np.arccos(np.clip(np.cos(a) * np.cos(b), -1.0, 1.0))))

    # RASAero II's angle-of-attack box takes whole degrees only (it drops a decimal point: 8.94 becomes 894)
    angles = list(range(0, int(np.ceil(max(alpha_t(c) for c in conds) - 1e-9)) + 1))
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
    scale = np.pi * big * big / sref  # RASAero's coefficients are on the body's cross-section area
    rows = []
    for c in conds:
        at = alpha_t(c)
        alt = round(c["Altitude"], 3)
        j = min(int(at), angles[-1] - 1)
        t = at - j
        lo_, hi_ = at_mach(tables[(alt, j)], c["Mach"]), at_mach(tables[(alt, j + 1)], c["Mach"])
        cn = (1 - t) * lo_[0] + t * hi_[0]
        ca = (1 - t) * lo_[1] + t * hi_[1]
        # the moment (CN x_cp) is linear between the runs.
        x_cp = ((1 - t) * lo_[0] * lo_[2] + t * hi_[0] * hi_[2]) / cn if abs(cn) > 1e-9 else lo_[3]
        cn, ca = cn * scale, ca * scale
        a, b = np.radians(c["AoA"]), np.radians(c["Beta"])
        v = np.array([-np.sin(b), np.sin(a) * np.cos(b)])  # the crossflow in (y, z)
        h = np.hypot(*v)
        f = np.array([ca, *(cn * v / h if h > 1e-12 else (0.0, 0.0))])  # CA along +x (aft): drag
        m = np.cross([x_cp + xs[0] - xr, axis_y - yr, axis_z - zr], f) / [bref, cref, bref]  # x_cp from the nose tip
        rows.append([c["Mach"], c["Altitude"], c["AoA"], c["Beta"], *f.tolist(), *m.tolist(),
                     sref, bref, cref, xr, yr, zr])
    write_csv(os.path.join(work, case + ".ntvsp.coefficients.csv"), HEADER_ALT, rows)
    log(f"{len(rows)} conditions written; RASAero II coefficients scaled by {scale:.4f} (body area over Sref)")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
