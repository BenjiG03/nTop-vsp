"""Body (nose) lift of an OpenVSP model's fuselages from their section-area distribution, driven from nTop.

    python -m vsp_ntop.run_body_lift SETTINGS_FILE

The chosen fuselages (``components``: names separated by ``|``, or empty for every fuselage) are sampled on their
OpenVSP surface (CompPnt01) in ``slices`` rings of 64 points; OpenVSP's fuselage cross-sections are planes normal to
x, so each ring gives the section area S(x) and its centroid (the tessellation-based Planar Slice was too coarse at
the nose: 80% off there on a rocket). Several fuselages add on one x grid. On that area
distribution, at every (Mach, AoA, Beta) of the flight condition within Mach Min to Mach Max, with alpha_t the total
angle of attack (acos(cos a cos b)) and the force along the crossflow:

* potential lift (``potential = slender``): slender-body theory, the normal force per length
  dCN/dx = (dS/dx / Sref) sin(2 alpha_t) cos(alpha_t / 2) (Munk; the large-angle form of Jorgensen, NASA TR R-474,
  equation 2.12), so lift sits where the area grows (the nose), none on a constant section, negative on a boattail;
* viscous crossflow (``crossflow = jorgensen``): Allen and Perkins, dCN/dx = eta Cdn (2 r(x) / Sref) sin^2(alpha_t),
  r the radius of the circle of equal area, Cdn of a circular cylinder against the crossflow Mach number M sin(alpha_t)
  (TR R-474 Figure 1) and eta against fineness and crossflow Mach (Figures 4 and 6), as in hypanel's body_lift.

Writes ``<case>.ntvsp.coefficients.csv`` (Mach, AoA, Beta, CA, CY, CN, Cl, Cm, Cn about x aft, y right, z up, with the
flight condition's reference: Cl and Cn over Sref bref, Cm over Sref cref) for Coefficients Contribution (VSP), and
``<case>.ntvsp.body_sections.csv`` (x, S, centroid). Standard library only.
"""
from __future__ import annotations

import bisect
import math
import os

from .common import DriverError, Log, Settings, main_wrapper, openvsp_dir, run_vspscript, vsp_path, write_csv
from .flight import Flight, _rows

# NASA TR R-474 (Jorgensen, 1977), read from its figures (the tables of hypanel's body_lift module).
CDN_MACH = (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 2.4, 2.8, 3.2, 4.0,
            5.0, 8.0)
CDN_VALUE = (1.20, 1.20, 1.22, 1.27, 1.36, 1.50, 1.57, 1.50, 1.42, 1.55, 2.08, 1.96, 1.79, 1.67, 1.55, 1.47, 1.40, 1.34,
             1.32, 1.31, 1.29, 1.27, 1.27)
ETA_LENGTH = (1.0, 2.0, 4.0, 8.0, 10.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0)
ETA_VALUE = (0.53, 0.57, 0.62, 0.66, 0.68, 0.70, 0.73, 0.76, 0.78, 0.80, 0.82, 0.83, 0.84)
ETA_MACH = (0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 2.0)
ETA_MACH_VALUE = (0.69, 0.70, 0.72, 0.78, 0.82, 0.88, 0.99, 0.98, 0.77, 0.91, 0.94, 0.98, 0.98, 1.0)
HEADER = ("Mach", "AoA", "Beta", "CA", "CY", "CN", "Cl", "Cm", "Cn", "Sref", "Bref", "Cref", "Xref", "Yref", "Zref")


def interp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    j = bisect.bisect_right(xs, x)
    t = (x - xs[j - 1]) / (xs[j] - xs[j - 1])
    return ys[j - 1] + t * (ys[j] - ys[j - 1])


def crossflow_drag(mn):
    return interp(mn, CDN_MACH, CDN_VALUE)


def crossflow_factor(fineness, mn):
    """Allen and Perkins' eta: Figure 4 against fineness, its rise with crossflow Mach (Figure 6, for fineness 10 to
    12) rescaled to start from the Figure 4 value (hypanel's assumption)."""
    low = interp(fineness, ETA_LENGTH, ETA_VALUE)
    measured = interp(mn, ETA_MACH, ETA_MACH_VALUE)
    return low + (measured - ETA_MACH_VALUE[0]) * (1.0 - low) / (1.0 - ETA_MACH_VALUE[0])


def ring_script(model, names, n_u, n_w):
    """Rings of points on each chosen fuselage's OpenVSP surface at n_u + 1 values of u (OpenVSP's fuselage
    cross-sections are planes normal to x), printed as ``NTVSP_RING component a x y z ...``."""
    pick = ([f'    {{ array<string> @ids = FindGeomsWithName( "{nm}" ); if ( ids.size() == 0 ) '
             f'{{ Print( "ERR no component named {nm}" ); return; }} chosen.insertLast( ids[0] ); }}' for nm in names]
            if names else
            ['    array<string> @all = FindGeoms();',
             '    for ( int i = 0; i < int( all.size() ); i++ ) { if ( GetGeomTypeName( all[i] ) == "Fuselage" ) '
             'chosen.insertLast( all[i] ); }'])
    return "\n".join([
        "void main()", "{", f'    ReadVSPFile( "{vsp_path(model)}" );', "    Update();", "    array<string> chosen;"]
        + pick + [
        '    if ( chosen.size() == 0 ) { Print( "ERR the model has no fuselage" ); return; }',
        "    for ( int g = 0; g < int( chosen.size() ); g++ ) {",
        f"      for ( int a = 0; a <= {int(n_u)}; a++ ) {{",
        '        string line = "NTVSP_RING " + g + " " + a;',
        f"        for ( int b = 0; b < {int(n_w)}; b++ ) {{",
        f"          vec3d q = CompPnt01( chosen[g], 0, a / {float(n_u)!r}, b / {float(n_w)!r} );",
        '          line += " " + formatFloat( q.x(), "", 0, 8 ) + " " + formatFloat( q.y(), "", 0, 8 ) + " " + '
        'formatFloat( q.z(), "", 0, 8 ); }',
        "        Print( line ); } }",
        '    Print( "NTVSP_OK" );', "}"]) + "\n"


def rings(text):
    """Per component: [(x, area, (y, z) centroid)] along the body, caps (rings at the same x) dropped."""
    comps = {}
    for line in text.splitlines():
        p = line.split()
        if len(p) < 9 or p[0] != "NTVSP_RING":
            continue
        v = [float(t) for t in p[3:]]
        pts = [v[i:i + 3] for i in range(0, len(v) - 2, 3)]
        x = sum(q[0] for q in pts) / len(pts)
        a2 = sum(pts[i][1] * pts[i - 1][2] - pts[i - 1][1] * pts[i][2] for i in range(len(pts)))
        area = 0.5 * abs(a2)
        if area > 1e-30:
            cy = sum((pts[i][1] + pts[i - 1][1]) * (pts[i][1] * pts[i - 1][2] - pts[i - 1][1] * pts[i][2])
                     for i in range(len(pts))) / (3.0 * a2)
            cz = sum((pts[i][2] + pts[i - 1][2]) * (pts[i][1] * pts[i - 1][2] - pts[i - 1][1] * pts[i][2])
                     for i in range(len(pts))) / (3.0 * a2)
        else:
            cy, cz = (sum(q[1] for q in pts) / len(pts), sum(q[2] for q in pts) / len(pts))
        out = comps.setdefault(p[1], [])
        if out and x <= out[-1][0] + 1e-9:
            continue  # an end cap: more rings at the last section's x
        out.append((x, area, (cy, cz)))
    return list(comps.values())


def combine(parts, n):
    """Several fuselages on one x grid: areas add, centroids are area-weighted."""
    x0, x1 = min(p[0][0] for p in parts), max(p[-1][0] for p in parts)
    xs = [x0 + (x1 - x0) * i / n for i in range(n + 1)]
    areas, centres = [], []
    for x in xs:
        tot = cy = cz = 0.0
        for p in parts:
            if not p[0][0] <= x <= p[-1][0]:
                continue
            px = [q[0] for q in p]
            a = interp(x, px, [q[1] for q in p])
            tot += a
            cy += a * interp(x, px, [q[2][0] for q in p])
            cz += a * interp(x, px, [q[2][1] for q in p])
        areas.append(tot)
        centres.append((cy / tot, cz / tot) if tot > 0 else (0.0, 0.0))
    return xs, areas, centres


def body_loads(xs, areas, centres, conds, ref, potential, crossflow):
    """Coefficient rows (HEADER) for each condition."""
    sref, bref, cref, (xr, yr, zr) = ref
    x0, x1 = xs[0], xs[-1]
    radius = [math.sqrt(max(a, 0.0) / math.pi) for a in areas]
    fineness = (x1 - x0) / max(2.0 * max(radius), 1e-12)
    rows = []
    for c in conds:
        a, b = math.radians(c["AoA"]), math.radians(c["Beta"])
        vy, vz = -math.sin(b), math.sin(a) * math.cos(b)
        s = math.hypot(vy, vz)
        at = math.atan2(s, math.cos(a) * math.cos(b))
        ey, ez = (vy / s, vz / s) if s > 1e-12 else (0.0, 0.0)
        mn = c["Mach"] * math.sin(at)
        k_v = crossflow_factor(fineness, mn) * crossflow_drag(mn) if crossflow else 0.0
        cn = my = mz = 0.0  # normal force along the crossflow, moments of it about the reference point
        mx = 0.0
        for i in range(len(xs) - 1):
            dx = xs[i + 1] - xs[i]
            if dx <= 0.0:
                continue
            xm = 0.5 * (xs[i] + xs[i + 1])
            ym = 0.5 * (centres[i][0] + centres[i + 1][0])
            zm = 0.5 * (centres[i][1] + centres[i + 1][1])
            f = 0.0
            if potential:
                f += (areas[i + 1] - areas[i]) / sref * math.sin(2.0 * at) * math.cos(0.5 * at)
            if crossflow:
                f += k_v * (radius[i] + radius[i + 1]) * dx / sref * math.sin(at) ** 2
            fy, fz = f * ey, f * ez
            cn += f
            rx, ry, rz = xm - xr, ym - yr, zm - zr
            mx += ry * fz - rz * fy
            my += -rx * fz
            mz += rx * fy
        rows.append([c["Mach"], c["AoA"], c["Beta"], 0.0, cn * ey, cn * ez, mx / bref, my / cref, mz / bref,
                     sref, bref, cref, xr, yr, zr])
    return rows, fineness


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = openvsp_dir(s)
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    model = s.text("model")
    if not model.lower().endswith(".vsp3"):
        model += ".vsp3"
    if not os.path.isfile(model):
        raise DriverError(f"OpenVSP model {model} not found")
    names = [n.strip() for n in s.text("components", "").split("|") if n.strip()]
    potential = s.text("potential", "slender") == "slender"
    crossflow = s.text("crossflow", "jorgensen") == "jorgensen"
    if not (potential or crossflow):
        raise DriverError("both the potential lift and the viscous crossflow are off")
    flight = Flight(s.text("flight"))
    rows = _rows(s.text("flight").strip() + ".ntvsp.inviscid.csv")
    lo, hi = s.real("mach_min", 0.0), s.real("mach_max", 1e9)
    conds = [r for r in rows if lo <= r["Mach"] < hi]
    n = s.integer("slices", 200)
    text = run_vspscript(vsp_dir, ring_script(model, names, n, 64), work, case + ".rings", log, s.real("timeout_s", 3600))
    parts = rings(text)
    if not parts or max(max(q[1] for q in p) for p in parts) <= 0.0:
        raise DriverError("the chosen fuselages have no section area")
    xs, areas, centres = combine(parts, n)
    ends = (xs[0], xs[-1])
    out, fineness = body_loads(xs, areas, centres, conds, flight.ref, potential, crossflow)
    write_csv(os.path.join(work, case + ".ntvsp.body_sections.csv"), ["x", "Area", "Yc", "Zc"],
              [[x, a, c[0], c[1]] for x, a, c in zip(xs, areas, centres)])
    write_csv(os.path.join(work, case + ".ntvsp.coefficients.csv"), HEADER, out)
    log(f"body lift of {', '.join(names) or 'every fuselage'}: {len(xs)} sections from x = {ends[0]:.4f} to "
        f"{ends[1]:.4f} m, largest area {max(areas):.6g} m^2, fineness {fineness:.2f}; "
        f"{'slender-body potential lift' if potential else ''}{' + ' if potential and crossflow else ''}"
        f"{'viscous crossflow (Jorgensen)' if crossflow else ''}; {len(conds)} conditions")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
