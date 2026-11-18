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


def _solve(a, rhs):
    """Gaussian elimination with partial pivoting: a (n x n) times columns = rhs columns (list of lists)."""
    n = len(a)
    m = [row[:] + [r[i] for r in rhs] for i, row in enumerate(a)]
    k = len(rhs)
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        piv = m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / piv
            if f:
                rr, rc = m[r], m[c]
                for j in range(c, n + k):
                    rr[j] -= f * rc[j]
    out = [[0.0] * n for _ in range(k)]
    for j in range(k):
        for r in range(n - 1, -1, -1):
            s = m[r][n + j] - sum(m[r][c] * out[j][c] for c in range(r + 1, n))
            out[j][r] = s / m[r][r]
    return out


def apparent_mass(pts):
    """The section's apparent-mass tensor (``_panel_mass``), Richardson-extrapolated from the ring's own edges and
    each edge split in two: the constant-source panels converge to first order (a 64-point circle: +2.2% on 64
    panels, +1.1% on 128, -0.02% extrapolated)."""
    fine = [(p[0] + (q[0] - p[0]) * j / 2, p[1] + (q[1] - p[1]) * j / 2)
            for p, q in zip(pts, pts[1:] + pts[:1]) for j in range(2)]
    m1, m2 = _panel_mass(pts), _panel_mass(fine)
    return [[2.0 * m2[r][c] - m1[r][c] for c in range(2)] for r in range(2)]


def _panel_mass(pts):
    """The 2D apparent-mass tensor of a closed section per unit fluid density, [[m_yy, m_yz], [m_zy, m_zz]] (m^2), from
    ``pts`` [(y, z)] around it: the exterior potential flow of the section moving with unit velocity along y and along
    z (constant-strength source panels, ``d phi / d n = n . e`` on the contour), m_jk = -(contour integral of
    phi_j n_k ds), n out of the body. A circle gives its area times the identity, an ellipse with semi-axes a (along
    y) and b gives m_zz = pi a^2 and m_yy = pi b^2, a flat plate of width w moving across itself pi w^2 / 4."""
    pts = list(pts)
    area2 = sum(pts[i - 1][0] * pts[i][1] - pts[i][0] * pts[i - 1][1] for i in range(len(pts)))
    if area2 < 0.0:
        pts = pts[::-1]
    pan = []
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        L = math.hypot(b[0] - a[0], b[1] - a[1])
        if L <= 1e-15:
            continue
        t = ((b[0] - a[0]) / L, (b[1] - a[1]) / L)
        pan.append((a, L, t, (t[1], -t[0]), (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))))  # n out of a CCW loop
    n = len(pan)
    if n < 3:
        return [[0.0, 0.0], [0.0, 0.0]]

    def F(u, z):
        r2 = u * u + z * z
        return (u * math.log(r2) if r2 > 0.0 else 0.0) - 2.0 * u + (2.0 * z * math.atan(u / z) if z != 0.0 else 0.0)

    A = [[0.0] * n for _ in range(n)]
    P = [[0.0] * n for _ in range(n)]
    for i, (_, Li, _, ni, ci) in enumerate(pan):
        for k, (ak, Lk, tk, nk, _) in enumerate(pan):
            dx, dz = ci[0] - ak[0], ci[1] - ak[1]
            x, z = dx * tk[0] + dz * tk[1], dx * nk[0] + dz * nk[1]
            if i == k:
                A[i][k] = 0.5
                P[i][k] = (Lk * math.log(0.25 * Lk * Lk) - 2.0 * Lk) / (4.0 * math.pi)
                continue
            ux = (math.log(x * x + z * z) - math.log((x - Lk) ** 2 + z * z)) / (4.0 * math.pi)
            w = (math.atan((Lk - x) / z) + math.atan(x / z)) / (2.0 * math.pi) if abs(z) > 1e-15 else 0.0
            vy, vz = ux * tk[0] + w * nk[0], ux * tk[1] + w * nk[1]
            A[i][k] = vy * ni[0] + vz * ni[1]
            P[i][k] = (F(Lk - x, z) - F(-x, z)) / (4.0 * math.pi)
    sig = _solve(A, [[p[3][0] for p in pan], [p[3][1] for p in pan]])
    m = [[0.0, 0.0], [0.0, 0.0]]
    for j in range(2):
        phi = [sum(P[i][k] * sig[j][k] for k in range(n)) for i in range(n)]
        for kk in range(2):
            m[j][kk] = -sum(phi[i] * pan[i][3][kk] * pan[i][1] for i in range(n))
    off = 0.5 * (m[0][1] + m[1][0])
    return [[m[0][0], off], [off, m[1][1]]]


def rings(text, shape=False):
    """Per component: [(x, area, (y, z) centroid, apparent mass)] along the body, caps (rings at the same x) dropped;
    the apparent-mass tensor (``apparent_mass``) only with ``shape``, else the area times the identity."""
    comps = {}
    raw = {}
    for line in text.splitlines():
        p = line.split()
        if len(p) < 9 or p[0] != "NTVSP_RING":
            continue
        v = [float(t) for t in p[3:]]
        pts = [v[i:i + 3] for i in range(0, len(v) - 2, 3)]
        raw.setdefault(p[1], []).append((sum(q[0] for q in pts) / len(pts), pts))
    comps = []
    for rl in raw.values():
        # OpenVSP's end caps are several rings in one plane, sweeping across the cap: partial outlines, not sections.
        # A cap group at the front becomes a point (its area enters as the step to the first section, as for a blunt
        # nose); one further aft (a flat base) is left out, so the body does not seem to shrink before its base.
        groups = []
        for x, pts in rl:
            if groups and abs(x - groups[-1][0][0]) <= 1e-9:
                groups[-1].append((x, pts))
            else:
                groups.append([(x, pts)])
        out = []
        for gi, grp in enumerate(groups):
            if len(grp) > 1:
                if gi == 0:
                    pts = [q for _, ring in grp for q in ring]
                    centre = (sum(q[1] for q in pts) / len(pts), sum(q[2] for q in pts) / len(pts))
                    out.append((grp[0][0], 0.0, centre, [[0.0, 0.0], [0.0, 0.0]]))
                continue
            x, pts = grp[0]
            a2 = sum(pts[i][1] * pts[i - 1][2] - pts[i - 1][1] * pts[i][2] for i in range(len(pts)))
            area = 0.5 * abs(a2)
            if area > 1e-30:
                cy = sum((pts[i][1] + pts[i - 1][1]) * (pts[i][1] * pts[i - 1][2] - pts[i - 1][1] * pts[i][2])
                         for i in range(len(pts))) / (3.0 * a2)
                cz = sum((pts[i][2] + pts[i - 1][2]) * (pts[i][1] * pts[i - 1][2] - pts[i - 1][1] * pts[i][2])
                         for i in range(len(pts))) / (3.0 * a2)
            else:
                cy, cz = (sum(q[1] for q in pts) / len(pts), sum(q[2] for q in pts) / len(pts))
            mass = (apparent_mass([(q[1], q[2]) for q in pts]) if shape and area > 1e-30
                    else [[area, 0.0], [0.0, area]])
            out.append((x, area, (cy, cz), mass))
        comps.append(out)
    return comps


def combine(parts, n):
    """Several fuselages on one x grid: areas and apparent-mass tensors add (no interference between fuselages),
    centroids are area-weighted."""
    x0, x1 = min(p[0][0] for p in parts), max(p[-1][0] for p in parts)
    xs = [x0 + (x1 - x0) * i / n for i in range(n + 1)]
    areas, centres, masses = [], [], []
    for x in xs:
        tot = cy = cz = 0.0
        mm = [0.0, 0.0, 0.0]
        for p in parts:
            if not p[0][0] <= x <= p[-1][0]:
                continue
            px = [q[0] for q in p]
            a = interp(x, px, [q[1] for q in p])
            tot += a
            cy += a * interp(x, px, [q[2][0] for q in p])
            cz += a * interp(x, px, [q[2][1] for q in p])
            for j, (r, c) in enumerate(((0, 0), (0, 1), (1, 1))):
                mm[j] += interp(x, px, [q[3][r][c] for q in p])
        areas.append(tot)
        centres.append((cy / tot, cz / tot) if tot > 0 else (0.0, 0.0))
        masses.append(((mm[0], mm[1]), (mm[1], mm[2])))
    return xs, areas, centres, masses


def body_loads(xs, areas, centres, conds, ref, potential, crossflow, masses=None):
    """Coefficient rows (HEADER) for each condition. With ``masses`` (apparent-mass tensors per station), the potential
    force per length follows d(M e)/dx, e the crossflow direction, instead of (dS/dx) e: the same for round sections,
    and a force off the crossflow where the section is not symmetric about it."""
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
        sy = sz = 0.0  # the force's y and z components (off the crossflow too, with section shapes)
        mx = 0.0
        for i in range(len(xs) - 1):
            dx = xs[i + 1] - xs[i]
            if dx <= 0.0:
                continue
            xm = 0.5 * (xs[i] + xs[i + 1])
            ym = 0.5 * (centres[i][0] + centres[i + 1][0])
            zm = 0.5 * (centres[i][1] + centres[i + 1][1])
            f = 0.0
            fy = fz = 0.0
            if potential:
                k_p = math.sin(2.0 * at) * math.cos(0.5 * at) / sref
                if masses is None:
                    f += (areas[i + 1] - areas[i]) * k_p
                    fy, fz = f * ey, f * ez
                else:
                    (a0, b0), (_, c0) = masses[i]
                    (a1, b1), (_, c1) = masses[i + 1]
                    fy = ((a1 - a0) * ey + (b1 - b0) * ez) * k_p
                    fz = ((b1 - b0) * ey + (c1 - c0) * ez) * k_p
                    f = fy * ey + fz * ez
            if crossflow:
                g = k_v * (radius[i] + radius[i + 1]) * dx / sref * math.sin(at) ** 2
                f += g
                fy, fz = fy + g * ey, fz + g * ez
            cn += f
            sy, sz = sy + fy, sz + fz
            rx, ry, rz = xm - xr, ym - yr, zm - zr
            mx += ry * fz - rz * fy
            my += -rx * fz
            mz += rx * fy
        fy_t, fz_t = (cn * ey, cn * ez) if masses is None else (sy, sz)
        rows.append([c["Mach"], c["AoA"], c["Beta"], 0.0, fy_t, fz_t, mx / bref, my / cref, mz / bref,
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
    mode = s.text("potential", "slender")
    if mode not in ("slender", "shape", "off"):
        raise DriverError(f"potential '{mode}': expected slender, shape or off")
    potential, shape = mode != "off", mode == "shape"
    crossflow = s.text("crossflow", "jorgensen") == "jorgensen"
    if not (potential or crossflow):
        raise DriverError("both the potential lift and the viscous crossflow are off")
    flight = Flight(s.text("flight"))
    rows = _rows(s.text("flight").strip() + ".ntvsp.inviscid.csv")
    lo, hi = s.real("mach_min", 0.0), s.real("mach_max", 1e9)
    conds = [r for r in rows if lo <= r["Mach"] < hi]
    n = s.integer("slices", 200)
    text = run_vspscript(vsp_dir, ring_script(model, names, n, 64), work, case + ".rings", log, s.real("timeout_s", 3600))
    parts = rings(text, shape)
    if not parts or max(max(q[1] for q in p) for p in parts) <= 0.0:
        raise DriverError("the chosen fuselages have no section area")
    xs, areas, centres, masses = combine(parts, n)
    ends = (xs[0], xs[-1])
    out, fineness = body_loads(xs, areas, centres, conds, flight.ref, potential, crossflow, masses if shape else None)
    write_csv(os.path.join(work, case + ".ntvsp.body_sections.csv"), ["x", "Area", "Yc", "Zc", "M_yy", "M_yz", "M_zz"],
              [[x, a, c[0], c[1], m[0][0], m[0][1], m[1][1]] for x, a, c, m in zip(xs, areas, centres, masses)])
    if shape:
        big = max(range(len(xs)), key=lambda i: areas[i])
        log(f"section apparent mass at the largest section (x = {xs[big]:.4f} m): M_yy {masses[big][0][0]:.6g}, M_zz "
            f"{masses[big][1][1]:.6g}, M_yz {masses[big][0][1]:.3g} m^2 against its area {areas[big]:.6g} m^2")
    write_csv(os.path.join(work, case + ".ntvsp.coefficients.csv"), HEADER, out)
    log(f"body lift of {', '.join(names) or 'every fuselage'}: {len(xs)} sections from x = {ends[0]:.4f} to "
        f"{ends[1]:.4f} m, largest area {max(areas):.6g} m^2, fineness {fineness:.2f}; "
        f"{('slender-body potential lift' + (' (section apparent mass)' if shape else '')) if potential else ''}{' + ' if potential and crossflow else ''}"
        f"{'viscous crossflow (Jorgensen)' if crossflow else ''}; {len(conds)} conditions")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
