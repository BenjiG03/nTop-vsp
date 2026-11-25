"""Body (nose) lift of an OpenVSP model's fuselages from their section-area distribution, driven from nTop.

    python -m vsp_ntop.run_body_lift SETTINGS_FILE
"""
from __future__ import annotations

import os

import numpy as np

from .common import DriverError, Log, Settings, main_wrapper, openvsp_dir, run_vspscript, vsp_path, write_csv
from .flight import Flight, _rows
from .geometry import shoelace

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


def crossflow_drag(mn):
    return np.interp(mn, CDN_MACH, CDN_VALUE)


def crossflow_factor(fineness, mn):
    """Allen and Perkins' eta: Figure 4 against fineness, its rise with crossflow Mach (Figure 6, for fineness 10 to
    12) rescaled to start from the Figure 4 value (hypanel's assumption)."""
    low = np.interp(fineness, ETA_LENGTH, ETA_VALUE)
    measured = np.interp(mn, ETA_MACH, ETA_MACH_VALUE)
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


def apparent_mass(pts):
    """The section's apparent-mass tensor (``_panel_mass``), Richardson-extrapolated from the ring's own edges and
    each edge split in two: the constant-source panels converge to first order (a 64-point circle: +2.2% on 64
    panels, +1.1% on 128, -0.02% extrapolated)."""
    pts = np.asarray(pts, float)
    fine = np.stack([pts, pts + 0.5 * (np.roll(pts, -1, axis=0) - pts)], axis=1).reshape(-1, 2)
    return 2.0 * _panel_mass(fine) - _panel_mass(pts)


def _panel_mass(pts):
    """The 2D apparent-mass tensor of a closed section per unit fluid density, [[m_yy, m_yz], [m_zy, m_zz]] (m^2),
    from ``pts`` (y, z) around it: the exterior potential flow of the section moving with unit velocity along y and
    along z (constant-strength source panels, ``d phi / d n = n . e`` on the contour), m_jk = -(contour integral of
    phi_j n_k ds), n out of the body."""
    pts = np.asarray(pts, float)
    if shoelace(pts)[0] < 0.0:
        pts = pts[::-1]
    a, d = pts, np.roll(pts, -1, axis=0) - pts
    L = np.hypot(d[:, 0], d[:, 1])
    keep = L > 1e-15
    a, d, L = a[keep], d[keep], L[keep]
    if len(L) < 3:
        return np.zeros((2, 2))
    t = d / L[:, None]
    nrm = np.column_stack([t[:, 1], -t[:, 0]])  # out of a counter-clockwise loop
    c = a + 0.5 * d
    rel = c[:, None, :] - a[None, :, :]  # collocation point i relative to panel k's start
    x, z = (rel * t[None]).sum(axis=2), (rel * nrm[None]).sum(axis=2)
    Lk = L[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        ux = (np.log(x * x + z * z) - np.log((x - Lk) ** 2 + z * z)) / (4.0 * np.pi)
        w = np.where(np.abs(z) > 1e-15, (np.arctan((Lk - x) / z) + np.arctan(x / z)) / (2.0 * np.pi), 0.0)

        def F(u, z):
            r2 = u * u + z * z
            return (np.where(r2 > 0.0, u * np.log(r2), 0.0) - 2.0 * u
                    + np.where(z != 0.0, 2.0 * z * np.arctan(u / z), 0.0))
        P = (F(Lk - x, z) - F(-x, z)) / (4.0 * np.pi)
    A = ux * (nrm @ t.T) + w * (nrm @ nrm.T)  # normal velocity at i of panel k's unit source
    diag = np.arange(len(L))
    A[diag, diag] = 0.5
    P[diag, diag] = (L * np.log(0.25 * L * L) - 2.0 * L) / (4.0 * np.pi)
    phi = P @ np.linalg.solve(A, nrm)  # columns: motion along y, along z
    m = -(phi * L[:, None]).T @ nrm
    off = 0.5 * (m[0, 1] + m[1, 0])
    return np.array([[m[0, 0], off], [off, m[1, 1]]])


def rings(text, shape=False):
    """Per component: [(x, area, (y, z) centroid, apparent mass)] along the body, caps (rings at the same x) dropped;
    the apparent-mass tensor (``apparent_mass``) only with ``shape``, else the area times the identity."""
    raw = {}
    for line in text.splitlines():
        p = line.split()
        if len(p) < 9 or p[0] != "NTVSP_RING":
            continue
        v = np.array(p[3:3 + 3 * ((len(p) - 3) // 3)], float).reshape(-1, 3)
        raw.setdefault(p[1], []).append((v[:, 0].mean(), v[:, 1:]))
    comps = []
    for rl in raw.values():
        # OpenVSP's end caps are several rings in one plane, sweeping across the cap
        groups = []
        for x, yz in rl:
            if groups and abs(x - groups[-1][0][0]) <= 1e-9:
                groups[-1].append((x, yz))
            else:
                groups.append([(x, yz)])
        out = []
        for gi, grp in enumerate(groups):
            if len(grp) > 1:
                if gi == 0:
                    out.append((grp[0][0], 0.0, np.vstack([yz for _, yz in grp]).mean(axis=0), np.zeros((2, 2))))
                continue
            x, yz = grp[0]
            a, c = shoelace(yz)
            area = abs(a)
            if area <= 1e-30:
                c = yz.mean(axis=0)
            mass = apparent_mass(yz) if shape and area > 1e-30 else area * np.eye(2)
            out.append((x, area, c, mass))
        comps.append(out)
    return comps


def combine(parts, n):
    """Several fuselages on one x grid: areas and apparent-mass tensors add (no interference between fuselages),
    centroids are area-weighted."""
    x0, x1 = min(p[0][0] for p in parts), max(p[-1][0] for p in parts)
    xs = x0 + (x1 - x0) * np.arange(n + 1) / n
    areas, first = np.zeros(n + 1), np.zeros((n + 1, 2))
    masses = np.zeros((n + 1, 2, 2))
    for p in parts:
        px = np.array([q[0] for q in p])
        on = (px[0] <= xs) & (xs <= px[-1])
        a = np.interp(xs, px, [q[1] for q in p])
        areas += np.where(on, a, 0.0)
        for k in range(2):
            first[:, k] += np.where(on, a * np.interp(xs, px, [q[2][k] for q in p]), 0.0)
        for r, c in ((0, 0), (0, 1), (1, 1)):
            masses[:, r, c] += np.where(on, np.interp(xs, px, [q[3][r][c] for q in p]), 0.0)
    masses[:, 1, 0] = masses[:, 0, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        centres = np.where(areas[:, None] > 0, first / areas[:, None], 0.0)
    return xs, areas, centres, masses


def body_loads(xs, areas, centres, conds, ref, potential, crossflow, masses=None):
    """Coefficient rows (HEADER) for each condition."""
    sref, bref, cref, xyz_ref = ref
    xs, areas, centres = np.asarray(xs), np.asarray(areas), np.asarray(centres)
    radius = np.sqrt(np.maximum(areas, 0.0) / np.pi)
    fineness = (xs[-1] - xs[0]) / max(2.0 * radius.max(), 1e-12)
    dx = np.diff(xs)
    ok = dx > 0.0
    arm = np.column_stack([0.5 * (xs[1:] + xs[:-1]), 0.5 * (centres[1:] + centres[:-1])]) - xyz_ref
    rows = []
    for c in conds:
        a, b = np.radians(c["AoA"]), np.radians(c["Beta"])
        vy, vz = -np.sin(b), np.sin(a) * np.cos(b)
        s = np.hypot(vy, vz)
        at = np.arctan2(s, np.cos(a) * np.cos(b))
        e = np.array([vy, vz]) / s if s > 1e-12 else np.zeros(2)  # the crossflow direction in (y, z)
        mn = c["Mach"] * np.sin(at)
        fyz = np.zeros((len(dx), 2))  # each section's force (y, z)
        if potential:
            k_p = np.sin(2.0 * at) * np.cos(0.5 * at) / sref
            if masses is None:
                fyz += np.outer(np.diff(areas) * k_p, e)
            else:
                fyz += np.diff(masses, axis=0) @ e * k_p
        if crossflow:
            k_v = crossflow_factor(fineness, mn) * crossflow_drag(mn)
            fyz += np.outer(k_v * (radius[:-1] + radius[1:]) * dx / sref * np.sin(at) ** 2, e)
        fyz[~ok] = 0.0
        force = np.column_stack([np.zeros(len(dx)), fyz])
        mx, my, mz = np.cross(arm, force).sum(axis=0)
        fy_t, fz_t = (fyz @ e).sum() * e if masses is None else fyz.sum(axis=0)
        rows.append([c["Mach"], c["AoA"], c["Beta"], 0.0, fy_t, fz_t, mx / bref, my / cref, mz / bref,
                     sref, bref, cref, *xyz_ref])
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
        big = int(np.argmax(areas))
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
