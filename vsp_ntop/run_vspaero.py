"""VSPAERO analysis of an OpenVSP model, driven from nTop.

    python -m vsp_ntop.run_vspaero SETTINGS_FILE
"""
from __future__ import annotations

import csv
import os
import re

import numpy as np

from . import adb
from . import pressure
from . import vspgeom
from .common import (openvsp_dir, DriverError, Log, Settings, fmt, main_wrapper, run, run_vspscript, vsp_path,
                     write_csv, write_point_map)
from .flight import Flight
from .geometry import ls_slope, wind_axes
from .vspmodel import MIRRORS

# Solver settings written to <case>.vspaero, in OpenVSP's order
SOLVER_KEYS = [
    ("FreezeMultiPoleAtIteration", "freeze_multipole_at_iteration", int),
    ("FreezeWakeAtIteration", "freeze_wake_at_iteration", int),
    ("FreezeWakeRootVortices", "freeze_wake_root_vortices", int),
    ("ImplicitWake", "implicit_wake", int),
    ("ImplicitWakeStartIteration", "implicit_wake_start_iteration", int),
]
SOLVER_KEYS_2 = [
    ("WakeRelax", "wake_relax", float),
    ("ForwardGMRESConvergenceFactor", "forward_gmres_factor", float),
    ("AdjointGMRESConvergenceFactor", "adjoint_gmres_factor", float),
    ("NonLinearConvergenceFactor", "nonlinear_factor", float),
    ("CoreSizeFactor", "core_size_factor", float),
    ("FarAway", "far_away", float),
    ("UpdateMatrixPreconditioner", "update_matrix_preconditioner", int),
    ("UseWakeNodeMatrixPreconditioner", "wake_node_matrix_preconditioner", int),
]


def geometry_script(model, vsp3, find_body_wakes):
    return f"""void main()
{{
    ReadVSPFile( "{vsp_path(model)}" );
    if ( FindGeoms().size() == 0 ) {{ Print( "ERR the model has no components" ); return; }}
    SetVSP3FileName( "{vsp_path(vsp3)}" );
    Update();
    int thick = GetGeomSet( "Thick" ).size() > 0 ? 3 : -1;
    int thin = GetGeomSet( "Thin" ).size() > 0 ? 4 : -1;
    if ( thick < 0 && thin < 0 ) {{ Print( "ERR no component is in the Thick or Thin set" ); return; }}
    Print( "NTVSP_SETS " + thick + " " + thin );
    SetIntAnalysisInput( "VSPAEROComputeGeometry", "GeomSet", array<int> = {{thick}} );
    SetIntAnalysisInput( "VSPAEROComputeGeometry", "ThinGeomSet", array<int> = {{thin}} );
    SetIntAnalysisInput( "VSPAEROComputeGeometry", "Symmetry", array<int> = {{0}} );
    SetIntAnalysisInput( "VSPAEROComputeGeometry", "FindBodyWakesFlag", array<int> = {{{1 if find_body_wakes else 0}}} );
    string rid = ExecAnalysis( "VSPAEROComputeGeometry" );
    PrintResults( rid );
    Print( "NTVSP_OK" );
}}
"""


def vspaero_file(s: Settings, machs, alphas, betas, recrefs):
    fixed_wake = s.flag("fixed_wake")
    lines = [
        ("Sref", fmt(s.real("sref"))), ("Cref", fmt(s.real("cref"))), ("Bref", fmt(s.real("bref"))),
        ("X_cg", fmt(s.real("xcg"))), ("Y_cg", fmt(s.real("ycg"))), ("Z_cg", fmt(s.real("zcg"))),
        ("Mach", ", ".join(fmt(v) for v in machs)), ("AoA", ", ".join(fmt(v) for v in alphas)),
        ("Beta", ", ".join(fmt(v) for v in betas)), ("ReCref", ", ".join(fmt(v) for v in recrefs)),
        ("Vinf", fmt(s.real("vinf"))),
    ]
    if s.flag("manual_vref"):
        lines += [("Vref", fmt(s.real("vref"))), ("Mref", fmt(s.real("machref")))]
    lines += [
        ("Rho", fmt(s.real("rho"))),
        # Stall model: not exposed (out of scope); OpenVSP's defaults.
        ("StallModel", "0"), ("Clo2D", "0"), ("CLMax2D", "1"),
        ("Symmetry", "0"),  # X-Z symmetry: out of scope
    ]
    lines += [(k, fmt(t(s.real(n)))) for k, n, t in SOLVER_KEYS]
    lines += [
        ("FarDist", fmt(s.real("far_dist")) if s.flag("far_dist_on") else "-1"),
        ("NumWakeNodes", fmt(s.integer("wake_nodes"))),
        ("WakeIters", "0" if fixed_wake else fmt(s.integer("wake_iterations"))),
    ]
    lines += [(k, fmt(t(s.real(n)))) for k, n, t in SOLVER_KEYS_2]
    lines += [
        ("Write2DFEMFile", "0"),  # out of scope
        ("WriteTecplotFile", "1"),  # used internally to check the surface pressures read from the .adb
        ("NumberOfControlGroups", "0"),  # control surfaces: out of scope
        ("NumberofSurveyPoints", "0"),
        ("QuadTreeBufferLevels", fmt(s.integer("quadtree_buffer_levels"))),
        ("NumberOfQuadTrees", "0"),  # slicer: out of scope
        ("NumberOfInlets", "0"), ("NumberOfNozzles", "0"),
        ("VSP_StabilityType", "0"),  # stability type: out of scope
    ]
    return "".join(f"{k} = {v} \n" for k, v in lines)


def read_polar(path):
    """VSPAERO's .polar: two banner lines, a header line, then one row per case."""
    with open(path) as f:
        rows = [l.split() for l in f.read().split("\n")]
    header = None
    data = []
    for r in rows:
        if not r:
            continue
        if r[0] == "Beta" and header is None:
            header = r
        elif header is not None and len(r) >= len(header):
            # NaN rows are kept (as nan) for the convergence check.
            values = [_float(v) for v in r[:len(header)]]
            if not np.isnan(values[:3]).all():
                data.append(values)
    if header is None:
        raise DriverError(f"{path}: no header line")
    return header, data


GMRES_RECORD = re.compile(r"Wake Iter:\s*(\d+)\s*/\s*(\d+)\s*\.\.\.\s*GMRES Iter:\s*(\d+)\s*\.\.\.\s*Red:\s*(\S+)\s*/\s*(\S+)"
                          r"\s*\.\.\.\s*Max:\s*(\S+)\s*/\s*(\S+)")
SOLVING = re.compile(r"Solving\.\.\.\s*Mach:\s*(\S+)\s*\.\.\.\s*Alpha:\s*(\S+)\s*\.\.\.\s*Beta:\s*(\S+)")
# Floors of the wake-change measure
WAKE_FLOORS = (("CLtot", 0.05), ("CDtot", 0.005), ("CMytot", 0.01))


def _float(x):
    try:
        return float(x)
    except ValueError:
        return np.nan


def _same(a, b):
    # VSPAERO solves Mach 0 as 0.001.
    return all(abs(max(x, 0.001) - max(y, 0.001)) < 1e-4 if i == 0 else abs(x - y) < 1e-4
               for i, (x, y) in enumerate(zip(a, b)))


def solver_records(log_path):
    """Per solve (Mach, alpha, beta) in VSPAERO's log: the last GMRES record of each wake iteration, as
    (wake iteration, GMRES iterations, residual reduction, its target, max residual, its target), logs base 10."""
    text = open(log_path, encoding="utf-8", errors="replace").read()
    marks = list(SOLVING.finditer(text))
    out = []
    for k, m in enumerate(marks):
        seg = text[m.end():marks[k + 1].start() if k + 1 < len(marks) else len(text)]
        last = {}
        for r in GMRES_RECORD.findall(seg):
            last[int(r[0])] = (int(r[0]), int(r[2])) + tuple(_float(v) for v in r[3:7])
        out.append((tuple(_float(v) for v in m.groups()), [last[w] for w in sorted(last)]))
    return out


def wake_changes(history_path):
    """Per solver case of VSPAERO's .history: ((Mach, alpha, beta), number of wake iterations, change of the forces
    over the last wake iteration: max over CLtot, CDtot, CMytot of |change| / max(|value|, floor))."""
    cases, header = [], None
    if not os.path.isfile(history_path):
        return cases
    for line in open(history_path, encoding="utf-8", errors="replace"):
        w = line.split()
        if not w:
            continue
        if w[0] == "Iter":
            header = w
            cases.append([])
        elif header and cases and w[0].isdigit() and len(w) >= len(header) - 2:
            cases[-1].append(dict(zip(header, (_float(v) for v in w))))
    out = []
    for rows in cases:
        if not rows:
            continue
        cond = (rows[-1]["Mach"], rows[-1]["AoA"], rows[-1]["Beta"])
        change = 0.0
        if len(rows) > 1:
            a, b = rows[-2], rows[-1]
            change = max(abs(b[k] - a[k]) / max(abs(b[k]), floor) for k, floor in WAKE_FLOORS)
        out.append((cond, len(rows), change))
    return out


def convergence(work, case, header, data, s: Settings, log):
    """Checks every case: VSPAERO's linear solve (GMRES) of the last wake iteration reached both its residual targets
    (VSPAERO stops it after 1500 iterations otherwise), every value is finite, and with a relaxed wake the forces
    changed by less than ``wake_tolerance`` over the last wake iteration."""
    solves = solver_records(os.path.join(work, case + ".vspaero.log"))
    wakes = wake_changes(os.path.join(work, case + ".history"))
    tol = s.real("wake_tolerance", 0.1)
    i_m, i_a, i_b = header.index("Mach"), header.index("AoA"), header.index("Beta")
    # Force and moment coefficients only
    coefs = [i for i, h in enumerate(header) if h.startswith("C")]
    rows, flags, failures = [], [], []
    for n, row in enumerate(data, start=1):
        cond = (row[i_m], row[i_a], row[i_b])
        rec = next((r for c, r in solves if r and _same(c, cond)), None)
        wake = next(((it, ch) for c, it, ch in wakes if _same(c, cond)), (0, np.nan))
        why = []
        if not np.isfinite([row[i] for i in coefs]).all():
            why.append("non-finite results")
        if rec is None:
            why.append("no solver record in the log")
            final = (0, 0) + (np.nan,) * 4
        else:
            final = rec[-1]
            _, its, red, red_t, mx, mx_t = final
            if not (red <= red_t and mx <= mx_t):
                why.append(f"linear solve stopped at {its} iterations with residual reduction {red:.2f} (target "
                           f"{red_t:.2f}) and max residual {mx:.2f} (target {mx_t:.2f}), logs base 10")
        if wake[0] > 1 and not wake[1] <= tol:
            why.append(f"forces changed {100 * wake[1]:.1f}% over the last wake iteration (tolerance {100 * tol:g}%)")
        ok = not why
        flags.append(int(ok))
        rows.append([n, *cond, final[1], final[2], final[3], final[4], final[5], wake[0], wake[1], int(ok),
                     "; ".join(why)])
        if not ok:
            failures.append(f"case {n} (Mach {cond[0]:g}, alpha {cond[1]:g}, beta {cond[2]:g}): {why[0]}")
    write_csv(os.path.join(work, case + ".ntvsp.convergence.csv"),
              ["Case", "Mach", "AoA", "Beta", "GMRES_Iterations", "Residual_Reduction", "Reduction_Target",
               "Max_Residual", "Max_Residual_Target", "Wake_Iterations", "Wake_Change", "Converged", "Reason"], rows)
    log(f"convergence: {sum(flags)} of {len(flags)} cases converged")
    msg = None
    if failures:
        msg = (f"VSPAERO did not converge in {len(failures)} of {len(flags)} cases: " + "; ".join(failures[:3])
               + ("; ..." if len(failures) > 3 else "") + f". See {case}.ntvsp.convergence.csv and {case}.vspaero.log.")
        log(("error: " if s.integer("convergence_check", 1) == 1 else "warning: ") + msg)
        if s.integer("convergence_check", 1) != 1:
            msg = None
    return flags, msg


def clean(name):
    return name.replace("/", "_over_").replace("Re_over_1e6", "Re_1e6")


PRESSURE_POINTS = ("off", "centroids", "nodes", "both")


def surface_outputs(work, case, g: vspgeom.VSPGeom, s: Settings, header, data, log, flight=None):
    """The panels' points and areas (``<case>.ntvsp.faces.{1,4}.csv``, thick components only) for the viscous
    correction, and, with ``export_pressure`` (off, centroids, nodes or both), surface pressure point maps per flight
    condition case (``pressure.py``): p = p_inf + q Cp on thick components."""
    tris = vspgeom.level0(g)
    thin = vspgeom.thin_parts(os.path.join(work, case + ".vkey"))
    faces = np.array([t for t, _ in tris], int)
    is_thin = np.array([g.parts[i] in thin for _, i in tris], bool)
    nf = len(faces)
    corners = g.nodes[faces]  # (face, corner, xyz)
    centroids = corners.mean(axis=1)
    areas = 0.5 * np.linalg.norm(np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]), axis=1)
    thick = ~is_thin
    if thick.any():
        write_point_map(os.path.join(work, case + ".ntvsp.faces.1.csv"), centroids[thick], areas[thick])
        p, q, r = corners[thick].transpose(1, 0, 2)
        pq, qr, rp = 0.5 * (p + q), 0.5 * (q + r), 0.5 * (r + p)
        four = np.stack([(p + pq + rp) / 3.0, (pq + q + qr) / 3.0, (rp + qr + r) / 3.0, (pq + qr + rp) / 3.0], axis=1)
        write_point_map(os.path.join(work, case + ".ntvsp.faces.4.csv"), four.reshape(-1, 3),
                        np.repeat(areas[thick] / 4.0, 4))
    mode = s.text("export_pressure", "off").strip().lower()
    if mode not in PRESSURE_POINTS:
        raise DriverError(f"export_pressure = {mode}: expected one of {', '.join(PRESSURE_POINTS)}")
    base = os.path.join(work, case)
    removed = pressure.clear(base)
    if mode == "off":
        if removed:
            log(f"removed {removed} pressure files of an earlier run (Export Pressures is off)")
        return
    if flight is None:
        raise DriverError("pressure export needs the flight condition")
    cases = adb.read_cases(os.path.join(work, case + ".adb.cases"))
    try:
        blocks, n_tri = adb.read_cp(os.path.join(work, case + ".adb"), cases, len(tris))
        err = adb.check_against_tecplot(blocks, os.path.join(work, case + ".tec.dat"))
    except ValueError as e:
        raise DriverError(f"surface pressures: {e}")
    log(f"surface pressures: {len(cases)} conditions, {n_tri} loops with the wakes (check against Tecplot output: "
        f"{err:.1e})")
    i_mach, i_aoa, i_beta = header.index("Mach"), header.index("AoA"), header.index("Beta")
    i_re = header.index("Re/1e6")
    kinds = {"pressure": np.flatnonzero(thick), "pressure_jump": np.flatnonzero(is_thin)}

    def solution(row):
        # VSPAERO solves Mach 0 as Mach 0.001 and lists it so in .adb.cases.
        cond = (max(row[i_mach], 0.001), row[i_aoa], row[i_beta])
        match = [k for k, c in enumerate(cases) if all(abs(a - b) < 1e-5 for a, b in zip(c, cond))]
        if not match:
            raise DriverError(f"no surface solution for Mach {cond[0]}, alpha {cond[1]}, beta {cond[2]}")
        return np.asarray(blocks[match[0]][:nf])

    written = []
    for fc in flight.conditions:
        rows = [r for r in data if abs(r[i_mach] - fc["Mach"]) < 1e-6 and abs(r[i_aoa] - fc["AoA"]) < 1e-6
                and abs(r[i_beta] - fc["Beta"]) < 1e-6]
        if not rows:
            continue  # outside the Mach numbers this run solved (Mach Min / Mach Max)
        row = min(rows, key=lambda r: abs(np.log(max(r[i_re] * 1e6, 1.0)) - np.log(max(fc["Re_cref"], 1.0))))
        cp = solution(row)
        q, p_inf = fc["Dynamic_Pressure"], fc["Pressure"]
        for kind, sel in kinds.items():
            if not len(sel):
                continue
            offset = p_inf if kind == "pressure" else 0.0  # the jump is q Delta Cp
            pts, vals = [], []
            if mode in ("centroids", "both"):
                pts.append(centroids[sel])
                vals.append(offset + q * cp[sel])
            if mode in ("nodes", "both"):  # area-weighted average of the faces around each node
                at = faces[sel].ravel()
                used = np.unique(at)
                num = np.bincount(at, np.repeat(areas[sel] * cp[sel], 3), len(g.nodes))
                den = np.bincount(at, np.repeat(areas[sel], 3), len(g.nodes))
                pts.append(g.nodes[used])
                vals.append(offset + q * num[used] / den[used])
            pressure.write(base, fc["Case"], np.vstack(pts), np.concatenate(vals), kind)
        written.append(int(fc["Case"]))
    which = " and ".join(k for k, sel in kinds.items() if len(sel))
    where = {"centroids": "centroids", "nodes": "nodes", "both": "centroids and nodes"}[mode]
    log(f"{which} (Pa) at the panels' {where} of flight condition cases {', '.join(map(str, written))}")


def pressure_centers(work, case, header, data, s: Settings, log):
    """Centre of pressure of every case and aerodynamic centre of every (Mach, beta, Re_cref) alpha sweep."""
    i = {h: header.index(h) for h in ("Mach", "AoA", "Beta", "Re/1e6", "CFiz", "CMiy")}
    xr, yr, zr, cref = s.real("xcg"), s.real("ycg"), s.real("zcg"), s.real("cref")
    rows = []
    for k, row in enumerate(data, start=1):
        cfz, cm = row[i["CFiz"]], row[i["CMiy"]]
        ok = abs(cfz) >= 1e-6
        rows.append([k, row[i["Mach"]], row[i["AoA"]], row[i["Beta"]], xr - cm * cref / cfz if ok else xr, yr, zr,
                     int(ok), cfz, cm])
    write_csv(os.path.join(work, case + ".ntvsp.center_of_pressure.csv"),
              ["Case", "Mach", "AoA", "Beta", "X", "Y", "Z", "Defined", "CFiz", "CMiy"], rows)
    groups = {}
    for row in data:
        groups.setdefault((row[i["Mach"]], row[i["Beta"]], row[i["Re/1e6"]]), []).append(row)
    ac = []
    for (m, b, re), rows_g in groups.items():
        slope = ls_slope([r[i["CFiz"]] for r in rows_g], [r[i["CMiy"]] for r in rows_g])
        if slope is None:
            continue
        ac.append([m, b, re * 1e6, xr - slope * cref, yr, zr, slope])
    write_csv(os.path.join(work, case + ".ntvsp.aerodynamic_center.csv"),
              ["Mach", "Beta", "Re_cref", "X", "Y", "Z", "dCMy_dCFz"], ac)
    if ac:
        log("aerodynamic centre x: " + ", ".join(f"Mach {r[0]:g} beta {r[1]:g}: {r[3]:.3f} m" for r in ac))


# --- Force Accounting: Exposed Strip

LOD_FORCES = {  # polar column: .lod strip column (0-based) of the same quantity
    "CLo": 14, "CDo": 15, "CSo": 16, "CLi": 17, "CDi": 18, "CSi": 19, "CLtot": 11, "CDtot": 12, "CStot": 13,
    "CFox": 23, "CFoy": 24, "CFoz": 25, "CFix": 26, "CFiy": 27, "CFiz": 28, "CFxtot": 20, "CFytot": 21, "CFztot": 22}
LOD_MOMENTS = {  # polar column: (.lod column, reference length key): per strip area and chord
    "CMox": (32, "bref"), "CMoy": (33, "cref"), "CMoz": (34, "bref"), "CMix": (35, "bref"), "CMiy": (36, "cref"),
    "CMiz": (37, "bref"), "CMxtot": (29, "bref"), "CMytot": (30, "cref"), "CMztot": (31, "bref")}


def read_lod(path):
    """VSPAERO's spanwise loads: {(Mach, AoA, Beta): strips}, each strip the row's numbers (Xavg, Yavg, Zavg at 3-5,
    chord 8, area 9, then the coefficients). VSPAERO writes one set per (Mach, AoA, Beta), whatever the number of
    Reynolds numbers."""
    cases, cur, key = {}, None, {}
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] in ("Mach_", "AoA_"):
                key[p[0]] = float(p[1])
            elif p[0] == "Beta_":
                cur = []
                cases[(round(key["Mach_"], 6), round(key["AoA_"], 6), round(float(p[1]), 6))] = cur
            elif cur is not None and p[0].isdigit() and len(p) > 38:
                try:
                    cur.append([float(x) for x in p[:39]])
                except ValueError:
                    pass
    return {k: np.array(v, float).reshape(-1, 39) for k, v in cases.items()}


def exposed_wings(model_vsp3):
    """The model's wings with Force Accounting Exposed Strip (from OpenVSP Model's ``.ntvsp.wings.csv``), each as a
    list of (span axis, plane normal, root point, root span) for it and its mirror copies."""
    path = model_vsp3[:-5] + ".ntvsp.wings.csv"
    if not os.path.isfile(path):
        return []
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["Accounting"] not in ("exposed", "carryover"):
                continue
            s, n, q = (np.array([float(r[k + c]) for c in "xyz"]) for k in ("S", "N", "Q"))
            out.append((r["Name"], [(m * s, m * n, m * q, float(r["Root_Span"])) for m in MIRRORS[r["Symmetry"]]]))
    return out


def inboard(strips, wings):
    """Which strips lie on an Exposed Strip wing inboard of its root: between the origin plane and the root, in the
    wing's own plane."""
    p = strips[:, 3:6]
    out = np.zeros(len(strips), bool)
    for _, copies in wings:
        for s, n, q, root in copies:
            along, off = p @ s, np.abs((p - q) @ n)
            out |= (-1e-6 <= along) & (along < root * (1.0 - 1e-6)) & (off <= 0.05 * max(root, 1e-6) + 1e-6)
    return out


def exposed_accounting(header, data, lod_path, wings, s: Settings, log):
    """Take the inboard strips' loads off every case: VSPAERO solves the whole (extended) wing, the forces count only
    where the wing is exposed."""
    cases = read_lod(lod_path)
    sref, ref = s.real("sref"), {"bref": s.real("bref"), "cref": s.real("cref")}
    col = {h: i for i, h in enumerate(header)}
    shares = []
    for row in data:
        k = tuple(round(row[col[c]], 6) for c in ("Mach", "AoA", "Beta"))
        strips = cases.get(k)
        if strips is None:
            raise DriverError(f"Force Accounting: no spanwise loads for Mach {k[0]:g}, alpha {k[1]:g}, beta {k[2]:g} "
                              f"in {lod_path}")
        inner = strips[inboard(strips, wings)]
        lift_all = strips[:, 17] @ strips[:, 9]
        if abs(lift_all) > 1e-3 * strips[:, 9].sum():  # skip conditions without lift
            shares.append(inner[:, 17] @ inner[:, 9] / lift_all)
        for name, j in LOD_FORCES.items():
            if name in col:
                row[col[name]] -= float(inner[:, j] @ inner[:, 9]) / sref
        for name, (j, length) in LOD_MOMENTS.items():
            if name in col:
                row[col[name]] -= float(inner[:, j] @ (inner[:, 9] * inner[:, 8])) / (sref * ref[length])
        if "L/D" in col and abs(row[col["CDtot"]]) > 1e-12:
            row[col["L/D"]] = row[col["CLtot"]] / row[col["CDtot"]]
    first = next(iter(cases.values()), np.empty((0, 39)))
    n_in = int(inboard(first, wings).sum())
    log(f"Force Accounting Exposed Strip ({', '.join(w for w, _ in wings)}): {n_in} inboard strips left out; they "
        + (f"carried {100 * min(shares):.0f} to {100 * max(shares):.0f}% of the inviscid lift" if shares else ""))


def read_carryover(model_vsp3):
    """Each copy of the wings with Force Accounting Exposed Strip + Carry-over, as OpenVSP Model measured it
    (``.ntvsp.carryover.csv``)."""
    path = model_vsp3[:-5] + ".ntvsp.carryover.csv"
    if not os.path.isfile(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return [{k: (v if k in ("Name", "Body") else float(v)) for k, v in r.items()} for r in csv.DictReader(f)]


def carryover_accounting(header, data, lod_path, fins, s: Settings, log):
    """Add the lift each fin carries onto the body (TR 1307 K_B(W), ``carryover.py``) to every case: K_B(W) times
    the fin's own force along its normal (its exposed strips in VSPAERO's spanwise loads), along that normal, at the
    planar model's centre of the carried load."""
    from . import carryover as CO
    cases = read_lod(lod_path)
    sref, ref = s.real("sref"), {"bref": s.real("bref"), "cref": s.real("cref")}
    xyz = np.array([s.real("xcg"), s.real("ycg"), s.real("zcg")])
    col = {h: i for i, h in enumerate(header)}
    factors, notes = {}, {}
    for row in data:
        mach, a, b = (row[col[c]] for c in ("Mach", "AoA", "Beta"))
        strips = cases[tuple(round(v, 6) for v in (mach, a, b))]
        p = strips[:, 3:6]
        force, moment = np.zeros(3), np.zeros(3)
        for f in fins:
            sv, nv, q = (np.array([f[k + c] for c in "xyz"]) for k in ("S", "N", "Q"))
            along, off = p @ sv, np.abs((p - q) @ nv)
            on = (f["Root"] * (1.0 - 1e-6) <= along) & (along <= f["Semispan"] * 1.02) & (off <= 0.05 * f["Semispan"])
            fin_n = (strips[on, 20:23] @ nv) @ strips[on, 9] / sref
            key = (f["Name"], f["Copy"], round(mach, 6))
            if key not in factors:
                factors[key] = CO.body_factor(mach, f["Root"], f["Semispan"], f["Root_Chord"], f["Tip_Chord"],
                                              f["LE_Slope"], f["Planform"], f["Body_Width"], f["Body_Length"])
            k, xi, eta, method = factors[key]
            notes[(f["Name"], round(mach, 6))] = (k, method)
            d_force = k * fin_n * nv
            force += d_force
            moment += np.cross(q + [xi, 0.0, 0.0] - eta * sv - xyz, d_force)
        v, lift, side = wind_axes(a, b)
        add = {"L": lift @ force, "D": v @ force, "S": side @ force, "Fx": force[0], "Fy": force[1], "Fz": force[2],
               "Mx": moment[0] / ref["bref"], "My": moment[1] / ref["cref"], "Mz": moment[2] / ref["bref"]}
        for q_, v in add.items():
            for name in (f"C{q_}i", f"C{q_}tot") if len(q_) == 1 else (f"C{q_[0]}i{q_[1]}", f"C{q_}tot"):
                if name in col:
                    row[col[name]] += float(v)
        if "L/D" in col and abs(row[col["CDtot"]]) > 1e-12:
            row[col["L/D"]] = row[col["CLtot"]] / row[col["CDtot"]]
    for (name, mach), (k, method) in sorted(notes.items()):
        log(f"carry-over (TR 1307): {name} at Mach {mach:g}: K_B(W) {k:.3f} ({method})")


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = openvsp_dir(s)
    timeout = s.real("timeout_s", 7200)
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    for stale in os.listdir(work):
        if stale.startswith(case + ".ntvsp.") and stale.endswith(".csv"):
            os.remove(os.path.join(work, stale))
    flight = Flight(s.text("flight"))
    # The reference
    sref, bref, cref, (xr, yr, zr) = flight.ref
    first = flight.points[0]
    s.raw.update(sref=repr(sref), bref=repr(bref), cref=repr(cref), xcg=repr(xr), ycg=repr(yr), zcg=repr(zr),
                 vinf=repr(first["Velocity"]), rho=repr(first["Density"]), manual_vref="0")
    log("flight condition: " + flight.describe())
    model = s.text("model")
    if not model.lower().endswith(".vsp3"):
        model += ".vsp3"
    if not os.path.isfile(model):
        raise DriverError(f"OpenVSP model {model} not found")

    # 1. Panel geometry from the model's Thick and Thin sets.
    log("OpenVSP: read the model and write the VSPAERO panel file")
    text = run_vspscript(vsp_dir, geometry_script(model, os.path.join(work, case + ".vsp3"), s.flag("find_body_wakes")),
                         work, case + ".geometry", log, timeout)
    sets = next((l.split()[1:] for l in text.splitlines() if l.strip().startswith("NTVSP_SETS")), ["?", "?"])
    log(f"Thick set {'used' if sets[0] == '3' else 'empty'}, Thin set {'used' if sets[1] == '4' else 'empty'}")
    geom = os.path.join(work, case + ".vspgeom")
    if not os.path.isfile(geom):
        raise DriverError("OpenVSP wrote no .vspgeom file")
    g = vspgeom.read(geom)
    log(f"panels: {len(g.faces)} faces, {len(g.nodes)} nodes, {len(g.wakes)} trailing-edge wakes")
    if not g.wakes:
        log("warning: no trailing-edge wakes (no Wing component); VSPAERO sheds no vorticity and lift will be small")

    # 2. Solver input.
    lo, hi = s.real("mach_min", 0.0), s.real("mach_max", 1e9)
    points = [p for p in flight.points if lo - 1e-9 <= p["Mach"] < hi - 1e-9]
    if not points:
        raise DriverError(f"the flight condition has no Mach number from {lo:g} up to {hi:g} (Mach Min and Max)")
    machs = [m for m in flight.machs() if any(p["Mach"] == m for p in points)]
    alphas, betas = flight.angles("AoA"), flight.angles("Beta")
    recrefs = sorted({round(p["Re_cref"], 3) for p in points})
    with open(os.path.join(work, case + ".vspaero"), "w", newline="\n") as f:
        f.write(vspaero_file(s, machs, alphas, betas, recrefs))
    n_cases = len(machs) * len(alphas) * len(betas) * len(recrefs)
    log(f"VSPAERO: {len(machs)} Mach x {len(alphas)} alpha x {len(betas)} beta x {len(recrefs)} Re_cref = {n_cases} cases")

    # 3. Solve.
    for old in (".polar", ".adb", ".adb.cases", ".tec.dat", ".history", ".lod"):
        p = os.path.join(work, case + old)
        if os.path.exists(p):
            os.remove(p)
    exe = os.path.join(vsp_dir, "vspaero.exe")
    cmd = [exe, "-omp", str(max(1, s.integer("ncpu")))]
    if s.flag("ground_effect_on"):
        cmd += ["-groundheight", fmt(s.real("ground_height"))]
    restart = s.integer("restart_mode", 0)  # 0 none, 1 save, 2 restart from saved
    if restart == 1:
        cmd += ["-save"]
    elif restart == 2:
        cmd += ["-restart"]
    cmd += [case]
    run(cmd, log, work, os.path.join(work, case + ".vspaero.log"), timeout)

    # 4. Results.
    polar = os.path.join(work, case + ".polar")
    if not os.path.isfile(polar):
        raise DriverError(f"VSPAERO wrote no .polar file; see {case}.vspaero.log")
    header, data = read_polar(polar)
    if len(data) != n_cases:
        raise DriverError(f"VSPAERO's polar has {len(data)} rows for {n_cases} cases; see {case}.vspaero.log")
    flags, failed = convergence(work, case, header, data, s, log)
    wings = exposed_wings(model)
    if wings:
        write_csv(os.path.join(work, case + ".ntvsp.polar_whole.csv"),
                  ["Case"] + [clean(h) for h in header] + ["Converged"],
                  [[i] + list(row) + [f] for i, (row, f) in enumerate(zip(data, flags), start=1)])
        exposed_accounting(header, data, os.path.join(work, case + ".lod"), wings, s, log)
        fins = read_carryover(model)
        if fins:
            carryover_accounting(header, data, os.path.join(work, case + ".lod"), fins, s, log)
    write_csv(os.path.join(work, case + ".ntvsp.polar.csv"), ["Case"] + [clean(h) for h in header] + ["Converged"],
              [[i] + row + [f] for i, (row, f) in enumerate(zip(data, flags), start=1)])
    if failed:
        raise DriverError(failed)
    pressure_centers(work, case, header, data, s, log)
    surface_outputs(work, case, g, s, header, data, log, flight)
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
