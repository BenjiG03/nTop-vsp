"""VSPAERO analysis of a tessellated surface (panel method or vortex lattice), driven from nTop.

    python run_vspaero.py SETTINGS_FILE

Steps (each logged to ``<case>.ntvsp.log``):

1. OpenVSP imports the surface mesh (``<case>.stl``, metres) as a MeshGeom and writes VSPAERO's panel file
   (``VSPAEROComputeGeometry``: ``<case>.vspgeom``, ``<case>.vkey``). With ``method = vlm`` the surface is then
   marked thin in the ``.vkey`` file, so VSPAERO solves it as a vortex lattice (a camber sheet) instead of panels.
2. Trailing-edge wakes are added to the panel file (``vspgeom.py``) along the mesh edges joining the trailing-edge
   vertices the user lists: VSP marks none on a tessellation. A thin surface's file is then rewritten in VSPAERO's
   polygon layout (``vspgeom.polygon_layout``), so VSPAERO solves on the triangles as given.
3. ``<case>.vspaero`` is written with every setting, in the form OpenVSP itself writes it.
4. ``vspaero.exe`` solves all Mach x alpha x beta x Re_cref combinations in one run.
5. Results for nTop: ``<case>.ntvsp.polar.csv`` (one row per case), per-case surface pressure point maps (the
   pressure jump across the sheet on a thin surface), the panels behind them, and the wake lines.

Settings (``key = value``; lengths in metres, angles in degrees, SI elsewhere) are listed in ``KEYS`` below.
"""
from __future__ import annotations

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import adb  # noqa: E402
import vspgeom  # noqa: E402
from common import (DriverError, Log, Settings, fmt, linspace, main_wrapper, run, run_vspscript,  # noqa: E402
                    vsp_path, write_csv, write_point_map)

# Solver settings written to <case>.vspaero, in OpenVSP's order: (file key, settings key, kind).
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


def geometry_script(stl, vsp3, find_body_wakes):
    return f"""void main()
{{
    string gid = ImportFile("{vsp_path(stl)}", IMPORT_STL, "");
    if (gid == "") {{ Print("ERR could not import the surface mesh"); return; }}
    Update();
    SetVSP3FileName("{vsp_path(vsp3)}");
    SetIntAnalysisInput("VSPAEROComputeGeometry", "GeomSet", array<int> = {{0}});
    SetIntAnalysisInput("VSPAEROComputeGeometry", "ThinGeomSet", array<int> = {{-1}});
    SetIntAnalysisInput("VSPAEROComputeGeometry", "Symmetry", array<int> = {{0}});
    SetIntAnalysisInput("VSPAEROComputeGeometry", "FindBodyWakesFlag", array<int> = {{{1 if find_body_wakes else 0}}});
    string rid = ExecAnalysis("VSPAEROComputeGeometry");
    PrintResults(rid);
    Print("NTVSP_OK");
}}
"""


def read_te_points(s: Settings):
    """Trailing-edge vertices (m): the mesh vertices that lie on trailing edges, in any order."""
    xs, ys, zs = s.reals("te_x"), s.reals("te_y"), s.reals("te_z")
    if not (len(xs) == len(ys) == len(zs)):
        raise DriverError("trailing-edge vertex coordinates have different lengths")
    return list(zip(xs, ys, zs))


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
        elif header is not None:
            if any("nan" in v.lower() for v in r[:len(header)]):
                raise DriverError(f"VSPAERO did not converge (NaN results) for Mach {r[1]}, alpha {r[2]}, beta "
                                  f"{r[0]}; see the .vspaero.log")
            try:
                data.append([float(v) for v in r[:len(header)]])
            except ValueError:
                pass
    if header is None:
        raise DriverError(f"{path}: no header line")
    return header, data


def clean(name):
    return name.replace("/", "_over_").replace("Re_over_1e6", "Re_1e6")


def surface_outputs(work, case, g: vspgeom.VSPGeom, order, method, s: Settings, header, data, log):
    """Per-case Cp and pressure difference point maps at triangle centroids and nodes.

    ``order``: the face of ``g`` behind each of VSPAERO's level-0 triangles. On a thin surface (``method`` vlm) the
    values are the pressure jump across the sheet, Cp on the side the face normals point to minus Cp on
    the other side. Also writes ``<case>.ntvsp.panels.txt``
    (the triangles in that order) for VSPAERO Thin Surface Projection (VSP).
    """
    cases = adb.read_cases(os.path.join(work, case + ".adb.cases"))
    try:
        blocks, n_tri = adb.read_cp(os.path.join(work, case + ".adb"), cases, len(order))
        err = adb.check_against_tecplot(blocks, os.path.join(work, case + ".tec.dat"))
    except ValueError as e:
        raise DriverError(f"surface pressures: {e}")
    log(f"surface pressures: {len(cases)} conditions, {n_tri} loops with the wakes (check against Tecplot output: "
        f"{err:.1e})")
    faces = [g.faces[i] for i in order]
    if any(len(f) != 3 for f in faces):
        raise DriverError("the panel file has non-triangular faces")
    nf = len(faces)
    centroids = [tuple(sum(g.nodes[i][c] for i in f) / 3.0 for c in range(3)) for f in faces]
    areas = []
    for f in faces:
        a, b, c = (g.nodes[i] for i in f)
        u, v = [b[k] - a[k] for k in range(3)], [c[k] - a[k] for k in range(3)]
        cx = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        areas.append(0.5 * math.sqrt(cx[0] ** 2 + cx[1] ** 2 + cx[2] ** 2))
    q = 0.5 * s.real("rho") * s.real("vinf") ** 2
    used = sorted({i for f in faces for i in f})
    index = {v: k for k, v in enumerate(used)}
    with open(os.path.join(work, case + ".ntvsp.panels.txt"), "w", newline="\n") as f:
        f.write(f"{method} {len(used)} {nf}\n")
        f.writelines(f"{g.nodes[i][0]!r} {g.nodes[i][1]!r} {g.nodes[i][2]!r}\n" for i in used)
        f.writelines(" ".join(str(index[i]) for i in t) + "\n" for t in faces)
    i_mach, i_aoa, i_beta = header.index("Mach"), header.index("AoA"), header.index("Beta")
    for row_number, row in enumerate(data, start=1):
        # VSPAERO solves Mach 0 as Mach 0.001 and lists it so in .adb.cases; the case list has 7 decimals.
        cond = (max(row[i_mach], 0.001), row[i_aoa], row[i_beta])
        match = [k for k, c in enumerate(cases) if all(abs(a - b) < 1e-5 for a, b in zip(c, cond))]
        if not match:
            raise DriverError(f"case {row_number}: no surface solution for Mach {cond[0]}, alpha {cond[1]}, beta {cond[2]}")
        cp = blocks[match[0]][:nf]
        node_sum, node_w = {}, {}
        for f, a, v in zip(faces, areas, cp):
            for i in f:
                node_sum[i] = node_sum.get(i, 0.0) + a * v
                node_w[i] = node_w.get(i, 0.0) + a
        node_cp = [node_sum[i] / node_w[i] for i in used]
        node_pts = [g.nodes[i] for i in used]
        for quantity, scale in (("cp", 1.0), ("dp", q)):
            write_point_map(os.path.join(work, f"{case}.ntvsp.{quantity}.centroids.{row_number}.csv"), centroids,
                            [v * scale for v in cp])
            write_point_map(os.path.join(work, f"{case}.ntvsp.{quantity}.nodes.{row_number}.csv"), node_pts,
                            [v * scale for v in node_cp])
            both_pts = centroids + node_pts
            both = [v * scale for v in cp] + [v * scale for v in node_cp]
            write_point_map(os.path.join(work, f"{case}.ntvsp.{quantity}.both.{row_number}.csv"), both_pts, both)


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = os.path.abspath(s.text("vsp_dir"))
    timeout = s.real("timeout_s", 7200)
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    for stale in os.listdir(work):
        if stale.startswith(case + ".ntvsp.") and stale.endswith(".csv"):
            os.remove(os.path.join(work, stale))
    stl = os.path.join(work, s.text("mesh_file", case + ".stl"))
    if not os.path.isfile(stl):
        raise DriverError(f"surface mesh {stl} not found")

    # 1. Panel geometry from the tessellation.
    log("OpenVSP: import the surface mesh and write the VSPAERO panel file")
    run_vspscript(vsp_dir, geometry_script(stl, os.path.join(work, case + ".vsp3"), s.flag("find_body_wakes")),
                  work, case + ".geometry", log, timeout)
    geom = os.path.join(work, case + ".vspgeom")
    if not os.path.isfile(geom):
        raise DriverError("OpenVSP wrote no .vspgeom file")
    method = s.text("method", "panel")
    if method not in ("panel", "vlm"):
        raise DriverError(f"method = {method}: expected panel or vlm")
    n_free = vspgeom.free_edges(vspgeom.read(geom))
    if method == "vlm":
        vspgeom.set_thick(os.path.join(work, case + ".vkey"), False)
        log(f"method: vortex lattice (the mesh is a thin camber surface; {n_free} free edges)")
    else:
        log("method: panel (the mesh is a thick, closed surface)")
        if n_free:
            log(f"warning: the mesh has {n_free} free edges, so it is not closed; the panel method needs a closed "
                "surface (use the vortex lattice method for a thin sheet)")

    # 2. Trailing-edge wakes.
    try:
        n_wakes, n_edges, polylines = vspgeom.add_wakes(geom, geom, read_te_points(s), s.real("te_tolerance"))
    except ValueError as e:
        raise DriverError(f"trailing edges: {e}")
    log(f"wakes: {n_wakes} trailing edges, {n_edges} edges")
    if n_wakes == 0:
        log("warning: no trailing-edge wakes; VSPAERO will shed no vorticity and the surface carries almost no lift")
    tri_geom = vspgeom.read(geom)
    order = list(range(len(tri_geom.faces)))
    if method == "vlm":
        # VSPAERO's own pairing of triangles can drop trailing-wake edges behind swept trailing edges (see
        # vspgeom.polygon_layout); hand it the triangles as its solution grid instead, as OpenVSP does for wings.
        try:
            order = vspgeom.polygon_layout(geom, geom)
        except ValueError as e:
            raise DriverError(f"thin surface: {e}")
        log("panel file in VSPAERO's polygon layout: VSPAERO solves on the mesh triangles as given")
    wake_pts, wake_ids = [], []
    for k, line in enumerate(polylines, start=1):
        wake_pts += line
        wake_ids += [float(k)] * len(line)
    write_point_map(os.path.join(work, case + ".ntvsp.wakes.csv"), wake_pts or [(0.0, 0.0, 0.0)], wake_ids or [0.0])

    # 3. Solver input.
    machs = linspace(s.real("mach_start"), s.real("mach_end"), s.integer("mach_npts"))
    alphas = linspace(s.real("alpha_start"), s.real("alpha_end"), s.integer("alpha_npts"))
    betas = linspace(s.real("beta_start"), s.real("beta_end"), s.integer("beta_npts"))
    recrefs = linspace(s.real("recref_start"), s.real("recref_end"), s.integer("recref_npts"))
    with open(os.path.join(work, case + ".vspaero"), "w", newline="\n") as f:
        f.write(vspaero_file(s, machs, alphas, betas, recrefs))
    n_cases = len(machs) * len(alphas) * len(betas) * len(recrefs)
    log(f"VSPAERO: {len(machs)} Mach x {len(alphas)} alpha x {len(betas)} beta x {len(recrefs)} Re_cref = {n_cases} cases")

    # 4. Solve.
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

    # 5. Results.
    polar = os.path.join(work, case + ".polar")
    if not os.path.isfile(polar):
        raise DriverError(f"VSPAERO wrote no .polar file; see {case}.vspaero.log")
    header, data = read_polar(polar)
    if len(data) != n_cases:
        raise DriverError(f"VSPAERO's polar has {len(data)} rows for {n_cases} cases; see {case}.vspaero.log")
    write_csv(os.path.join(work, case + ".ntvsp.polar.csv"), ["Case"] + [clean(h) for h in header],
              [[i] + row for i, row in enumerate(data, start=1)])
    surface_outputs(work, case, tri_geom, order, method, s, header, data, log)
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
