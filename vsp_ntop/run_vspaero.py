"""VSPAERO analysis of an OpenVSP model, driven from nTop.

    python -m vsp_ntop.run_vspaero SETTINGS_FILE

Steps (each logged to ``<case>.ntvsp.log``):

1. OpenVSP reads the model (``model`` = a ``.vsp3`` made by OpenVSP Model) and writes VSPAERO's panel file
   (``VSPAEROComputeGeometry``: ``<case>.vspgeom``, ``<case>.vkey``). Components in the model's "Thick" set become
   panel-method surfaces, those in the "Thin" set vortex-lattice camber surfaces. OpenVSP marks the trailing edges
   of its wings as wakes itself; with ``find_body_wakes`` it also marks wakes behind fuselage bases.
2. ``<case>.vspaero`` is written with every setting, in the form OpenVSP itself writes it.
3. ``vspaero.exe`` solves all Mach x alpha x beta x Re_cref combinations in one run.
4. Results for nTop: ``<case>.ntvsp.polar.csv`` (one row per case), the centre of pressure of each case and the
   aerodynamic centre of each alpha sweep (``pressure_centers``), per-case surface pressure point maps (on thin
   components, the pressure jump across the sheet) and the panels' points and areas (``.ntvsp.faces.{1,4}.csv``).

Settings (``key = value``; lengths in metres, angles in degrees, SI elsewhere) are listed in the blocks.
"""
from __future__ import annotations

import math
import os


from . import adb
from . import vspgeom
from .common import (openvsp_dir, DriverError, Log, Settings, fmt, linspace, main_wrapper, run, run_vspscript,
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


def surface_outputs(work, case, g: vspgeom.VSPGeom, s: Settings, header, data, log):
    """Per-case Cp and pressure difference point maps at the panel triangles' centroids and nodes, and the panels'
    points and areas (``<case>.ntvsp.faces.{1,4}.csv``, thick components only) for the viscous correction.

    On thin (vortex lattice) components the values are the pressure jump across the sheet: Cp on the side the
    panel normals point to minus Cp on the other side."""
    tris = vspgeom.level0(g)
    cases = adb.read_cases(os.path.join(work, case + ".adb.cases"))
    try:
        blocks, n_tri = adb.read_cp(os.path.join(work, case + ".adb"), cases, len(tris))
        err = adb.check_against_tecplot(blocks, os.path.join(work, case + ".tec.dat"))
    except ValueError as e:
        raise DriverError(f"surface pressures: {e}")
    log(f"surface pressures: {len(cases)} conditions, {n_tri} loops with the wakes (check against Tecplot output: "
        f"{err:.1e})")
    thin = vspgeom.thin_parts(os.path.join(work, case + ".vkey"))
    faces = [t for t, _ in tris]
    is_thin = [g.parts[i] in thin for _, i in tris]
    nf = len(faces)
    centroids = [tuple(sum(g.nodes[i][c] for i in f) / 3.0 for c in range(3)) for f in faces]
    areas = []
    for f in faces:
        a, b, c = (g.nodes[i] for i in f)
        u, v = [b[k] - a[k] for k in range(3)], [c[k] - a[k] for k in range(3)]
        cx = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        areas.append(0.5 * math.sqrt(cx[0] ** 2 + cx[1] ** 2 + cx[2] ** 2))
    one_pts, one_areas, four_pts, four_areas = [], [], [], []
    for f, a, c, t in zip(faces, areas, centroids, is_thin):
        if t:
            continue
        one_pts.append(c)
        one_areas.append(a)
        p, q_, r = (g.nodes[i] for i in f)
        pq, qr, rp = ([(u[k] + v[k]) / 2 for k in range(3)] for u, v in ((p, q_), (q_, r), (r, p)))
        for tri in ((p, pq, rp), (pq, q_, qr), (rp, qr, r), (pq, qr, rp)):
            four_pts.append(tuple(sum(v[k] for v in tri) / 3.0 for k in range(3)))
            four_areas.append(a / 4.0)
    if one_pts:
        write_point_map(os.path.join(work, case + ".ntvsp.faces.1.csv"), one_pts, one_areas)
        write_point_map(os.path.join(work, case + ".ntvsp.faces.4.csv"), four_pts, four_areas)
    q = 0.5 * s.real("rho") * s.real("vinf") ** 2
    used = sorted({i for f in faces for i in f})
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
            write_point_map(os.path.join(work, f"{case}.ntvsp.{quantity}.both.{row_number}.csv"), centroids + node_pts,
                            [v * scale for v in cp] + [v * scale for v in node_cp])


def pressure_centers(work, case, header, data, s: Settings, log):
    """Centre of pressure of every case and aerodynamic centre of every (Mach, beta, Re_cref) alpha sweep.

    VSPAERO's axes are x aft, z up, and CMy (nose up) is about the moment reference point over Sref cref. The
    resultant of the normal force CFz and the moment crosses the line through the reference point along x at
    x_cp = x_ref - CMy cref / CFz (inviscid coefficients, CFiz and CMiy). The aerodynamic centre is where the moment
    does not change with lift: x_ac = x_ref - (dCMy/dCFz) cref, the slope fitted over the alpha sweep. Both are
    given as points (x, y_ref, z_ref) in metres. A case with no normal force (|CFz| < 1e-6) has no centre of
    pressure: its point is the reference point and Defined is 0.
    """
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
        f = [r[i["CFiz"]] for r in rows_g]
        c = [r[i["CMiy"]] for r in rows_g]
        n = len(f)
        mf, mc = sum(f) / n, sum(c) / n
        sff = sum((x - mf) ** 2 for x in f)
        if n < 2 or sff < 1e-12:
            continue
        slope = sum((x - mf) * (y - mc) for x, y in zip(f, c)) / sff
        ac.append([m, b, re * 1e6, xr - slope * cref, yr, zr, slope])
    write_csv(os.path.join(work, case + ".ntvsp.aerodynamic_center.csv"),
              ["Mach", "Beta", "Re_cref", "X", "Y", "Z", "dCMy_dCFz"], ac)
    if ac:
        log("aerodynamic centre x: " + ", ".join(f"Mach {r[0]:g} beta {r[1]:g}: {r[3]:.3f} m" for r in ac))


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
    machs = linspace(s.real("mach_start"), s.real("mach_end"), s.integer("mach_npts"))
    alphas = linspace(s.real("alpha_start"), s.real("alpha_end"), s.integer("alpha_npts"))
    betas = linspace(s.real("beta_start"), s.real("beta_end"), s.integer("beta_npts"))
    recrefs = linspace(s.real("recref_start"), s.real("recref_end"), s.integer("recref_npts"))
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
    write_csv(os.path.join(work, case + ".ntvsp.polar.csv"), ["Case"] + [clean(h) for h in header],
              [[i] + row for i, row in enumerate(data, start=1)])
    pressure_centers(work, case, header, data, s, log)
    surface_outputs(work, case, g, s, header, data, log)
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
