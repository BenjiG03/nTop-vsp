"""VSPAERO custom blocks: settings, analysis, and result readers."""
from __future__ import annotations

from ntopfile import Graph, InputSpec

import pressure_blocks as pb
from kit import (default_text, LENGTH, F_choice, F_int, F_length, F_real, F_switch, Kit, analysis, d_choice, d_int,
                 d_text, make_block, model_inputs, run_inputs)

# --- Solver

SOLVER_FIELDS = [
    F_int("Processors", "ncpu", 4, "Threads for VSPAERO (-omp)."),
    F_switch("Wake", "fixed_wake", False, "Relaxed: the wake is aligned with the local flow over the wake iterations. "
             "Fixed: the wake stays along the freestream (no iterations).", labels=("Relaxed", "Fixed")),
    F_int("Wake Iterations", "wake_iterations", 5, "Number of wake relaxation iterations (3 or more)."),
    F_int("Wake Nodes", "wake_nodes", 16, "Number of nodes along each trailing wake line. VSPAERO's run time grows "
          "steeply with it: on the X-59 at Mach 1.6, 16 gave the same coefficients as 64 (within 0.2%) in a twentieth "
          "of the time."),
    F_switch("Far Field Distance", "far_dist_on", False, "Automatic: VSPAERO picks the wake length. Manual: use the "
             "distance below.", labels=("Automatic", "Manual")),
    F_length("Manual Far Field Distance", "far_dist", 100.0, "Distance the wake is adapted to, when Far Field Distance is Manual."),
    F_switch("Ground Effect", "ground_effect_on", False, "Image method for flight near the ground (-groundheight)."),
    F_length("Ground Height", "ground_height", 10.0, "Height of the moment reference point above the ground."),
    F_int("Freeze Multipole At Iteration", "freeze_multipole_at_iteration", 10000, "Stop updating the multipole expansion "
          "after this iteration."),
    F_int("Freeze Wake At Iteration", "freeze_wake_at_iteration", 10000, "Stop moving the wake after this iteration."),
    F_switch("Freeze Wake Root Vortices", "freeze_wake_root_vortices", False, "Hold the wake root vortices fixed."),
    F_switch("Implicit Wake", "implicit_wake", False, "Solve the wake shape implicitly with the surface solution."),
    F_int("Implicit Wake Start Iteration", "implicit_wake_start_iteration", 0, "Iteration at which the implicit wake starts."),
    F_real("Wake Relaxation", "wake_relax", 1.0, description="Under-relaxation of the wake update (0 to 1)."),
    F_real("Forward GMRES Convergence Factor", "forward_gmres_factor", 1.0, description="Scales the default residual "
           "reduction of the forward GMRES solve."),
    F_real("Adjoint GMRES Convergence Factor", "adjoint_gmres_factor", 1.0, description="Scales the default residual "
           "reduction of the adjoint GMRES solve."),
    F_real("Nonlinear Convergence Factor", "nonlinear_factor", 1.0, description="Scales the default residual reduction "
           "of the overall nonlinear solve."),
    F_real("Core Size Factor", "core_size_factor", 1.0, description="Scales the wake vortex core size model."),
    F_real("Multipole Far Away", "far_away", 5.0, description="Distance ratio beyond which the multipole expansion is used."),
    F_switch("Update Matrix Preconditioner", "update_matrix_preconditioner", False, "Rebuild the matrix preconditioners "
             "every Newton iteration."),
    F_switch("Wake Node Matrix Preconditioner", "wake_node_matrix_preconditioner", False, "Use the wake node matrix "
             "preconditioner."),
    F_int("Quad Tree Buffer Levels", "quadtree_buffer_levels", 0, "Buffer levels of VSPAERO's quad tree (more levels, "
          "more cells)."),
    F_choice("Restart File", "restart_mode", ("None", "Save (-save)", "Restart from saved (-restart)"), (0, 1, 2), 0,
             "VSPAERO's restart file: write one, or start from the one saved under this Case Name."),
    F_choice("Convergence Check", "convergence_check", ("Stop on unconverged cases", "Warn in the log only"), (1, 0), 0,
             "After the run every case is checked: VSPAERO's linear solve (GMRES) of the last wake iteration reached "
             "both its residual targets, every force and moment coefficient is finite, and with a relaxed wake the "
             "forces changed less than the Wake Change Tolerance over the last wake iteration. Stop: an unconverged "
             "case stops the notebook with the reason. Either way `<Case Name>.ntvsp.convergence.csv` lists every "
             "case, and the polar's Converged column flags them."),
    F_real("Wake Change Tolerance", "wake_tolerance", 0.1, description="Largest change of CLtot, CDtot and CMytot "
           "over the last wake iteration, relative to their values (at least 0.05, 0.005 and 0.01), for a relaxed-wake "
           "case to count as converged. 0.1 is 10%: converged runs stay within about 6%, diverging ones change by "
           "tens of percent or more."),
]

def settings_block(key, display, description, fields, comments):
    specs = [f.spec for f in fields]
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    for c in comments:
        g.comment(c)
    lines = k.field_lines(fields, lambda f: inp(f.spec.name))
    lst = g.block("core.list<text>", "list<text>", lines, top=False)
    out = g.block("concatenate<list<text>,text>", "text", [lst, k.t("")], "Settings Text")
    return make_block(key, display, description, specs, "text", g, out)


def solver_block():
    return settings_block(
        "vspaero_solver", "VSPAERO Solver Settings (VSP)",
        "Solver options of **VSPAERO Analysis (VSP)**: threads, wake model and iterations, far field, ground effect, "
        "multipole and wake freezing, implicit wake, convergence factors, preconditioners, quad tree buffer levels "
        "and VSPAERO's restart file. Defaults are OpenVSP's, except Wake Nodes (16, see its description).",
        SOLVER_FIELDS,
        ["**VSPAERO solver settings.** The advanced options of OpenVSP's VSPAERO screen and the vspaero command line "
         "(-omp, -groundheight, -save, -restart), as `key = value` lines for the analysis block."])


# --- Analysis

PRESSURE_EXPORT = ("Off", "Panel centroids", "Nodes", "Centroids and nodes")
PRESSURE_EXPORT_KEYS = ("off", "centroids", "nodes", "both")
PRESSURE_EXPORT_DOC = (
    "Off: no pressure files. Otherwise, per flight condition case, at the chosen points: Thick-set components' "
    "static pressure p = p_inf + q Cp (Pa) for **Surface Pressure (VSP)**, Thin-set components' jump q Delta Cp (Pa) "
    "for **Surface Pressure Jump (VSP)**. Centroids carry VSPAERO's own values, nodes area-weighted averages. Each "
    "run first removes the pressure files of an earlier run.")


def analysis_block(solver=SOLVER_FIELDS):
    own = (
        InputSpec("OpenVSP Model", "text", "Output of OpenVSP Model (VSP): the model whose Thick set VSPAERO solves "
                  "with its panel method and whose Thin set with its vortex lattice."),
        InputSpec("Flight Condition", "text", "Output of Flight Condition (VSP). VSPAERO runs its distinct Mach numbers, "
                  "angles of attack and sideslips (every combination) and the Reynolds numbers of its flight points, "
                  "with its reference area, span, chord and moment point."),
        InputSpec("Mach Min", "real", "Lowest Mach number of the flight condition that VSPAERO runs.",
                  default=lambda g: g.real(0.0)),
        InputSpec("Mach Max", "real", "VSPAERO runs the flight condition's Mach numbers below this one (for example 3, "
                  "leaving higher Mach numbers to Modified Newtonian).", default=lambda g: g.real(1000.0)),
        InputSpec("Solver Settings", "text", "From VSPAERO Solver Settings (VSP).", default=d_text(default_text(solver))),
        InputSpec("Body Wakes", "choice", "OpenVSP's FindBodyWakesFlag: also shed wakes from the bases of thick "
                  "fuselages. Wing trailing edges always shed wakes.", default=d_choice(("Off", "On"))),
        pb.export_input(PRESSURE_EXPORT_DOC, PRESSURE_EXPORT),
    )
    specs = own + run_inputs("vspaero")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**VSPAERO analysis of an OpenVSP model.** 1.")
    lines = [
        k.line("model", inp("OpenVSP Model")),
        k.line("flight", inp("Flight Condition")),
        k.line("mach_min", k.num(inp("Mach Min"))),
        k.line("mach_max", k.num(inp("Mach Max"))),
        k.line("find_body_wakes", k.pick(inp("Body Wakes"), ("0", "1"))),
        k.line("export_pressure", k.pick(inp("Export Pressures"), PRESSURE_EXPORT_KEYS)),
        k.cat(["\n", inp("Solver Settings"), "\n"]),
    ]
    out = analysis(g, k, inp, "run_vspaero", lambda base: [], lines, "VSPAERO")
    return make_block(
        "vspaero_analysis", "VSPAERO Analysis (VSP)",
        "Runs OpenVSP's **VSPAERO** on an OpenVSP model (from **OpenVSP Model (VSP)**) and returns the path of its "
        "results (connect it to the VSPAERO reader blocks). Thick-set components are solved with the panel method, "
        "Thin-set components with the vortex lattice method; wakes leave every wing trailing edge. The conditions, "
        "reference area, span, chord and moment point come from **Flight Condition (VSP)**, the solver options from "
        "**VSPAERO Solver Settings (VSP)**. Every case is checked for convergence (linear-solve residuals, finite coefficients, wake "
        "settled); by default an unconverged case stops the notebook with the reason, and "
        "`<Case Name>.ntvsp.convergence.csv` lists every case. Files are written to the Working Directory under the "
        "Case Name; the run's log is `<Case Name>.ntvsp.log`.\n\nNot exposed (out of scope): X-Z symmetry, 2D FEM and Tecplot output, reference "
        "quantities and CG from the model, stall model, propellers and actuator disks, stability analyses, slicer, "
        "control surfaces.", specs, "text", g, out)


def prop_of(ref, c):
    from kit import prop
    return prop(ref, c)




# --- Readers

POLAR_COLUMNS = ("CLtot", "CDtot", "CStot", "CMxtot", "CMytot", "CMztot", "CFxtot", "CFytot", "CFztot", "L_over_D", "E",
                 "CLo", "CLi", "CDo", "CDi", "CSo", "CSi", "CMox", "CMoy", "CMoz", "CMix", "CMiy", "CMiz",
                 "CFox", "CFoy", "CFoz", "CFix", "CFiy", "CFiz", "CLwtot", "CDwtot", "CSwtot", "CLiw", "CDiw", "CSiw",
                 "CFwxtot", "CFwytot", "CFwztot", "CFiwx", "CFiwy", "CFiwz", "LoDw", "Ew", "Mach", "AoA", "Beta",
                 "Re_1e6", "Case", "Converged")


def table_column(g, k, path, header):
    table = g.block("import_table<file_path,bool>", "table", [path, g.boolean(True)], "Import Results Table")
    col = g.block("extract_column<table,text>", "list<text>", [table, header], top=False)
    return g.block("scalar_from_text<text>", "real", [col], "Values")


def coefficient_block():
    specs = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
        InputSpec("Coefficient", "choice", "Column of VSPAERO's polar. tot: total; o: skin friction (CDo...); i: "
                  "inviscid surface pressure; w: from the wake (Trefftz plane). Mach, AoA, Beta, Re_1e6 and Case give "
                  "each case's conditions; Converged is 1 for a case that passed the convergence check, 0 otherwise.",
                  default=d_choice(POLAR_COLUMNS)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**VSPAERO coefficients.** The driver writes VSPAERO's `.polar` as `<Case Name>.ntvsp.polar.csv`, one "
              "row per case (Re_cref outermost, then Mach, beta, alpha).")
    path = k.cat([inp("Results"), ".ntvsp.polar.csv"])
    out = table_column(g, k, path, k.pick(inp("Coefficient"), POLAR_COLUMNS))
    return make_block("vspaero_coefficient", "VSPAERO Coefficient (VSP)",
                      "One force or moment coefficient (or a case condition) of every case of a VSPAERO run, as a list "
                      "in case order.", specs, "list<real>", g, out)


FACE_POINTS = ("4 per face (each face split in four)", "1 per face (face centroids)")


def mesh_faces_block():
    specs = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
        InputSpec("Points per Face", "choice", "Face centroids, or the centroids of each face split in four at its "
                  "edge midpoints (a finer rule; moments converge better, at a little more time).",
                  default=d_choice(FACE_POINTS)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Mesh faces.** The driver writes the panel triangles of the Thick-set components as surface points "
              "with the area each stands for: `x, y, z, area` (metres, square metres).")
    n = k.pick(inp("Points per Face"), ("4", "1"))
    path = k.cat([inp("Results"), ".ntvsp.faces.", n, ".csv"], name="Faces File")
    scale = g.leaf("real_field", {"expression": "1", "units": {}})
    out = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                  [path, g.unit_length("m"), scale], "Mesh Faces")
    return make_block("vspaero_mesh_faces", "VSPAERO Mesh Faces (VSP)",
                      "The thick (panel) surfaces of a VSPAERO run as points weighted by area (Scalar Point Map, value: "
                      "m^2), for the Mesh Faces input of the Viscous Correction on Mesh blocks.",
                      specs, "real_point_map", g, out)


CENTERS = ("Center of pressure (one point per condition)", "Aerodynamic center (one point per alpha sweep)")


def center_block():
    specs = (
        InputSpec("Aero Database", "text", "Output of Aero Database (VSP): the centers come from its totals (every "
                  "contribution, about the flight condition's moment point)."),
        InputSpec("Quantity", "choice", "Center of pressure: where the resultant of each condition's normal force "
                                        "and pitching moment crosses the line through the moment reference point "
                                        "along x (it runs away as the normal force goes to zero).", default=d_choice(CENTERS)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Center of pressure.** The Aero Database driver computes, from the total normal force CN (z up) and "
              "pitching moment Cm about the moment reference point (over Sref and cref, x aft): x_cp = x_ref - Cm "
              "cref / CN for every condition, and x_ac = x_ref - (dCm/dCN) cref for every alpha sweep at one Mach "
              "number, altitude and sideslip (a straight-line fit).")
    path = k.cat([inp("Aero Database"), k.pick(inp("Quantity"), (".ntvsp.center_of_pressure.csv",
                                                            ".ntvsp.aerodynamic_center.csv"))], name="Table File")
    table = g.block("import_table<file_path,bool>", "table", [path, g.boolean(True)], "Import Table")
    metre = g.real(1.0, LENGTH)
    xyz = []
    for c in "XYZ":
        col = g.block("extract_column<table,text>", "list<text>", [table, k.t(c)], top=False)
        val = g.block("scalar_from_text<text>", "real", [col], top=False)
        xyz.append(g.block("multiply<real,real>", "real", [val, metre], top=False))
    out = g.block("point<real,real,real>", "list<point>", xyz, "Centers")
    return make_block(
        "center_of_pressure", "Center of Pressure (VSP)",
        "Center of pressure of every condition of an Aero Database, or the aerodynamic center (neutral point) of "
        "every alpha sweep, as a list of points: absolute positions in the model's frame with length units (not "
        "fractions of the reference chord), on the line through the moment reference point along x.", specs, "list<point>", g, out)
