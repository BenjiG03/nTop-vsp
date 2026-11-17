"""VSPAERO custom blocks: settings, analysis, and result readers."""
from __future__ import annotations

from ntopfile import Graph, InputSpec

from kit import (default_text, AREA, DENSITY, LENGTH, VELOCITY, F_angle, F_choice, F_int, F_length,
                 F_real, F_switch, Kit, analysis, d_choice, d_int, d_point_mm,
                 d_real, d_text, make_block, model_inputs, run_inputs)

# --- Flow condition ---------------------------------------------------------------------------------------------

FLOW_FIELDS = [
    F_angle("Alpha Start", "alpha_start", 0.0, "First angle of attack of the sweep."),
    F_angle("Alpha End", "alpha_end", 10.0, "Last angle of attack of the sweep."),
    F_int("Alpha Points", "alpha_npts", 3, "Number of angles of attack, evenly spaced from Start to End (1: Start only)."),
    F_angle("Beta Start", "beta_start", 0.0, "First sideslip angle."),
    F_angle("Beta End", "beta_end", 0.0, "Last sideslip angle."),
    F_int("Beta Points", "beta_npts", 1, "Number of sideslip angles."),
    F_real("Mach Start", "mach_start", 0.0, description="First freestream Mach number (Prandtl-Glauert/Karman-Tsien "
           "compressibility in VSPAERO; subsonic)."),
    F_real("Mach End", "mach_end", 0.0, description="Last freestream Mach number."),
    F_int("Mach Points", "mach_npts", 1, "Number of Mach numbers."),
    F_real("Re_cref Start", "recref_start", 1.0e7, description="First Reynolds number based on the reference chord "
           "(used for VSPAERO's skin-friction drag estimate)."),
    F_real("Re_cref End", "recref_end", 2.0e7, description="Last Reynolds number based on the reference chord."),
    F_int("Re_cref Points", "recref_npts", 1, "Number of Reynolds numbers."),
    F_real("Freestream Velocity", "vinf", 100.0, VELOCITY, "Vinf: scales dimensional loads and pressures."),
    F_real("Freestream Density", "rho", 1.225, DENSITY, "Rho: scales dimensional loads and pressures."),
    F_switch("Reference Velocity", "manual_vref", False, "Automatic: Vref = Vinf. Manual: use the values below.",
             labels=("Automatic (Vinf)", "Manual")),
    F_real("Manual Reference Velocity", "vref", 100.0, VELOCITY, "Vref when Reference Velocity is Manual."),
    F_real("Manual Reference Mach", "machref", 0.3, description="Mref when Reference Velocity is Manual."),
]

# --- Solver ----------------------------------------------------------------------------------------------------------

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


def flow_block():
    return settings_block(
        "vspaero_flow", "VSPAERO Flow Condition (VSP)",
        "Flight conditions for **VSPAERO Analysis (VSP)**: the alpha, beta, Mach and Re_cref sweeps (VSPAERO solves "
        "every combination in one run, like the Flow Condition section of OpenVSP's VSPAERO screen), the freestream "
        "velocity and density, and the reference velocity. Connect the output to the analysis block's **Flow "
        "Condition**.", FLOW_FIELDS,
        ["**VSPAERO flow condition.** Each input becomes one `key = value` line of text, in SI units with angles in "
         "degrees. Sweeps are Start, End and Points, evenly spaced, as in OpenVSP. The analysis block adds these lines "
         "to its settings file."])


def solver_block():
    return settings_block(
        "vspaero_solver", "VSPAERO Solver Settings (VSP)",
        "Solver options of **VSPAERO Analysis (VSP)**: threads, wake model and iterations, far field, ground effect, "
        "multipole and wake freezing, implicit wake, convergence factors, preconditioners, quad tree buffer levels "
        "and VSPAERO's restart file. Defaults are OpenVSP's, except Wake Nodes (16, see its description).",
        SOLVER_FIELDS,
        ["**VSPAERO solver settings.** The advanced options of OpenVSP's VSPAERO screen and the vspaero command line "
         "(-omp, -groundheight, -save, -restart), as `key = value` lines for the analysis block."])


# --- Analysis --------------------------------------------------------------------------------------------------------

def analysis_block(flow=FLOW_FIELDS, solver=SOLVER_FIELDS):
    own = (
        InputSpec("OpenVSP Model", "text", "Output of OpenVSP Model (VSP): the model whose Thick set VSPAERO solves "
                  "with its panel method and whose Thin set with its vortex lattice."),
        InputSpec("Reference Area", "real", "Sref.", AREA, d_real(1.0, AREA)),
        InputSpec("Reference Span", "real", "bref.", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Reference Chord", "real", "cref.", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Moment Reference Point", "point", "Moments are taken about this point (VSPAERO's X_cg, Y_cg, Z_cg).",
                  LENGTH, d_point_mm(0, 0, 0)),
        InputSpec("Flow Condition", "text", "From VSPAERO Flow Condition (VSP).", default=d_text(default_text(flow))),
        InputSpec("Solver Settings", "text", "From VSPAERO Solver Settings (VSP).", default=d_text(default_text(solver))),
        InputSpec("Body Wakes", "choice", "OpenVSP's FindBodyWakesFlag: also shed wakes from the bases of thick "
                  "fuselages. Wing trailing edges always shed wakes.", default=d_choice(("Off", "On"))),
    )
    specs = own + run_inputs("vspaero")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**VSPAERO analysis of an OpenVSP model.** 1. Write the settings file. 2. The Python driver has OpenVSP "
              "read the model and write VSPAERO's panel file: components in the Thick set become panel-method surfaces, "
              "those in the Thin set vortex-lattice camber surfaces, and OpenVSP marks every wing trailing edge as a "
              "wake. 3. It writes the `.vspaero` input and runs `vspaero.exe` for every Mach x alpha x beta x Re_cref "
              "case. 4. It checks that every case converged (Solver Settings' Convergence Check) and stops with the "
              "reason if one did not. 5. The output is the path of the results, for the VSPAERO reader blocks.")
    lines = [
        k.line("model", inp("OpenVSP Model")),
        k.line("sref", k.num(inp("Reference Area"), 1.0, AREA)),
        k.line("bref", k.num(inp("Reference Span"), 1.0, LENGTH)),
        k.line("cref", k.num(inp("Reference Chord"), 1.0, LENGTH)),
    ] + [k.line(f"{c}cg", k.num(prop_of(inp("Moment Reference Point"), c), 1.0, LENGTH)) for c in "xyz"] + [
        k.line("find_body_wakes", k.pick(inp("Body Wakes"), ("0", "1"))),
        k.cat(["\n", inp("Flow Condition"), "\n", inp("Solver Settings"), "\n"]),
    ]
    out = analysis(g, k, inp, "run_vspaero", lambda base: [], lines, "VSPAERO")
    return make_block(
        "vspaero_analysis", "VSPAERO Analysis (VSP)",
        "Runs OpenVSP's **VSPAERO** on an OpenVSP model (from **OpenVSP Model (VSP)**) and returns the path of its "
        "results (connect it to the VSPAERO reader blocks). Thick-set components are solved with the panel method, "
        "Thin-set components with the vortex lattice method; wakes leave every wing trailing edge. Settings come from "
        "**VSPAERO Flow Condition** and **VSPAERO Solver Settings** (VSP). Reference area, span, chord and moment "
        "point are inputs. Every case is checked for convergence (linear-solve residuals, finite coefficients, wake "
        "settled); by default an unconverged case stops the notebook with the reason, and "
        "`<Case Name>.ntvsp.convergence.csv` lists every case. Files are written to the Working Directory under the "
        "Case Name; the run's log is `<Case Name>.ntvsp.log`.\n\nNot exposed (out of scope): X-Z symmetry, 2D FEM and Tecplot output, reference "
        "quantities and CG from the model, stall model, propellers and actuator disks, stability analyses, slicer, "
        "control surfaces.", specs, "text", g, out)


def prop_of(ref, c):
    from kit import prop
    return prop(ref, c)




# --- Readers ---------------------------------------------------------------------------------------------------------

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
    g.comment("**VSPAERO coefficients.** The driver writes VSPAERO's `.polar` as `<Case Name>.ntvsp.polar.csv`, one row "
              "per case (Re_cref outermost, then Mach, beta, alpha). Import Table reads it and the chosen column comes "
              "out as a list of reals, one per case.")
    path = k.cat([inp("Results"), ".ntvsp.polar.csv"])
    out = table_column(g, k, path, k.pick(inp("Coefficient"), POLAR_COLUMNS))
    return make_block("vspaero_coefficient", "VSPAERO Coefficient (VSP)",
                      "One force or moment coefficient (or a case condition) of every case of a VSPAERO run, as a list "
                      "in case order.", specs, "list<real>", g, out)


SURFACE_QUANTITIES = ("Pressure coefficient Cp", "Pressure difference p - p_inf")
SURFACE_POINTS = ("Panel centroids", "Nodes", "Centroids and nodes")


def surface_block():
    specs = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
        InputSpec("Case Number", "integer", "Case, counted from 1 in the order of VSPAERO Coefficient's lists.",
                  default=d_int(1)),
        InputSpec("Quantity", "choice", "Cp, or p - p_inf = 0.5 Rho Vinf^2 Cp in pascals (Flow Condition's velocity and "
                  "density). On Thin-set components: the jump across the sheet instead, Cp on the side the panel "
                  "normals point to minus Cp on the other side, or the same jump in pascals.",
                  default=d_choice(SURFACE_QUANTITIES)),
        InputSpec("Sample Points", "choice", "Panel centroids carry VSPAERO's own values; node values are area-weighted "
                  "averages of the panels around each node; both together give the densest colouring.",
                  default=d_choice(SURFACE_POINTS, 2)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Surface pressure as a Scalar Point Map.** The driver reads every case's panel pressures from "
              "VSPAERO's `.adb` file and writes them as `x, y, z, value` files (metres). This block picks the file "
              "for the case, quantity and sample points and imports it. Show the point map to see the contours; "
              "Field from Point Map turns it into a field.")
    qty = k.pick(inp("Quantity"), ("cp", "dp"))
    pts = k.pick(inp("Sample Points"), ("centroids", "nodes", "both"))
    path = k.cat([inp("Results"), ".ntvsp.", qty, ".", pts, ".", k.num(inp("Case Number")), ".csv"], name="Point Map File")
    scale = g.leaf("real_field", {"expression": "1", "units": {}})
    out = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                  [path, g.unit_length("m"), scale], "Surface Pressure")
    return make_block("vspaero_surface", "VSPAERO Surface Pressure (VSP)",
                      "Surface pressure of one VSPAERO case as a **Scalar Point Map** for contour plots in nTop: Cp or "
                      "p - p_inf (Pa), at panel centroids, nodes, or both. On Thin-set components these are the "
                      "pressure jump across the sheet.", specs, "real_point_map", g, out)


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
    g.comment("**Mesh faces.** The driver writes the panel triangles of the Thick-set components as surface points with "
              "the area each stands for: `x, y, z, area` (metres, square metres). Feed this point map to the Mesh Faces "
              "input of ntop_newtonian's Viscous Correction on Mesh blocks, which then sum the skin friction over the "
              "points instead of integrating it over the body: seconds instead of minutes.")
    n = k.pick(inp("Points per Face"), ("4", "1"))
    path = k.cat([inp("Results"), ".ntvsp.faces.", n, ".csv"], name="Faces File")
    scale = g.leaf("real_field", {"expression": "1", "units": {}})
    out = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                  [path, g.unit_length("m"), scale], "Mesh Faces")
    return make_block("vspaero_mesh_faces", "VSPAERO Mesh Faces (VSP)",
                      "The thick (panel) surfaces of a VSPAERO run as points weighted by area (Scalar Point Map, value: "
                      "m^2), for the Mesh Faces input of the Viscous Correction on Mesh blocks.",
                      specs, "real_point_map", g, out)


CENTERS = ("Center of pressure (one point per case)", "Aerodynamic center (one point per Mach, beta, Re_cref)")


def center_block():
    specs = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
        InputSpec("Quantity", "choice", "Center of pressure: where the resultant of each case's normal force and "
                  "pitching moment crosses the line through the moment reference point along x (it runs away as the "
                  "normal force goes to zero). Aerodynamic center: the point about which the pitching moment does not "
                  "change with lift (the neutral point), from the alpha sweep at each Mach number, beta and Re_cref.",
                  default=d_choice(CENTERS)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Center of pressure.** The driver computes, from VSPAERO's inviscid normal force CFz and pitching "
              "moment CMy about the moment reference point (over Sref and cref, x aft): x_cp = x_ref - CMy cref / CFz "
              "for every case, and x_ac = x_ref - (dCMy/dCFz) cref for every alpha sweep (a straight-line fit). The "
              "points lie on the line through the reference point along x; a case with no normal force gets the "
              "reference point (column Defined = 0 in the file).")
    path = k.cat([inp("Results"), k.pick(inp("Quantity"), (".ntvsp.center_of_pressure.csv",
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
        "vspaero_center_of_pressure", "VSPAERO Center of Pressure (VSP)",
        "Center of pressure of every VSPAERO case, or the aerodynamic center (neutral point) of every alpha sweep, as "
        "a list of points: absolute positions in the model's frame with length units (not fractions of the "
        "reference chord), on the line through the moment reference point along x. Compare the aerodynamic center "
        "with the CG for the static margin.", specs, "list<point>", g, out)
