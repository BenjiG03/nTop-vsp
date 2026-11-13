"""VSPAERO custom blocks: settings, analysis, and result readers."""
from __future__ import annotations

from ntopfile import Graph, InputSpec

from kit import (default_text, AREA, DENSITY, LENGTH, VELOCITY, F_angle, F_choice, F_int, F_length,
                 F_real, F_switch, Field, Kit, analysis, d_choice, d_int, d_point_mm,
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
    F_int("Wake Nodes", "wake_nodes", 64, "Number of nodes along each trailing wake line."),
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
]

# --- Wakes -----------------------------------------------------------------------------------------------------------

WAKE_FIELDS = [
    Field(InputSpec("Trailing Edge Vertices", "list<point>",
                    "The vertices of the Surface Mesh that lie on trailing edges, in any order. Every mesh edge whose two "
                    "ends are both listed becomes part of a wake; connected edges form one wake. List every vertex "
                    "along each trailing edge. An empty list sheds no wakes.",
                    LENGTH, lambda g: g.block("core.list<point>", "list<point>", [g.point(0, 0, 0)], top=False)),
          "te", "points", default=((0.0, 0.0, 0.0),)),
    F_length("Vertex Tolerance", "te_tolerance", 1e-4, "How far a listed point may be from its mesh vertex (the "
             "nearest vertex within this distance is taken; a point with none is an error)."),
    F_switch("Find Body Wakes", "find_body_wakes", False, "OpenVSP's FindBodyWakesFlag (an option of its "
             "experimental mode): OpenVSP also marks wakes behind the blunt bases of its own body components. On an "
             "imported mesh it found none in testing (the HL-20 meshes), so leave it off and list the trailing-edge "
             "vertices instead."),
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
        "and VSPAERO's restart file. Defaults are OpenVSP's, except Wake Nodes (64, as OpenVSP's GUI writes).",
        SOLVER_FIELDS,
        ["**VSPAERO solver settings.** The advanced options of OpenVSP's VSPAERO screen and the vspaero command line "
         "(-omp, -groundheight, -save, -restart), as `key = value` lines for the analysis block."])


def wake_block():
    return settings_block(
        "vspaero_wakes", "VSPAERO Wake Settings (VSP)",
        "Where VSPAERO sheds wakes from the surface mesh. OpenVSP marks wakes only on its own wings, never on a "
        "tessellation, and without wakes VSPAERO (panel or vortex lattice) carries almost no lift. You list the mesh "
        "vertices on the trailing edges; the mesh edges joining them shed the wakes. Check the result with **VSPAERO "
        "Wake Lines (VSP)**.", WAKE_FIELDS,
        ["**Wake settings.** VSPAERO sheds vorticity only from edges listed as wakes in its panel file. The driver "
         "adds them after OpenVSP converts the mesh.",
         "**Trailing-edge vertices.** Each listed point is matched to the mesh vertex within Vertex Tolerance. Every "
         "mesh edge whose two ends are listed becomes a wake edge, and connected wake edges form one wake. The run "
         "stops with an error if a point is not a mesh vertex, if a listed vertex has no listed neighbour, if three "
         "wake edges meet at a vertex, or if all the vertices of one face are listed."])


# --- Analysis --------------------------------------------------------------------------------------------------------

METHODS = ("Panel (thick surface)", "Vortex lattice (thin surface)")
METHOD_KEYS = ("panel", "vlm")

def analysis_block(flow=FLOW_FIELDS, solver=SOLVER_FIELDS, wakes=WAKE_FIELDS):
    own = (
        InputSpec("Surface Mesh", "mesh", "Surface tessellation (a Mesh), used as given: VSPAERO's panels are its "
                  "triangles, so mesh quality is the user's to set. Panel method: the closed, outward-facing surface of "
                  "the vehicle. Vortex lattice: the thin camber surface (an open sheet, e.g. a planform)."),
        InputSpec("Method", "choice", "Panel: VSPAERO's panel method on a thick, closed surface (surface pressures on "
                  "both sides). Vortex lattice: VSPAERO's vortex lattice method on a thin surface (one sheet; the "
                  "result is the pressure jump across it). OpenVSP's thin geometry set builds thin surfaces only from "
                  "its own wings, so the driver passes the mesh through the thick set and marks it thin for VSPAERO.",
                  default=d_choice(METHODS)),
        InputSpec("Reference Area", "real", "Sref.", AREA, d_real(1.0, AREA)),
        InputSpec("Reference Span", "real", "bref.", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Reference Chord", "real", "cref.", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Moment Reference Point", "point", "Moments are taken about this point (VSPAERO's X_cg, Y_cg, Z_cg).",
                  LENGTH, d_point_mm(0, 0, 0)),
        InputSpec("Flow Condition", "text", "From VSPAERO Flow Condition (VSP).", default=d_text(default_text(flow))),
        InputSpec("Solver Settings", "text", "From VSPAERO Solver Settings (VSP).", default=d_text(default_text(solver))),
        InputSpec("Wake Settings", "text", "From VSPAERO Wake Settings (VSP).", default=d_text(default_text(wakes))),
    )
    specs = own + run_inputs("vspaero")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**VSPAERO analysis of a surface mesh.** 1. Export the mesh (STL, metres). 2. Write the settings file. "
              "3. The Python driver has OpenVSP convert the mesh into VSPAERO panels, adds trailing-edge wakes, writes "
              "the `.vspaero` input and runs `vspaero.exe` for every Mach x alpha x beta x Re_cref case. 4. The "
              "output is the path of the results, for the VSPAERO reader blocks.")
    g.comment("**Method.** Panel: the mesh is a closed, thick surface. Vortex lattice: the mesh is a thin sheet; the "
              "driver marks it thin in OpenVSP's `.vkey` file and writes the panel file in the polygon layout OpenVSP "
              "uses for its own wings, so VSPAERO solves on the triangles as given (its own regrouping of triangles "
              "can lose wake edges behind swept trailing edges).")

    def exports(base):
        g.comment("**Geometry.** Export Mesh writes the Surface Mesh as `<Case Name>.stl` in metres; every length "
                  "sent to VSPAERO is in metres too, so the reference quantities match the panels.")
        path = k.cat([base, ".stl"])
        return [g.block("export_mesh<file_path,mesh,unit_length_enum>", "mesh_file_data",
                        [path, inp("Surface Mesh"), g.unit_length("m")], "Export Surface Mesh")]

    lines = [
        k.line("method", k.pick(inp("Method"), METHOD_KEYS)),
        k.line("sref", k.num(inp("Reference Area"), 1.0, AREA)),
        k.line("bref", k.num(inp("Reference Span"), 1.0, LENGTH)),
        k.line("cref", k.num(inp("Reference Chord"), 1.0, LENGTH)),
    ] + [k.line(f"{c}cg", k.num(prop_of(inp("Moment Reference Point"), c), 1.0, LENGTH)) for c in "xyz"] + [
        k.cat(["\n", inp("Flow Condition"), "\n", inp("Solver Settings"), "\n", inp("Wake Settings"), "\n"]),
    ]
    out = analysis(g, k, inp, "run_vspaero.py", exports, lines, "VSPAERO")
    return make_block(
        "vspaero_analysis", "VSPAERO Analysis (VSP)",
        "Runs OpenVSP's **VSPAERO** on a surface mesh, with the panel method (thick, closed surface) or the vortex "
        "lattice method (thin sheet), and returns the path of its results (connect it to the VSPAERO reader blocks). "
        "Settings come from **VSPAERO Flow Condition**, **VSPAERO Solver Settings** and **VSPAERO Wake Settings** "
        "(VSP).\n\nThe mesh is used as given. Panel method: closed and outward-facing. Vortex lattice: one open sheet "
        "(a camber surface or planform); its surface results are the pressure jump across the sheet, and **VSPAERO "
        "Thin Surface Projection (VSP)** shows them on a thick mesh. Either way, wakes leave the trailing-edge "
        "vertices listed in Wake Settings. Reference area, span, chord and moment point are inputs (they are not "
        "taken from the geometry). Files are written to the Working Directory under the Case Name; the run's log is "
        "`<Case Name>.ntvsp.log`.\n\nNot exposed (out of scope): X-Z symmetry, 2D FEM and Tecplot output, reference "
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
                 "Re_1e6", "Case")


def table_column(g, k, path, header):
    table = g.block("import_table<file_path,bool>", "table", [path, g.boolean(True)], "Import Results Table")
    col = g.block("extract_column<table,text>", "list<text>", [table, header], top=False)
    return g.block("scalar_from_text<text>", "real", [col], "Values")


def coefficient_block():
    specs = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
        InputSpec("Coefficient", "choice", "Column of VSPAERO's polar. tot: total; o: skin friction (CDo...); i: "
                  "inviscid surface pressure; w: from the wake (Trefftz plane). Mach, AoA, Beta, Re_1e6 and Case give "
                  "each case's conditions.", default=d_choice(POLAR_COLUMNS)),
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
                  "density). Vortex lattice runs: the jump across the sheet instead, Cp on the side the triangle "
                  "normals point to minus Cp on the other side (negative where the sheet lifts toward its normals), "
                  "or the same jump in pascals.", default=d_choice(SURFACE_QUANTITIES)),
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
                      "p - p_inf (Pa), at panel centroids, nodes, or both. For a vortex lattice run these are the "
                      "pressure jump across the sheet; **VSPAERO Thin Surface Projection (VSP)** maps them onto a "
                      "thick mesh.", specs, "real_point_map", g, out)


def wake_lines_block():
    specs = (InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Wake lines.** The nodes of every trailing edge the driver gave VSPAERO, as a Scalar Point Map whose "
              "value is the wake number. Use it to check that wakes leave the intended edges.")
    path = k.cat([inp("Results"), ".ntvsp.wakes.csv"])
    scale = g.leaf("real_field", {"expression": "1", "units": {}})
    out = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                  [path, g.unit_length("m"), scale], "Wake Lines")
    return make_block("vspaero_wake_lines", "VSPAERO Wake Lines (VSP)",
                      "The trailing edges VSPAERO shed wakes from, as a Scalar Point Map (value: wake number).",
                      specs, "real_point_map", g, out)


# --- Thin surface projection -------------------------------------------------------------------------------------

PROJECTION_QUANTITIES = ("Cp (half the jump, by side)", "p - p_inf (half the jump, by side)", "Pressure jump (Cp)")
PROJECTION_KEYS = ("cp", "dp", "jump")


def projection_block():
    own = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP) run with the vortex lattice method."),
        InputSpec("Thick Mesh", "mesh", "The vehicle's thick surface mesh (closed, outward-facing) to colour. It need "
                  "not match the thin sheet: each of its points takes the value of the closest point of the sheet."),
        InputSpec("Case Number", "integer", "Case, counted from 1 in the order of VSPAERO Coefficient's lists.",
                  default=d_int(1)),
        InputSpec("Maximum Distance", "real", "Points of the thick mesh farther than this from the sheet are left out "
                  "(for parts the sheet does not model, such as fins).", LENGTH, d_real(1000.0, LENGTH)),
        InputSpec("Quantity", "choice", "Cp: each point gets half the jump, with the sign of its side of the sheet "
                  "(+ where its outward normal points the way the sheet's normals do, - otherwise): the lifting part of "
                  "the pressure. p - p_inf: the same in pascals (Flow Condition's velocity and density). Pressure jump: "
                  "the jump itself, unsigned by side.", default=d_choice(PROJECTION_QUANTITIES)),
        InputSpec("Sample Points", "choice", "Nodes, triangle centroids, or both, of the thick mesh.",
                  default=d_choice(SURFACE_POINTS, 2)),
    )
    specs = own + run_inputs("hl20_projection")
    g = Graph()
    k = Kit(g)
    model = model_inputs(specs)

    def inp(name):  # no OpenVSP program runs here; helper commands run in the driver folder
        return model("Driver Directory" if name == "OpenVSP Directory" else name)

    g.comment("**Thin to thick projection.** A vortex lattice gives only the pressure jump across its sheet. The "
              "driver `run_project.py` finds, for every node and centroid of the thick mesh, the closest point of the "
              "sheet, interpolates the jump there, and gives the point half of it with the sign of its side: Cp = "
              "+jump/2 where the thick surface faces the way the sheet normals point, -jump/2 on the other side. That "
              "is the lifting part of the pressure only; a thin solution has no thickness effect.")

    def exports(base):
        g.comment("**Thick mesh.** Exported as `<Case Name>.stl` in metres for the driver.")
        return [g.block("export_mesh<file_path,mesh,unit_length_enum>", "mesh_file_data",
                        [k.cat([base, ".stl"]), inp("Thick Mesh"), g.unit_length("m")], "Export Thick Mesh")]

    lines = [k.line("vspaero_results", inp("Results")), k.line("case_number", k.num(inp("Case Number"))),
             k.line("max_distance", k.num(inp("Maximum Distance"), 1.0, LENGTH))]
    base = analysis(g, k, inp, "run_project.py", exports, lines, "Thin Surface Projection")
    g.comment("**Point map.** The driver writes every quantity at nodes, centroids and both; the block imports the "
              "chosen one as a Scalar Point Map (x, y, z in metres).")
    path = k.cat([base, ".ntvsp.proj.", k.pick(inp("Quantity"), PROJECTION_KEYS), ".",
                  k.pick(inp("Sample Points"), ("centroids", "nodes", "both")), ".csv"], name="Point Map File")
    scale = g.leaf("real_field", {"expression": "1", "units": {}})
    out = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                  [path, g.unit_length("m"), scale], "Projected Pressure")
    return make_block(
        "vspaero_thin_projection", "VSPAERO Thin Surface Projection (VSP)",
        "Shows a vortex lattice (thin surface) VSPAERO solution on the vehicle's thick mesh, for visualization: a "
        "**Scalar Point Map** at the thick mesh's points. Each point takes the pressure jump at the closest point of "
        "the thin sheet; Cp and p - p_inf give it half the jump with the sign of its side (the lifting part of the "
        "pressure), Pressure jump gives the jump itself. Points beyond Maximum Distance from the sheet are left out. "
        "Coefficients still come from **VSPAERO Coefficient (VSP)**.", specs, "real_point_map", g, out)
