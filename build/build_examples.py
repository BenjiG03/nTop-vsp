"""Write (and optionally run) the HL-20 example notebooks.

    sh build/py.sh build/build_examples.py [--run NAME ...]

Each notebook embeds the custom blocks it uses, so it opens in nTop as it is. ``--run`` runs it with ntopcl
(``--trustnotebook``, since the blocks use Run Command) and prints its outputs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import drag_blocks as d  # noqa: E402
import vspaero_blocks as v  # noqa: E402
from kit import ANGLE, AREA, DENSITY, LENGTH, ROOT, VELOCITY  # noqa: E402
from ntopfile import Graph, Ref, save_host  # noqa: E402
from ntrun import NTOPCL  # noqa: E402

# The viscous correction is ntop_newtonian's (nTop blocks only); its builder sits beside this repository.
NEWTONIAN_DIR = os.path.join(os.path.dirname(ROOT), "ntop_newtonian")

EXAMPLES = os.path.join(ROOT, "examples")
RUNS = ROOT + "/runs"
NTOP_MESH = ROOT + "/h20_ntop_mesh.stl"  # exported from nTop in millimetres
VSP_MODEL = ROOT + "/hl20.vsp3"

# Trailing-edge vertices of h20_ntop_mesh.stl (millimetres, as in the file): the vertices along the trailing edges
# of the HL-20's tip fins and vertical tail. Each fin and the tail end in a blunt base; these vertices run along one
# corner of each base (see the README).
TE_VERTICES = [
    # tip fin at +y
    (20971.9082, 7426.3037, 5219.8213), (20924.2949, 7424.5347, 5061.6953), (20814.5254, 7259.1743, 4848.5864),
    (20702.4785, 7091.3218, 4632.2314), (20581.5371, 6909.0264, 4397.3018), (20460.7500, 6727.5806, 4163.4434),
    (20340.3730, 6546.7534, 3930.3840), (20219.9473, 6367.3438, 3699.0898), (20124.2539, 6222.0962, 3511.9458),
    (20028.8359, 6078.7593, 3327.2043), (19912.9199, 5904.6362, 3102.7834), (19797.8770, 5731.8105, 2880.0352),
    (19683.3906, 5560.0225, 2658.6172), (19566.4375, 5384.7144, 2432.6548), (19444.7129, 5203.7227, 2199.3096),
    (19336.9258, 5039.8564, 1988.1854), (19228.2090, 4876.2520, 1777.3334), (19024.3086, 4559.1123, 1685.8671),
    # tip fin at -y
    (19031.4453, -4574.0435, 1706.3049), (19144.4629, -4733.9282, 1907.3702), (19279.4551, -4930.1904, 2133.6677),
    (19466.7930, -5234.4966, 2239.0676), (19577.7480, -5402.7520, 2455.8625), (19690.9766, -5571.3501, 2673.2202),
    (19803.8945, -5740.8516, 2891.6880), (19910.7871, -5901.4316, 3098.6523), (20017.6875, -6062.0088, 3305.6145),
    (20124.5820, -6222.5874, 3512.5789), (20232.6699, -6384.9595, 3721.8525), (20349.6758, -6560.7261, 3948.3928),
    (20467.2305, -6737.3164, 4175.9927), (20585.3145, -6914.7085, 4404.6250), (20704.1484, -7093.9189, 4635.5732),
    (20820.6367, -7268.5864, 4860.7080), (20928.2363, -7430.6177, 5069.5288), (20971.2969, -7413.2900, 5230.3086),
    # vertical tail
    (22839.4688, -118.0499, 3694.1597), (22657.9941, -116.5693, 3606.5386), (22299.3145, -117.4662, 3433.4231),
    (21940.6328, -117.0558, 3260.3110), (21662.9434, -116.7380, 3126.2871), (21388.2852, -116.4237, 2993.7266),
    (21113.6270, -116.1094, 2861.1660), (20845.5938, -115.8027, 2731.8042), (20575.9453, -115.4941, 2601.6621),
    (20303.8809, -115.1828, 2470.3538), (20028.3242, -114.8675, 2337.3604), (19747.9395, -114.5467, 2202.0374),
    (19461.2637, -114.2187, 2063.6765), (19162.7051, -113.8771, 1919.5811), (18953.6777, -113.6379, 1818.6980),
    (18877.2461, 0.0000, 1781.7574), (18862.1055, 123.3809, 1774.3940),
]


# The thin HL-20 for the vortex lattice: a planform sheet at z = 0 (millimetres), not the same model as NTOP_MESH.
VLM_MESH = ROOT + "/h20_vlm_mesh_from ntop.stl"

# Trailing-edge vertices of VLM_MESH (millimetres): the two wing trailing edges (x about 16.4 to 16.9 m) and the
# body's base edge (x = 20.64 m), boundary vertices of the sheet.
VLM_TE_VERTICES = [
    (16439.3613, -3841.5801, 0.0000), (16614.8105, -3596.3240, 0.0000), (16742.0059, -3322.9771, 0.0000),
    (16822.1855, -3031.9111, 0.0000), (16856.9727, -2732.9224, 0.0000), (16871.7324, -2431.1794, 0.0000),
    (16859.3633, -2129.1995, 0.0000), (16829.9824, -1828.6278, 0.0000), (16772.3750, -1532.1338, 0.0000),
    (16716.3008, -1235.1611, 0.0000), (16652.6250, -939.7012, 0.0000), (16652.6250, 939.7012, 0.0000),
    (16716.3008, 1235.1613, 0.0000), (16772.3750, 1532.1340, 0.0000), (16829.9824, 1828.6281, 0.0000),
    (16859.3633, 2129.2000, 0.0000), (16871.7324, 2431.1802, 0.0000), (16856.9727, 2732.9231, 0.0000),
    (16822.1855, 3031.9119, 0.0000), (16742.0039, 3322.9778, 0.0000), (16614.8105, 3596.3247, 0.0000),
    (16439.3613, 3841.5811, 0.0000), (20638.2422, -1457.4578, 0.0000), (20638.2422, -1165.9663, 0.0000),
    (20638.2422, -874.4747, 0.0000), (20638.2422, -582.9832, 0.0000), (20638.2422, -291.4916, 0.0000),
    (20638.2422, -0.0000, 0.0000), (20638.2422, 291.4916, 0.0000), (20638.2422, 582.9832, 0.0000),
    (20638.2422, 874.4747, 0.0000), (20638.2422, 1165.9663, 0.0000), (20638.2422, 1457.4578, 0.0000),
]


def call(g: Graph, block, given, name):
    """An instance of a custom block: ``given`` maps input names to instances (or callables of g); defaults elsewhere."""
    inputs = []
    for spec in block.inputs:
        if spec.name in given:
            value = given[spec.name]
            inputs.append(value(g) if callable(value) else value)
        else:
            inputs.append(spec.default(g))
    return g.block(block.signature, block.output_type, inputs, name)


def choose(labels, label):
    labels = list(labels)
    return lambda g: g.choice(labels, labels.index(label))


def import_mesh(g):
    path = g.file_path(NTOP_MESH)
    return g.block("import_mesh<file_path,unit_length_enum>[1.1.0]", "mesh", [path, g.unit_length("mm")], "HL-20 Mesh")


# --- VSPAERO ---------------------------------------------------------------------------------------------------------

def vspaero_example():
    flow, solver, wakes = v.flow_block(), v.solver_block(), v.wake_block()
    ana, coef, surf, lines = v.analysis_block(), v.coefficient_block(), v.surface_block(), v.wake_lines_block()
    blocks = [flow, solver, wakes, ana, coef, surf, lines]
    g = Graph()
    g.comment("# HL-20 with VSPAERO\nThis notebook runs OpenVSP's VSPAERO panel method on the HL-20 lifting body, "
              "meshed in nTop, and shows the surface pressure as a Scalar Point Map. The conditions match the "
              "OpenVSP run that came with the model (`hl20.vspaero`): Mach 0, alpha 1, 5.5 and 10 degrees, "
              "Sref = 100, cref = bref = 1, moments about the origin.\n\nThe blocks call Python and OpenVSP with Run "
              "Command: allow them when nTop asks. A run takes a few minutes; its log is `runs/hl20_vspaero.ntvsp.log`.")
    g.comment("**1. Geometry.** The surface mesh nTop made of the HL-20 (Import Mesh, millimetres). VSPAERO uses its "
              "triangles as panels, so the mesh is the user's to get right: closed, outward-facing, and fine where "
              "the pressure changes quickly.")
    mesh = import_mesh(g)
    g.comment("**2. Flow condition.** VSPAERO sweeps every combination of the alpha, beta, Mach and Re_cref lists in "
              "one run. Here alpha goes from 1 to 10 degrees in 3 points at Mach 0.")
    fc = call(g, flow, {"Alpha Start": lambda g: g.real(math.radians(1.0), ANGLE),
                        "Alpha End": lambda g: g.real(math.radians(10.0), ANGLE),
                        "Alpha Points": lambda g: g.integer(3)}, "Flow Condition")
    g.comment("**3. Solver settings.** OpenVSP's defaults with 8 threads.")
    ss = call(g, solver, {"Processors": lambda g: g.integer(8)}, "Solver Settings")
    g.comment("**4. Wakes.** VSPAERO needs to know which mesh edges shed wakes, and OpenVSP marks none on a mesh. "
              "Trailing Edge Vertices lists the mesh vertices along the trailing edges of the two tip fins and the "
              "vertical tail (coordinates of vertices of the imported mesh); the mesh edges joining them become the "
              "wakes. The HL-20's fins and tail end in a blunt base, and these vertices run along one corner of it: "
              "the result depends on which corner (see the README).")
    pts = [g.point(x, y, z) for x, y, z in TE_VERTICES]
    te_points = g.block("core.list<point>", "list<point>", pts, "Trailing Edge Vertices")
    ws = call(g, wakes, {"Trailing Edge Vertices": te_points}, "Wake Settings")
    g.comment("**5. Analysis.** Reference quantities are inputs (here OpenVSP's defaults, as in `hl20.vspaero`). The "
              "output is the path of the results, which the reader blocks below take.")
    res = call(g, ana, {"Surface Mesh": mesh, "Reference Area": lambda g: g.real(100.0, AREA),
                        "Reference Span": lambda g: g.real(1.0, LENGTH), "Reference Chord": lambda g: g.real(1.0, LENGTH),
                        "Flow Condition": fc, "Solver Settings": ss, "Wake Settings": ws,
                        "Working Directory": lambda g: g.file_path(RUNS), "Case Name": lambda g: g.text("hl20_vspaero")},
               "HL-20 VSPAERO")
    g.comment("**6. Coefficients.** One list per coefficient, one value per case (here the three angles of attack).")
    outs = {}
    for name in ("AoA", "CLtot", "CDtot", "CMytot"):
        outs[name] = call(g, coef, {"Results": res, "Coefficient": choose(v.POLAR_COLUMNS, name)}, name)
    g.comment("**7. Contours.** Surface Cp of case 2 (alpha 5.5 degrees) as a Scalar Point Map, at panel centroids and "
              "nodes. Show it and set the colour range in its HUD. Field from Point Map turns it into a field.")
    cp = call(g, surf, {"Results": res, "Case Number": lambda g: g.integer(2)}, "Surface Cp (alpha 5.5)")
    field = g.block("field_from_point_map<real_point_map,point_cloud_interpolation_enum,point_cloud_extrapolation_enum>",
                    "real_field", [cp, g.enum("point_cloud_interpolation_enum", 0), g.enum("point_cloud_extrapolation_enum", 0)],
                    "Cp Field")
    g.comment("**8. Wake check.** The trailing-edge nodes the wakes left from, coloured by wake number.")
    wl = call(g, lines, {"Results": res}, "Wake Lines")
    group = g.block("core.group<any>", "group", [outs["AoA"], outs["CLtot"], outs["CDtot"], outs["CMytot"], cp, wl], "Outputs")
    return g, blocks, group


def vlm_example():
    flow, solver, wakes = v.flow_block(), v.solver_block(), v.wake_block()
    ana, coef, surf, lines, proj = (v.analysis_block(), v.coefficient_block(), v.surface_block(), v.wake_lines_block(),
                                    v.projection_block())
    blocks = [flow, solver, wakes, ana, coef, surf, lines, proj]
    g = Graph()
    g.comment("# HL-20 with VSPAERO's vortex lattice\nThis notebook runs VSPAERO's vortex lattice method on a thin "
              "model of the HL-20: a planform sheet at z = 0, meshed in nTop. It is not the same geometry as the thick "
              "mesh of the panel example, so the numbers differ from it. The conditions match that example: Mach 0, "
              "alpha 1, 5.5 and 10 degrees, Sref = 100, cref = bref = 1, moments about the origin. The thin solution "
              "is then shown on the thick mesh.\n\nThe blocks call Python and OpenVSP with Run Command: allow them "
              "when nTop asks. Logs: `runs/hl20_vlm.ntvsp.log` and `runs/hl20_vlm_projection.ntvsp.log`.")
    g.comment("**1. Geometry.** The thin sheet (Import Mesh, millimetres): one open surface whose triangles become "
              "the vortex lattice. The thick mesh is only for the projection in step 7.")
    sheet = g.block("import_mesh<file_path,unit_length_enum>[1.1.0]", "mesh",
                    [g.file_path(VLM_MESH), g.unit_length("mm")], "HL-20 Thin Mesh")
    thick = import_mesh(g)
    g.comment("**2. Flow condition and solver settings**, as in the panel example.")
    fc = call(g, flow, {"Alpha Start": lambda g: g.real(math.radians(1.0), ANGLE),
                        "Alpha End": lambda g: g.real(math.radians(10.0), ANGLE),
                        "Alpha Points": lambda g: g.integer(3)}, "Flow Condition")
    ss = call(g, solver, {"Processors": lambda g: g.integer(8)}, "Solver Settings")
    g.comment("**3. Wakes.** The vertices along the sheet's trailing edges: the two wing trailing edges and the body's "
              "base edge. The edges of the sheet joining them shed the wakes.")
    pts = [g.point(x, y, z) for x, y, z in VLM_TE_VERTICES]
    te_points = g.block("core.list<point>", "list<point>", pts, "Trailing Edge Vertices")
    ws = call(g, wakes, {"Trailing Edge Vertices": te_points}, "Wake Settings")
    g.comment("**4. Analysis** with Method = Vortex lattice.")
    res = call(g, ana, {"Surface Mesh": sheet, "Method": choose(v.METHODS, "Vortex lattice (thin surface)"),
                        "Reference Area": lambda g: g.real(100.0, AREA),
                        "Reference Span": lambda g: g.real(1.0, LENGTH), "Reference Chord": lambda g: g.real(1.0, LENGTH),
                        "Flow Condition": fc, "Solver Settings": ss, "Wake Settings": ws,
                        "Working Directory": lambda g: g.file_path(RUNS), "Case Name": lambda g: g.text("hl20_vlm")},
               "HL-20 VSPAERO VLM")
    g.comment("**5. Coefficients**, one value per angle of attack.")
    outs = {}
    for name in ("AoA", "CLtot", "CDtot", "CMytot"):
        outs[name] = call(g, coef, {"Results": res, "Coefficient": choose(v.POLAR_COLUMNS, name)}, name)
    g.comment("**6. Pressure jump on the sheet** for case 2 (alpha 5.5 degrees): Cp on the side the triangle normals "
              "point to (here +z, the upper side) minus Cp on the lower side, so negative where the sheet lifts.")
    jump = call(g, surf, {"Results": res, "Case Number": lambda g: g.integer(2)}, "Pressure Jump (alpha 5.5)")
    g.comment("**7. Thin solution on the thick mesh.** Each point of the thick HL-20 mesh takes half the jump at the "
              "closest point of the sheet, with the sign of its side: suction on the upper surface, pressure on the "
              "lower surface. This is the lifting part of the pressure only. The fins and the vertical tail are not "
              "in the thin model; Maximum Distance (3 m) leaves out their outer parts (the body lies within about 2.5 m of the sheet).")
    cp = call(g, proj, {"Results": res, "Thick Mesh": thick, "Case Number": lambda g: g.integer(2),
                        "Maximum Distance": lambda g: g.real(3.0, LENGTH),
                        "Working Directory": lambda g: g.file_path(RUNS),
                        "Case Name": lambda g: g.text("hl20_vlm_projection")}, "Cp on Thick Mesh (alpha 5.5)")
    g.comment("**8. Wake check.** The trailing-edge nodes the wakes left from, coloured by wake number.")
    wl = call(g, lines, {"Results": res}, "Wake Lines")
    group = g.block("core.group<any>", "group",
                    [outs["AoA"], outs["CLtot"], outs["CDtot"], outs["CMytot"], jump, cp, wl], "Outputs")
    return g, blocks, group

# --- Parasite drag ---------------------------------------------------------------------------------------------------

def parasite_example():
    comp, pd, result, values = d.component_block(), d.parasite_block(), d.parasite_result_block(), d.parasite_components_block()
    blocks = [comp, pd, result, values]
    g = Graph()
    g.comment("# HL-20 parasite drag\nOpenVSP's Parasite Drag tool on the HL-20's OpenVSP model (`hl20.vsp3`). The "
              "tool needs OpenVSP's own components, so its geometry input is the model file rather than an nTop "
              "mesh: OpenVSP takes each component's wetted area, reference length and thickness or fineness ratio "
              "from its parametric definition.")
    g.comment("**1. Component settings.** One Parasite Drag Component block per component to change; the others keep "
              "the settings saved in the model. Here the wings get 10% laminar flow and an interference factor of 1.1.")
    wings = call(g, comp, {"Component Name": lambda g: g.text("Wings"),
                           "Percent Laminar": lambda g: g.real(10.0), "Interference Factor": lambda g: g.real(1.1)},
                 "Wings Settings")
    comps = g.block("core.list<text>", "list<text>", [wings], "Components")
    g.comment("**2. Parasite drag.** Model dimensions in metres (as the HL-20 was brought into nTop), Sref = 100 m^2, "
              "Mach 0.5 at 6096 m (20,000 ft) in the US Standard Atmosphere 1976, Schlichting compressible "
              "turbulent friction.")
    res = call(g, pd, {"OpenVSP Model": lambda g: g.file_path(VSP_MODEL),
                       "Model Length Unit": choose(d.MODEL_UNITS, "m"),
                       "Reference Area": lambda g: g.real(100.0, AREA),
                       "Speed": lambda g: g.real(0.5), "Speed Unit": choose(d.SPEED_UNITS, "Mach"),
                       "Altitude": lambda g: g.real(6096.0, LENGTH), "Components": comps,
                       "Working Directory": lambda g: g.file_path(RUNS), "Case Name": lambda g: g.text("hl20_parasite")},
               "HL-20 Parasite Drag")
    g.comment("**3. Results.** The total CD0 and, per component, CD and wetted area (in OpenVSP's table order).")
    total = call(g, result, {"Results": res, "Quantity": choose(d.PARASITE_RESULTS, "Total_CD")}, "CD0")
    cd = call(g, values, {"Results": res, "Column": choose(d.PARASITE_COMPONENT_COLUMNS, "CD")}, "Component CD")
    swet = call(g, values, {"Results": res, "Column": choose(d.PARASITE_COMPONENT_COLUMNS, "Swet")}, "Component Swet")
    group = g.block("core.group<any>", "group", [total, cd, swet], "Outputs")
    return g, blocks, group


# --- Wave drag -------------------------------------------------------------------------------------------------------

def wave_example():
    wave, result, areas = d.wave_block(), d.wave_result_block(), d.wave_areas_block()
    blocks = [wave, result, areas]
    g = Graph()
    g.comment("# HL-20 wave drag\nOpenVSP's Wave Drag (supersonic area rule) on the HL-20 surface mesh made in nTop.")
    g.comment("**1. Geometry.** The nTop mesh of the HL-20 (millimetres); the block exports it in metres for OpenVSP.")
    mesh = import_mesh(g)
    g.comment("**2. Wave drag** at Mach 2 with OpenVSP's defaults: 20 slices, 10 rotation sections, x-z symmetry. "
              "Sref = 100 m^2.")
    res = call(g, wave, {"Surface Mesh": mesh, "Mach": lambda g: g.real(2.0), "Reference Area": lambda g: g.real(100.0, AREA),
                         "Working Directory": lambda g: g.file_path(RUNS), "Case Name": lambda g: g.text("hl20_wave")},
               "HL-20 Wave Drag")
    g.comment("**3. Results.** CD_wave, and the equivalent area distribution of the first rotation section (roll angle "
              "0): x/L stations and areas.")
    cdw = call(g, result, {"Results": res}, "CD Wave")
    x = call(g, areas, {"Results": res, "Value": choose(("X_Norm", "Area"), "X_Norm")}, "Area Stations x/L")
    a = call(g, areas, {"Results": res, "Value": choose(("X_Norm", "Area"), "Area")}, "Areas")
    group = g.block("core.group<any>", "group", [cdw, x, a], "Outputs")
    return g, blocks, group


# --- Drag buildup ---------------------------------------------------------------------------------------------------

# Flight condition of the buildup example: Mach 0.5 at 6096 m (20,000 ft) in the US Standard Atmosphere 1976, as
# the Parasite Drag example. Subsonic, because VSPAERO's panel method gave no usable supersonic result on the HL-20
# (see the README); so there is no wave drag and the viscous correction is ntop_newtonian's subsonic block.
BUILDUP_MACH = 0.5
BUILDUP_ALTITUDE = 6096.0  # m
BUILDUP_T = 288.15 - 0.0065 * BUILDUP_ALTITUDE  # K, troposphere
BUILDUP_P = 101325.0 * (BUILDUP_T / 288.15) ** (9.80665 / (0.0065 * 287.05287))  # Pa
GAS_R, GAMMA = 287.05, 1.4


def var(g, typ, leaf, name):
    """A named notebook variable holding leaf."""
    return g.block(f"core.var<{typ}>", typ, [leaf], name)


def buildup_conditions():
    """(Vinf m/s, rho kg/m^3, Re per metre) of the buildup flight condition, for VSPAERO's dimensional inputs."""
    vinf = BUILDUP_MACH * math.sqrt(GAMMA * GAS_R * BUILDUP_T)
    rho = BUILDUP_P / (GAS_R * BUILDUP_T)
    mu = 1.458e-6 * BUILDUP_T ** 1.5 / (BUILDUP_T + 110.4)
    return vinf, rho, rho * vinf / mu


def buildup_example():
    if NEWTONIAN_DIR not in sys.path:
        sys.path.append(NEWTONIAN_DIR)
    import viscous  # ntop_newtonian

    flow, solver, wakes = v.flow_block(), v.solver_block(), v.wake_block()
    ana, surf = v.analysis_block(), v.surface_block()
    friction, buildup = viscous.block("Coefficients", subsonic=True), d.buildup_block()
    blocks = [flow, solver, wakes, ana, surf, friction, buildup]
    vinf, rho, re_per_m = buildup_conditions()
    cref = 1.0
    g = Graph()
    g.comment("# HL-20 drag buildup\nCL, CD and CMy of the HL-20 at Mach 0.5 (6096 m, US Standard Atmosphere 1976) "
              "from two analyses that share one flight condition and one reference: VSPAERO's panel method (inviscid "
              "pressure) and the nTop-native subsonic viscous correction from `ntop_newtonian` (skin friction on "
              "VSPAERO's surface pressure). **Drag Buildup (VSP)** adds them and exports the table "
              "`runs/hl20_buildup.ntvsp.buildup.csv`. Subsonic, so there is no wave drag: Wave Drag Results is left "
              "empty (above Mach 1, connect a Wave Drag (VSP) run at the same Mach number and Sref, and use the "
              "supersonic viscous block).\n\nThe blocks call Python and OpenVSP with Run Command: allow them when "
              "nTop asks.")
    g.comment("**1. Shared flight condition and reference.** Every analysis takes these values from here, so they "
              "cannot drift apart. Moments are about the origin over Sref and cref (OpenVSP's defaults, as in "
              "`hl20.vspaero`).")
    mach = var(g, "real", g.real(BUILDUP_MACH), "Mach")
    a0 = var(g, "real", g.real(math.radians(1.0), ANGLE), "Alpha Start")
    a1 = var(g, "real", g.real(math.radians(10.0), ANGLE), "Alpha End")
    na = var(g, "integer", g.integer(3), "Alpha Points")
    sref = var(g, "real", g.real(100.0, AREA), "Reference Area")
    cref_i = var(g, "real", g.real(cref, LENGTH), "Reference Chord")
    center = var(g, "point", g.point(0.0, 0.0, 0.0), "Moment Reference Point")
    g.comment("**Freestream state** at 6096 m: static pressure and temperature for the viscous correction. VSPAERO's "
              f"dimensional inputs are set to match: Vinf = {vinf:.1f} m/s, rho = {rho:.5f} kg/m^3, Re_cref = "
              f"{re_per_m * cref:.4g}. They do not change VSPAERO's coefficients, and its own viscous estimate (the only "
              "thing Re_cref changes) is left out of the buildup.")
    p_inf = var(g, "real", g.real(BUILDUP_P, viscous.PRESSURE), "Freestream Pressure")
    t_inf = var(g, "real", g.real(BUILDUP_T, viscous.TEMPERATURE), "Freestream Temperature")

    g.comment("**2. Geometry.** VSPAERO uses the nTop mesh of the HL-20 (millimetres) as given. The viscous correction "
              "integrates over nTop's implicit HL-20 (`h20.implicit`, the body the mesh was made from): its surface "
              "integrals run about a third faster than on an implicit made from the mesh, with the same result "
              "(friction CD within 0.2%).")
    mesh = import_mesh(g)
    body = g.block("import_implicit<file_path>", "implicit", [g.file_path(ROOT + "/h20.implicit")], "HL-20 Body")

    g.comment("**3. VSPAERO** at the shared Mach number and angles, with the trailing-edge vertices of the VSPAERO "
              "example.")
    fc = call(g, flow, {"Mach Start": mach, "Mach End": mach, "Alpha Start": a0, "Alpha End": a1, "Alpha Points": na,
                        "Re_cref Start": lambda g: g.real(re_per_m * cref),
                        "Re_cref End": lambda g: g.real(re_per_m * cref),
                        "Freestream Velocity": lambda g: g.real(vinf, VELOCITY),
                        "Freestream Density": lambda g: g.real(rho, DENSITY)}, "Flow Condition")
    ss = call(g, solver, {"Processors": lambda g: g.integer(8)}, "Solver Settings")
    pts = [g.point(x, y, z) for x, y, z in TE_VERTICES]
    te_points = g.block("core.list<point>", "list<point>", pts, "Trailing Edge Vertices")
    ws = call(g, wakes, {"Trailing Edge Vertices": te_points}, "Wake Settings")
    aero = call(g, ana, {"Surface Mesh": mesh, "Reference Area": sref, "Reference Chord": cref_i,
                         "Reference Span": lambda g: g.real(1.0, LENGTH), "Moment Reference Point": center,
                         "Flow Condition": fc, "Solver Settings": ss, "Wake Settings": ws,
                         "Working Directory": lambda g: g.file_path(RUNS),
                         "Case Name": lambda g: g.text("hl20_buildup_vspaero")}, "HL-20 VSPAERO")

    g.comment("**4. Viscous correction.** VSPAERO's surface Cp of each case (at the mesh nodes) becomes a field (Field "
              "from Point Map), and ntop_newtonian's Viscous Correction Subsonic Coefficients block computes skin "
              "friction on the implicit body at the edge state that Cp gives, at that case's angle of attack (from the "
              "shared Alpha inputs). One call per case gives all eight friction coefficients; CD, CL and Cm are picked "
              "from each list. Reference Length is cref, so its Cm matches VSPAERO's CMy. Fully turbulent, adiabatic "
              "wall. Each call is six surface integrals over the HL-20: several minutes per case.")
    alphas = g.block("sequence_from_bounds<real,real,integer>", "list<real>", [a0, a1, na], "Angles of Attack")
    shared = {"Body": body, "Mach": mach, "Reference Area": sref, "Reference Length": cref_i, "Moment Center": center,
              "Freestream Pressure": p_inf, "Freestream Temperature": t_inf}
    per_case = []
    for i in range(3):
        cp_map = call(g, surf, {"Results": aero, "Case Number": lambda g, n=i + 1: g.integer(n),
                                "Sample Points": choose(v.SURFACE_POINTS, "Nodes")}, f"Surface Cp, Case {i + 1}")
        cp_field = g.block("field_from_point_map<real_point_map,point_cloud_interpolation_enum,point_cloud_extrapolation_enum>",
                           "real_field", [cp_map, g.enum("point_cloud_interpolation_enum", 0),
                                          g.enum("point_cloud_extrapolation_enum", 0)], f"Cp Field, Case {i + 1}")
        per_case.append(call(g, friction, {**shared, "Pressure Coefficient": cp_field,
                                           "Angle of Attack": Ref(alphas, (i,))}, f"Viscous, Case {i + 1}"))
    visc = {}
    for label in ("CD", "CL", "Cm"):
        k = viscous.COEFFICIENTS.index(label)
        visc[label] = g.block("core.list<real>", "list<real>", [Ref(c, (k,)) for c in per_case], f"Viscous {label}")

    g.comment("**5. Buildup.** CL = CLi + CL_viscous, CD = CDi + CD_viscous, CMy = CMiy + Cm_viscous per case (no "
              "wave drag below Mach 1). The block writes the table to `runs/hl20_buildup.ntvsp.buildup.csv` and "
              "outputs it.")
    table = call(g, buildup, {"VSPAERO Results": aero, "Viscous CD": visc["CD"],
                              "Viscous CL": visc["CL"], "Viscous CMy": visc["Cm"],
                              "Working Directory": lambda g: g.file_path(RUNS),
                              "Case Name": lambda g: g.text("hl20_buildup")}, "HL-20 Drag Buildup")
    group = g.block("core.group<any>", "group", [table, visc["CD"], visc["CL"], visc["Cm"]], "Outputs")
    return g, blocks, group


EXAMPLE_BUILDERS = {
    "HL-20 VSPAERO": vspaero_example,
    "HL-20 VSPAERO VLM": vlm_example,
    "HL-20 Parasite Drag": parasite_example,
    "HL-20 Wave Drag": wave_example,
    "HL-20 Drag Buildup": buildup_example,
}


def run_notebook(path, timeout=3600):
    import subprocess
    res = os.path.splitext(path)[0] + ".json"
    if os.path.exists(res):
        os.remove(res)
    t0 = time.time()
    log_path = os.path.splitext(path)[0] + ".log"
    # ntopcl writes its output file and then may never exit after notebooks that use Run Command inside custom
    # blocks, so wait for the output file to settle and stop it.
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        p = subprocess.Popen([NTOPCL, "--trustnotebook", "-v", "2", "-o", res, path], stdout=logf,
                             stderr=subprocess.STDOUT)
        last = -1
        while p.poll() is None and time.time() - t0 < timeout:
            time.sleep(2)
            if os.path.exists(res):
                size = os.path.getsize(res)
                if size == last and size > 0:
                    break
                last = size
        if p.poll() is None:
            p.kill()
            p.wait()
    log = open(log_path, encoding="utf-8", errors="replace").read()
    out = json.load(open(res)) if os.path.exists(res) else None
    return out, log, time.time() - t0


def summarize(out):
    if not out:
        return "no output"
    lines = []
    for c in out[0].get("components", []):
        val = c.get("value")
        if c["type"] in ("real", "list<real>"):
            vals = val.get("value") if isinstance(val, dict) and "value" in val else [val]
            nums = [x.get("val") for x in vals if isinstance(x, dict)]
            lines.append(f"  {c.get('name', c['type'])}: {nums}")
        else:
            lines.append(f"  {c.get('name', '')}: {c['type']}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", nargs="*", default=None, help="names to run (default: none; empty: all)")
    args = parser.parse_args()
    os.makedirs(EXAMPLES, exist_ok=True)
    paths = {}
    for name, builder in EXAMPLE_BUILDERS.items():
        g, blocks, group = builder()
        paths[name] = save_host(g, blocks, os.path.join(EXAMPLES, name + ".ntop"), output=group, output_type="group",
                                name=name)
        print(paths[name])
    if args.run is not None:
        for name in args.run or list(EXAMPLE_BUILDERS):
            out, log, dt = run_notebook(paths[name])
            errs = [l for l in log.splitlines() if "[E]" in l and "version of your file" not in l]
            print(f"== {name}: {dt:.0f} s, {'ok' if out else 'FAILED'}")
            print(summarize(out))
            if errs:
                print("\n".join(errs[:15]))


if __name__ == "__main__":
    main()
