"""OpenVSP CAD Export (VSP): an OpenVSP model (from OpenVSP Model) written as STEP or IGES by OpenVSP itself."""
from __future__ import annotations

from ntopfile import Graph, InputSpec

from kit import F_choice, F_int, F_real, F_switch, Kit, analysis, make_block, model_inputs, run_inputs

METHODS = ("Trimmed (intersected, one solid)", "Stitched (one solid per component)", "Untrimmed Surfaces")
FORMATS = ("STEP", "IGES")

CAD_FIELDS = (
    F_choice("Method", "method", METHODS, ("trimmed", "stitched", "surfaces"), 0,
             "How OpenVSP writes its surfaces. Trimmed: OpenVSP's Surface Intersection intersects the components "
             "with each other, trims away the parts inside other components and writes trimmed surfaces: the closed "
             "ones become one BREP solid (STEP with Representation BREP), open ones (a duct's wall) stay surfaces. "
             "Stitched: each component's patches are split at feature lines and stitched into its own watertight "
             "BREP (an open shell for an open surface); components are not intersected with each other. Untrimmed "
             "Surfaces: each component's surface patches, no topology."),
    F_choice("Format", "format", FORMATS, ("step", "iges"), 0, "STEP (.stp) or IGES (.igs)."),
    F_choice("Representation", "representation", ("BREP Solid", "Shell"), ("brep", "shell"), 0,
             "Trimmed STEP only: a BREP solid (MANIFOLD_SOLID_BREP) for closed surfaces, or a shell-based surface "
             "model."),
    F_real("Tolerance", "tolerance", 1e-6, None, "Tolerance written to the STEP file, which trimmed curves are also "
           "built to (OpenVSP's STEP Tol), in metres. The file is in metres (OpenVSP Model's unit); CAD programs "
           "convert on import."),
    F_switch("Demote to Cubic", "to_cubic", False, "Untrimmed and stitched: demote higher-order surfaces to cubic "
             "(OpenVSP's To Cubic), within 1e-6."),
    F_switch("Merge Points", "merge_points", False, "STEP: merge coincident points (OpenVSP's Merge Points)."),
    F_switch("Merge LE/TE Patches", "merge_le_te", False, "Untrimmed and stitched: merge wing leading- and "
             "trailing-edge patches (OpenVSP's Merge LE/TE)."),
    F_switch("Label Names", "label_names", True, "Name the CAD surfaces after the OpenVSP components."),
    F_switch("Planar End Caps", "planar_caps", False, "STEP only. Replace each flat end cap (a blunt fuselage end, a "
             "wing's root and tip) with one plane face bounded by its rim. OpenVSP builds these caps from patches "
             "collapsed to a point (the cap's centre, a wing cap's leading and trailing edges), which nTop's CAD import "
             "reports as degenerate geometry; the plane face has the same shape. A pointed tip keeps its collapsed "
             "point. OpenVSP's own file is kept as <Case Name>.openvsp.stp."),
    F_choice("Geometry Set", "geom_set", ("All", "Thick", "Thin"), ("all", "thick", "thin"), 0,
             "Which components: all, or those in the Thick or Thin set (the component blocks' Set)."),
)


def cad_block():
    own = (InputSpec("OpenVSP Model", "text", "Output of OpenVSP Model (VSP), or the path of any .vsp3 model "
                     "(in metres)."),) + tuple(f.spec for f in CAD_FIELDS)
    specs = own + run_inputs("cad")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**OpenVSP to CAD.** The driver `run_cad.py` opens the model in OpenVSP and writes it with OpenVSP's "
              "own CAD export: trimmed and intersected (Surface Intersection), split and stitched (one BREP per "
              "component), or untrimmed surfaces.")
    lines = k.field_lines(CAD_FIELDS, lambda f: inp(f.spec.name)) + [k.line("model", inp("OpenVSP Model"))]
    base = analysis(g, k, inp, "run_cad", lambda b: [], lines, "CAD Export")
    out = k.cat([base, k.pick(inp("Format"), (".stp", ".igs"))], name="CAD File")
    return make_block(
        "openvsp_cad", "OpenVSP CAD Export (VSP)",
        "Writes an OpenVSP model (from **OpenVSP Model (VSP)**, so from nTop implicit bodies) as **STEP or IGES** "
        "with OpenVSP's CAD export and returns the file's path (`<Working Directory>/<Case Name>.stp` or `.igs`). "
        "Trimmed: the components intersected into one BREP solid.", specs, "text", g, out)


# --- Lofted CAD: B-spline solids fitted to the nTop slices, without OpenVSP's skinning

LOFT_FIELDS = (
    F_int("Points Around", "points_around", 48,
          "Control points round each body (a multiple of 4, at least 8). More follow tighter section shapes (a "
          "flattened nose, a canopy fillet) more closely; 48 kept the test bodies within 0.01 mm of their exact "
          "surfaces."),
    F_int("Faces Around", "faces_around", 4,
          "CAD faces round each body: 4 (default) splits it at the right, bottom, left and top; 8, 12, 16 ... cut "
          "each quarter evenly into more faces. A multiple of 4 that divides Points Around. The surface is the same "
          "either way (cut exactly, faces join smoothly); smaller faces can suit downstream CAD work."),
    F_real("Slices per Control Point", "rows_per_ctrl", 2.0, None,
           "nTop slices per control point along the body (at least 1). The surface is fitted through every slice; "
           "fewer control points (a larger number) give a lighter file and a smoother, less exact surface."),
    F_choice("Sharp Edges", "sharp_edges", ("Auto", "None"), ("auto", "none"), 0,
             "Auto: where the slices turn a corner at the body's sides, top or bottom (a chine: the contours turn by "
             "more than 30 degrees there along most of the body), the surface keeps the edge sharp. None: smooth all "
             "round. Edges elsewhere on the section are rounded over about one control-point spacing."),
)


def loft_block():
    own = (InputSpec("Components", "list<text>", "Outputs of OpenVSP Fuselage (VSP) blocks, in a List: their nTop "
                     "slices are what is fitted (their OpenVSP settings are not used). Slice Spacing 0.0025 (twice "
                     "as fine as the default) was used in the tests."),) + tuple(f.spec for f in LOFT_FIELDS)
    specs = own + run_inputs("loft")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Lofted CAD.** The driver `run_loft.py` fits a smooth bicubic B-spline surface through all of each "
              "body's nTop slices (least squares, no OpenVSP), cuts it into Faces Around faces round the body, "
              "closes a pointed end at a single point and a blunt end with a flat face, and writes one closed solid "
              "per body to `<Case Name>.stp` (millimetres).")
    lines = k.field_lines(LOFT_FIELDS, lambda f: inp(f.spec.name)) + [k.cat([k.join(inp("Components"), "\n"), "\n"])]
    base = analysis(g, k, inp, "run_loft", lambda b: [], lines, "Lofted CAD", needs_openvsp=False)
    out = k.cat([base, ".stp"], name="STEP File")
    return make_block(
        "lofted_cad", "Lofted CAD Export (VSP)",
        "Writes implicit bodies as **STEP solids fitted directly to their nTop slices**, with no OpenVSP skinning, "
        "and returns the file's path (`<Working Directory>/<Case Name>.stp`). Give it the **OpenVSP Fuselage (VSP)** "
        "blocks of the bodies. Each body becomes one closed solid: a smooth B-spline surface through every slice, "
        "cut into faces round the body (4 by default), a single point at a pointed end and a flat face at a blunt "
        "one.", specs, "text", g, out)
