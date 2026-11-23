"""OpenVSP model blocks: nTop implicit bodies to OpenVSP Fuselage and Wing components.

A component block slices its implicit bodies in nTop (Slice Body, exported as CLI contours; a list of bodies is
sliced body by body through nTop's list processing) and outputs one line of text describing them; OpenVSP Model
collects those lines, and its driver (``run_model.py``) makes one OpenVSP component per body, fits OpenVSP
sections to its contours and builds the ``.vsp3`` model the analysis blocks take.
"""
from __future__ import annotations

from ntopfile import Graph, InputSpec, Ref

from kit import (LENGTH, Kit, analysis, d_choice, d_int, d_real, d_text, make_block, model_inputs, prop, run_inputs)

SETS = ("Thick (panel method)", "Thin (vortex lattice)", "None (Wave and Parasite Drag only)")
SET_KEYS = ("thick", "thin", "none")
SYMMETRY = ("None", "x-z plane (left and right)", "x-y plane (top and bottom)", "x-z and x-y planes")
SYMMETRY_KEYS = ("none", "xz", "xy", "xzxy")
ACCOUNTING = ("Whole Wing", "Exposed Strip", "Exposed Strip + Carry-over (TR 1307)")
ACCOUNTING_KEYS = ("whole", "exposed", "carryover")


def _r(g, op, a, b):
    return g.block(f"{op}<real,real>", "real", [a, b], top=False)


def _sqrt(g, a):
    return g.block("sqrt<real>", "real", [a], top=False)


def _sq(g, a):
    return _r(g, "multiply", a, a)


def _sum(g, items):
    out = items[0]
    for x in items[1:]:
        out = _r(g, "add", out, x)
    return out


def _box(g, body):
    """Bounding box of the body: min and max coordinate lists, centre, size per axis and diagonal."""
    mn = [prop(prop(prop(body, "bounding box"), "min point"), c) for c in "xyz"]
    mx = [prop(prop(prop(body, "bounding box"), "max point"), c) for c in "xyz"]
    half = g.real(0.5)
    centre = [_r(g, "multiply", _r(g, "add", a, b), half) for a, b in zip(mn, mx)]
    size = [_r(g, "subtract", b, a) for a, b in zip(mn, mx)]
    diag = _sqrt(g, _sum(g, [_sq(g, s) for s in size]))
    return mn, mx, centre, size, diag


def _slice(g, k, inp, bodies, origin, ex, ey, height, tolerance, what):
    """Slice each body in its frame (origin, ex, ey) and export its contours to a temporary CLI file; returns the
    exported file paths (a list, one per body)."""
    frame = g.block("frame<point,vector,vector>", "frame", [origin, ex, ey], "Slicing Frames")
    stack = g.block("slice_implicit<implicit,real,real,bool,frame,implicit>[3.0.0]", "slice_stack",
                    [bodies, height, tolerance, g.boolean(False), frame, Ref(0)], "Sections")
    temp = g.block("temporary_directory<list_interface>[5.21.0]", "temporary_directory", [inp("Names")],
                   "Temporary Folder")
    index = g.block("sequence<real,real,integer>", "list<real>", [g.real(1.0), g.real(1.0), prop(bodies, "size")],
                    top=False)
    path = k.ecat([prop(temp, "path"), "/", k.join(inp("Names"), "_"), ".", k.num(index), "." + what + ".cli"])
    written = g.block("export_slices_to_cli<file_path,slice_stack,ascii_or_binary_enum>", "cli_slice_data",
                      [path, stack, g.enum("ascii_or_binary_enum", 0)], "Export Sections")
    return k.wait(written, path, name="Section Files")


def _per_body(k, values, ref=1.0, units=None):
    """A per-body list of reals as ``|``-separated text."""
    return k.join(k.num(values, ref, units), "|")


BODIES_DOC = ("A list of bodies: each becomes its own OpenVSP component, with the settings below (sliced one by one "
              "in nTop).")
NAMES_DOC = ("Component names in OpenVSP, one per body in order; one name for several bodies numbers them (Name_1, "
             "Name_2, ...). Names repeated across components get a number too.")


def _num(k, v, ref=1.0, units=None):
    return k.num(v, ref, units)


# --- Fuselage ---------------------------------------------------------------------------------------------------------

def fuselage_block():
    specs = (
        InputSpec("Bodies", "list<implicit>", "The bodies (implicits), each with its axis along +x: any shapes whose "
                  "cross-sections normal to x are single closed curves (fuselages, pods, closed nacelles, lifting "
                  "bodies). " + BODIES_DOC + " Union touching bodies in nTop first (one fuselage), rather than giving "
                  "OpenVSP overlapping ones."),
        InputSpec("Names", "list<text>", NAMES_DOC, default=lambda g: g.block("core.list<text>", "list<text>",
                                                                             [g.text("Fuselage")], top=False)),
        InputSpec("Set", "choice", "How VSPAERO sees the body. None (default): not at all; VSPAERO solves the lifting "
                  "surfaces, and the body still counts in Wave Drag and Parasite Drag. Thick: the panel method on the "
                  "closed surface, subsonic only (VSPAERO's panel method gave wrong drag above Mach 1 in every test, "
                  "on OpenVSP's own models too). Thin: a vortex lattice on its camber plates. See the README for the "
                  "tests behind the default.", default=d_choice(SETS, 2)),
        InputSpec("Slice Spacing", "real", "Distance between nTop's section cuts, as a fraction of the body's length.",
                  default=d_real(0.005)),
        InputSpec("Contour Tolerance", "real", "Point spacing of the section contours, as a fraction of the body's "
                  "bounding-box diagonal.", default=d_real(1e-4)),
        InputSpec("Sections", "integer", "OpenVSP cross-sections (including the nose and tail).", default=d_int(24)),
        InputSpec("Nose Clustering", "real", "Section spacing near the nose: 1 even, smaller is finer there "
                  "(OpenVSP's clustering function).", default=d_real(0.5)),
        InputSpec("Tail Clustering", "real", "Section spacing near the tail, as Nose Clustering.", default=d_real(0.5)),
        InputSpec("Tess W", "integer", "Tessellation points around each cross-section (OpenVSP Tess_W).",
                  default=d_int(33)),
        InputSpec("Panel Aspect Ratio", "real", "Panel length along the body over its width around it; sets the "
                  "points along the body (OpenVSP SectTess_U of each section). VSPAERO's panel method goes wrong when "
                  "body panels are much longer than wide.", default=d_real(1.5)),
        InputSpec("Nose Tess Clustering", "real", "OpenVSP FwdCluster of the first section (tessellation spacing at "
                  "the nose; 1 even).", default=d_real(1.0)),
        InputSpec("Tail Tess Clustering", "real", "OpenVSP AftCluster of the last section.", default=d_real(1.0)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    body = inp("Bodies")
    g.comment("**OpenVSP Fuselages from implicit bodies.** nTop cuts each body with planes normal to x (Slice Body) "
              "and exports the section contours. OpenVSP Model fits an OpenVSP Fuselage to them: cross-sections at "
              "clustered stations, each a closed curve through its contour; a pointed nose or tail becomes a point "
              "section, a blunt end a flat cap.")
    mn, mx, centre, size, diag = _box(g, body)
    g.comment("**Slicing.** The cuts start just ahead of the body's bounding box, every Slice Spacing x length; "
              "the contour points are Contour Tolerance x the box diagonal apart.")
    reach = _r(g, "multiply", diag, g.real(0.505))
    origin = g.block("point<real,real,real>", "point", [_r(g, "subtract", centre[0], reach), centre[1], centre[2]],
                     "Slicing Origin")
    height = _r(g, "multiply", size[0], inp("Slice Spacing"))
    tol = _r(g, "multiply", diag, inp("Contour Tolerance"))
    cli = _slice(g, k, inp, body, origin, g.vector(0, 1, 0), g.vector(0, 0, 1), height, tol, "fuselage")
    g.comment("**Component line** for OpenVSP Model: the names, the contour files and slicing origins (metres) of "
              "the bodies (`|`-separated lists, one item per body), and the OpenVSP settings.")
    parts = ["component = kind: fuselage; names: ", k.join(inp("Names"), "|"), "; files: ", k.join(cli, "|"),
             "; set: ", k.pick(inp("Set"), SET_KEYS)]
    for key, v in zip(("ox", "oy", "oz"), (prop(origin, "x"), prop(origin, "y"), prop(origin, "z"))):
        parts += [f"; {key}: ", _per_body(k, v, 1.0, LENGTH)]
    # the bounding box's x range places blunt ends between the outermost slice and the next (empty) one
    parts += ["; bx0: ", _per_body(k, mn[0], 1.0, LENGTH), "; bx1: ", _per_body(k, mx[0], 1.0, LENGTH)]
    for key, name in (("sections", "Sections"), ("nose_cluster", "Nose Clustering"), ("tail_cluster", "Tail Clustering"),
                      ("tess_w", "Tess W"), ("panel_aspect", "Panel Aspect Ratio"),
                      ("fwd_tess_cluster", "Nose Tess Clustering"), ("aft_tess_cluster", "Tail Tess Clustering")):
        parts += [f"; {key}: ", k.num(inp(name))]
    out = k.cat(parts, name="Component")
    return make_block(
        "openvsp_fuselage", "OpenVSP Fuselage (VSP)",
        "OpenVSP **Fuselages** fitted to a list of implicit bodies (one component each), as one component line for "
        "**OpenVSP Model (VSP)**. nTop slices each body normal to x; the fit places Sections cross-sections along it "
        "and gives each the shape of the body's contour there. Set puts them in OpenVSP's Thick or Thin set for "
        "VSPAERO.", specs, "text", g, out)


# --- Wing ---------------------------------------------------------------------------------------------------------------

def wing_block():
    specs = (
        InputSpec("Bodies", "list<implicit>", "The lifting surfaces (implicits): wings, tails, canards or fins, each "
                  "with airfoil sections in planes normal to its span axis and its chord along +x. " + BODIES_DOC),
        InputSpec("Names", "list<text>", NAMES_DOC, default=lambda g: g.block("core.list<text>", "list<text>",
                                                                             [g.text("Wing")], top=False)),
        InputSpec("Set", "choice", "Thin (default): VSPAERO's vortex lattice on the camber surface, which held from "
                  "subsonic to supersonic Mach in testing. Thick: the panel method on the closed surface (subsonic "
                  "only). None: left out of VSPAERO, still in Wave Drag and Parasite Drag.", default=d_choice(SETS, 1)),
        InputSpec("Span Axes", "list<vector>", "Spanwise direction of each body, in the order of Bodies (as many as "
                  "there are bodies), normal to x: +y for a wing or tailplane (root to tip), +z for a fin, in between "
                  "for a canted surface.", default=lambda g: g.block("core.list<vector>", "list<vector>",
                                                                       [g.vector(0, 1, 0)], top=False)),
        InputSpec("Symmetry", "choice", "OpenVSP's planar symmetry. x-z plane: each surface is mirrored left to right "
                  "(a wing, a tailplane); x-y plane: top to bottom (a dorsal fin mirrored into a ventral one); both: "
                  "four copies (cruciform or X tails). With symmetry, only the part of each body on the +span side of "
                  "the plane through the origin normal to its span axis is used, so give the whole surface or that "
                  "half. None: each whole body is one surface.", default=d_choice(SYMMETRY, 1)),
        InputSpec("Extension", "real", "How far to extend each wing inboard of its root (where the given body starts), "
                  "along its span axis: a new root section continues the first section's leading and trailing edges "
                  "and thickness, with the root airfoil. VSPAERO then solves panels inside the fuselage, which "
                  "removes the free root edge (and its tip-like loss of lift) that a wing ending at the body surface "
                  "has. 0: no extension. It must stop short of the plane through the origin normal to the span axis "
                  "(the body's axis for fins on a body centred on x): wings that meet or cross there, such as four "
                  "fins reaching the axis, made VSPAERO return no finite result. Stopping a little short, about 1 to 2% "
                  "of the root's distance, keeps the root vortices of opposite wings nearly cancelling. The OpenVSP "
                  "wing, and so Wave and Parasite Drag, include the extension.", LENGTH, d_real(0.0, LENGTH)),
        InputSpec("Force Accounting", "choice", "Which part of each wing's VSPAERO loads counts. Whole Wing: every "
                  "strip, the Extension included (the extension's lift then stands for the body's carry-over lift, "
                  "which it over-predicts on a slender body). Exposed Strip: only the strips outboard of the wing "
                  "body's own root (where the given implicit starts: give the wing trimmed at the fuselage surface); "
                  "the driver takes the inboard strips' forces and moments off every case, from VSPAERO's spanwise "
                  "loads. Exposed Strip + Carry-over (TR 1307): as Exposed Strip, plus the lift each wing carries "
                  "onto the body, K_B(W) times the wing's own normal force, from NACA TR 1307 (supersonic planar "
                  "model on the body's measured width and length behind the root, or slender-body theory where its "
                  "selection rule asks); the body must be an OpenVSP Fuselage in the model, in any Set.",
                  default=d_choice(ACCOUNTING, 0)),
        InputSpec("Slice Spacing", "real", "Distance between nTop's section cuts, as a fraction of the body's size along "
                  "Span Axis.", default=d_real(0.01)),
        InputSpec("Contour Tolerance", "real", "Point spacing of the section contours, as a fraction of the body's "
                  "bounding-box diagonal.", default=d_real(1e-4)),
        InputSpec("Sections", "integer", "Most OpenVSP wing sections: stations are added where the planform, twist or "
                  "thickness stop varying linearly, up to this many.", default=d_int(8)),
        InputSpec("Section Tolerance", "real", "Largest departure from linear variation between stations, as a "
                  "fraction of the root chord, before another section is added.", default=d_real(0.01)),
        InputSpec("Airfoil Points", "integer", "Points on each surface of each airfoil (cosine spaced).",
                  default=d_int(41)),
        InputSpec("Trailing Edge", "choice", "Sharp: a blunt trailing edge is closed (each surface sheared onto the "
                  "chord line), so VSPAERO's wake leaves one edge. As sliced: kept blunt.",
                  default=d_choice(("Sharp", "As sliced"), 0)),
        InputSpec("Tess U", "integer", "Spanwise tessellation points over the surface (shared between sections by "
                  "span).", default=d_int(20)),
        InputSpec("Tess W", "integer", "Chordwise tessellation points around each airfoil (OpenVSP Tess_W).",
                  default=d_int(33)),
        InputSpec("LE Clustering", "real", "OpenVSP LECluster: chordwise panel spacing at the leading edge (smaller is "
                  "finer).", default=d_real(0.25)),
        InputSpec("TE Clustering", "real", "OpenVSP TECluster at the trailing edge.", default=d_real(0.25)),
        InputSpec("Root Clustering", "real", "OpenVSP InCluster of the first section (spanwise spacing at the root).",
                  default=d_real(1.0)),
        InputSpec("Tip Clustering", "real", "OpenVSP OutCluster of the last section (spanwise spacing at the tip).",
                  default=d_real(1.0)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    body = inp("Bodies")
    g.comment("**OpenVSP Wings from implicit bodies.** nTop cuts each body with planes normal to its span axis and "
              "exports the airfoil contours. OpenVSP Model fits an OpenVSP Wing: at each station the chord line runs "
              "from the leading edge to the trailing edge, giving chord, twist and the airfoil, and stations become "
              "sections where those stop varying linearly.")
    mn, mx, centre, size, diag = _box(g, body)
    s = inp("Span Axes")
    sc = [prop(s, c) for c in "xyz"]
    norm = _sqrt(g, _sum(g, [_sq(g, c) for c in sc]))
    unit = [_r(g, "divide", c, norm) for c in sc]
    g.comment("**Slicing.** Cuts every Slice Spacing x the body's extent along the span axis, starting below the "
              "body (or, with Symmetry, at the plane through the origin normal to the span axis). The frame's x axis "
              "is x cross span and its y axis is x, so the contours are airfoils with x along the chord.")
    extent = _sum(g, [_r(g, "multiply", _sqrt(g, _sq(g, u)), d) for u, d in zip(unit, size)])
    reach = _r(g, "multiply", diag, g.real(0.505))
    along = _sum(g, [_r(g, "multiply", c, u) for c, u in zip(centre, unit)])
    # back = reach without symmetry, along (the plane through the origin) with it: a 0/1 weighting, since reach and
    # along are lists (one per body) and a list cannot be an item of a list.
    mirrored = g.block("core.select_by_choice<choice,list_interface>", "real",
                       [inp("Symmetry"), g.block("core.list<real>", "list<real>",
                                                 [g.real(0.0), g.real(1.0), g.real(1.0), g.real(1.0)], top=False)],
                       top=False)
    back = _r(g, "subtract", reach, _r(g, "multiply", mirrored, _r(g, "subtract", reach, along)))
    origin = g.block("point<real,real,real>", "point",
                     [_r(g, "subtract", c, _r(g, "multiply", back, u)) for c, u in zip(centre, unit)], "Slicing Origin")
    span_vec = g.block("vector<real,real,real>", "vector", unit, top=False)
    ex = g.block("cross_product<vector,vector>", "vector", [g.vector(1, 0, 0), span_vec], top=False)
    height = _r(g, "multiply", extent, inp("Slice Spacing"))
    tol = _r(g, "multiply", diag, inp("Contour Tolerance"))
    cli = _slice(g, k, inp, body, origin, ex, g.vector(1, 0, 0), height, tol, "wing")
    g.comment("**Component line** for OpenVSP Model: the names, contour files, slicing origins (metres) and span "
              "axes of the bodies (`|`-separated lists, one item per body), and the OpenVSP settings.")
    parts = ["component = kind: wing; names: ", k.join(inp("Names"), "|"), "; files: ", k.join(cli, "|"),
             "; set: ", k.pick(inp("Set"), SET_KEYS), "; symmetry: ", k.pick(inp("Symmetry"), SYMMETRY_KEYS),
             "; sharp_te: ", k.pick(inp("Trailing Edge"), ("1", "0")),
             "; extend: ", k.num(inp("Extension"), 1.0, LENGTH),
             "; accounting: ", k.pick(inp("Force Accounting"), ACCOUNTING_KEYS)]
    for key, v in zip(("ox", "oy", "oz"), (prop(origin, "x"), prop(origin, "y"), prop(origin, "z"))):
        parts += [f"; {key}: ", _per_body(k, v, 1.0, LENGTH)]
    for key, v in zip(("sx", "sy", "sz"), unit):
        parts += [f"; {key}: ", _per_body(k, v)]
    for key, name in (("sections", "Sections"), ("section_tolerance", "Section Tolerance"),
                      ("airfoil_points", "Airfoil Points"), ("tess_u", "Tess U"), ("tess_w", "Tess W"),
                      ("le_cluster", "LE Clustering"), ("te_cluster", "TE Clustering"),
                      ("root_cluster", "Root Clustering"), ("tip_cluster", "Tip Clustering")):
        parts += [f"; {key}: ", k.num(inp(name))]
    out = k.cat(parts, name="Component")
    return make_block(
        "openvsp_wing", "OpenVSP Wing (VSP)",
        "OpenVSP **Wings** fitted to a list of implicit lifting surfaces (one component each, with its own span "
        "axis), as one component line for **OpenVSP Model (VSP)**. nTop slices each body normal to its span axis; the "
        "fit finds each station's chord line, twist and airfoil and keeps the stations needed to follow the planform. "
        "Symmetry mirrors them about the x-z plane, the x-y plane or both. OpenVSP marks their trailing edges as wakes "
        "for VSPAERO, so no wake input is needed.", specs, "text", g, out)


# --- Duct ---------------------------------------------------------------------------------------------------------------

WALLS = ("Outer", "Inner")
WALL_KEYS = ("outer", "inner")
DUCT_INSET = 0.002  # end cuts this fraction of the way inside the body (a cut on an end face is degenerate)


def _vadd(g, a, b, s=None):
    """a + s b, component lists (s may be a list)."""
    return [_r(g, "add", x, _r(g, "multiply", y, s) if s is not None else y) for x, y in zip(a, b)]


def _dot(g, a, b):
    return _sum(g, [_r(g, "multiply", x, y) for x, y in zip(a, b)])


def _unit(g, a):
    n = _sqrt(g, _dot(g, a, a))
    return [_r(g, "divide", x, n) for x in a]


def _plane(g, plane):
    o = prop(plane, "origin")
    n = prop(plane, "normal")
    return [prop(o, c) for c in "xyz"], [prop(n, c) for c in "xyz"]


def _vector(g, comps):
    return g.block("vector<real,real,real>", "vector", comps, top=False)


def duct_block():
    specs = (
        InputSpec("Body", "implicit", "The duct as a closed body: a flow passage, or a tube with a wall. Its two ends "
                  "are planar and lie in Inlet Plane and Outlet Plane (which need not be parallel). The OpenVSP "
                  "surface is the body's side wall between them; the end faces are left out."),
        InputSpec("Name", "text", "Component name in OpenVSP.", default=d_text("Duct")),
        InputSpec("Inlet Plane", "plane", "The plane the body's inlet end lies in. The duct may run along any direction except y (the cuts' in-plane axes come from n x y)."),
        InputSpec("Outlet Plane", "plane", "The plane the body's outlet end lies in."),
        InputSpec("Wall", "choice", "Which wall becomes the surface: Outer (the body's outside, the only one of a "
                  "solid flow passage) or Inner (the hole of a tube: its inside wall).", default=d_choice(WALLS, 0)),
        InputSpec("Set", "choice", "How VSPAERO sees the duct (as OpenVSP Fuselage). Default None: only Wave Drag, "
                  "Parasite Drag and CAD export use it.", default=d_choice(SETS, 2)),
        InputSpec("Sections", "integer", "Cross-sections, the two ends included. nTop cuts the body with planes that "
                  "turn gradually from the inlet plane to the outlet plane; each cut becomes one OpenVSP section.",
                  default=d_int(24)),
        InputSpec("Inlet Clustering", "real", "Section spacing at the inlet: 1 even, smaller is finer there "
                  "(OpenVSP's clustering function).", default=d_real(1.0)),
        InputSpec("Outlet Clustering", "real", "Section spacing at the outlet, as Inlet Clustering.",
                  default=d_real(1.0)),
        InputSpec("Contour Tolerance", "real", "Point spacing of the section contours, as a fraction of the body's "
                  "bounding-box diagonal.", default=d_real(1e-4)),
        InputSpec("Tess W", "integer", "Tessellation points around each cross-section (OpenVSP Tess_W).",
                  default=d_int(33)),
        InputSpec("Panel Aspect Ratio", "real", "Panel length along the duct over its width around it; sets the "
                  "points between sections (OpenVSP SectTess_U).", default=d_real(1.5)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    body = inp("Body")
    g.comment("**OpenVSP duct from an implicit body.** The side wall of a body whose two planar ends lie in the inlet "
              "and outlet planes, as an OpenVSP Fuselage open at both ends (no caps). nTop cuts the body with "
              "Sections planes that turn gradually from the inlet plane to the outlet plane, so every section is a "
              "closed curve and the ends land exactly in their planes even when the planes are not parallel.")
    mn, mx, centre, size, diag = _box(g, body)
    o_in, n_in = _plane(g, inp("Inlet Plane"))
    o_out, n_out = _plane(g, inp("Outlet Plane"))
    g.comment("**End points and normals.** The centre of the body's bounding box projected onto each plane; each "
              "plane's normal turned to point from the inlet towards the outlet.")
    neg = g.real(-1.0)
    p_in = _vadd(g, centre, n_in, _r(g, "multiply", _dot(g, _vadd(g, centre, o_in, neg), n_in), neg))
    p_out = _vadd(g, centre, n_out, _r(g, "multiply", _dot(g, _vadd(g, centre, o_out, neg), n_out), neg))
    d = _vadd(g, p_out, p_in, neg)

    def towards(n):
        dn = _dot(g, d, n)
        sign = _r(g, "divide", dn, _sqrt(g, _sq(g, dn)))
        return [_r(g, "multiply", c, sign) for c in n]

    n_in, n_out = towards(n_in), towards(n_out)
    g.comment("**Cutting planes.** s from 0 to 1 by OpenVSP's clustering function (Inlet and Outlet Clustering), the "
              "end cuts a little inside the body; plane s passes through (1 - s) inlet point + s outlet point with "
              "the normal (1 - s) inlet normal + s outlet normal.")
    t = g.block("sequence_from_bounds<real,real,integer>", "list<real>", [g.real(0.0), g.real(1.0), inp("Sections")],
                top=False)
    mt = _r(g, "subtract", g.real(1.0), t)
    s = _sum(g, [_r(g, "multiply", _r(g, "multiply", _sq(g, mt), t), inp("Inlet Clustering")),
                 _r(g, "multiply", _r(g, "multiply", mt, _sq(g, t)), _r(g, "subtract", g.real(3.0),
                                                                           inp("Outlet Clustering"))),
                 _r(g, "multiply", _sq(g, t), t)])
    s = _r(g, "add", g.real(DUCT_INSET), _r(g, "multiply", s, g.real(1.0 - 2.0 * DUCT_INSET)))
    origin_c = _vadd(g, p_in, d, s)
    normal_c = _unit(g, _vadd(g, n_in, _vadd(g, n_out, n_in, neg), s))
    # in-plane axes of every cut: y axis = n x (world y), x axis = y axis x n. nTop writes a cut on a plane that is
    # not normal to a world axis in exactly these axes, whatever frame it is given (found on this block's test duct:
    # the axes of the frame passed in were ignored), so the frame is built to match them. A duct along y has none.
    nvec = _vector(g, normal_c)
    ey = g.block("cross_product<vector,vector>", "vector", [nvec, g.vector(0, 1, 0)], top=False)
    ey_c = _unit(g, [prop(ey, c) for c in "xyz"])
    eyv = _vector(g, ey_c)
    ex = g.block("cross_product<vector,vector>", "vector", [eyv, nvec], top=False)
    ex_c = [prop(ex, c) for c in "xyz"]
    exv = _vector(g, ex_c)
    origin = g.block("point<real,real,real>", "point", origin_c, "Cut Origins")
    frame = g.block("frame<point,vector,vector>", "frame", [origin, exv, eyv], "Cutting Frames")
    g.comment("**Cuts.** One Slice Body per plane, with a slice spacing larger than the body so that each holds only "
              "its cut through the frame origin, exported to a temporary CLI file each.")
    stack = g.block("slice_implicit<implicit,real,real,bool,frame,implicit>[3.0.0]", "slice_stack",
                    [body, _r(g, "multiply", diag, g.real(2.0)), _r(g, "multiply", diag, inp("Contour Tolerance")),
                     g.boolean(False), frame, Ref(0)], "Cuts")
    temp = g.block("temporary_directory<list_interface>[5.21.0]", "temporary_directory", [inp("Name")],
                   "Temporary Folder")
    index = g.block("sequence<real,real,integer>", "list<real>", [g.real(1.0), g.real(1.0), inp("Sections")],
                    top=False)
    path = k.ecat([prop(temp, "path"), "/", inp("Name"), ".", k.num(index), ".duct.cli"])
    written = g.block("export_slices_to_cli<file_path,slice_stack,ascii_or_binary_enum>", "cli_slice_data",
                      [path, stack, g.enum("ascii_or_binary_enum", 0)], "Export Cuts")
    cli = k.wait(written, path, name="Cut Files")
    g.comment("**Component line** for OpenVSP Model: the cut files and their frames (`|`-separated, one per section, "
              "metres), the end planes, and the OpenVSP settings.")
    parts = ["component = kind: duct; name: ", inp("Name"), "; cuts: ", k.join(cli, "|"),
             "; wall: ", k.pick(inp("Wall"), WALL_KEYS), "; set: ", k.pick(inp("Set"), SET_KEYS)]
    for key, comps, units in (("p", origin_c, LENGTH), ("u", ex_c, None), ("v", ey_c, None)):
        for c, v in zip("xyz", comps):
            parts += [f"; {key}{c}: ", _per_body(k, v, 1.0, units)]
    for key, comps, units in (("inlet_p", p_in, LENGTH), ("inlet_n", n_in, None), ("outlet_p", p_out, LENGTH),
                              ("outlet_n", n_out, None)):
        parts += [f"; {key}: ", k.cat([_num(k, comps[0], 1.0, units), ",", _num(k, comps[1], 1.0, units), ",",
                                       _num(k, comps[2], 1.0, units)])]
    for key, name in (("tess_w", "Tess W"), ("panel_aspect", "Panel Aspect Ratio")):
        parts += [f"; {key}: ", k.num(inp(name))]
    out = k.cat(parts, name="Component")
    return make_block(
        "openvsp_duct", "OpenVSP Duct (VSP)",
        "An OpenVSP **duct**: the side wall of an implicit body whose ends are planar and lie in the given inlet and "
        "outlet planes (not necessarily parallel), as an OpenVSP Fuselage open at both ends, for **OpenVSP Model "
        "(VSP)**. nTop cuts the body on planes turning from the inlet plane to the outlet plane; each cut's outer "
        "contour (or, for a tube, its inner one) becomes a section in that plane.", specs, "text", g, out)


# --- Model ----------------------------------------------------------------------------------------------------------------

def model_block():
    own = (InputSpec("Components", "list<text>", "Outputs of OpenVSP Fuselage, OpenVSP Wing and OpenVSP Duct blocks, in a List.",
                     default=lambda g: g.block("core.list<text>", "list<text>", [g.text("")], top=False)),)
    specs = own + run_inputs("model")
    g = Graph()
    k = Kit(g)
    model = model_inputs(specs)
    inp = model
    g.comment("**OpenVSP model from nTop bodies.** The driver `run_model.py` reads each component's section contours, "
              "fits OpenVSP sections to them, builds the components in OpenVSP (in its Thick and Thin sets), saves "
              "`<Case Name>.vsp3`, exports OpenVSP's tessellation (`<Case Name>.stl`) and compares the OpenVSP "
              "surfaces with the nTop contours (`<Case Name>.ntvsp.model.csv`).")
    lines = [k.cat([k.join(inp("Components"), "\n"), "\n"])]
    out = analysis(g, k, inp, "run_model", lambda base: [], lines, "OpenVSP Model")
    return make_block(
        "openvsp_model", "OpenVSP Model (VSP)",
        "Builds an OpenVSP model (`.vsp3`) from OpenVSP Fuselage and OpenVSP Wing component lines and returns its path "
        "base for **VSPAERO Analysis**, **Wave Drag** and **Parasite Drag** (VSP). Each component is fitted to its nTop "
        "body's section contours; **OpenVSP Model Check (VSP)** reports how closely, and **OpenVSP Model Mesh (VSP)** "
        "shows the result in nTop.", specs, "text", g, out)


MODEL_COLUMNS = ("Max_Deviation", "Mean_Deviation", "Max_Deviation_Fraction", "Sections", "Slices", "Size",
                 "Check_Points")


def model_check_block():
    specs = (InputSpec("Model", "text", "Output of OpenVSP Model (VSP)."),
             InputSpec("Column", "choice", "Per component, in the order given: largest and mean distance (metres) from "
                       "the OpenVSP surface to the nTop section contours, the largest as a fraction of the component's "
                       "length (fuselage) or root chord (wing), sections, slices used, that size (metres) and the number "
                       "of surface points checked.", default=d_choice(MODEL_COLUMNS)))
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Fit check.** One row per component in `<Case Name>.ntvsp.model.csv`; the chosen column comes out as a "
              "list.")
    table = g.block("import_table<file_path,bool>", "table", [k.cat([inp("Model"), ".ntvsp.model.csv"]), g.boolean(True)],
                    "Model Table")
    col = g.block("extract_column<table,text>", "list<text>", [table, k.pick(inp("Column"), MODEL_COLUMNS)], top=False)
    out = g.block("scalar_from_text<text>", "real", [col], "Values")
    return make_block("openvsp_model_check", "OpenVSP Model Check (VSP)",
                      "How closely each OpenVSP component follows its nTop body (distances in metres), one value per "
                      "component.", specs, "list<real>", g, out)


def model_mesh_block():
    specs = (InputSpec("Model", "text", "Output of OpenVSP Model (VSP)."),)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**OpenVSP tessellation** of the model (`<Case Name>.stl`, metres), to compare with the nTop bodies.")
    out = g.block("import_mesh<file_path,unit_length_enum>[1.1.0]", "mesh", [k.cat([inp("Model"), ".stl"]),
                                                                            g.unit_length("m")], "Model Mesh")
    return make_block("openvsp_model_mesh", "OpenVSP Model Mesh (VSP)",
                      "The OpenVSP model's tessellation as a mesh, to view against the nTop bodies.", specs, "mesh", g, out)
