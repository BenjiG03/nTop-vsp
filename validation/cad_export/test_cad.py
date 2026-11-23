"""nTop implicit bodies -> OpenVSP model -> CAD (STEP/IGES), checked against the implicit bodies.

    sh build/py.sh validation/cad_export/test_cad.py [build] [check]

Bodies (mm): the hobby rocket's body (ogive nose, flat base) as an OpenVSP Fuselage, its +y fin (flat plate with
bevelled leading and trailing edges) as an OpenVSP Wing, and a duct: a waisted tube (inner radius 60 mm narrowing to
51 mm mid-length, 6 mm wall) cut by an inlet plane turned 15 degrees about y and an outlet plane turned 25 degrees
about z, as two OpenVSP Ducts (its outer and its inner wall).

build: one notebook makes the components, the OpenVSP model, its CAD exports (trimmed STEP, stitched STEP,
untrimmed STEP, trimmed IGES) and nTop meshes of the bodies (STL).
check: each CAD file is read back with OpenCASCADE (gmsh): solids, shells, faces, area, volume. Deviation both ways:
points meshed on the CAD surfaces are evaluated in nTop's own field of the body they belong to (CAD -> implicit),
and nTop mesh vertices are projected onto the CAD surfaces (implicit -> CAD). Runs: runs/cad_export.
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))
import rocket as rk  # noqa: E402
import cad_blocks as cb  # noqa: E402
import model_blocks as mb  # noqa: E402
from ntopfile import Graph  # noqa: E402

WORK = os.path.join(rk.ROOT, "runs", "cad_export").replace("\\", "/")
LENGTH, ANGLE = rk.LENGTH, rk.ANGLE
# duct (mm)
DZ = 300.0                     # axis at y = 0, z = DZ (clear of the rocket)
RI, T, LD = 60.0, 6.0, 300.0   # inner radius at the ends, wall, length along the axis between the end points
A_IN, A_OUT = 15.0, 25.0       # inlet plane turned about y, outlet plane about z (degrees)
P_IN, P_OUT = (0.0, 0.0, DZ), (LD, 0.0, DZ)
N_IN = (math.cos(math.radians(A_IN)), 0.0, -math.sin(math.radians(A_IN)))
N_OUT = (math.cos(math.radians(A_OUT)), math.sin(math.radians(A_OUT)), 0.0)
MESH_TOL = {"Body": 2.0, "Fin": 0.5, "Duct": 1.0}  # nTop mesh tolerance (mm) for the implicit -> CAD check
# (method, format, planar end caps)
CASES = (("trimmed", "step", False), ("stitched", "step", False), ("surfaces", "step", False),
         ("trimmed", "iges", False), ("trimmed", "step", True), ("stitched", "step", True))


def case_name(method, fmt, caps):
    return f"parts_{method}_{fmt}" + ("_caps" if caps else "")


def r_inner(x):
    return RI * (1.0 - 0.15 * math.sin(math.pi * min(max(x, 0.0), LD) / LD) ** 2)


def duct_body(g):
    """The tube, longer than the duct, cut by the two half spaces (big boxes turned about the end points)."""
    xs = [-120.0 + (LD + 240.0) * i / 160 for i in range(161)]
    inner = [(x, 0.0, DZ + r_inner(x)) for x in xs]
    outer = [(x, 0.0, DZ + r_inner(x) + T) for x in reversed(xs)]
    prof = rk._profile(g, inner + outer, "Tube Wall Profile")
    axis = g.block("axis<point,vector>", "axis", [g.point(0, 0, DZ), g.vector(1, 0, 0)], "Tube Axis", top=False)
    tube = g.block("revolve<new_profile,axis,real>[5.20.0]", "implicit", [prof, axis, g.real(rk.FULL_TURN, ANGLE)],
                   "Tube")
    big = 1000.0
    keep_in = g.block("box_from_corners<point,point>", "implicit",
                      [g.point(0, -big, DZ - big), g.point(2 * big, big, DZ + big)], "Inlet Half Space", top=False)
    keep_in = g.block("rotate<spatial3d,point,vector,real>[1.1.0]", "implicit",
                      [keep_in, g.point(*P_IN), g.vector(0, 1, 0), g.real(math.radians(A_IN), ANGLE)],
                      "Inlet Half Space Turned", top=False)
    keep_out = g.block("box_from_corners<point,point>", "implicit",
                       [g.point(LD - 2 * big, -big, DZ - big), g.point(LD, big, DZ + big)], "Outlet Half Space",
                       top=False)
    keep_out = g.block("rotate<spatial3d,point,vector,real>[1.1.0]", "implicit",
                       [keep_out, g.point(*P_OUT), g.vector(0, 0, 1), g.real(math.radians(A_OUT), ANGLE)],
                       "Outlet Half Space Turned", top=False)
    return rk.boolean(g, "intersect", [tube, keep_in, keep_out], "Duct Body")


def plane(g, p, n, name):
    return g.block("plane_from_normal<point,vector>[1.1.0]", "plane", [g.point(*p), g.vector(*n)], name)


def lst(g, typ, items):
    return g.block(f"core.list<{typ}>", f"list<{typ}>", items, top=False)


def bodies(g):
    return {"Body": rk.body(g), "Fin": rk.fin(g), "Duct": duct_body(g)}


def build():
    os.makedirs(WORK, exist_ok=True)
    fus, wing, duct, model, cad = (mb.fuselage_block(), mb.wing_block(), mb.duct_block(), mb.model_block(),
                                   cb.cad_block())
    g = Graph()
    b = bodies(g)
    comps = [
        rk.call(g, fus, {"Bodies": lst(g, "implicit", [b["Body"]]), "Names": lst(g, "text", [g.text("Body")])},
                "Body Component"),
        rk.call(g, wing, {"Bodies": lst(g, "implicit", [b["Fin"]]), "Names": lst(g, "text", [g.text("Fin")]),
                          "Span Axes": lst(g, "vector", [g.vector(0, 1, 0)]), "Symmetry": g.choice(list(mb.SYMMETRY), 0)},
                "Fin Component"),
    ]
    pin, pout = plane(g, P_IN, N_IN, "Inlet Plane"), plane(g, P_OUT, N_OUT, "Outlet Plane")
    for wall in (0, 1):
        name = ("DuctOuter", "DuctInner")[wall]
        comps.append(rk.call(g, duct, {"Body": b["Duct"], "Name": g.text(name), "Inlet Plane": pin,
                                       "Outlet Plane": pout, "Wall": g.choice(list(mb.WALLS), wall)},
                             name + " Component"))
    m = rk.call(g, model, {"Components": lst(g, "text", comps), "Working Directory": g.file_path(WORK),
                           "Case Name": g.text("parts")}, "Model")
    outs = []
    for method, fmt, caps in CASES:
        outs.append(rk.call(g, cad, {"OpenVSP Model": m, "Working Directory": g.file_path(WORK),
                                     "Case Name": g.text(case_name(method, fmt, caps)),
                                     "Method": g.choice(list(cb.METHODS), ("trimmed", "stitched", "surfaces").index(method)),
                                     "Format": g.choice(list(cb.FORMATS), ("step", "iges").index(fmt)),
                                     "Planar End Caps": g.choice(["Off", "On"], 1 if caps else 0)},
                            f"CAD {case_name(method, fmt, caps)}"))
    for name, body in b.items():
        mesh = rk.mesh(g, body, MESH_TOL[name], name + " Mesh")
        outs.append(g.block("export_mesh<file_path,mesh,unit_length_enum>", "mesh_file_data",
                            [g.file_path(f"{WORK}/{name}.ntop_mesh.stl"), mesh, g.unit_length("mm")], name + " STL"))
    out = g.block("core.group<any>", "group", outs, "Outputs")
    res, log, dt = rk.run(g, [fus, wing, duct, model, cad], out, "group", os.path.join(WORK, "parts.ntop"))
    print(f"build: {dt:.0f} s, output {'yes' if res else 'no'}")
    for line in log.splitlines():
        if "[E]" in line or "rror" in line:
            print("  ", line[:300])
    if os.path.exists(os.path.join(WORK, "parts.ntvsp.log")):
        print(open(os.path.join(WORK, "parts.ntvsp.log")).read()[-3000:])


if __name__ == "__main__":
    steps = sys.argv[1:] or ["build", "check"]
    if "build" in steps:
        build()
    if "check" in steps:
        import check_cad
        check_cad.main()
