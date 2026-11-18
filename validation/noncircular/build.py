"""Run the nTop_vsp blocks on the trapezoid-section rocket at the SU2 conditions, in one nTop notebook.

    sh build/py.sh validation/noncircular/build.py

OpenVSP Fuselage (Set None) and four OpenVSP Wing fins (Thin, Extend Into Body, Exposed Strip + Carry-over (TR 1307),
each fin its own body and span axis: the roots differ, 25 mm at the sides, 24 mm top and bottom) make the OpenVSP
model. On it: VSPAERO (fixed wake), Body Lift (slender body + crossflow), RASAero II (an axisymmetric code: its driver
fits a round body to the model) and Modified Newtonian on Mesh on the body alone. Flight condition: Mach 2, 3, 5 at sea
level, alpha -4, 0, 4, 8, beta 0; Sref = pi 54^2 / 4 mm^2, reference lengths 54 mm, moments about x = 0.7 m (as the
round rocket's study). Outputs go to runs/noncircular.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))

import geometry as geo  # noqa: E402
import rocket as rk  # noqa: E402
import build_models as bm  # noqa: E402
import aero_blocks as a  # noqa: E402
import drag_blocks as d  # noqa: E402
import model_blocks as mb  # noqa: E402
import vspaero_blocks as v  # noqa: E402
from ntopfile import Graph  # noqa: E402

WORK = os.path.join(rk.ROOT, "runs", "noncircular").replace("\\", "/")
MACHS = (2.0, 3.0, 5.0)
ALPHAS = (-4.0, 0.0, 4.0, 8.0)
BETAS = (0.0,)
SPAN_AXES = ((0, 1, 0), (0, 0, 1), (0, -1, 0), (0, 0, -1))
ACCOUNTING = os.environ.get("NC_ACCOUNTING", "carryover")


def build(work=WORK):
    fus, wing, model = mb.fuselage_block(), mb.wing_block(), mb.model_block()
    flight, solver, ana = a.flight_block(), v.solver_block(), v.analysis_block()
    body_b, ras_b, newt = d.body_lift_block(), d.rasaero_block(), a.newtonian_mesh_block()
    blocks = [fus, wing, model, flight, solver, ana, body_b, ras_b, newt]
    g = Graph()
    body = geo.body(g)
    fins = geo.fins(g)
    run = {"Working Directory": g.file_path(work)}
    lst = bm.lst
    comps = lst(g, "text", [
        rk.call(g, fus, {"Bodies": lst(g, "implicit", [body]), "Names": lst(g, "text", [g.text("Body")]),
                         "Set": g.choice(list(mb.SETS), 2)}, "Body Component"),
        rk.call(g, wing, {"Bodies": lst(g, "implicit", fins), "Names": lst(g, "text", [g.text("Fin")]),
                          "Set": g.choice(list(mb.SETS), 1),
                          "Span Axes": lst(g, "vector", [g.vector(*s) for s in SPAN_AXES]),
                          "Symmetry": g.choice(list(mb.SYMMETRY), 0), "Extend Into Body": g.boolean(True),
                          "Force Accounting": g.choice(list(mb.ACCOUNTING), mb.ACCOUNTING_KEYS.index(ACCOUNTING))},
                "Fin Components"),
    ], "Components")
    vsp_model = rk.call(g, model, {"Components": comps, **run, "Case Name": g.text("model")}, "Model")
    fc = rk.call(g, flight, bm.flight_inputs(g, work, "flight", MACHS, ALPHAS, BETAS), "Flight Condition")
    ss = rk.call(g, solver, {"Processors": g.integer(16), "Wake": g.choice(["Relaxed", "Fixed"], 1)}, "Solver Settings")
    outs = [vsp_model,
            rk.call(g, ana, {"OpenVSP Model": vsp_model, "Flight Condition": fc, "Solver Settings": ss, **run,
                             "Case Name": g.text("vspaero")}, "VSPAERO"),
            rk.call(g, body_b, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                                "Case Name": g.text("body_lift")}, "Body Lift"),
            rk.call(g, ras_b, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                               "Case Name": g.text("rasaero")}, "RASAero II"),
            rk.call(g, newt, {"Body": body, "Flight Condition": fc, "Mesh Tolerance": g.real(1e-3, rk.LENGTH), **run,
                              "Case Name": g.text("newtonian_body")}, "Newtonian Body")]
    return g, blocks, g.block("core.group<any>", "group", outs, "Results")


def build_plate(work=WORK):
    """VSPAERO with the body as a horizontal vortex-lattice plate: the body cut by a 2 mm slab in the fins' plane
    (z = 0), an OpenVSP Wing along +y with x-z symmetry, so it is as wide as the body at z = 0 and stops short of the
    side fins' roots (0.7 mm). The fins run Thin from the body surface (not extended: they would overlap the plate),
    with Exposed Strip + Carry-over so the carry-over term can be added or left off afterwards."""
    fus, wing, model = mb.fuselage_block(), mb.wing_block(), mb.model_block()
    flight, solver, ana = a.flight_block(), v.solver_block(), v.analysis_block()
    blocks = [fus, wing, model, flight, solver, ana]
    g = Graph()
    body = geo.body(g)
    slab = g.block("box<point,real,real,real>", "box", [g.point(600, 0, 0), g.real(1.3, rk.LENGTH),
                                                       g.real(0.2, rk.LENGTH), g.real(0.002, rk.LENGTH)], "Slab",
                   top=False)
    plate = rk.boolean(g, "intersect", [body, slab], "Body Plate", top=False)
    run = {"Working Directory": g.file_path(work)}
    lst = bm.lst
    comps = lst(g, "text", [
        rk.call(g, fus, {"Bodies": lst(g, "implicit", [body]), "Names": lst(g, "text", [g.text("Body")]),
                         "Set": g.choice(list(mb.SETS), 2)}, "Body Component"),
        rk.call(g, wing, {"Bodies": lst(g, "implicit", [plate]), "Names": lst(g, "text", [g.text("Plate")]),
                          "Set": g.choice(list(mb.SETS), 1), "Span Axes": lst(g, "vector", [g.vector(0, 1, 0)]),
                          "Symmetry": g.choice(list(mb.SYMMETRY), 1), "Slice Spacing": g.real(0.02)},
                "Plate Component"),
        rk.call(g, wing, {"Bodies": lst(g, "implicit", geo.fins(g)), "Names": lst(g, "text", [g.text("Fin")]),
                          "Set": g.choice(list(mb.SETS), 1),
                          "Span Axes": lst(g, "vector", [g.vector(*s) for s in SPAN_AXES]),
                          "Symmetry": g.choice(list(mb.SYMMETRY), 0), "Extend Into Body": g.boolean(False),
                          "Force Accounting": g.choice(list(mb.ACCOUNTING), 2)}, "Fin Components"),
    ], "Components")
    vsp_model = rk.call(g, model, {"Components": comps, **run, "Case Name": g.text("plate_model")}, "Model")
    fc = rk.call(g, flight, bm.flight_inputs(g, work, "flight", MACHS, ALPHAS, BETAS), "Flight Condition")
    ss = rk.call(g, solver, {"Processors": g.integer(16), "Wake": g.choice(["Relaxed", "Fixed"], 1)}, "Solver Settings")
    outs = [vsp_model, rk.call(g, ana, {"OpenVSP Model": vsp_model, "Flight Condition": fc, "Solver Settings": ss, **run,
                                        "Case Name": g.text("vspaero_plate")}, "VSPAERO")]
    return g, blocks, g.block("core.group<any>", "group", outs, "Results")


def main():
    os.makedirs(WORK, exist_ok=True)
    plate = "plate" in sys.argv[1:]
    g, blocks, out = build_plate() if plate else build()
    res, log, dt = rk.run(g, blocks, out, "group", os.path.join(WORK, "plate.ntop" if plate else "noncircular.ntop"))
    print(f"{dt:.0f} s, output {'yes' if res else 'no'}", flush=True)
    for line in log.splitlines():
        if "[E]" in line or "failed" in line.lower():
            print("   ", line[:300])


if __name__ == "__main__":
    main()
