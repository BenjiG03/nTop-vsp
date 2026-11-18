"""OpenVSP models of the hobby rocket from its nTop implicits, with a small VSPAERO, Wave Drag and Parasite Drag check run
through the blocks (their settings files are the templates of the full sweeps in ``sweep.py``).

    sh build/py.sh validation/rasaero/build_models.py [variant ...]

Variants (fins: OpenVSP Wing blocks, one per fin, span axes +y, +z, -y, -z; body: OpenVSP Fuselage):
  thin_exposed      fins Thin from the body surface, body in no VSPAERO set
  thin_axis         fins Thin from the axis (through the body: carry-over lift), body in no set
  thin_exposed_body fins Thin from the body surface, body Thick
  thick_exposed     fins Thick from the body surface, body in no set
  thick_body        fins Thick, body Thick
"""
from __future__ import annotations

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import rocket as rk  # noqa: E402
import aero_blocks as a  # noqa: E402
import drag_blocks as d  # noqa: E402
import model_blocks as mb  # noqa: E402
import vspaero_blocks as v  # noqa: E402
from ntopfile import Graph  # noqa: E402

DEG = math.radians(1.0)
VARIANTS = {  # fin set, fin inner fraction of the radius, body set (, extend into body, force accounting)
    "thin_extend_exposed": ("thin", 1.0, "none", True, 1),
    "thin_exposed": ("thin", 1.0, "none"),
    "thin_axis": ("thin", 0.0, "none"),
    "thin_exposed_body": ("thin", 1.0, "thick"),
    "thick_exposed": ("thick", 1.0, "none"),
    "thick_body": ("thick", 1.0, "thick"),
}
CHECK_MACHS = (2.0, 4.0)
CHECK_ALPHAS = (0.0, 8.0)
CHECK_BETAS = (0.0, 4.0)
X_REF_MM = 700.0  # moment reference point on the axis, about 58% of the length (a typical CG)


def lst(g, typ, items, name=None):
    return g.block(f"core.list<{typ}>", f"list<{typ}>", items, name, top=name is not None)


def flight_inputs(g, work, case, machs, alphas, betas):
    r = rk.R
    return {"Mach": lst(g, "real", [g.real(m) for m in machs]),
            "Altitude": lst(g, "real", [g.real(0.0, rk.LENGTH)]),
            "Angle of Attack": lst(g, "real", [g.real(x * DEG, rk.ANGLE) for x in alphas]),
            "Sideslip": lst(g, "real", [g.real(x * DEG, rk.ANGLE) for x in betas]),
            "Reference Area": g.real(r.reference_area_mm2 * 1e-6, rk.AREA),
            "Reference Span": g.real(r.diameter * 1e-3, rk.LENGTH),
            "Reference Chord": g.real(r.diameter * 1e-3, rk.LENGTH),
            "Moment Reference Point": g.point(X_REF_MM, 0, 0),
            "Working Directory": g.file_path(work), "Case Name": g.text(case)}


def build(variant, work):
    fin_set, inner, body_set, *wing_opts = VARIANTS[variant]
    extend, accounting = (wing_opts + [False, 0])[:2]
    fus, wing, model = mb.fuselage_block(), mb.wing_block(), mb.model_block()
    flight, solver, ana, wave, parasite = (a.flight_block(), v.solver_block(), v.analysis_block(), d.wave_block(),
                                           d.parasite_block())
    body_b, ras_b = d.body_lift_block(), d.rasaero_block()
    blocks = [fus, wing, model, flight, solver, ana, wave, parasite, body_b, ras_b]
    g = Graph()
    body = rk.body(g)
    fins = rk.fins(g, inner_fraction=inner)
    run = {"Working Directory": g.file_path(work)}
    sets = list(mb.SETS)
    key = {"thick": 0, "thin": 1, "none": 2}
    comps = lst(g, "text", [
        rk.call(g, fus, {"Bodies": lst(g, "implicit", [body]), "Names": lst(g, "text", [g.text("Body")]),
                         "Set": g.choice(sets, key[body_set])}, "Body Component"),
        rk.call(g, wing, {"Bodies": lst(g, "implicit", fins), "Names": lst(g, "text", [g.text("Fin")]),
                          "Set": g.choice(sets, key[fin_set]),
                          "Span Axes": lst(g, "vector", [g.vector(0, 1, 0), g.vector(0, 0, 1), g.vector(0, -1, 0),
                                                          g.vector(0, 0, -1)]),
                          "Symmetry": g.choice(list(mb.SYMMETRY), 0), "Extend Into Body": g.boolean(extend),
                          "Force Accounting": g.choice(list(mb.ACCOUNTING), accounting)}, "Fin Components"),
    ], "Components")
    vsp_model = rk.call(g, model, {"Components": comps, **run, "Case Name": g.text(f"{variant}_model")}, "Model")
    fc = rk.call(g, flight, flight_inputs(g, work, f"{variant}_check_flight", CHECK_MACHS, CHECK_ALPHAS, CHECK_BETAS),
                 "Flight Condition")
    ss = rk.call(g, solver, {"Processors": g.integer(16), "Wake": g.choice(["Relaxed", "Fixed"], 1)}, "Solver Settings")
    outs = [vsp_model]
    if variant in ("thin_exposed", "thin_extend_exposed"):  # the check runs: templates of the sweeps
        outs.append(rk.call(g, ana, {"OpenVSP Model": vsp_model, "Flight Condition": fc, "Solver Settings": ss, **run,
                                     "Case Name": g.text(f"{variant}_check_vspaero")}, "VSPAERO"))
    if variant == "thin_extend_exposed":
        outs.append(rk.call(g, body_b, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                                      "Case Name": g.text(f"{variant}_check_body_lift")}, "Body Lift"))
        outs.append(rk.call(g, ras_b, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                                     "Case Name": g.text(f"{variant}_check_rasaero")}, "RASAero II"))
    if variant == "thin_exposed":
        outs.append(rk.call(g, wave, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                                      "Case Name": g.text(f"{variant}_check_wave")}, "Wave Drag"))
        outs.append(rk.call(g, parasite, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                                          "Case Name": g.text(f"{variant}_check_parasite")}, "Parasite Drag"))
    out = g.block("core.group<any>", "group", outs, "Results")
    return g, blocks, out


def main():
    names = sys.argv[1:] or list(VARIANTS)
    work = os.path.join(rk.RUNS, "models").replace("\\", "/")
    os.makedirs(work, exist_ok=True)
    for name in names:
        g, blocks, out = build(name, work)
        res, log, dt = rk.run(g, blocks, out, "group", os.path.join(work, f"{name}.ntop"))
        errs = [l for l in log.splitlines() if "[E]" in l or "failed" in l.lower()]
        print(f"{name}: {dt:.0f} s, output {'yes' if res else 'no'}", flush=True)
        for e in errs[:10]:
            print("   ", e[:300])


if __name__ == "__main__":
    main()
