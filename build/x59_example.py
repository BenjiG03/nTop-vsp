"""Write (and optionally run) the X-59 aero database example notebook.

    sh build/py.sh build/x59_example.py [--run]

The notebook is the X-59 hackathon starter's parametric geometry with blocks appended: each component body becomes
an OpenVSP component, OpenVSP Model builds the model, VSPAERO, Wave Drag and Parasite Drag run on it over a shared
flight envelope, and Aero Database combines them into a six-degree-of-freedom table. Wing Span and Wing Inner Sweep
are notebook inputs and are recorded in every row, so ``examples/x59_doe.py`` can run the notebook over a design of
experiments and collect one table (``runs/x59_doe.csv``).
"""
from __future__ import annotations

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import aero_blocks as a  # noqa: E402
import buildup_blocks as bb  # noqa: E402
import drag_blocks as d  # noqa: E402
import model_blocks as mb  # noqa: E402
import vspaero_blocks as v  # noqa: E402
import x59  # noqa: E402
from kit import ANGLE, AREA, LENGTH, ROOT  # noqa: E402
from ntopfile import InputSpec  # noqa: E402

EXAMPLES = ROOT + "/examples"
RUNS = ROOT + "/runs"
NOTEBOOK = EXAMPLES + "/X-59 Aero Database.ntop"
DOE_TABLE = RUNS + "/x59_doe.csv"

# Flight condition shared by every analysis
MACHS = (0.6, 0.9, 1.2, 1.5, 1.8)
ALTITUDE_FT = 55000.0         # the starter's cruise altitude
ALPHAS = (-2.0, 0.0, 2.0, 4.0, 6.0, 8.0)   # degrees
BETAS = (0.0, 4.0)
# Reference: the OpenVSP wing fitted to the starter's baseline wing (area and span to the centreline, mean
# aerodynamic chord), and the starter's CG station.
SREF, BREF, CREF = 61.5, 9.04, 8.72
CG = (18269.0, 0.0, 0.0)      # mm


def call(g, block, given, name):
    inputs = []
    for spec in block.inputs:
        if spec.name in given:
            value = given[spec.name]
            inputs.append(value(g) if callable(value) else value)
        else:
            inputs.append(spec.default(g))
    return g.block(block.signature, block.output_type, inputs, name)


def build(path=NOTEBOOK):
    fus, wing, model, check, mesh = (mb.fuselage_block(), mb.wing_block(), mb.model_block(), mb.model_check_block(),
                                     mb.model_mesh_block())
    flight, solver, ana, coef = a.flight_block(), v.solver_block(), v.analysis_block(), v.coefficient_block()
    wave, parasite = d.wave_block(), d.parasite_block()
    c_vspaero, c_wave, c_parasite = (bb.vspaero_contribution_block(), bb.drag_contribution_block("wave"),
                                     bb.drag_contribution_block("parasite"))
    database, db_table, center, plot = bb.database_block(), bb.database_table_block(), v.center_block(), bb.plot_block()
    blocks = [fus, wing, model, check, mesh, flight, solver, ana, coef, wave, parasite, c_vspaero, c_wave, c_parasite,
              database, db_table, center, plot]
    bodies = ["Fuselage Assembly", "Main Wing Installed", "Canards Installed", "Main HStab", "Top T Stab", "VStab",
              "Wing Span", "Wing Inner Sweep"]
    dv = x59.Derived()
    dv.keep(bodies)
    dv.expose("Wing Span", InputSpec("Wing Span", "real", "Main wing span, tip to tip (the starter's variable).", LENGTH,
                                     lambda g: g.real(9.0, LENGTH)))
    dv.expose("Wing Inner Sweep", InputSpec("Wing Inner Sweep", "real", "Quarter-chord sweep of the inner wing panel "
                                            "(the starter's variable).", ANGLE,
                                            lambda g: g.real(math.radians(68.3), ANGLE)))
    g = dv.graph()
    g.comment("# X-59 aero database\nThe X-59 hackathon starter's parametric geometry (above, unchanged; " + x59.CREDIT +
              ") turned into an OpenVSP model and analysed on one Flight Condition: VSPAERO (vortex lattice), Wave Drag "
              "and Parasite Drag, combined by **Aero Database (VSP)** into CL, CD, CY, Cl, Cm and Cn at every Mach x "
              "alpha x beta condition.\n\n**Wing Span** and **Wing Inner Sweep** are inputs of this notebook and are "
              "recorded in every row: `examples/x59_doe.py` runs it over a design of experiments and collects one "
              "table. The blocks call Python and OpenVSP with Run Command: allow them when nTop asks.")
    g.comment("**1. Components.** Each body of the starter becomes an OpenVSP component; a component block takes a list "
              "of bodies. The fuselage assembly (fuselage and nacelle, unioned in the starter) is an OpenVSP Fuselage "
              "in no VSPAERO set: VSPAERO solves the lifting surfaces with its vortex lattice, while Wave Drag and "
              "Parasite Drag take the whole model (see the README for why). The wing, canards and tails are one Wing "
              "block (four bodies, each with its span axis, mirrored about the x-z plane), the fin another (no "
              "symmetry, span axis +z), all in the Thin set. The fin is cut 65 mm below the T-tail's plane, because thin "
              "surfaces that cross near a tip make VSPAERO fail above Mach 1.")
    cut = g.block("box<point,real,real,real>", "box", [g.point(25000, 0, 3100 + 5000), g.real(20.0, LENGTH),
                                                        g.real(20.0, LENGTH), g.real(10.0, LENGTH)], "Above T-Tail")
    cut_list = g.block("core.list<implicit>", "list<implicit>", [cut], top=False)
    fin = g.block("boolean_subtract<blend_enum,real_field,implicit,list<implicit>>[5.44.0]", "implicit",
                  [g.enum("blend_enum", 0), g.real(0.0, LENGTH), dv.id_of("VStab"), cut_list], "Fin Below T-Tail")

    def texts(*items):
        return lambda g: g.block("core.list<text>", "list<text>", [g.text(t) for t in items], top=False)

    def implicits(name, items):
        return g.block("core.list<implicit>", "list<implicit>", items, name)

    def vectors(*items):
        return lambda g: g.block("core.list<vector>", "list<vector>", [g.vector(*x) for x in items], top=False)

    lifting = implicits("Lifting Surfaces", [dv.id_of(n) for n in ("Main Wing Installed", "Canards Installed",
                                                                   "Main HStab", "Top T Stab")])
    comps = [
        call(g, fus, {"Bodies": implicits("Fuselage Bodies", [dv.id_of("Fuselage Assembly")]),
                      "Names": texts("Fuselage")}, "Fuselage Component"),
        call(g, wing, {"Bodies": lifting, "Names": texts("Wing", "Canard", "HStab", "TStab"),
                       "Span Axes": vectors((0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0))}, "Wing Components"),
        call(g, wing, {"Bodies": implicits("Fin Bodies", [fin]), "Names": texts("VStab"),
                       "Span Axes": vectors((0, 0, 1)), "Symmetry": lambda g: g.choice(list(mb.SYMMETRY), 0)},
             "Fin Component"),
    ]
    comp_list = g.block("core.list<text>", "list<text>", comps, "Components")
    g.comment("**2. OpenVSP model.** Fitted to the bodies' section contours; Model Check gives each component's largest "
              "distance from its body (metres) and Model Mesh shows OpenVSP's surfaces.")
    run = {"Working Directory": lambda g: g.file_path(RUNS)}
    vsp_model = call(g, model, {"Components": comp_list, **run, "Case Name": lambda g: g.text("x59_model")}, "X-59 Model")
    call(g, check, {"Model": vsp_model}, "Largest Deviation")
    call(g, mesh, {"Model": vsp_model}, "OpenVSP Surfaces")
    g.comment(f"**3. Flight condition**, shared by every analysis and the database: Mach {', '.join('%g' % m for m in MACHS)} "
              f"at {ALTITUDE_FT:,.0f} ft (the US Standard Atmosphere 1976 gives the rest), alpha "
              f"{', '.join('%g' % x for x in ALPHAS)} deg, beta {', '.join('%g' % x for x in BETAS)} deg; Sref, bref and "
              "cref of the baseline wing, moments about the starter's CG station.")
    deg = math.radians(1.0)
    fc = call(g, flight, {"Mach": lambda g: g.block("core.list<real>", "list<real>", [g.real(m) for m in MACHS], "Mach"),
                          "Altitude": lambda g: g.block("core.list<real>", "list<real>",
                                                        [g.real(ALTITUDE_FT * 0.3048, LENGTH)], "Altitude"),
                          "Angle of Attack": lambda g: g.block("core.list<real>", "list<real>",
                                                               [g.real(x * deg, ANGLE) for x in ALPHAS], "Angle of Attack"),
                          "Sideslip": lambda g: g.block("core.list<real>", "list<real>",
                                                        [g.real(x * deg, ANGLE) for x in BETAS], "Sideslip"),
                          "Reference Area": lambda g: g.real(SREF, AREA), "Reference Span": lambda g: g.real(BREF, LENGTH),
                          "Reference Chord": lambda g: g.real(CREF, LENGTH),
                          "Moment Reference Point": lambda g: g.point(*CG),
                          **run, "Case Name": lambda g: g.text("x59_flight")}, "Flight Condition")
    g.comment("**4. VSPAERO** (vortex lattice on the Thin set), every condition in one run, with a fixed wake for speed "
              "(a relaxed wake changed CL by under 0.1% on this model at Mach 1.6).")
    ss = call(g, solver, {"Processors": lambda g: g.integer(16),
                          "Wake": lambda g: g.choice(["Relaxed", "Fixed"], 1)}, "Solver Settings")
    aero = call(g, ana, {"OpenVSP Model": vsp_model, "Flight Condition": fc, "Solver Settings": ss,
                         **run, "Case Name": lambda g: g.text("x59_vspaero")}, "VSPAERO")
    g.comment("**5. Wave Drag** at the supersonic Mach numbers of the flight condition (zero-lift area rule over the "
              "whole model).")
    wv = call(g, wave, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                        "Case Name": lambda g: g.text("x59_wave")}, "Wave Drag")
    g.comment("**6. Parasite Drag** of every component at each flight point (Mach number and altitude).")
    pd = call(g, parasite, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                            "Case Name": lambda g: g.text("x59_parasite")}, "Parasite Drag")
    g.comment("**7. Aero database.** Each analysis becomes a contribution: VSPAERO's inviscid forces and moments, wave "
              "drag (supersonic) and parasite drag. Aero Database sums them at every condition of the flight condition "
              "into CL, CD, CY, Cl, Cm, Cn, with the two design variables in the first columns; each run's rows are "
              "also appended to `runs/x59_doe.csv`. Other sources (Modified Newtonian above some Mach number, a "
              "Viscous Correction in place of Parasite Drag) go into the same list through Coefficients Contribution "
              "(VSP).")
    span_m = g.block("divide<real,real>", "real", [dv.id_of("Wing Span"), g.real(1.0, LENGTH)], top=False)
    sweep_deg = g.block("divide<real,real>", "real", [dv.id_of("Wing Inner Sweep"), g.real(math.radians(1.0), ANGLE)],
                        top=False)
    names = g.block("core.list<text>", "list<text>", [g.text("Wing_Span_m"), g.text("Wing_Inner_Sweep_deg")], top=False)
    vals = g.block("core.list<real>", "list<real>", [span_m, sweep_deg], top=False)
    contributions = g.block("core.list<text>", "list<text>", [
        call(g, c_vspaero, {"Results": aero}, "VSPAERO Contribution"),
        call(g, c_wave, {"Results": wv}, "Wave Drag Contribution"),
        call(g, c_parasite, {"Results": pd}, "Parasite Drag Contribution")], "Contributions")
    db = call(g, database, {"Flight Condition": fc, "Contributions": contributions, "Variable Names": names,
                            "Variable Values": vals, "Append To": lambda g: g.text(DOE_TABLE), **run,
                            "Case Name": lambda g: g.text("x59_database")}, "Aero Database")
    table = call(g, db_table, {"Results": db}, "Database Table")
    cl = call(g, coef, {"Results": aero, "Coefficient": lambda g: g.choice(list(v.POLAR_COLUMNS), v.POLAR_COLUMNS.index("CLi"))},
              "CLi")
    g.comment("**8. Aerodynamic center** (neutral point) at each Mach number and beta, from the database's pitching "
              "moment against normal force: positions, to compare with the CG station (18.27 m) for the static margin.")
    ac = call(g, center, {"Aero Database": db, "Quantity": lambda g: g.choice(list(v.CENTERS), 1)}, "Aerodynamic Center")
    g.comment("**9. Drag polar**, CL against CD at zero sideslip, one curve per Mach number. It opens in a window and is "
              "saved as `runs/x59_polar.png`.")
    img = call(g, plot, {"Results": lambda g: g.block("core.list<text>", "list<text>", [db], top=False),
                         "Filter": lambda g: g.choice(list(bb.FILTERS), 1), **run,
                         "Case Name": lambda g: g.text("x59_polar")}, "Drag Polar")
    out = g.block("core.group<any>", "group", [table, cl, ac, img], "Outputs")
    os.makedirs(EXAMPLES, exist_ok=True)
    dv.save(path, g, blocks=blocks, output=out, output_type="group", display="X-59 Aero Database",
            description="X-59 starter geometry to an OpenVSP model and a six-degree-of-freedom aero database.")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    path = build()
    print(path)
    if args.run:
        from ntrun import run_notebook
        out, log, dt = run_notebook(path, timeout=7200)
        print(f"{dt:.0f} s, output {'ok' if out else 'missing'}")
        print("\n".join(l for l in log.splitlines() if "[E]" in l)[:3000])


if __name__ == "__main__":
    main()
