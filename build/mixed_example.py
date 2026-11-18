"""Write (and optionally run) the mixed-analysis example: one aero database from VSPAERO, Wave Drag, Modified Newtonian
and the Viscous Correction, each over its own Mach range.

    sh build/py.sh build/mixed_example.py [--run]

The X-59 starter's geometry again (its parametric graph, unchanged), on one Flight Condition of three flight points
at three altitudes, Mach 0.8, 1.6 and 4: VSPAERO below Mach 2.5 (its vortex lattice on the lifting surfaces),
Wave Drag up to Mach 2.5, Modified Newtonian on the whole aircraft from Mach 2.5, and the Viscous Correction's skin
friction at every Mach number in place of Parasite Drag. Aero Database sums them, Center of Pressure gives the
aerodynamic center and Polar Plot draws CL against CD.
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
from x59_example import BREF, CG, CREF, SREF, call  # noqa: E402

EXAMPLES = ROOT + "/examples"
RUNS = ROOT + "/runs"
NOTEBOOK = EXAMPLES + "/X-59 Mixed Analyses.ntop"
POINTS = ((0.8, 30000.0), (1.6, 50000.0), (4.0, 80000.0))   # (Mach, altitude ft)
ALPHAS = (0.0, 4.0)
HANDOVER = 2.5   # Mach: VSPAERO and Wave Drag below, Modified Newtonian from here


def build(path=NOTEBOOK):
    fus, wing, model = mb.fuselage_block(), mb.wing_block(), mb.model_block()
    flight, solver, ana = a.flight_block(), v.solver_block(), v.analysis_block()
    wave, newton, visc = d.wave_block(), a.newtonian_block(), a.viscous_block()
    c_vspaero, c_wave, c_coef = (bb.vspaero_contribution_block(), bb.drag_contribution_block("wave"),
                                 bb.coefficients_contribution_block())
    database, db_table, center, plot = bb.database_block(), bb.database_table_block(), v.center_block(), bb.plot_block()
    blocks = [fus, wing, model, flight, solver, ana, wave, newton, visc, c_vspaero, c_wave, c_coef, database, db_table, center,
              plot]
    names = ["Fuselage Assembly", "Main Wing Installed", "Canards Installed", "Main HStab", "Top T Stab", "VStab"]
    dv = x59.Derived()
    dv.keep(names)
    g = dv.graph()
    g.comment("# X-59: mixing analyses in one aero database\nThe X-59 hackathon starter's geometry (above, unchanged; "
              + x59.CREDIT + ") on one **Flight Condition** of three flight points: Mach 0.8 at 30,000 ft, 1.6 at "
              "50,000 ft and 4 at 80,000 ft, each at alpha 0 and 4 degrees. Each analysis takes the flight condition "
              "and uses what it needs; the contributions say over which Mach numbers each counts:\n\n"
              "* **VSPAERO** (vortex lattice on the lifting surfaces) below Mach 2.5 (its own Mach Max);\n"
              "* **Wave Drag** (area rule) at the supersonic Mach numbers up to 2.5 (the contribution's Mach Max);\n"
              "* **Modified Newtonian** on the whole aircraft from Mach 2.5 (its own Mach Min);\n"
              "* **Viscous Correction** (skin friction from the standard atmosphere at each altitude) at every Mach "
              "number, in place of Parasite Drag.\n\nThe blocks call Python and OpenVSP with Run Command: allow them "
              "when nTop asks.")
    g.comment("**1. OpenVSP model** for VSPAERO and Wave Drag: the fuselage in no VSPAERO set, the wing, canards and "
              "horizontal tails as one Wing block (mirrored about the x-z plane), the fin (cut below the T-tail) as "
              "another.")
    cut = g.block("box<point,real,real,real>", "box", [g.point(25000, 0, 3100 + 5000), g.real(20.0, LENGTH),
                                                        g.real(20.0, LENGTH), g.real(10.0, LENGTH)], "Above T-Tail")
    fin = g.block("boolean_subtract<blend_enum,real_field,implicit,list<implicit>>[5.44.0]", "implicit",
                  [g.enum("blend_enum", 0), g.real(0.0, LENGTH), dv.id_of("VStab"),
                   g.block("core.list<implicit>", "list<implicit>", [cut], top=False)], "Fin Below T-Tail")

    def lst(typ, items, name=None):
        return g.block(f"core.list<{typ}>", f"list<{typ}>", items, name, top=name is not None)

    comps = lst("text", [
        call(g, fus, {"Bodies": lst("implicit", [dv.id_of("Fuselage Assembly")]),
                      "Names": lst("text", [g.text("Fuselage")])}, "Fuselage Component"),
        call(g, wing, {"Bodies": lst("implicit", [dv.id_of(n) for n in names[1:5]], "Lifting Surfaces"),
                       "Names": lst("text", [g.text(t) for t in ("Wing", "Canard", "HStab", "TStab")]),
                       "Span Axes": lst("vector", [g.vector(0, 1, 0)] * 4)}, "Wing Components"),
        call(g, wing, {"Bodies": lst("implicit", [fin]), "Names": lst("text", [g.text("VStab")]),
                       "Span Axes": lst("vector", [g.vector(0, 0, 1)]),
                       "Symmetry": lambda g: g.choice(list(mb.SYMMETRY), 0)}, "Fin Component")], "Components")
    run = {"Working Directory": lambda g: g.file_path(RUNS)}
    vsp_model = call(g, model, {"Components": comps, **run, "Case Name": lambda g: g.text("mixed_model")}, "X-59 Model")
    g.comment("**2. Flight condition**: Mach numbers paired with altitudes (US Standard Atmosphere 1976), every pair at "
              "each angle of attack; the reference area, span and chord of the baseline wing and the CG station.")
    deg = math.radians(1.0)
    fc = call(g, flight, {"Mach": lambda g: lst("real", [g.real(m) for m, _ in POINTS], "Mach"),
                          "Altitude": lambda g: lst("real", [g.real(h * 0.3048, LENGTH) for _, h in POINTS], "Altitude"),
                          "Angle of Attack": lambda g: lst("real", [g.real(x * deg, ANGLE) for x in ALPHAS],
                                                           "Angle of Attack"),
                          "Sideslip": lambda g: lst("real", [g.real(0.0, ANGLE)], "Sideslip"),
                          "Reference Area": lambda g: g.real(SREF, AREA), "Reference Span": lambda g: g.real(BREF, LENGTH),
                          "Reference Chord": lambda g: g.real(CREF, LENGTH),
                          "Moment Reference Point": lambda g: g.point(*CG),
                          **run, "Case Name": lambda g: g.text("mixed_flight")}, "Flight Condition")
    g.comment(f"**3. Analyses.** VSPAERO only below Mach {HANDOVER:g} (its Mach Max), with a fixed wake (a relaxed one "
              "was still swinging the pitching moment by 40% over its last iteration at Mach 0.8, and VSPAERO's "
              "convergence check stopped the run); Wave Drag at the supersonic Mach "
              f"numbers; Modified Newtonian and the Viscous Correction on the union of the starter's bodies, Newtonian "
              f"from Mach {HANDOVER:g} (its Mach Min), friction at every Mach number with the freestream as edge state.")
    ss = call(g, solver, {"Processors": lambda g: g.integer(16),
                          "Wake": lambda g: g.choice(["Relaxed", "Fixed"], 1)}, "Solver Settings")
    aero = call(g, ana, {"OpenVSP Model": vsp_model, "Flight Condition": fc, "Solver Settings": ss,
                         "Mach Max": lambda g: g.real(HANDOVER), **run,
                         "Case Name": lambda g: g.text("mixed_vspaero")}, "VSPAERO")
    wv = call(g, wave, {"OpenVSP Model": vsp_model, "Flight Condition": fc, **run,
                        "Case Name": lambda g: g.text("mixed_wave")}, "Wave Drag")
    whole = g.block("boolean_union<blend_enum,real_field,list<implicit>>[5.44.0]", "implicit",
                    [g.enum("blend_enum", 0), g.real(0.0, LENGTH),
                     lst("implicit", [dv.id_of(n) for n in names])], "Whole Aircraft")
    nw = call(g, newton, {"Body": whole, "Flight Condition": fc, "Mach Min": lambda g: g.real(HANDOVER)},
              "Modified Newtonian")
    vc = call(g, visc, {"Body": whole, "Flight Condition": fc,
                        "Inviscid Pressure": lambda g: g.choice(list(a.PRESSURE_SOURCES), 1)}, "Viscous Correction")
    g.comment(f"**4. Contributions and database.** Wave Drag's contribution stops at Mach {HANDOVER:g} (Mach Max), where "
              "Modified Newtonian's pressure (which includes the wave drag) takes over. Each contribution's share has "
              "its own columns in the database.")
    contributions = lst("text", [
        call(g, c_vspaero, {"Results": aero}, "VSPAERO Contribution"),
        call(g, c_wave, {"Results": wv, "Mach Max": lambda g: g.real(HANDOVER)}, "Wave Drag Contribution"),
        call(g, c_coef, {"Results": nw, "Name": lambda g: g.text("Newtonian")}, "Newtonian Contribution"),
        call(g, c_coef, {"Results": vc, "Name": lambda g: g.text("Friction")}, "Friction Contribution")], "Contributions")
    db = call(g, database, {"Flight Condition": fc, "Contributions": contributions, **run,
                            "Case Name": lambda g: g.text("mixed_database")}, "Aero Database")
    table = call(g, db_table, {"Results": db}, "Database Table")
    ac = call(g, center, {"Aero Database": db, "Quantity": lambda g: g.choice(list(v.CENTERS), 1)}, "Aerodynamic Center")
    g.comment("**5. Drag polar**, CL against CD, one curve per flight point; it opens in a window and is saved as "
              "`runs/mixed_polar.png`.")
    img = call(g, plot, {"Results": lambda g: lst("text", [db]), **run,
                         "Case Name": lambda g: g.text("mixed_polar")}, "Drag Polar")
    out = g.block("core.group<any>", "group", [table, ac, img], "Outputs")
    os.makedirs(EXAMPLES, exist_ok=True)
    dv.save(path, g, blocks=blocks, output=out, output_type="group", display="X-59 Mixed Analyses",
            description="X-59 starter geometry: VSPAERO, Wave Drag, Modified Newtonian and Viscous Correction in one "
                        "aero database.")
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
