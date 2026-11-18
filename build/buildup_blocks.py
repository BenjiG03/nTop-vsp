"""Aero buildup blocks: contributions from any analysis, the Aero Database that sums them, its readers and plots.

A **contribution** is one line of text naming a source of forces and moments and how to use it: a VSPAERO run, a
Wave Drag run, a Parasite Drag run, or a coefficient table made in nTop (Modified Newtonian, Viscous Correction,
Aero Coefficients). **Aero Database (VSP)** takes a Flight Condition (its conditions and reference) and a list of
contributions; it takes each at the database's conditions (interpolating where needed), moves it to the reference
and moment point, keeps it only inside its Mach range, and sums them (``vsp_ntop.run_database``).
"""
from __future__ import annotations

from ntopfile import Graph, InputSpec

from kit import (FLIGHT_INPUT, Kit, analysis, d_choice, d_real, d_text, make_block, model_inputs, run_inputs)

USE = ("All forces and moments", "Drag only (force along the freestream)", "All but drag")
USE_KEYS = ("all", "drag", "nodrag")


def range_inputs(mach_min=0.0, what="it"):
    return (
        InputSpec("Mach Min", "real", f"Lowest Mach number at which {what} applies.", default=d_real(mach_min)),
        InputSpec("Mach Max", "real", f"{what[0].upper() + what[1:]} applies below this Mach number (not at it), so "
                  "contributions that hand over at one Mach number do not overlap.", default=d_real(1000.0)),
        InputSpec("Scale", "real", "Factor on every force and moment of this contribution (1: as computed).",
                  default=d_real(1.0)),
    )


USE_INPUT = InputSpec("Use", "choice", "All: every force and moment. Drag only: just the force component along the "
                      "freestream, no moments (to take only drag from this source). All but drag: everything else "
                      "(to take lift, side force and moments from this source and drag from others).",
                      default=d_choice(USE))


def range_text(k, inp, with_use=True):
    parts = ["; mach_min: ", k.num(inp("Mach Min")), "; mach_max: ", k.num(inp("Mach Max")), "; scale: ",
             k.num(inp("Scale"))]
    if with_use:
        parts += ["; use: ", k.pick(inp("Use"), USE_KEYS)]
    return parts


RANGE_DOC = ("**Mach range.** The database uses this contribution from Mach Min up to (not including) Mach Max, so "
             "one source can hand over to another: for example VSPAERO below Mach 3 and Modified Newtonian from Mach 3.")

# --- contributions ---------------------------------------------------------------------------------------------------

VSPAERO_TERMS = ("Inviscid (surface pressure: CFi, CMi)", "Total (with VSPAERO's own skin-friction estimate)")


def vspaero_contribution_block():
    specs = (
        InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
        InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text("VSPAERO")),
        InputSpec("Terms", "choice", "Inviscid: VSPAERO's surface-pressure forces and moments, for databases whose "
                  "friction comes from Parasite Drag or a viscous correction. Total: with VSPAERO's own skin-friction "
                  "estimate (the o columns).", default=d_choice(VSPAERO_TERMS)),
        USE_INPUT,
    ) + range_inputs(0.0)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**VSPAERO contribution.** One line of text for **Aero Database (VSP)**: the run's results path and how "
              "to use them. The database reads the body-axis force and moment coefficients of every case from the "
              "run's polar, and its reference and moment point from the run's settings.\n\n" + RANGE_DOC)
    out = k.cat(["kind: vspaero; name: ", inp("Name"), "; results: ", inp("Results"), "; terms: ",
                 k.pick(inp("Terms"), ("inviscid", "total"))] + range_text(k, inp), name="Contribution")
    return make_block("contribution_vspaero", "VSPAERO Contribution (VSP)",
                      "A VSPAERO run as a contribution to **Aero Database (VSP)**: its forces and moments over its Mach x "
                      "alpha x beta cases, inviscid or with VSPAERO's own friction, within a Mach range.",
                      specs, "text", g, out)


def drag_contribution_block(kind):
    wave = kind == "wave"
    title = "Wave Drag" if wave else "Parasite Drag"
    specs = (
        InputSpec("Results", "text", f"Output of {title} (VSP)."),
        InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text(title)),
    ) + range_inputs(1.0 if wave else 0.0)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment(f"**{title} contribution.** One line of text for **Aero Database (VSP)**. {title} is a drag coefficient "
              "per Mach number: the database interpolates it linearly in Mach and applies it along the freestream at "
              "every angle, through the database's moment point (no moments)."
              + (" Mach Min defaults to 1: wave drag is zero subsonic." if wave else "") + "\n\n" + RANGE_DOC)
    out = k.cat([f"kind: {kind}; name: ", inp("Name"), "; results: ", inp("Results")] + range_text(k, inp, False),
                name="Contribution")
    return make_block(f"contribution_{kind}", f"{title} Contribution (VSP)",
                      f"A {title} run as a contribution to **Aero Database (VSP)**: a drag coefficient per Mach number, "
                      "applied along the freestream at every angle, within a Mach range.", specs, "text", g, out)


def coefficients_contribution_block():
    specs = (
        InputSpec("Results", "text", "Output of Modified Newtonian (VSP), Viscous Correction (VSP) or Aero Coefficients "
                  "(VSP): a coefficient table, one row per flight condition."),
        InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text("Newtonian")),
        USE_INPUT,
    ) + range_inputs(0.0)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Coefficients contribution.** A coefficient table computed in nTop as one line of text for **Aero "
              "Database (VSP)**. The database takes it where it has values: over its own Mach range (an analysis "
              "computed between its own Mach Min and Max) and, for tables by altitude (Viscous Correction, Aero "
              "Coefficients), at its own altitudes; between its Mach numbers and angles it interpolates.\n\n" + RANGE_DOC)
    out = k.cat(["kind: coefficients; name: ", inp("Name"), "; results: ", inp("Results")] + range_text(k, inp),
                name="Contribution")
    return make_block("contribution_coefficients", "Coefficients Contribution (VSP)",
                      "Coefficients computed in nTop (Modified Newtonian, Viscous Correction, Aero Coefficients) as a "
                      "contribution to **Aero Database (VSP)**, within a Mach range.", specs, "text", g, out)


# --- database ----------------------------------------------------------------------------------------------------------

def _list_default(values, units=None):
    return lambda g: g.block("core.list<real>", "list<real>", [g.real(v, units) for v in values], top=False)


def database_block():
    own = (
        InputSpec("Flight Condition", "text", "Output of Flight Condition (VSP): the database's conditions (every "
                  "flight point x sideslip x angle of attack) and its reference area, span, chord and moment point. "
                  "Normally the one the analyses ran on."),
        InputSpec("Contributions", "list<text>", "The contributions to sum: outputs of VSPAERO, Wave Drag, Parasite "
                  "Drag and Coefficients Contribution (VSP), in any number and order."),
        InputSpec("Variable Names", "list<text>", "Names of design variables to record with every row (e.g. geometry "
                  "inputs of a design of experiments). Empty list item: none.",
                  default=lambda g: g.block("core.list<text>", "list<text>", [g.text("")], top=False)),
        InputSpec("Variable Values", "list<real>", "Their values, in the same order (as plain numbers: convert lengths "
                  "to one unit first).", default=_list_default([0.0])),
        InputSpec("Append To", "text", "A CSV file that collects the rows of every run (created with a header when "
                  "new): a design-of-experiments table. Empty: none.", default=d_text("")),
    )
    specs = own + run_inputs("database")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Aero database.** A Python driver sums the contributions at every condition of the flight condition. "
              "Each contribution applies only inside its Mach range (and, for coefficient tables, where it has values); "
              "it is taken at the condition or interpolated (multilinear in Mach, alpha and beta; Wave and Parasite "
              "Drag in Mach, at the condition's altitude), turned into body-axis forces and moments, moved to the "
              "flight condition's moment point (M + r x F) and normalised by its reference. The totals give CL, CD, "
              "CY (wind axes, VSPAERO's convention) and Cl, Cm, Cn (aircraft body-axis signs: right wing down, nose "
              "up, nose right).")
    g.comment("**Outputs.** `<Case Name>.ntvsp.database.csv`: design variables, Mach, Altitude, AoA, Beta, CL, CD, "
              "CY, Cl, Cm, "
              "Cn, L_over_D, CA, CN, then each contribution's own six coefficients (columns named after it). Also the "
              "centre of pressure and aerodynamic centre files read by **Center of Pressure (VSP)**. With Append To, "
              "the rows are added to that file as well.")
    contrib = k.cat(["contribution = ", k.join(inp("Contributions"), "\ncontribution = "), "\n"])
    lines = [
        contrib,
        k.line("flight", inp("Flight Condition")),
        k.line("append_to", inp("Append To")),
        k.line("variables", k.join(inp("Variable Names"))),
        k.line("values", k.join(k.num(inp("Variable Values")))),
    ]
    out = analysis(g, k, inp, "run_database", lambda ready: [], lines, "Aero Database", needs_openvsp=False)
    return make_block(
        "aero_database", "Aero Database (VSP)",
        "A six-degree-of-freedom aero database on a Flight Condition from any mix of contributions: VSPAERO, Wave Drag, "
        "Parasite Drag, and coefficients computed in nTop (Modified Newtonian, Viscous Correction, Aero "
        "Coefficients). Each applies within its own Mach range and is converted to the flight condition's reference "
        "and moment point; the sum gives CL, CD, CY, Cl, Cm and Cn at every condition (Mach, altitude, angle of "
        "attack, sideslip), with each contribution's share in its own columns. Outputs the results path for the database readers, Center of Pressure and Polar Plot (VSP); Append "
        "To collects runs with different design variables into one design-of-experiments file.",
        specs, "text", g, out)


DATABASE_COLUMNS = ("CL", "CD", "CY", "Cl", "Cm", "Cn", "L_over_D", "CA", "CN", "Mach", "Altitude", "AoA", "Beta")


def database_coefficient_block():
    specs = (
        InputSpec("Results", "text", "Output of Aero Database (VSP)."),
        InputSpec("Coefficient", "choice", "A total coefficient (CL, CD, CY in wind axes; Cl, Cm, Cn in body-axis "
                  "signs; CA, CN body-axis forces) or a condition (Mach, Altitude in metres, AoA and Beta in degrees).",
                  default=d_choice(DATABASE_COLUMNS)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Database column.** Import Table reads `<Case Name>.ntvsp.database.csv`; the chosen column comes out "
              "as a list of reals, one per flight condition.")
    table = g.block("import_table<file_path,bool>", "table", [k.cat([inp("Results"), ".ntvsp.database.csv"]),
                                                              g.boolean(True)], "Import Database")
    col = g.block("extract_column<table,text>", "list<text>", [table, k.pick(inp("Coefficient"), DATABASE_COLUMNS)],
                  top=False)
    out = g.block("scalar_from_text<text>", "real", [col], "Values")
    return make_block("aero_database_coefficient", "Aero Database Coefficient (VSP)",
                      "One column of an Aero Database (VSP), as a list in row order.", specs, "list<real>", g, out)


def database_table_block():
    specs = (InputSpec("Results", "text", "Output of Aero Database (VSP)."),)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Database table.** `<Case Name>.ntvsp.database.csv` through Import Table, every column including "
              "each contribution's own coefficients.")
    out = g.block("import_table<file_path,bool>", "table", [k.cat([inp("Results"), ".ntvsp.database.csv"]),
                                                            g.boolean(True)], "Database Table")
    return make_block("aero_database_table", "Aero Database Table (VSP)",
                      "The whole Aero Database (VSP) as a table.", specs, "table", g, out)


# --- plots -------------------------------------------------------------------------------------------------------------

PLOT_COLUMNS = ("CL", "CD", "CY", "Cl", "Cm", "Cn", "L_over_D", "CA", "CN", "AoA", "Beta", "Mach", "Altitude")
ALONG = ("Angle of attack", "Sideslip", "Mach number")
ALONG_KEYS = ("AoA", "Beta", "Mach")
FILTERS = ("All conditions", "Zero sideslip only", "Zero angle of attack only")
FILTER_KEYS = ("all", "beta0", "alpha0")


def plot_block():
    own = (
        InputSpec("Results", "list<text>", "Outputs of Aero Database (VSP): one or more databases to compare."),
        InputSpec("Labels", "list<text>", "A legend label per result, in order. Empty: the case names.",
                  default=lambda g: g.block("core.list<text>", "list<text>", [g.text("")], top=False)),
        InputSpec("X Axis", "choice", "Quantity along the horizontal axis.", default=d_choice(PLOT_COLUMNS, 1)),
        InputSpec("Y Axis", "choice", "Quantity along the vertical axis.", default=d_choice(PLOT_COLUMNS, 0)),
        InputSpec("Curves Along", "choice", "The condition that varies along each curve; every combination of the "
                  "other conditions (Mach, altitude, angle of attack, sideslip) and every result is its own curve.",
                  default=d_choice(ALONG)),
        InputSpec("Filter", "choice", "Keep only some rows before plotting.", default=d_choice(FILTERS)),
        InputSpec("Title", "text", "Plot title. Empty: Y against X.", default=d_text("")),
        InputSpec("Show Window", "choice", "On: the plot also opens in a window (matplotlib's viewer, which zooms, "
                  "pans and saves) as soon as it is made; nTop does not wait for it to close. Off: only the image file.",
                  default=d_choice(("On", "Off"))),
    )
    specs = own + run_inputs("polar")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Polar plot.** A Python driver plots the chosen columns of each database with matplotlib, one curve "
              "per combination of the conditions that are not along the curve. Curves of one result are shaded from "
              "light to dark by condition, results are told apart by line style (or by colour when each result has "
              "one curve). The image is `<Case Name>.png` in the Working Directory; Show Window opens it in a window "
              "too. Needs matplotlib in the Python of `NTOP_PYTHON_PATH`.")
    lines = [
        k.line("results", k.join(inp("Results"))),
        k.line("labels", k.join(inp("Labels"))),
        k.line("x", k.pick(inp("X Axis"), PLOT_COLUMNS)),
        k.line("y", k.pick(inp("Y Axis"), PLOT_COLUMNS)),
        k.line("along", k.pick(inp("Curves Along"), ALONG_KEYS)),
        k.line("filter", k.pick(inp("Filter"), FILTER_KEYS)),
        k.line("title", inp("Title")),
        k.line("show", k.pick(inp("Show Window"), ("1", "0"))),
    ]
    base = analysis(g, k, inp, "run_plot", lambda ready: [], lines, "Polar Plot", needs_openvsp=False)
    out = k.cat([base, ".png"], name="Image")
    return make_block("polar_plot", "Polar Plot (VSP)",
                      "Plots one quantity of one or more Aero Databases against another (CL against CD, Cm against "
                      "alpha, ...) with matplotlib, picked from dropdowns. Saves `<Case Name>.png` and, by default, "
                      "opens the plot in a window. Outputs the image path.", specs, "text", g, out)
