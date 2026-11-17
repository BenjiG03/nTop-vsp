"""Aero buildup blocks: contributions from any analysis, the Aero Database that sums them, its readers and plots.

A **contribution** is one line of text naming a source of forces and moments and how to use it: a VSPAERO run, a
Wave Drag run, a Parasite Drag run, or coefficient rows made in nTop (Modified Newtonian, Viscous Correction, or
any six numbers through Aero Coefficient Row). **Aero Database (VSP)** takes a list of contributions, interpolates
each onto the database's flight conditions, moves it to one reference and moment point, keeps it only inside its
Mach range, and sums them (``vsp_ntop.run_database``).
"""
from __future__ import annotations

import math

from ntopfile import Graph, InputSpec

from kit import (ANGLE, AREA, LENGTH, Kit, analysis, d_choice, d_point_mm, d_real, d_text, make_block, model_inputs,
                 prop, run_inputs)

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


def rows_contribution_block():
    specs = (
        InputSpec("Coefficient Rows", "list<text>", "One row per flight condition: the output of Aero Coefficient Row "
                  "(VSP), or of ntop_newtonian's Database Row blocks (Modified Newtonian, Viscous Correction), swept "
                  "over a list of conditions."),
        InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text("Newtonian")),
        USE_INPUT,
    ) + range_inputs(0.0)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Coefficient rows contribution.** The rows (each a flight condition, six coefficients, a reference and "
              "a moment point) become one line of text for **Aero Database (VSP)**. At a database condition that a "
              "row has, the database uses that row; between rows it interpolates, which needs every combination of "
              "the neighbouring Mach numbers, angles of attack and sideslips. Lists swept in nTop are paired item by "
              "item, so either compute the rows at the database's own conditions (or make the rows the database's "
              "conditions), or give paired lists that form a grid.\n\n" + RANGE_DOC)
    out = k.cat(["kind: rows; name: ", inp("Name"), "; rows: ", k.join(inp("Coefficient Rows"), "|")]
                + range_text(k, inp), name="Contribution")
    return make_block("contribution_rows", "Coefficient Rows Contribution (VSP)",
                      "Coefficients computed in nTop (Modified Newtonian, Viscous Correction, or any six numbers) as a "
                      "contribution to **Aero Database (VSP)**, within a Mach range.", specs, "text", g, out)


ROW_AXES = ("Body: CA, CY, CN; moments about x aft, y right, z up (VSPAERO, Modified Newtonian, Viscous Correction)",
            "Wind: CD, CY, CL; moments in aircraft body-axis signs (Cl right wing down, Cm nose up, Cn nose right)")


def row_block():
    zero = d_real(0.0)
    specs = (
        InputSpec("Axes", "choice", "How the six coefficients are given. Body: axial, side and normal force along x "
                  "(aft), y (right), z (up), and right-handed moments about those axes, the convention of VSPAERO and "
                  "of ntop_newtonian's blocks. Wind: drag, side force and lift, with rolling, pitching and yawing "
                  "moments in aircraft body-axis signs (as Aero Database reports them).", default=d_choice(ROW_AXES)),
        InputSpec("Mach", "real", "Mach number of this condition.", default=d_real(0.5)),
        InputSpec("Angle of Attack", "real", "Angle of attack of this condition.", ANGLE, d_real(0.0, ANGLE)),
        InputSpec("Sideslip", "real", "Sideslip angle of this condition.", ANGLE, d_real(0.0, ANGLE)),
        InputSpec("CA or CD", "real", "Axial force coefficient (Body) or drag coefficient (Wind).", default=zero),
        InputSpec("CY", "real", "Side force coefficient, positive to the right.", default=zero),
        InputSpec("CN or CL", "real", "Normal force coefficient, positive up (Body), or lift coefficient (Wind).",
                  default=zero),
        InputSpec("Cl", "real", "Rolling moment coefficient (over area and span).", default=zero),
        InputSpec("Cm", "real", "Pitching moment coefficient, nose up positive (over area and chord).", default=zero),
        InputSpec("Cn", "real", "Yawing moment coefficient (over area and span).", default=zero),
        InputSpec("Reference Area", "real", "Area the coefficients are normalised by.", AREA, d_real(1.0, AREA)),
        InputSpec("Reference Span", "real", "Length the rolling and yawing moments are normalised by.", LENGTH,
                  d_real(1.0, LENGTH)),
        InputSpec("Reference Chord", "real", "Length the pitching moment is normalised by.", LENGTH,
                  d_real(1.0, LENGTH)),
        InputSpec("Moment Center", "point", "Point the moments are taken about.", LENGTH, d_point_mm(0, 0, 0)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Coefficient row.** One flight condition's six coefficients, with their reference and moment point, as "
              "one line of text: `axes, Mach, alpha, beta, c1..c6, Sref, bref, cref, x, y, z` (degrees and metres). "
              "Give any input a list and nTop makes one row per item, which is how a sweep of another block's "
              "coefficients becomes rows for **Coefficient Rows Contribution (VSP)**.")
    deg = (math.radians(1.0), ANGLE)
    parts = [k.pick(inp("Axes"), ("body", "wind")), k.num(inp("Mach")), k.num(inp("Angle of Attack"), *deg),
             k.num(inp("Sideslip"), *deg)]
    parts += [k.num(inp(n)) for n in ("CA or CD", "CY", "CN or CL", "Cl", "Cm", "Cn")]
    parts += [k.num(inp("Reference Area"), 1.0, AREA), k.num(inp("Reference Span"), 1.0, LENGTH),
              k.num(inp("Reference Chord"), 1.0, LENGTH)]
    parts += [k.num(prop(inp("Moment Center"), c), 1.0, LENGTH) for c in "xyz"]
    out = k.cat(parts, name="Row", sep=",")
    return make_block("aero_coefficient_row", "Aero Coefficient Row (VSP)",
                      "Six force and moment coefficients at one flight condition, with their reference and moment "
                      "point, as a row for **Coefficient Rows Contribution (VSP)**. Give the inputs lists (for example "
                      "the outputs of Modified Newtonian Coefficient blocks swept over angles) to get one row per "
                      "condition.", specs, "text", g, out)


# --- database ----------------------------------------------------------------------------------------------------------

CONDITIONS = ("Those of the first contribution", "Every combination of the lists below")
REFERENCES = ("The first contribution's", "The values below")


def _list_default(values, units=None):
    return lambda g: g.block("core.list<real>", "list<real>", [g.real(v, units) for v in values], top=False)


def database_block():
    own = (
        InputSpec("Contributions", "list<text>", "The contributions to sum: outputs of VSPAERO, Wave Drag, Parasite "
                  "Drag and Coefficient Rows Contribution (VSP), in any number and order."),
        InputSpec("Flight Conditions", "choice", "The database's conditions: the cases of the first contribution (a "
                  "VSPAERO run or coefficient rows), or every combination of the Mach, Angle of Attack and Sideslip "
                  "lists below. Each contribution is interpolated onto them.", default=d_choice(CONDITIONS)),
        InputSpec("Mach", "list<real>", "Mach numbers, when Flight Conditions is the lists.", default=_list_default([0.5])),
        InputSpec("Angle of Attack", "list<real>", "Angles of attack, when Flight Conditions is the lists.", ANGLE,
                  _list_default([0.0], ANGLE)),
        InputSpec("Sideslip", "list<real>", "Sideslip angles, when Flight Conditions is the lists.", ANGLE,
                  _list_default([0.0], ANGLE)),
        InputSpec("Reference", "choice", "Reference area, span, chord and moment point of the database's coefficients: "
                  "the first contribution's (a VSPAERO run or coefficient rows), or the values below. Every "
                  "contribution is converted to it.", default=d_choice(REFERENCES)),
        InputSpec("Reference Area", "real", "Sref, when Reference is the values below.", AREA, d_real(1.0, AREA)),
        InputSpec("Reference Span", "real", "bref (rolling and yawing moments).", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Reference Chord", "real", "cref (pitching moment).", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Moment Reference Point", "point", "Moments are about this point.", LENGTH, d_point_mm(0, 0, 0)),
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
    g.comment("**Aero database.** A Python driver sums the contributions at every flight condition. Each contribution "
              "applies only inside its Mach range; it is interpolated onto the condition (multilinear in Mach, alpha "
              "and beta; Wave and Parasite Drag in Mach only), turned into body-axis forces and moments, moved to the "
              "database's moment point (M + r x F) and normalised by the database's reference. The totals give CL, CD, "
              "CY (wind axes, VSPAERO's convention) and Cl, Cm, Cn (aircraft body-axis signs: right wing down, nose "
              "up, nose right).")
    g.comment("**Outputs.** `<Case Name>.ntvsp.database.csv`: design variables, Mach, AoA, Beta, CL, CD, CY, Cl, Cm, "
              "Cn, L_over_D, CA, CN, then each contribution's own six coefficients (columns named after it). Also the "
              "centre of pressure and aerodynamic centre files read by **Center of Pressure (VSP)**. With Append To, "
              "the rows are added to that file as well.")
    deg = (math.radians(1.0), ANGLE)
    contrib = k.cat(["contribution = ", k.join(inp("Contributions"), "\ncontribution = "), "\n"])
    point = inp("Moment Reference Point")
    lines = [
        contrib,
        k.line("conditions_from_first", k.pick(inp("Flight Conditions"), ("1", "0"))),
        k.line("mach", k.join(k.num(inp("Mach")))),
        k.line("alpha", k.join(k.num(inp("Angle of Attack"), *deg))),
        k.line("beta", k.join(k.num(inp("Sideslip"), *deg))),
        k.line("reference_from_first", k.pick(inp("Reference"), ("1", "0"))),
        k.line("sref", k.num(inp("Reference Area"), 1.0, AREA)),
        k.line("bref", k.num(inp("Reference Span"), 1.0, LENGTH)),
        k.line("cref", k.num(inp("Reference Chord"), 1.0, LENGTH)),
    ] + [k.line(f"{c}ref", k.num(prop(point, c), 1.0, LENGTH)) for c in "xyz"] + [
        k.line("append_to", inp("Append To")),
        k.line("variables", k.join(inp("Variable Names"))),
        k.line("values", k.join(k.num(inp("Variable Values")))),
    ]
    out = analysis(g, k, inp, "run_database", lambda ready: [], lines, "Aero Database", needs_openvsp=False)
    return make_block(
        "aero_database", "Aero Database (VSP)",
        "A six-degree-of-freedom aero database from any mix of contributions: VSPAERO, Wave Drag, Parasite Drag, and "
        "coefficients computed in nTop (Modified Newtonian, Viscous Correction, anything else as Aero Coefficient "
        "Rows). Each applies within its own Mach range and is converted to one reference and moment point; the sum "
        "gives CL, CD, CY, Cl, Cm and Cn at every flight condition, with each contribution's share in its own "
        "columns. Outputs the results path for the database readers, Center of Pressure and Polar Plot (VSP); Append "
        "To collects runs with different design variables into one design-of-experiments file.",
        specs, "text", g, out)


DATABASE_COLUMNS = ("CL", "CD", "CY", "Cl", "Cm", "Cn", "L_over_D", "CA", "CN", "Mach", "AoA", "Beta")


def database_coefficient_block():
    specs = (
        InputSpec("Results", "text", "Output of Aero Database (VSP)."),
        InputSpec("Coefficient", "choice", "A total coefficient (CL, CD, CY in wind axes; Cl, Cm, Cn in body-axis "
                  "signs; CA, CN body-axis forces) or a condition (Mach, AoA and Beta in degrees).",
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

PLOT_COLUMNS = ("CL", "CD", "CY", "Cl", "Cm", "Cn", "L_over_D", "CA", "CN", "AoA", "Beta", "Mach")
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
                  "other two conditions (and every result) is its own curve.", default=d_choice(ALONG)),
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
