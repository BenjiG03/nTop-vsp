"""Aero buildup blocks: contributions from any analysis, the Aero Database that sums them, its readers and plots.

A **contribution** tells the database where one source of forces and moments is and how to use it. Every
contribution is made by one block, **Aero Contribution (VSP)**, whose inputs all have defaults; the per-source
blocks (VSPAERO, Wave Drag, Parasite Drag, Coefficients, nTop Fluids Contribution) wrap it with that source's
defaults and only the inputs that matter for it. Its output is a few ``contribution.<key> = <value>`` lines, one
setting per line, so names and paths may hold any character but a line break. **Aero Database (VSP)** takes a
Flight Condition (its conditions and reference) and a list of contributions; it takes each at the database's
conditions (interpolating where needed), moves it to the reference and moment point, keeps it only inside its Mach
range, and sums them (``vsp_ntop.run_database``).
"""
from __future__ import annotations

from ntopfile import Graph, InputSpec

from kit import (FLIGHT_INPUT, Kit, analysis, d_choice, d_real, d_text, make_block, model_inputs, prop, run_inputs)

USE = ("All forces and moments", "Drag only (force along the freestream)", "All but drag")
USE_KEYS = ("all", "drag", "nodrag")
SOURCES = ("VSPAERO", "Wave Drag", "Parasite Drag",
           "Coefficient table (Modified Newtonian, Viscous Correction, nTop Fluids, Aero Coefficients)")
SOURCE_KEYS = ("vspaero", "wave", "parasite", "coefficients")
SOURCE_NAMES = ("VSPAERO", "Wave Drag", "Parasite Drag", "Coefficients")
VSPAERO_TERMS = ("Inviscid (surface pressure: CFi, CMi)", "Total (with VSPAERO's own skin-friction estimate)")

USE_INPUT = InputSpec("Use", "choice", "All: every force and moment. Drag only: just the force component along the "
                      "freestream, no moments (to take only drag from this source). All but drag: everything else "
                      "(to take lift, side force and moments from this source and drag from others).",
                      default=d_choice(USE))
TERMS_INPUT = InputSpec("Terms", "choice", "VSPAERO only. Inviscid: VSPAERO's surface-pressure forces and moments, for "
                        "databases whose friction comes from Parasite Drag or a viscous correction. Total: with "
                        "VSPAERO's own skin-friction estimate (the o columns).", default=d_choice(VSPAERO_TERMS))
RANGE_DOC = ("**Mach range.** The database uses a contribution where its analysis has values and from Mach Min up to "
             "(not including) Mach Max, so one source can hand over to another: for example VSPAERO below Mach 3 and "
             "Modified Newtonian from Mach 3.")


def range_inputs(mach_min=0.0, mach_max=1000.0):
    return (
        InputSpec("Mach Min", "real", "Lowest Mach number at which the database uses this contribution.",
                  default=d_real(mach_min)),
        InputSpec("Mach Max", "real", "The database uses this contribution below this Mach number (not at it), so "
                  "contributions that hand over at one Mach number do not overlap.", default=d_real(mach_max)),
        InputSpec("Scale", "real", "Factor on every force and moment of this contribution (1: as computed).",
                  default=d_real(1.0)),
    )


# --- the common contribution block ------------------------------------------------------------------------------------

def contribution_block():
    specs = (
        InputSpec("Source", "choice", "Which kind of analysis Results comes from.", default=d_choice(SOURCES)),
        InputSpec("Results", "text", "Output of the analysis: VSPAERO Analysis, Wave Drag, Parasite Drag, or a block "
                  "that writes a coefficient table (Modified Newtonian, Viscous Correction, nTop Fluids, Aero "
                  "Coefficients).", default=d_text("")),
        InputSpec("Name", "text", "Name of this contribution in the database's columns. Empty: the source's name.",
                  default=d_text("")),
        TERMS_INPUT,
        USE_INPUT,
    ) + range_inputs()
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Aero contribution.** The one block every contribution goes through: the per-source contribution "
              "blocks call it with their own defaults. It checks that Results is connected and writes the settings "
              "the database reads, one per line: `contribution.kind`, `name`, `results`, `terms`, `use`, `mach_min`, "
              "`mach_max`, `scale`.\n\n" + RANGE_DOC)
    results = k.cat([inp("Results")])
    connected = g.block("not_equal<text,text>", "bool", [results, k.t("")], top=False)
    check = g.block("core.validate<bool,text,message_type_enum>", "text",
                    [connected, k.t("Aero Contribution: connect an analysis's output to Results"),
                     g.enum("message_type_enum", 1)], "Results Connected")
    name = k.choose(g.block("equals<text,text>", "bool", [inp("Name"), k.t("")], top=False),
                    k.pick(inp("Source"), SOURCE_NAMES), inp("Name"))
    lines = k.cat(["contribution.kind = ", k.pick(inp("Source"), SOURCE_KEYS),
                   "\ncontribution.name = ", name,
                   "\ncontribution.results = ", results,
                   "\ncontribution.terms = ", k.pick(inp("Terms"), ("inviscid", "total")),
                   "\ncontribution.use = ", k.pick(inp("Use"), USE_KEYS),
                   "\ncontribution.mach_min = ", k.num(inp("Mach Min")),
                   "\ncontribution.mach_max = ", k.num(inp("Mach Max")),
                   "\ncontribution.scale = ", k.num(inp("Scale"))])
    out = k.wait(check, lines, name="Contribution")
    return make_block("aero_contribution", "Aero Contribution (VSP)",
                      "One source of forces and moments for **Aero Database (VSP)**: which analysis (Source), its "
                      "output (Results), a name for its columns, which forces to use, the Mach range it applies over "
                      "and a scale factor. Every input has a default; the per-source contribution blocks call this "
                      "one.", specs, "text", g, out)


def _wrapper(key, display, description, source, own, comment, base):
    """A contribution block for one source: its own inputs (by name, among Results, Name, Terms, Use, Mach Min,
    Mach Max, Scale), the rest fixed at the common block's defaults, calling the common block."""
    names = [s.name for s in own]
    g = Graph()
    k = Kit(g)
    inp = model_inputs(own)
    g.comment(comment + "\n\n" + RANGE_DOC)
    values = []
    for spec in base.inputs:
        if spec.name == "Source":
            values.append(g.choice(list(SOURCES), SOURCE_KEYS.index(source)))
        elif spec.name in names:
            values.append(inp(spec.name))
        else:
            values.append(spec.default(g))
    out = g.block(base.signature, base.output_type, values, "Contribution")
    return make_block(key, display, description, own, "text", g, out, requires=(base,))


def vspaero_contribution_block(base=None):
    base = base or contribution_block()
    own = (InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP)."),
           InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text("VSPAERO")),
           TERMS_INPUT, USE_INPUT) + range_inputs()
    return _wrapper("contribution_vspaero", "VSPAERO Contribution (VSP)",
                    "A VSPAERO run as a contribution to **Aero Database (VSP)**: its forces and moments over its cases, "
                    "inviscid or with VSPAERO's own friction, within a Mach range.", "vspaero", own,
                    "**VSPAERO contribution**, through Aero Contribution (VSP). The database reads the body-axis force "
                    "and moment coefficients of every case from the run's polar (at each condition the run's Reynolds "
                    "number nearest the condition's), and the reference and moment point from the run's flight "
                    "condition.", base)


def drag_contribution_block(kind, base=None):
    base = base or contribution_block()
    wave = kind == "wave"
    title = "Wave Drag" if wave else "Parasite Drag"
    own = (InputSpec("Results", "text", f"Output of {title} (VSP)."),
           InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text(title))) + \
        range_inputs(1.0 if wave else 0.0)
    return _wrapper(f"contribution_{kind}", f"{title} Contribution (VSP)",
                    f"A {title} run as a contribution to **Aero Database (VSP)**: a drag coefficient per Mach number"
                    + ("" if wave else " (and altitude)") + ", applied along the freestream at every angle, within a "
                    "Mach range.", kind, own,
                    f"**{title} contribution**, through Aero Contribution (VSP). {title} is a drag coefficient per Mach "
                    "number" + ("" if wave else " at each flight point's altitude") + ": the database interpolates it "
                    "linearly in Mach and applies it along the freestream at every angle, through the moment point (no "
                    "moments)." + (" Mach Min defaults to 1: wave drag is zero subsonic." if wave else ""), base)


def coefficients_contribution_block(base=None):
    base = base or contribution_block()
    own = (InputSpec("Results", "text", "Output of Modified Newtonian (VSP), Viscous Correction (VSP), nTop Fluids "
                     "(VSP) or Aero Coefficients (VSP): a coefficient table, one row per flight condition."),
           InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text("Newtonian")),
           USE_INPUT) + range_inputs()
    return _wrapper("contribution_coefficients", "Coefficients Contribution (VSP)",
                    "Coefficients computed in nTop (Modified Newtonian, Viscous Correction, nTop Fluids, Aero "
                    "Coefficients) as a contribution to **Aero Database (VSP)**, within a Mach range.", "coefficients",
                    own, "**Coefficients contribution**, through Aero Contribution (VSP). The database takes the table "
                    "where it has values: over its own Mach range and, for tables by altitude (Viscous Correction, nTop "
                    "Fluids, Aero Coefficients), at its own altitudes; between its Mach numbers and angles it "
                    "interpolates.", base)


# --- database ----------------------------------------------------------------------------------------------------------

def _list_default(values, units=None):
    return lambda g: g.block("core.list<real>", "list<real>", [g.real(v, units) for v in values], top=False)


def database_block():
    own = (
        InputSpec("Flight Condition", "text", "Output of Flight Condition (VSP): the database's conditions (every "
                  "flight point x sideslip x angle of attack) and its reference area, span, chord and moment point. "
                  "Normally the one the analyses ran on."),
        InputSpec("Contributions", "list<text>", "The contributions to sum: outputs of Aero Contribution (VSP) or of "
                  "the per-source contribution blocks (VSPAERO, Wave Drag, Parasite Drag, Coefficients, nTop Fluids), "
                  "in any number and order."),
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
    contrib = k.cat([k.join(inp("Contributions"), "\n"), "\n"])
    lines = [
        contrib,
        k.line("contribution_count", k.num(prop(inp("Contributions"), "size"))),
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


def rasaero_contribution_block(base=None):
    base = base or contribution_block()
    own = (InputSpec("Results", "text", "Output of RASAero II (VSP): its coefficient table, one row per flight condition."),
           InputSpec("Name", "text", "Name of this contribution in the database's columns.", default=d_text("RASAero")),
           USE_INPUT) + range_inputs()
    return _wrapper("contribution_rasaero", "RASAero Contribution (VSP)",
                    "A RASAero II run as a contribution to **Aero Database (VSP)**: its normal force, axial force "
                    "(drag included: friction, wave, base) and the moments of its centre of pressure, within a Mach "
                    "range. Use All but drag to take only its lift and moments, Drag only for its drag.", "coefficients",
                    own, "**RASAero II contribution**, through Aero Contribution (VSP). RASAero II's table is by Mach "
                    "number, altitude, angle of attack and sideslip; the database takes it where it has values.", base)
