"""Parasite Drag and Wave Drag custom blocks and their result readers."""
from __future__ import annotations

from ntopfile import Graph, InputSpec

from kit import (F_choice, F_int, F_length, F_real, F_switch, F_text, Kit, analysis, d_choice,
                 d_int, make_block, model_inputs, run_inputs)

# --- Parasite drag -------------------------------------------------------------------------------------------------

BODY_FF = ("Manual", "Schemensky Fuselage", "Schemensky Nacelle", "Hoerner Streamlined Body", "Torenbeek", "Shevell",
           "Covert", "Jenkinson Fuselage", "Jenkinson Wing Nacelle", "Jenkinson Aft Fuselage Nacelle")
WING_FF = ("Manual", "EDET Conventional Airfoil", "EDET Advanced Airfoil", "Hoerner", "Covert", "Shevell", "Kroo",
           "Torenbeek", "DATCOM", "Schemensky 6 Series Airfoil", "Schemensky 4 Series Airfoil", "Jenkinson Wing",
           "Jenkinson Tail", "Schemensky Supercritical Airfoil")
LAM_CF = ("Blasius", "Blasius with Heat Transfer")
TURB_CF = [("Explicit Fit of Spalding", 0), ("Explicit Fit of Spalding-Chi", 1), ("Explicit Fit of Schoenherr", 2),
           ("Implicit Schoenherr", 4), ("Implicit Karman-Schoenherr", 5), ("Power Law Blasius", 6),
           ("Power Law Prandtl Low Re", 7), ("Power Law Prandtl Medium Re", 8), ("Power Law Prandtl High Re", 9),
           ("Schlichting Compressible", 10), ("Schultz-Grunow Estimate of Schoenherr", 14),
           ("Roughness Schlichting Avg", 16), ("Roughness Schlichting Avg Compressible", 19)]
MODEL_UNITS = ("mm", "cm", "m", "in", "ft", "yd", "unitless")

COMPONENT_FIELDS = [
    F_text("Component Name", "name", "Wings", "Name of the component in the OpenVSP model (every component with "
           "this name gets these settings)."),
    F_choice("Body Form Factor Equation", "ff_body", BODY_FF, range(len(BODY_FF)), 3,
             "Form factor equation used when the component is body-like."),
    F_choice("Wing Form Factor Equation", "ff_wing", WING_FF, range(len(WING_FF)), 3,
             "Form factor equation used when the component is wing-like."),
    F_real("Manual Form Factor", "ff_user", 1.0, description="Form factor when the equation is Manual."),
    F_real("Percent Laminar", "percent_laminar", 0.0, description="Percentage of the component's length with laminar flow."),
    F_real("Interference Factor", "q", 1.0, description="Q, multiplies the component's drag."),
    F_length("Roughness Height", "roughness_m", 0.0, "Surface roughness height (for the roughness Cf equations)."),
    F_real("Te/Tw", "te_tw", 1.0, description="Ratio of freestream to wall temperature (heat-transfer Cf equations)."),
    F_real("Taw/Tw", "taw_tw", 1.0, description="Ratio of adiabatic wall to wall temperature (heat-transfer Cf equations)."),
    F_int("Incorporated Generation", "incorporated_gen", 0, "Ancestor generation whose row this component's wetted "
          "area is added to (0: its own row)."),
]


def component_block():
    specs = [f.spec for f in COMPONENT_FIELDS]
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**One component's Parasite Drag settings**, as one line of text: the name, then the form factor "
              "equations, manual form factor, % laminar, Q, roughness (metres), Te/Tw, Taw/Tw and incorporated "
              "generation, separated by `|`. Put several of these in a List for Parasite Drag (VSP)'s Components.")
    parts = []
    for i, f in enumerate(COMPONENT_FIELDS):
        v = inp(f.spec.name)
        if f.kind == "text":
            parts.append(v)
        elif f.kind == "choice":
            parts.append(k.pick(v, f.values))
        elif f.kind == "int":
            parts.append(k.num(v))
        else:
            parts.append(k.num(v, *f.reference))
        if i < len(COMPONENT_FIELDS) - 1:
            parts.append(" | ")
    out = k.cat(["component = "] + parts, name="Component Settings")
    return make_block("parasite_component", "Parasite Drag Component (VSP)",
                      "Per-component settings of OpenVSP's Parasite Drag table (form factor equation, % laminar, "
                      "interference factor Q, roughness, temperature ratios), for the component of the OpenVSP model "
                      "with this name. Components not listed keep the values stored in the model.",
                      specs, "text", g, out)


PARASITE_FIELDS = [
    F_choice("Model Length Unit", "length_unit", MODEL_UNITS, range(len(MODEL_UNITS)), 2,
             "Length unit of the OpenVSP model's dimensions (OpenVSP models carry no units)."),
    F_int("Geometry Set", "geom_set", 0, "OpenVSP set of components to include (0: All, 1: Shown, 2: Not Shown, "
          "3 and up: the model's user sets)."),
    F_choice("Laminar Cf Equation", "lam_cf", LAM_CF, (0, 1), 0, "Laminar skin-friction equation."),
    F_choice("Turbulent Cf Equation", "turb_cf", [n for n, _ in TURB_CF], [v for _, v in TURB_CF], 9,
             "Turbulent skin-friction equation."),
    F_switch("Export Sub-Components", "export_subcomp", False, "Also list each component's sub-surfaces in OpenVSP's "
             "own CSV table (<Case Name>.ParasiteDragBuildUp.csv)."),
]


def parasite_block():
    own = (
        InputSpec("OpenVSP Model", "text", "Output of OpenVSP Model (VSP), or the path of any .vsp3 model. Parasite Drag "
                  "works on OpenVSP's own components (wings, fuselages...): it needs their wetted areas, reference "
                  "lengths and thickness or fineness ratios."),
        InputSpec("Flight Condition", "text", "Output of Flight Condition (VSP): its flight points' Mach numbers and altitudes (one OpenVSP run per flight point), "
                  "and its reference area."),
    ) + tuple(f.spec for f in PARASITE_FIELDS) + (
        InputSpec("Components", "list<text>", "Per-component settings from Parasite Drag Component (VSP) blocks, in a "
                  "List. Leave the default (one empty line) to use the settings stored in the model.",
                  default=lambda g: g.block("core.list<text>", "list<text>", [g.text("")], top=False)),
    )
    specs = own + run_inputs("parasite")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**OpenVSP Parasite Drag.** A component drag buildup: each component's skin friction (laminar and "
              "turbulent Cf equations at its Reynolds number), times its form factor and interference factor, "
              "times its wetted area over Sref, at each flight point of the flight condition (its Mach number at its "
              "altitude in the US Standard Atmosphere 1976, which OpenVSP evaluates itself).")
    lines = k.field_lines(PARASITE_FIELDS, lambda f: inp(f.spec.name)) + [k.line("flight", inp("Flight Condition"))]
    lines.append(k.line("model", inp("OpenVSP Model")))
    lines.append(k.cat([k.join(inp("Components"), "\n"), "\n"]))
    out = analysis(g, k, inp, "run_parasite", lambda base: [], lines, "Parasite Drag")
    return make_block(
        "parasite_drag", "Parasite Drag (VSP)",
        "Runs OpenVSP's **Parasite Drag** tool on an OpenVSP model at every flight point (Mach number and altitude) "
        "of a Flight Condition, with its reference area, and returns the path of its results (connect it to Parasite "
        "Drag Result, Parasite Drag Component Values or Parasite Drag Contribution). The other settings of OpenVSP's "
        "Parasite Drag screen are inputs: geometry set, model length unit, laminar and turbulent Cf equations, "
        "sub-component export, and each component's form factor equation, % laminar, Q, roughness and temperature "
        "ratios. OpenVSP takes the form-factor inputs (reference length, thickness or fineness ratio) from the "
        "model's components.\n\nNot exposed: atmospheres other than the US Standard Atmosphere 1976 (no temperature "
        "offset), excrescences, reference area from a wing.",
        specs, "text", g, out)


PARASITE_RESULTS = ("Total_CD", "Geom_CD", "Total_f", "Geom_f", "Re_per_length", "Mach", "Altitude_m", "Vinf",
                    "Temperature", "Pressure", "Density", "Sref", "Num_Components", "Point")
PARASITE_COMPONENT_COLUMNS = ("CD", "Percent_Total_CD", "f", "Cf", "FF", "Q", "Swet", "Lref", "t_c_or_fineness", "Re",
                              "Percent_Laminar", "Roughness", "Te_Tw", "Taw_Tw", "Surface")


def _reader(key, display, description, file_suffix, labels, out_type, comment, choice_name="Quantity",
            choice_doc="Which value."):
    specs = (
        InputSpec("Results", "text", "Output of the analysis block."),
        InputSpec(choice_name, "choice", choice_doc, default=d_choice(labels)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment(comment)
    path = k.cat([inp("Results"), file_suffix])
    table = g.block("import_table<file_path,bool>", "table", [path, g.boolean(True)], "Import Results Table")
    col = g.block("extract_column<table,text>", "list<text>", [table, k.pick(inp(choice_name), labels)], top=False)
    out = g.block("scalar_from_text<text>", "real", [col], "Values")
    return make_block(key, display, description, specs, out_type, g, out)


def parasite_result_block():
    return _reader("parasite_result", "Parasite Drag Result (VSP)",
                   "A total of a Parasite Drag run, one value per flight point: CD (total and geometry), drag area f, "
                   "Reynolds number per unit length, and the flight condition OpenVSP used (Mach, altitude in metres, "
                   "Vinf, temperature in kelvin, pressure in pascals, density), Sref, components and flight point.",
                   ".ntvsp.parasite.csv", PARASITE_RESULTS, "list<real>",
                   "**Parasite Drag totals.** The driver writes one row per flight point to `<Case Name>.ntvsp.parasite.csv`; "
                   "Import Table reads it and the chosen column becomes a list.")


def parasite_components_block():
    specs = (
        InputSpec("Results", "text", "Output of Parasite Drag (VSP)."),
        InputSpec("Flight Point", "integer", "Which flight point of the flight condition, counted from 1.",
                  default=d_int(1)),
        InputSpec("Column", "choice", "Which column of the table.", default=d_choice(PARASITE_COMPONENT_COLUMNS)),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Parasite Drag table.** One row per component in `<Case Name>.ntvsp.parasite_components.<n>.csv` for "
              "flight point n; the chosen column comes out as a list.")
    path = k.cat([inp("Results"), ".ntvsp.parasite_components.", k.num(inp("Flight Point")), ".csv"])
    table = g.block("import_table<file_path,bool>", "table", [path, g.boolean(True)], "Import Results Table")
    col = g.block("extract_column<table,text>", "list<text>", [table, k.pick(inp("Column"), PARASITE_COMPONENT_COLUMNS)],
                  top=False)
    out = g.block("scalar_from_text<text>", "real", [col], "Values")
    return make_block("parasite_components", "Parasite Drag Component Values (VSP)",
                      "One column of the Parasite Drag table at one flight point, one value per component in "
                      "OpenVSP's order (CD, % of total, f, Cf, form factor, Q, wetted area, reference length, t/c or "
                      "fineness ratio, Reynolds number...).", specs, "list<real>", g, out)


# --- Wave drag --------------------------------------------------------------------------------------------------------

WAVE_FIELDS = [
    F_int("Number of Slices", "num_slices", 30, "Mach-plane cuts along the configuration for each rotation section. The "
          "result keeps changing with this when the area distribution has steps (a blunt base, an inlet face): check "
          "it at 20 and 40."),
    F_int("Number of Rotation Sections", "num_rot_sects", 30, "Roll angles at which the configuration is cut. Wings of "
          "high aspect ratio need many (the result falls towards linear theory as they increase): 30 is converged for "
          "slender aircraft such as the X-59, a straight wing of aspect ratio 12 needs over 120."),
    F_switch("Symmetry", "symmetry", True, "Use the configuration's x-z symmetry to cut only half the rotation sections."),
    F_choice("Geometry Set", "geom_set", ("All", "Thick", "Thin"), ("all", "thick", "thin"), 0,
             "Components of the model to include (by OpenVSP set)."),
]


def wave_block():
    own = (InputSpec("OpenVSP Model", "text", "Output of OpenVSP Model (VSP), or the path of any .vsp3 model."),
        InputSpec("Flight Condition", "text", "Output of Flight Condition (VSP): its flight points' Mach numbers above 1 (wave drag is zero subsonic), "
                  "and its reference area."),) + tuple(f.spec for f in WAVE_FIELDS)
    specs = own + run_inputs("wave")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**OpenVSP Wave Drag** at each supersonic Mach number of the flight condition. OpenVSP cuts the model with Mach planes "
              "(inclined at the Mach angle) at each rotation section, turns each set of cut areas into an equivalent "
              "area distribution, and averages the von Karman slender-body wave drag over the roll angles: the "
              "supersonic area rule, zero-lift linear theory. CD_wave = D / (q Sref).")
    lines = k.field_lines(WAVE_FIELDS, lambda f: inp(f.spec.name)) + [k.line("model", inp("OpenVSP Model")),
                                                                      k.line("flight", inp("Flight Condition"))]
    out = analysis(g, k, inp, "run_wave", lambda base: [], lines, "Wave Drag")
    return make_block(
        "wave_drag", "Wave Drag (VSP)",
        "Runs OpenVSP's **Wave Drag** (supersonic area rule) on an OpenVSP model at the supersonic Mach numbers of a "
        "Flight Condition, with its reference area, and returns the path of its results (connect it to Wave Drag "
        "Result, Wave Drag Area Distribution or Wave Drag Contribution). The other inputs are those of OpenVSP's Wave "
        "Drag screen: number of slices, number of rotation sections, symmetry and geometry set. It is linear theory at zero lift: see the README for the Mach range it holds over. "
        "Not exposed (out of scope): inflow and outflow faces.", specs, "text", g, out)


WAVE_RESULTS = ("CD_wave", "Mach", "Num_Slices", "Num_Rotation_Sections", "Sref")


def wave_result_block():
    return _reader("wave_result", "Wave Drag Result (VSP)",
                   "CD_wave of a Wave Drag run, one value per Mach number of the sweep (or the Mach numbers, slice and "
                   "section counts and Sref it used; Sref in square metres).", ".ntvsp.wave.csv", WAVE_RESULTS,
                   "list<real>", "**Wave Drag result.** One row per Mach number in `<Case Name>.ntvsp.wave.csv`; the "
                   "chosen column becomes a list.")


def wave_areas_block():
    specs = (
        InputSpec("Results", "text", "Output of Wave Drag (VSP)."),
        InputSpec("Mach Number", "integer", "Which Mach number of the sweep, counted from 1.", default=d_int(1)),
        InputSpec("Rotation Section", "integer", "Which roll angle, counted from 1.", default=d_int(1)),
        InputSpec("Value", "choice", "x/L station of each slice, or its equivalent area (square metres).",
                  default=d_choice(("X_Norm", "Area"))),
    )
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Area distribution.** The driver writes each Mach number's and rotation section's equivalent area "
              "distribution to `<Case Name>.ntvsp.wave_areas.<m>.<n>.csv` (roll angle, x/L, area). The chosen column "
              "comes out as a list, one value per slice.")
    path = k.cat([inp("Results"), ".ntvsp.wave_areas.", k.num(inp("Mach Number")), ".", k.num(inp("Rotation Section")),
                  ".csv"])
    table = g.block("import_table<file_path,bool>", "table", [path, g.boolean(True)], "Import Area Table")
    col = g.block("extract_column<table,text>", "list<text>", [table, k.pick(inp("Value"), ("X_Norm", "Area"))], top=False)
    out = g.block("scalar_from_text<text>", "real", [col], "Values")
    return make_block("wave_areas", "Wave Drag Area Distribution (VSP)",
                      "The equivalent area distribution of one rotation section at one Mach number of a Wave Drag run: "
                      "x/L stations or slice areas, as a list.", specs, "list<real>", g, out)
