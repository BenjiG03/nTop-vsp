"""Surface pressure readers shared by every analysis that exports pressures (``vsp_ntop/pressure.py``)."""
from __future__ import annotations

from kit import PRESSURE, Kit, d_choice, d_int, make_block, model_inputs
from ntopfile import Graph, InputSpec, Ref

EXPORT_DOC = ("Off: no pressure files. On: the static pressure p = p_inf + q Cp (Pa) per flight condition case, "
              "`<Case Name>.ntvsp.pressure.<Case>.csv`, for **Surface Pressure (VSP)**. Each run first removes "
              "the pressure files of an earlier run.")
EXPORT_CHOICES = ("Off", "On")
EXPORT_KEYS = ("0", "1")


def export_input(doc=EXPORT_DOC, choices=EXPORT_CHOICES):
    return InputSpec("Export Pressures", "choice", doc, default=d_choice(choices, 0))


CASE_DOC = ("Flight condition case, counted from 1 in the order of the Flight Condition's conditions (every flight "
            "point x sideslip x angle of attack, angle of attack fastest; the Case column of "
            "`<Flight Condition>.ntvsp.flight.csv`).")
RESULTS_DOC = ("Output of an analysis run with Export Pressures on: VSPAERO Analysis, Modified Newtonian on Mesh or "
               "Tangent Wedge on Mesh (VSP).")


def _import(g, k, path, name):
    scale = g.leaf("real_field", {"expression": "1", "units": {}})
    return g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                   [path, g.unit_length("m"), scale], name, top=name is not None)


def surface_pressure_block():
    specs = (InputSpec("Results", "text", RESULTS_DOC), InputSpec("Case Number", "integer", CASE_DOC, default=d_int(1)))
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Surface pressure** (Pa) of one case, from `<Results>.ntvsp.pressure.<Case>.csv`.")
    path = k.cat([inp("Results"), ".ntvsp.pressure.", k.num(inp("Case Number")), ".csv"], name="Pressure File")
    out = _import(g, k, path, "Surface Pressure")
    return make_block("surface_pressure", "Surface Pressure (VSP)",
                      "Static surface pressure (Pa) of one flight condition case of an analysis run with Export "
                      "Pressures on, as a **Scalar Point Map** for contour plots: VSPAERO (Thick-set components), "
                      "Modified Newtonian on Mesh, Tangent Wedge on Mesh.", specs, "real_point_map", g, out)


def surface_pressure_jump_block():
    specs = (InputSpec("Results", "text", "Output of VSPAERO Analysis (VSP) run with Export Pressures on, with "
                       "Thin-set components."),
             InputSpec("Case Number", "integer", CASE_DOC, default=d_int(1)),
             InputSpec("Reference Pressure", "real", "Required: a vortex-lattice sheet carries only the jump across "
                       "it, Delta p = q (Cp on its normals' side - Cp on the other), so the map shows Reference "
                       "Pressure + Delta p.",
                       PRESSURE))
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Pressure jump** across VSPAERO's thin surfaces for one case, "
              "`<Results>.ntvsp.pressure_jump.<Case>.csv` (Pa), shown as Reference Pressure + Delta p.")
    path = k.cat([inp("Results"), ".ntvsp.pressure_jump.", k.num(inp("Case Number")), ".csv"], name="Jump File")
    jump = _import(g, k, path, "Pressure Jump")
    reference = g.block("divide<real,real>", "real", [inp("Reference Pressure"), g.real(1.0, PRESSURE)], top=False)
    values = g.block("add<real,real>", "real", [Ref(jump, ("values",)), reference], top=False)
    out = g.block("point_map<list<point>,list<real>>", "real_point_map", [Ref(jump, ("points",)), values],
                  "Surface Pressure")
    return make_block("surface_pressure_jump", "Surface Pressure Jump (VSP)",
                      "The pressure jump across VSPAERO's Thin-set (vortex-lattice) surfaces for one flight condition "
                      "case, plus a required Reference Pressure, as a **Scalar Point Map** (Pa) for contour plots.",
                      specs, "real_point_map", g, out)
