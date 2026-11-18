"""The Flight Condition block, and the analyses computed in nTop itself on it: Modified Newtonian, the viscous
(skin-friction) correction and coefficients given directly.

Flight Condition (VSP) runs ``vsp_ntop.run_flight`` (US Standard Atmosphere 1976 through ``ambiance``) and outputs
its results base; every analysis takes it. The nTop-native analyses read its tables with Import Table, evaluate
their graphs once per condition (nTop's list processing over the condition lists), and write
``<flight base>_<Case Name>.ntvsp.coefficients.csv`` with Export Text: one row per condition, body-axis force and
moment coefficients (x aft, y right, z up) with the reference they are normalised by. Their output is that path
base, for Coefficients Contribution (VSP).
"""
from __future__ import annotations

import math

import newtonian
import viscous
from ntopfile import Graph, InputSpec

from kit import (ANGLE, AREA, FLIGHT_INPUT, LENGTH, Kit, analysis, d_choice, d_int, d_point_mm, d_real, d_text,
                 make_block, model_inputs, run_inputs)

DEG = math.radians(1.0)
FEET = 0.3048
GAMMA = 1.4  # the standard atmosphere's air
COEFFICIENT_HEADER = "Mach,AoA,Beta,CA,CY,CN,Cl,Cm,Cn,Sref,Bref,Cref,Xref,Yref,Zref"
COEFFICIENT_HEADER_ALT = "Mach,Altitude,AoA,Beta,CA,CY,CN,Cl,Cm,Cn,Sref,Bref,Cref,Xref,Yref,Zref"
WIND_HEADER_ALT = "Mach,Altitude,AoA,Beta,CD,CY,CL,Cl,Cm,Cn,Sref,Bref,Cref,Xref,Yref,Zref"


def _list_default(values, units=None):
    return lambda g: g.block("core.list<real>", "list<real>", [g.real(v, units) for v in values], top=False)


# --- Flight Condition ----------------------------------------------------------------------------------------------------

def flight_block():
    own = (
        InputSpec("Mach", "list<real>", "Mach numbers. Each is paired with the altitude at the same position in "
                  "Altitude (a flight point); one altitude serves every Mach number.", default=_list_default([0.5])),
        InputSpec("Altitude", "list<real>", "Geometric altitudes of the flight points, in feet or any length unit, "
                  "from -16,400 ft to 265,800 ft (the standard atmosphere's range).", LENGTH,
                  _list_default([0.0], LENGTH)),
        InputSpec("Angle of Attack", "list<real>", "Angles of attack; every flight point is analysed at each.", ANGLE,
                  _list_default([0.0, DEG * 2.0, DEG * 4.0], ANGLE)),
        InputSpec("Sideslip", "list<real>", "Sideslip angles; every flight point is analysed at each.", ANGLE,
                  _list_default([0.0], ANGLE)),
        InputSpec("Reference Area", "real", "Sref: every force and moment coefficient is over it.", AREA,
                  d_real(1.0, AREA)),
        InputSpec("Reference Span", "real", "bref: rolling and yawing moment coefficients are over Sref bref.", LENGTH,
                  d_real(1.0, LENGTH)),
        InputSpec("Reference Chord", "real", "cref: the pitching moment coefficient is over Sref cref; also the "
                  "length of VSPAERO's Reynolds number.", LENGTH, d_real(1.0, LENGTH)),
        InputSpec("Moment Reference Point", "point", "Moments are about this point (in the model's coordinates: x aft, "
                  "y right, z up).", LENGTH, d_point_mm(0, 0, 0)),
    )
    specs = own + run_inputs("flight")
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Flight condition.** The flight points are (Mach, altitude) pairs; the analyses run each at every "
              "sideslip and angle of attack. A Python driver evaluates the US Standard Atmosphere 1976 at each "
              "altitude (the `ambiance` package): temperature, pressure, density, speed of sound and viscosity, and "
              "from the Mach number the velocity, dynamic pressure and Reynolds number. The reference area, span, "
              "chord and moment point travel with the conditions to every analysis and to the Aero Database.")
    g.comment("**Files** (in the Working Directory, under the Case Name): `.ntvsp.points.csv` (per flight point), "
              "`.ntvsp.flight.csv` (per condition: point x sideslip x angle of attack, angle of attack fastest), "
              "`.ntvsp.inviscid.csv` (the distinct Mach, angle of attack and sideslip combinations, for analyses that "
              "do not depend on altitude).")
    point = inp("Moment Reference Point")
    lines = [
        k.line("mach", k.join(k.num(inp("Mach")))),
        k.line("altitude", k.join(k.num(inp("Altitude"), 1.0, LENGTH))),
        k.line("alpha", k.join(k.num(inp("Angle of Attack"), DEG, ANGLE))),
        k.line("beta", k.join(k.num(inp("Sideslip"), DEG, ANGLE))),
        k.line("sref", k.num(inp("Reference Area"), 1.0, AREA)),
        k.line("bref", k.num(inp("Reference Span"), 1.0, LENGTH)),
        k.line("cref", k.num(inp("Reference Chord"), 1.0, LENGTH)),
    ] + [k.line(f"{c}ref", k.num(_prop(point, c), 1.0, LENGTH)) for c in "xyz"]
    out = analysis(g, k, inp, "run_flight", lambda ready: [], lines, "Flight Condition", needs_openvsp=False)
    return make_block(
        "flight_condition", "Flight Condition (VSP)",
        "The flight conditions and reference shared by every analysis and by the Aero Database: Mach numbers paired "
        "with altitudes (the US Standard Atmosphere 1976 gives temperature, pressure, density, speed and Reynolds "
        "number), each at every angle of attack and sideslip, with the reference area, span, chord and moment "
        "point. Connect its output to the analyses (VSPAERO, Wave Drag, Parasite Drag, Modified Newtonian, Viscous "
        "Correction) and to Aero Database (VSP); each uses what it needs.", specs, "text", g, out)


def _prop(ref, c):
    from kit import prop
    return prop(ref, c)


# --- reading a flight condition in nTop ------------------------------------------------------------------------------------

class Conditions:
    """A flight condition's columns as lists, kept to the conditions with Mach Min <= Mach < Mach Max."""

    def __init__(self, g, k, flight, which, columns, mach_min, mach_max):
        self.table = k.table(k.cat([flight, f".ntvsp.{which}.csv"]), name="Flight Conditions")
        mach = k.column(self.table, "Mach")
        low = {c: k.keep(k.column(self.table, c), mach, ">=", mach_min) for c in columns}
        self.cols = {c: k.keep(v, low["Mach"], "<", mach_max) for c, v in low.items()}
        self.k = k
        self.g = g

    def __getitem__(self, c):
        return self.cols[c]

    def angle(self, c):
        return self.g.block("multiply<real,real>", "real", [self.cols[c], self.g.real(DEG, ANGLE)], top=False)

    def ref(self):
        k = self.k
        first = {c: k.first(k.column(self.table, c)) for c in ("Sref", "Bref", "Cref", "Xref", "Yref", "Zref")}
        area = k.first(k.column(self.table, "Sref"), AREA)
        chord = k.first(k.column(self.table, "Cref"), LENGTH)
        centre = self.g.block("point<real,real,real>", "point",
                              [k.first(k.column(self.table, c), LENGTH) for c in ("Xref", "Yref", "Zref")], "Moment Center")
        return first, area, chord, centre


def write_coefficients(g, k, flight, case, cond, coefficients, ref_first, with_altitude, header=None):
    """Export the rows; returns the results base (available once the file is written)."""
    cols = [cond["Mach"]] + ([cond["Altitude"]] if with_altitude else []) + [cond["AoA"], cond["Beta"]]
    parts = [k.num(c) for c in cols] + [k.num(coefficients[n]) for n in ("CA", "CY", "CN", "Cl", "Cm", "Cn")]
    # Newtonian normalises every moment by one reference length, the chord: say so in Bref and Cref.
    parts += [k.num(ref_first[c]) for c in ("Sref", "Cref", "Cref", "Xref", "Yref", "Zref")]
    rows = k.ecat(parts, sep=",")
    text = k.cat([header or (COEFFICIENT_HEADER_ALT if with_altitude else COEFFICIENT_HEADER), "\n",
                  k.join(rows, "\n"), "\n"], name="Coefficient Table")
    base = k.cat([flight, "_", case], name="Results Base")
    written = g.block("export_text<file_path,text>", "file_path", [k.cat([base, ".ntvsp.coefficients.csv"]), text],
                      "Write Coefficients")
    return k.wait(written, base, name="Results")


RANGE_DOC = ("Only the flight conditions with Mach Min <= Mach < Mach Max are computed (the others are left out of "
             "the table, and an Aero Database takes this analysis only where it has values).")


def _range_inputs(mach_min):
    return (InputSpec("Mach Min", "real", "Lowest Mach number computed.", default=d_real(mach_min)),
            InputSpec("Mach Max", "real", "Conditions at or above this Mach number are not computed.",
                      default=d_real(1000.0)))


def _axes_inputs():
    return (InputSpec("Aft Axis", "vector", "Body axis from nose to tail (the model's +x). At zero angles the flow "
                      "runs along it.", default=lambda g: g.vector(1, 0, 0)),
            InputSpec("Up Axis", "vector", "Approximately up in the body frame (the model's +z); sets the pitch plane.",
                      default=lambda g: g.vector(0, 0, 1)))


STAGNATION_INPUT = InputSpec("Stagnation Pressure", "choice", "Modified: Cp_max at the stagnation point, behind a "
                             "normal shock above Mach 1 (Rayleigh pitot) and isentropic below. Classical: Cp_max = 2. "
                             "User Defined: the Cp Max Field.", default=d_choice(newtonian.CP_MAX_MODES))
CP_MAX_FIELD_INPUT = InputSpec("Cp Max Field", "real_field", newtonian.CP_MAX_FIELD_DOC,
                               default=newtonian.default_cp_max_field)
TOLERANCE_INPUT = InputSpec("Integration Tolerance", "real", "Relative error of nTop's Surface Integral.",
                            default=d_real(1e-3))


def _newtonian_values(g, k, inp, cond, ref, extra=None):
    _, area, chord, centre = ref
    values = {"Body": inp("Body"), "Mach": cond["Mach"], "Angle of Attack": cond.angle("AoA"),
              "Sideslip": cond.angle("Beta"), "Gamma": g.real(GAMMA), "Aft Axis": inp("Aft Axis"),
              "Up Axis": inp("Up Axis"), "Reference Area": area, "Reference Length": chord, "Moment Center": centre,
              "Stagnation Pressure": inp("Stagnation Pressure"), "Cp Max Field": inp("Cp Max Field"),
              "Integration Tolerance": inp("Integration Tolerance")}
    values.update(extra or {})
    return values


# --- Modified Newtonian ---------------------------------------------------------------------------------------------------

def newtonian_block():
    specs = (InputSpec("Body", "implicit", "The body (an implicit). Its gradient gives the surface normal, so any "
                       "implicit works."), FLIGHT_INPUT) + _axes_inputs() + (STAGNATION_INPUT, CP_MAX_FIELD_INPUT, TOLERANCE_INPUT) + \
        _range_inputs(1.0) + (InputSpec("Case Name", "text", "Its results are `<Flight Condition>_<Case Name>`.",
                                        default=d_text("newtonian")),)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("# Modified Newtonian on a flight condition\nThe flight condition's distinct Mach numbers, angles of "
              "attack and sideslips (altitude does not matter here) within Mach Min to Mach Max are read as lists; "
              "every block below then runs once per condition (nTop's list processing). The reference area, the "
              "chord (one reference length for every moment) and the moment point come with the flight condition.")
    cond = Conditions(g, k, inp("Flight Condition"), "inviscid", ("Mach", "AoA", "Beta"), inp("Mach Min"), inp("Mach Max"))
    ref = cond.ref()
    sol = newtonian.build(g, _newtonian_values(g, k, inp, cond, ref).__getitem__)
    g.comment("**6. Table.** One row per condition: Mach, AoA, Beta (degrees), CA, CY, CN, Cl, Cm, Cn and the "
              "reference, written to `<Flight Condition>_<Case Name>.ntvsp.coefficients.csv` for Coefficients "
              "Contribution (VSP).")
    out = write_coefficients(g, k, inp("Flight Condition"), inp("Case Name"), cond, sol.coefficients, ref[0], False)
    return make_block(
        "modified_newtonian", "Modified Newtonian (VSP)",
        "Modified Newtonian (Lees) pressure forces and moments of an implicit body at every condition of a Flight "
        "Condition, computed entirely with nTop blocks. " + RANGE_DOC + " Outputs the results path for Coefficients "
        "Contribution (VSP).\n" + newtonian.METHOD_DOC, specs, "text", g, out)


def newtonian_pressure_block():
    specs = (InputSpec("Body", "implicit", "The body (an implicit)."), FLIGHT_INPUT,
             InputSpec("Condition Number", "integer", "Row of the flight condition's distinct (Mach, angle of attack, "
                       "sideslip) combinations, counted from 1 (`.ntvsp.inviscid.csv`).", default=d_int(1))) + \
        _axes_inputs() + (STAGNATION_INPUT, CP_MAX_FIELD_INPUT)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Modified Newtonian pressure** of one condition, as a scalar field to show on the body.")
    table = k.table(k.cat([inp("Flight Condition"), ".ntvsp.inviscid.csv"]), name="Flight Conditions")
    index = g.block("add<integer,integer>", "integer", [inp("Condition Number"), g.integer(-1)], top=False)

    def at(c, units=None):
        v = g.block("core.list_element<list_interface,integer>", "real", [k.column(table, c), index], top=False)
        return g.block("multiply<real,real>", "real", [v, g.real(1.0, units)], top=False) if units else v

    zero = g.real(0.0)
    values = {"Body": inp("Body"), "Mach": at("Mach"), "Angle of Attack": g.block("multiply<real,real>", "real",
                                                                                  [at("AoA"), g.real(DEG, ANGLE)], top=False),
              "Sideslip": g.block("multiply<real,real>", "real", [at("Beta"), g.real(DEG, ANGLE)], top=False),
              "Gamma": g.real(GAMMA), "Aft Axis": inp("Aft Axis"), "Up Axis": inp("Up Axis"),
              "Reference Area": g.real(1.0, AREA), "Reference Length": g.real(1.0, LENGTH),
              "Moment Center": g.point(0, 0, 0), "Stagnation Pressure": inp("Stagnation Pressure"),
              "Cp Max Field": inp("Cp Max Field"),
              "Integration Tolerance": zero}
    sol = newtonian.build(g, values.__getitem__, with_moments=False, integrate_forces=False)
    out = g.block("core.var<real_field>", "real_field", [sol.cp], "Cp")
    return make_block("modified_newtonian_pressure", "Modified Newtonian Pressure (VSP)",
                      "The Modified Newtonian pressure coefficient of one flight condition as a scalar field, to show on "
                      "the body (Field Viewer) or evaluate at points.", specs, "real_field", g, out)


# --- Modified Newtonian on a mesh ----------------------------------------------------------------------------------------

MESH_CP_MAX_KEYS = ("modified", "classical", "user")


TW_EXPANSION = ("Prandtl-Meyer", "Newtonian Shadow")
TW_EXPANSION_KEYS = ("prandtl_meyer", "shadow")
TW_SEPARATED = ("Freestream (Cp 0)", "Gabeaud Base Pressure (-1/M^2)")
TW_SEPARATED_KEYS = ("freestream", "gabeaud")


def newtonian_mesh_block():
    """Modified Newtonian of an implicit body with the force integrals done over nTop's mesh of it (a Python driver
    sums them), instead of nTop's Surface Integral: the same pressure, much faster."""
    return _mesh_block("newtonian")


def tangent_wedge_mesh_block():
    """Tangent wedge (with Prandtl-Meyer or Newtonian-shadow leeward faces) on nTop's mesh of an implicit body."""
    return _mesh_block("tangent_wedge")


def _mesh_block(method):
    wedge = method == "tangent_wedge"
    from kit import PYTHON_VARIABLE, environment
    from ops import Ops
    from ntopfile import NONE, Ref

    specs = (InputSpec("Body", "implicit", "The body (an implicit). Its gradient gives the surface normal of the "
                       "pressure; nTop meshes it for the integrals."), FLIGHT_INPUT,
             InputSpec("Mesh Tolerance", "real", "Largest distance between nTop's mesh of the Body (Mesh from Implicit "
                       "Body, sharp edges kept) and the body. Smaller is more accurate and slower; a few thousandths of "
                       "the body's length is a good start.", LENGTH, d_real(0.0005, LENGTH))) + _axes_inputs() + (
        STAGNATION_INPUT, CP_MAX_FIELD_INPUT)
    if wedge:
        specs += (
            InputSpec("Expansion", "choice", "Leeward faces (turned away from the flow). Prandtl-Meyer: the isentropic "
                      "expansion from the freestream through the face's angle, limited at vacuum. Newtonian Shadow: "
                      "Cp 0.", default=d_choice(TW_EXPANSION, 0)),
            InputSpec("Separated Faces", "choice", "With Prandtl-Meyer, faces turned more than 30 degrees away blend "
                      "(by 60 degrees) to a separated-flow pressure. Freestream: Cp 0, so a flat base adds no base "
                      "drag. Gabeaud Base Pressure: Cp = -1/M^2 (base drag).", default=d_choice(TW_SEPARATED, 0)))
    specs += _range_inputs(1.0 if not wedge else 1.0001) + run_inputs("tangent_wedge_mesh" if wedge else "newtonian_mesh")
    g = Graph()
    k = Kit(g)
    o = Ops(g)
    inp = model_inputs(specs)
    g.comment("# Modified Newtonian on a mesh\nThe pressure is the implicit body's, as in Modified Newtonian (VSP): "
              "Cp = Cp_max max(-n.V, 0)^2 with n the unit gradient of the Body. The force and moment integrals are sums "
              "over nTop's mesh of the Body, -(1/Sref) sum Cp S_f and -(1/(Sref cref)) sum Cp (r_f - c) x S_f over its "
              "faces (S_f the face's outward vector area, r_f its centroid), done by a Python driver at every "
              "(Mach, alpha, beta) of the flight condition within Mach Min to Mach Max: far faster than nTop's Surface "
              "Integral.")

    def exports(ready):
        g.comment("**1. Mesh.** Mesh from Implicit Body at the Mesh Tolerance (sharp edges kept), exported as an STL "
                  "in millimetres. The driver lists its faces' centroids and areas as a point map.")
        mesh = g.block("implicit_to_mesh<implicit,real,real,bool,bool>[2.5.0]", "mesh",
                       [inp("Body"), inp("Mesh Tolerance"), NONE, g.boolean(True), g.boolean(False)], "Body Mesh")
        stl = k.cat([ready, ".ntvsp.mesh.stl"], name="Mesh File")
        exported = g.block("export_mesh<file_path,mesh,unit_length_enum>", "mesh_file_data",
                           [stl, mesh, g.unit_length("mm")], "Export Mesh")
        python = environment(g, PYTHON_VARIABLE, "Python (Faces)")
        temp = environment(g, "TEMP", "Temporary Folder (Faces)")
        faces_csv = k.cat([ready, ".ntvsp.mesh_faces.csv"])
        command = k.wait(exported, k.cat(['"', python, '" -m vsp_ntop.run_newtonian_mesh faces "', stl, '" "',
                                          faces_csv, '"']))
        rc = k.run_command(command, temp, inp("Timeout"), name="List Mesh Faces")
        ok = k.check(rc, "Listing the mesh faces", name="Check Faces")
        g.comment("**2. The body at the faces.** nTop evaluates the Body's gradient (its normal) and the Cp Max Field "
                  "at every face centroid and exports them as point maps (`<Case Name>.ntvsp.mesh_gradient.csv`, "
                  "`<Case Name>.ntvsp.mesh_cp_max.csv`), in the faces' order.")
        faces = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                        [k.wait(ok, faces_csv), g.unit_length("mm"), o.field("1")], "Mesh Faces")
        points = Ref(faces, ("points",))
        gradient = g.block("evaluate_field_gradient<real_field,list<point>>", "list<vector>", [inp("Body"), points],
                           "Body Gradient at Faces")
        cp_max = g.block("evaluate_field<real_field,list<point>>", "list<real>", [inp("Cp Max Field"), points],
                         "Cp Max at Faces")
        gmap = g.block("point_map<list<point>,list<vector>>", "vector_point_map", [points, gradient], top=False)
        cmap = g.block("point_map<list<point>,list<real>>", "real_point_map", [points, cp_max], top=False)
        written = []
        for m, suffix, label in ((gmap, "gradient", "Gradient"), (cmap, "cp_max", "Cp Max")):
            kind = "vector_point_map" if m is gmap else "real_point_map"
            written.append(g.block(f"export_point_map<file_path,{kind},unit_length_enum>", "point_map_file",
                                   [k.cat([ready, f".ntvsp.mesh_{suffix}.csv"]), m, g.unit_length("mm")],
                                   f"Export {label}"))
        return written

    aft, up = inp("Aft Axis"), inp("Up Axis")
    lines = [k.line("flight", inp("Flight Condition")),
             k.line("mach_min", k.num(inp("Mach Min"))), k.line("mach_max", k.num(inp("Mach Max"))),
             k.line("stagnation", k.pick(inp("Stagnation Pressure"), MESH_CP_MAX_KEYS)),
             k.line("method", k.t(method))]
    if wedge:
        lines += [k.line("expansion", k.pick(inp("Expansion"), TW_EXPANSION_KEYS)),
                  k.line("separated", k.pick(inp("Separated Faces"), TW_SEPARATED_KEYS))]
    lines += [k.line(f"aft_{c}", k.num(_prop(aft, c))) for c in "xyz"]
    lines += [k.line(f"up_{c}", k.num(_prop(up, c))) for c in "xyz"]
    g.comment("**3. Sums.** The driver `vsp_ntop.run_newtonian_mesh` computes Cp on every face at every condition and "
              "sums the forces and moments; it writes `<Case Name>.ntvsp.coefficients.csv` (Mach, AoA, Beta, CA, CY, "
              "CN, Cl, Cm, Cn and the reference), the layout of Modified Newtonian (VSP).")
    title = "Tangent Wedge on Mesh" if wedge else "Modified Newtonian on Mesh"
    out = analysis(g, k, inp, "run_newtonian_mesh", exports, lines, title, needs_openvsp=False)
    if wedge:
        return make_block(
            "tangent_wedge_mesh", "Tangent Wedge on Mesh (VSP)",
            "Local-inclination pressure forces and moments of an implicit body at every supersonic condition of a "
            "Flight Condition, with the integrals over nTop's mesh of it (as Modified Newtonian on Mesh (VSP)). Each "
            "face's pressure comes from its impact angle delta (sin delta = -n.V, n the body's unit gradient): "
            "windward, the tangent wedge, the exact two-dimensional oblique shock that turns the flow through delta "
            "(weak solution), Cp = 2/(gamma M^2) (p2/p1 - 1); above the shock-detachment angle Cp rises with sin^2 "
            "delta from the detachment value to Cp_max at 90 degrees (Stagnation Pressure, or the Cp Max Field). "
            "Leeward, the Expansion choice. In the linear limit the windward and leeward faces of a thin surface give "
            "its 4 alpha / beta; unlike Newtonian it holds at moderate hypersonic similarity K = M sin(delta) (thin "
            "fins, wedges, flat plates), and tends to Newtonian as K grows. Faces are independent (no shock-expansion "
            "marching, no interference); gamma 1.4. Needs Mach > 1 (Mach Min). " + RANGE_DOC + " Outputs the "
            "results path for Coefficients Contribution (VSP).", specs, "text", g, out)
    return make_block(
        "modified_newtonian_mesh", "Modified Newtonian on Mesh (VSP)",
        "Modified Newtonian (Lees) pressure forces and moments of an implicit body at every condition of a Flight "
        "Condition, like Modified Newtonian (VSP), with the integrals done over nTop's mesh of the body instead of "
        "Surface Integrals: the pressure comes from the implicit (its gradient normal and the Cp_max choice, "
        "including a User Defined Cp Max Field), the areas, centroids and force directions from the mesh. Much "
        "faster on detailed bodies; its accuracy is the mesh's (Mesh Tolerance). " + RANGE_DOC + " Outputs the "
        "results path for Coefficients Contribution (VSP).\n" + newtonian.METHOD_DOC, specs, "text", g, out)


# --- Viscous correction -----------------------------------------------------------------------------------------------------

PRESSURE_SOURCES = ("Modified Newtonian (computed here)", "Freestream (Cp = 0: flat-plate friction)")


def viscous_block():
    specs = (InputSpec("Body", "implicit", "The body (an implicit)."), FLIGHT_INPUT,
             InputSpec("Inviscid Pressure", "choice", "Edge state of the boundary layer: from the Modified Newtonian "
                       "pressure of each condition (hypersonic practice), or the freestream (Cp = 0), which gives "
                       "flat-plate skin friction at the local running length (subsonic and low supersonic practice).",
                       default=d_choice(PRESSURE_SOURCES))) + _axes_inputs() + (STAGNATION_INPUT, CP_MAX_FIELD_INPUT) + tuple(
        s for s in viscous.INPUTS if s.name in ("Boundary Layer", "Transition Reynolds Number", "Wall",
                                                 "Wall Temperature")) + (TOLERANCE_INPUT,) + _range_inputs(0.0) + (
        InputSpec("Case Name", "text", "Its results are `<Flight Condition>_<Case Name>`.", default=d_text("viscous")),)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("# Viscous correction on a flight condition\nEvery condition of the flight condition within Mach Min to "
              "Mach Max (Mach, altitude, angle of attack, sideslip) is read as lists, with its standard-atmosphere "
              "pressure and temperature; every block below runs once per condition. The edge state's total pressure "
              "is behind a normal shock above Mach 1 and isentropic below it, condition by condition.")
    cond = Conditions(g, k, inp("Flight Condition"), "flight", ("Mach", "Altitude", "AoA", "Beta", "Pressure",
                                                                  "Temperature"), inp("Mach Min"), inp("Mach Max"))
    ref = cond.ref()
    values = _newtonian_values(g, k, inp, cond, ref)
    sol = newtonian.build(g, values.__getitem__, with_moments=False, integrate_forces=False)
    use = g.block("core.select_by_choice<choice,list_interface>", "real",
                  [inp("Inviscid Pressure"), g.block("core.list<real>", "list<real>", [g.real(1.0), g.real(0.0)],
                                                     top=False)], top=False)
    cp = g.block("multiply<real_field,real_field>", "real_field", [sol.cp, use], "Inviscid Cp")
    values.update({
        "Pressure Coefficient": cp,
        "Freestream Pressure": g.block("multiply<real,real>", "real", [cond["Pressure"], g.real(1.0, viscous.PRESSURE)],
                                       top=False),
        "Freestream Temperature": g.block("multiply<real,real>", "real", [cond["Temperature"],
                                                                         g.real(1.0, viscous.TEMPERATURE)], top=False),
        "Boundary Layer": inp("Boundary Layer"), "Transition Reynolds Number": inp("Transition Reynolds Number"),
        "Wall": inp("Wall"), "Wall Temperature": inp("Wall Temperature")})
    _, coefficients, _ = viscous.build(g, values.__getitem__)
    g.comment("**9. Table.** One row per condition (with its altitude): Mach, Altitude (m), AoA, Beta, CA ... Cn and "
              "the reference, written to `<Flight Condition>_<Case Name>.ntvsp.coefficients.csv` for Coefficients "
              "Contribution (VSP).")
    out = write_coefficients(g, k, inp("Flight Condition"), inp("Case Name"), cond, coefficients, ref[0], True)
    return make_block(
        "viscous_correction", "Viscous Correction (VSP)",
        "Skin-friction forces and moments of an implicit body at every condition of a Flight Condition (standard "
        "atmosphere pressure and temperature at each altitude), computed entirely with nTop blocks: the friction "
        "part of an aero database, in place of Parasite Drag. " + RANGE_DOC + " Outputs the results path for "
        "Coefficients Contribution (VSP).\n" + viscous.VISCOUS_DOC, specs, "text", g, out)


# --- Coefficients given directly ----------------------------------------------------------------------------------------------

AXES = ("Body: CA, CY, CN; moments about x aft, y right, z up",
        "Wind: CD, CY, CL; moments in aircraft body-axis signs (Cl right wing down, Cm nose up, Cn nose right)")


def coefficients_block():
    names = ("CA or CD", "CY", "CN or CL", "Cl", "Cm", "Cn")
    specs = (FLIGHT_INPUT,
             InputSpec("Axes", "choice", "How the six coefficients are given.", default=d_choice(AXES))) + tuple(
        InputSpec(n, "list<real>", "One value per condition of the flight condition, in its order (flight point, "
                  "then sideslip, then angle of attack), over its reference area, span (Cl, Cn) and chord (Cm), about "
                  "its moment point.", default=_list_default([0.0])) for n in names) + (
        InputSpec("Case Name", "text", "Its results are `<Flight Condition>_<Case Name>`.", default=d_text("coefficients")),)
    g = Graph()
    k = Kit(g)
    inp = model_inputs(specs)
    g.comment("**Coefficients from anywhere.** Six lists, one value per condition of the flight condition (test data, "
              "another tool, other nTop blocks), written as a coefficient table for Coefficients Contribution (VSP). "
              "Each list must have as many values as the flight condition has conditions.")
    cond = Conditions(g, k, inp("Flight Condition"), "flight", ("Mach", "Altitude", "AoA", "Beta"), g.real(-1.0),
                      g.real(1e9))
    first, _, _, _ = cond.ref()
    coefficients = dict(zip(("CA", "CY", "CN", "Cl", "Cm", "Cn"), (inp(n) for n in names)))
    header = k.pick(inp("Axes"), (COEFFICIENT_HEADER_ALT, WIND_HEADER_ALT))
    first = dict(first, Cref=first["Cref"])
    cols = [cond["Mach"], cond["Altitude"], cond["AoA"], cond["Beta"]]
    parts = [k.num(c) for c in cols] + [k.num(coefficients[n]) for n in ("CA", "CY", "CN", "Cl", "Cm", "Cn")]
    parts += [k.num(first[c]) for c in ("Sref", "Bref", "Cref", "Xref", "Yref", "Zref")]
    rows = k.ecat(parts, sep=",")
    text = k.cat([header, "\n", k.join(rows, "\n"), "\n"], name="Coefficient Table")
    base = k.cat([inp("Flight Condition"), "_", inp("Case Name")], name="Results Base")
    written = g.block("export_text<file_path,text>", "file_path", [k.cat([base, ".ntvsp.coefficients.csv"]), text],
                      "Write Coefficients")
    out = k.wait(written, base, name="Results")
    return make_block("aero_coefficients", "Aero Coefficients (VSP)",
                      "Force and moment coefficients from any source, one value per condition of a Flight Condition, as "
                      "a coefficient table for Coefficients Contribution (VSP).", specs, "text", g, out)
