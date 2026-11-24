"""Skin friction after the pressure solve, hypanel's ``viscous.correct``, as nTop blocks (from ntop_newtonian, now
part of nTop_vsp; lists of flight conditions as in ``newtonian``)."""
from __future__ import annotations

import math
from dataclasses import replace

import newtonian
from ntopfile import Graph, InputSpec, Ref
from ops import ANGLE, AREA, LENGTH, More

GAS_CONSTANT = 287.05
PRANDTL = 0.72
PRESSURE = {"length": -1, "mass": 1, "time": -2}
TEMPERATURE = {"temp": 1}
SEPARATION_START, SEPARATION_END = math.radians(30.0), math.radians(60.0)
BOUNDARY_LAYERS = ("Fully Turbulent", "Laminar", "Laminar to Turbulent")
WALLS = ("Adiabatic", "Fixed Temperature")
COEFFICIENTS = ("CA", "CY", "CN", "CD", "CL", "Cl", "Cm", "Cn")

_FLOW = {s.name: s for s in newtonian.INPUTS}
INPUTS = (
    InputSpec("Body", "implicit", "$markdown-enabled$The implicit body, the same one the pressure was solved on."),
    InputSpec("Pressure Coefficient", "real_field", "Cp of the inviscid solve, e.g. the Modified Newtonian Pressure block's output."),
    _FLOW["Mach"], _FLOW["Angle of Attack"], _FLOW["Sideslip"], _FLOW["Gamma"], _FLOW["Aft Axis"], _FLOW["Up Axis"],
    _FLOW["Reference Area"], _FLOW["Reference Length"], _FLOW["Moment Center"],
    InputSpec("Freestream Pressure", "real", "Static pressure of the freestream.", PRESSURE, lambda g: g.real(101325.0, PRESSURE)),
    InputSpec("Freestream Temperature", "real", "Static temperature of the freestream.", TEMPERATURE, lambda g: g.real(288.15, TEMPERATURE)),
    InputSpec("Boundary Layer", "choice", "Laminar to Turbulent switches at Transition Reynolds Number.",
              default=lambda g: g.choice(list(BOUNDARY_LAYERS))),
    InputSpec("Transition Reynolds Number", "real", "Running-length Reynolds number where the boundary layer turns turbulent.",
              default=lambda g: g.real(5e5)),
    InputSpec("Wall", "choice", "Adiabatic wall, or a fixed Wall Temperature.", default=lambda g: g.choice(list(WALLS))),
    InputSpec("Wall Temperature", "real", "Used when Wall is Fixed Temperature.", TEMPERATURE, lambda g: g.real(300.0, TEMPERATURE)),
    _FLOW["Integration Tolerance"],
)


def sutherland(o: More, t):
    """Air viscosity (Pa s) of a temperature field in kelvin."""
    return o.fdiv(o.fmul(o.fpow(t, o.field("1.5")), o.field("1.458e-6")), o.fadd(t, o.field("110.4")))


def build(g: Graph, inp, faces=None):
    """(skin friction coefficient field on q_inf, coefficient dict) for inputs ``inp(name)``."""
    o = More(g)
    g.comment("**1.")
    x, y, z, v, (ca, sa, cb, sb) = newtonian.freestream(o, inp)
    body, mach, gamma = inp("Body"), inp("Mach"), inp("Gamma")
    grad = o.grad(body)
    normal = o._b("divide<vector_field,real_field>", "vector_field", [grad, o.fsqrt(o.fdot(grad, grad))])
    v_f = o.vconst(v)
    sine = o.fmul(o.fdot(normal, v_f), o.field("-1"))  # sin(delta), positive windward
    delta = o.fasin(o.fclamp(sine, o.field("-1"), o.field("1")))

    g.comment("**2.")
    one = o.real(1.0)
    half_gm1 = o.mul(o.sub(gamma, one), o.real(0.5))
    p_inf = o.div(inp("Freestream Pressure"), o.real(1.0, PRESSURE))
    t_inf = o.div(inp("Freestream Temperature"), o.real(1.0, TEMPERATURE))
    m2 = o.sq(mach)
    w, m_sup, m_sub = o.supersonic(mach)
    pitot_sub = o.pow(o.add(one, o.mul(half_gm1, o.sq(m_sub))), o.div(gamma, o.sub(gamma, one)))
    pitot_sup = o.add(one, o.mul(o.mul(o.mul(gamma, o.sq(m_sup)), o.real(0.5)), newtonian.cp_max_modified(o, m_sup, gamma)))
    pitot = o.blend(w, pitot_sup, pitot_sub, name="Pitot Ratio")
    total_t = o.mul(t_inf, o.add(one, o.mul(half_gm1, m2)), name="Total Temperature")

    g.comment("**3.")
    p_ratio = o.fadd(o.field("1"), o.fmul(inp("Pressure Coefficient"), o.mul(o.mul(gamma, m2), o.real(0.5))))
    q = o.fclamp(o.fdiv(p_ratio, pitot), o.field("1e-12"), o.field("1"))
    exponent = o.neg(o.div(o.sub(gamma, one), gamma))
    edge_m2 = o.fmax(o.fdiv(o.fsub(o.fpow(q, exponent), o.field("1")), half_gm1), o.field("0"))
    energy = o.fmul(edge_m2, half_gm1)  # (gamma - 1) / 2 Me^2
    t_e = o.fdiv(o.fmul(o.field("1"), total_t), o.fadd(o.field("1"), energy), name="Edge Temperature")
    rho = o.fdiv(o.fmul(p_ratio, p_inf), o.fmul(t_e, o.real(GAS_CONSTANT)))
    speed2 = o.fmul(o.fmul(edge_m2, t_e), o.mul(gamma, o.real(GAS_CONSTANT)))
    mu_e = sutherland(o, t_e)

    g.comment("**4.")
    far = o.vscale(v, o.real(-100.0, LENGTH))
    front = g.block("closest_point_on_implicit<implicit,point,real>[1.2.0]", "closest_point_properties",
                    [body, o.point(o.comp(far, "x"), o.comp(far, "y"), o.comp(far, "z")), o.real(1e-8, LENGTH)], top=False)
    front_s = o.add(o.add(o.mul(Ref(front, ("point", "x")), o.comp(v, "x")), o.mul(Ref(front, ("point", "y")), o.comp(v, "y"))),
                    o.mul(Ref(front, ("point", "z")), o.comp(v, "z")))
    position = o.vfield(o.field("x", LENGTH), o.field("y", LENGTH), o.field("z", LENGTH))
    downstream = o.fmax(o.fdiv(o.fsub(o.fdot(position, v_f), front_s), o.real(1.0, LENGTH)), o.field("0"))
    cosine = o.fmax(o.fsqrt(o.fmax(o.fsub(o.field("1"), o.fsq(sine)), o.field("0"))), o.field("0.05"))
    running = o.fdiv(downstream, cosine, name="Running Length")  # metres
    reynolds = o.fdiv(o.fmul(o.fmul(rho, o.fsqrt(speed2)), running), mu_e, name="Running Length Reynolds Number")

    fixed = o.select(inp("Wall"), [o.real(0.0), o.real(1.0)])
    t_fixed = o.div(inp("Wall Temperature"), o.real(1.0, TEMPERATURE))

    def wall(recovery):
        adiabatic = o.fmul(t_e, o.fadd(o.field("1"), o.fmul(energy, o.field(repr(recovery)))))
        return adiabatic, o.fadd(o.fmul(adiabatic, o.sub(one, fixed)), o.fmul(o.field("1"), o.mul(t_fixed, fixed)))

    g.comment("**5.")
    t_aw, t_w = wall(PRANDTL ** 0.5)
    t_star = o.fadd(o.fadd(t_e, o.fmul(o.fsub(t_w, t_e), o.field("0.5"))), o.fmul(o.fsub(t_aw, t_e), o.field("0.22")))
    density_ratio = o.fdiv(t_e, t_star)
    reduced = o.fmax(o.fmul(o.fmul(reynolds, density_ratio), o.fdiv(mu_e, sutherland(o, t_star))), o.field("1"))
    cf_lam = o.fmul(density_ratio, o.fdiv(o.field("0.664"), o.fsqrt(reduced)), name="Laminar Cf")

    g.comment("**6.")
    t_aw_t, t_w_t = wall(PRANDTL ** (1.0 / 3.0))
    m = o.fmax(o.fmul(energy, o.field(repr(PRANDTL ** (1.0 / 3.0)))), o.field("1e-6"))
    ratio = o.fdiv(t_w_t, t_e)
    a2 = o.fdiv(m, ratio)
    b = o.fsub(o.fdiv(o.fadd(o.field("1"), m), ratio), o.field("1"))
    root = o.fsqrt(o.fadd(o.fsq(b), o.fmul(a2, o.field("4"))))
    stretch = o.fdiv(o.fsqrt(m), o.fadd(o.fasin(o.fdiv(o.fsub(o.fmul(a2, o.field("2")), b), root)), o.fasin(o.fdiv(b, root))))
    reduced_t = o.fdiv(o.fmul(o.fmul(reynolds, o.fdiv(mu_e, sutherland(o, t_w_t))), o.fsqrt(o.fdiv(t_e, t_w_t))), stretch)
    log_term = o.flog(o.fmul(o.fmax(reduced_t, o.field("1000")), o.field("0.06")))
    cf_turb = o.fdiv(o.fdiv(o.field("0.455"), o.fsq(log_term)), o.fsq(stretch), name="Turbulent Cf")

    transition = o.select(inp("Boundary Layer"), [o.real(0.0), o.real(1e300), inp("Transition Reynolds Number")])
    turbulent = o.fclamp(o.fadd(o.fsub(reynolds, transition), o.field("0.5")), o.field("0"), o.field("1"), name="Turbulent")
    cf_edge = o.fadd(o.fmul(turbulent, cf_turb), o.fmul(o.fsub(o.field("1"), turbulent), cf_lam))

    g.comment("**7.")
    t = o.fclamp(o.fdiv(o.fsub(o.fmul(delta, o.field("-1")), o.field(repr(SEPARATION_START))),
                        o.field(repr(SEPARATION_END - SEPARATION_START))), o.field("0"), o.field("1"))
    attached = o.fsub(o.field("1"), o.fmul(o.fsq(t), o.fsub(o.field("3"), o.fmul(t, o.field("2")))))
    q_inf = o.mul(o.mul(o.mul(gamma, p_inf), m2), o.real(0.5))
    cf = o.fdiv(o.fmul(o.fmul(o.fmul(cf_edge, rho), speed2), o.fmul(attached, o.field("0.5"))), q_inf, name="Skin Friction Coefficient")

    if faces is None:
        g.comment("**8.")
    else:
        g.comment("**8.")
    tangential = o.vfsub(v_f, o.vfscale(o.fdot(v_f, normal), normal))
    size = o.fmax(o.fsqrt(o.fdot(tangential, tangential)), o.field("1e-9"))
    direction = o._b("divide<vector_field,real_field>", "vector_field", [tangential, size])

    if faces is None:
        area, length, tol = inp("Reference Area"), inp("Reference Length"), inp("Integration Tolerance")

        def integral(f):
            return o.surface_integral(body, f, tol)
    else:
        area, length = inp("Reference Area"), inp("Reference Length")
        points, weights = replace(faces, prop=("points",)), replace(faces, prop=("values",))
        unit_area = o.real(1.0, AREA)

        def integral(f):
            at_points = o._b("evaluate_field<real_field,list<point>>", "list<real>", [f, points])
            total = o._b("add<list<real>>", "real", [o._b("multiply<real,real>", "list<real>", [at_points, weights])])
            return o.mul(total, unit_area)
    c = {}
    for label, e in (("CA", x), ("CY", y), ("CN", z)):
        c[label] = o.div(integral(o.fmul(cf, o.fdot(direction, o.vconst(e)))), area, name=label)
    c["CD"] = o.add(o.add(o.mul(o.mul(ca, cb), c["CA"]), o.mul(o.neg(sb), c["CY"])), o.mul(o.mul(sa, cb), c["CN"]), name="CD")
    c["CL"] = o.add(o.mul(o.neg(sa), c["CA"]), o.mul(ca, c["CN"]), name="CL")
    centre = inp("Moment Center")
    arm = o.vfsub(position, o.vfield(o.comp(centre, "x"), o.comp(centre, "y"), o.comp(centre, "z")))
    for label, e in (("Cl", x), ("Cm", y), ("Cn", z)):
        c[label] = o.div(integral(o.fmul(cf, o.fdot(direction, o.fcross(o.vconst(e), arm)))), o.mul(area, length),
                         name=label)
    return cf, c, {"reynolds": reynolds, "running": running, "turbulent": turbulent, "edge_temperature": t_e}


VISCOUS_DOC = """
**Method.** Skin friction after the pressure solve; the pressure is not changed. Connect the Modified Newtonian
Pressure block's output to **Pressure Coefficient** and use the same flight inputs.

* Edge state: `p/p_inf = 1 + gamma M^2 Cp / 2`; edge Mach isentropic from the pitot pressure; total temperature
  conserved. Air: R = 287.05 J/(kg K), Pr = 0.72, Sutherland viscosity.
* Laminar: Blasius with Eckert's reference temperature. Turbulent: White's closed form of van Driest II.
  Laminar to Turbulent switches suddenly at Transition Reynolds Number (running-length Reynolds number).
* Wall: adiabatic (recovery Pr^1/2 laminar, Pr^1/3 turbulent) or a fixed Wall Temperature.
* Shear acts along the freestream projected on the surface and fades out between 30 and 60 degrees of leeward
  inclination (separated flow).
* Running length: distance downstream of the body's most upstream point along the flow, over the cosine of the
  local inclination. Exact on flat plates and on cones at zero incidence; an approximation where surface
  streamlines curve.

**Totals.** Add these coefficients to the pressure coefficients. Axes and moment conventions are those of the
Modified Newtonian blocks.
"""


FACES_INPUT = InputSpec("Mesh Faces", "real_point_map", "$markdown-enabled$Points on the body's surface with the area each "
                        "stands for (m^2) as the value, e.g. **VSPAERO Mesh Faces (VSP)**, or Import Scalar Point Map "
                        "of an `x, y, z, area` file of mesh face centroids.")

FACES_DOC = """
**On mesh.** These coefficients are sums over the **Mesh Faces** points (Cf times the friction direction times each point's
area) instead of Surface Integrals over the body: one Evaluate Field per coefficient, seconds instead of minutes on
a full vehicle. The body still gives the surface normal at each point and the most upstream point, so the points
should lie on (or within a mesh tolerance of) the body. Accuracy is that of the point set: face centroids of a
mesh are a second-order rule, finest where faces are small; friction peaks at leading edges (laminar Cf goes as
1/sqrt(running length)), so coarse faces there under-resolve it.
"""
