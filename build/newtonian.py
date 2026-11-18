"""Modified Newtonian aerodynamics as a graph of nTop blocks (from ntop_newtonian, now part of nTop_vsp).

The graph takes lists for Mach, Angle of Attack and Sideslip (one item per flight condition) and nTop evaluates
every block once per item: the Modified Newtonian (VSP) block feeds it a Flight Condition's conditions.

Everything the solver computes is an nTop block: the freestream direction and Cp_max from scalar
blocks, the surface pressure as a field built from the body's own gradient, and forces and moments
as integrals by nTop's Surface Integral. Python only writes the graph.

    Cp = Cp_max sin^2(delta),   sin(delta) = max(-n . V, 0),   n = grad f / |grad f|

so surfaces facing away from the flow carry no pressure. There is no self-shadowing. Forces are

    C_F . e = -(1 / S_ref) oint Cp (n . e) dA

and moments about the moment centre c are -oint Cp n . (e x (r - c)) dA, by the scalar triple product.
"""
from __future__ import annotations

from dataclasses import dataclass

from ntopfile import Graph, InputSpec, Ref
from ops import ANGLE, AREA, LENGTH, Ops

CP_MAX_MODES = ("Modified (Rayleigh Pitot)", "Classical (Cp_max = 2)", "User Defined (Cp Max Field)")
CP_MAX_FIELD_DOC = ("Cp_max when Stagnation Pressure is User Defined: a scalar field, so it may vary over the body "
                    "(for example one value on the fins and another on the nose). A constant is a field too. "
                    "Ignored by the other choices.")


def default_cp_max_field(g):
    """The User Defined Cp_max default: the constant field 2 (classical Newtonian)."""
    return g.leaf("real_field", {"expression": "2", "units": {}})
COEFFICIENTS = ("CA", "CY", "CN", "CD", "CL", "Cl", "Cm", "Cn", "Cp_max")
OUTPUTS = ("Coefficients", "Pressure Coefficient")


def _real(value, units=None):
    return lambda g: g.real(value, units)


def _choice(labels):
    return lambda g: g.choice(list(labels))


INPUTS = (
    InputSpec("Body", "implicit", "$markdown-enabled$The implicit body. Its gradient sets the surface normal, so any implicit works; "
              "it need not be an exact distance field."),
    InputSpec("Mach", "real", "Freestream Mach number, above 1.", default=_real(8.0)),
    InputSpec("Angle of Attack", "real", "Positive with the flow coming from below the body's Up Axis.", ANGLE, _real(0.0, ANGLE)),
    InputSpec("Sideslip", "real", "Positive with the flow coming from the +y side (y = Up Axis x Aft Axis).", ANGLE, _real(0.0, ANGLE)),
    InputSpec("Gamma", "real", "Ratio of specific heats of the perfect gas.", default=_real(1.4)),
    InputSpec("Aft Axis", "vector", "Body axis pointing from nose to tail. At zero angles the flow runs along it.",
              default=lambda g: g.vector(1, 0, 0)),
    InputSpec("Up Axis", "vector", "Approximately up in the body frame; sets the pitch plane.", default=lambda g: g.vector(0, 0, 1)),
    InputSpec("Reference Area", "real", "Normalizes force and moment coefficients.", AREA, _real(1.0, AREA)),
    InputSpec("Reference Length", "real", "Normalizes moment coefficients.", LENGTH, _real(1.0, LENGTH)),
    InputSpec("Moment Center", "point", "Point the moments are taken about.", default=lambda g: g.point(0, 0, 0)),
    InputSpec("Stagnation Pressure", "choice",
              "Modified: Cp_max behind a normal shock (Rayleigh pitot formula). Classical: Cp_max = 2. User Defined: "
              "the Cp Max Field.", default=_choice(CP_MAX_MODES)),
    InputSpec("Cp Max Field", "real_field", CP_MAX_FIELD_DOC, default=default_cp_max_field),
    InputSpec("Integration Tolerance", "real", "Relative error of nTop's Surface Integral.", default=_real(1e-3)),
)


@dataclass
class Newtonian:
    """Instance ids of the solver's results in a graph."""

    cp: int  # real_field
    coefficients: dict  # name -> real instance
    cp_max: int
    direction: int  # vector
    frame: tuple  # (x, y, z) vectors
    normal: int  # vector_field, unit outward normal from the gradient


def freestream(o: Ops, inp):
    """(x, y, z, V, sines and cosines): aero axes and the unit freestream velocity."""
    x = o.unit(inp("Aft Axis"))
    up = inp("Up Axis")
    z = o.unit(o.vsub(up, o.vscale(x, o.dot(up, x))))
    y = o.cross(z, x)
    ca, sa = o.cos(inp("Angle of Attack")), o.sin(inp("Angle of Attack"))
    cb, sb = o.cos(inp("Sideslip")), o.sin(inp("Sideslip"))
    v = o.vadd(o.vadd(o.vscale(x, o.mul(ca, cb)), o.vscale(y, o.mul(sb, o.real(-1.0)))), o.vscale(z, o.mul(sa, cb)),
               name="Freestream Direction")
    return x, y, z, v, (ca, sa, cb, sb)


def cp_max_stagnation(o: Ops, mach, gamma):
    """Stagnation pressure coefficient at any Mach number: behind a normal shock above Mach 1 (Rayleigh pitot),
    isentropic below it, 2/(gamma M^2) ((1 + (gamma-1)/2 M^2)^(gamma/(gamma-1)) - 1). Both are 1.839 at Mach 1."""
    w, m_sup, m_sub = o.supersonic(mach)
    one = o.real(1.0)
    gm1 = o.sub(gamma, one)
    m2 = o.sq(m_sub)
    isentropic = o.mul(o.div(o.real(2.0), o.mul(gamma, m2)),
                       o.sub(o.pow(o.add(one, o.mul(o.mul(gm1, o.real(0.5)), m2)), o.div(gamma, gm1)), one))
    return o.blend(w, cp_max_modified(o, m_sup, gamma), isentropic)


def cp_max_modified(o: Ops, mach, gamma):
    """2 / (gamma M^2) (p02 / p_inf - 1), with p02 / p_inf from the Rayleigh pitot formula."""
    one, two = o.real(1.0), o.real(2.0)
    m2 = o.sq(mach)
    gp1 = o.add(gamma, one)
    gm1 = o.sub(gamma, one)
    ratio = o.div(o.mul(o.sq(gp1), m2), o.sub(o.mul(o.real(4.0), o.mul(gamma, m2)), o.mul(two, gm1)))
    p02 = o.mul(o.pow(ratio, o.div(gamma, gm1)), o.div(o.add(o.sub(one, gamma), o.mul(two, o.mul(gamma, m2))), gp1))
    return o.mul(o.div(two, o.mul(gamma, m2)), o.sub(p02, one))


def surface_fields(o: Ops, body, v_up):
    """(unit normal field, sin(delta)+ field) for an upstream direction ``v_up`` = -V."""
    grad = o.grad(body)
    magnitude = o.fsqrt(o.fdot(grad, grad))
    normal = o._b("divide<vector_field,real_field>", "vector_field", [grad, magnitude])
    sine = o.fmax(o.fdot(normal, o.vconst(v_up)), o.field("0"), name="Impact Sine")
    return normal, sine


def build(g: Graph, inp, with_moments=True, integrate_forces=True) -> Newtonian:
    """Add the solver to ``g``. ``inp(name)`` gives the instance (or model input) for each of INPUTS."""
    o = Ops(g)
    g.comment("**1. Flight axes and freestream.** x is the unit Aft Axis, z the Up Axis made normal to x, y = z x x. "
              "The freestream direction is V = cos(a)cos(b) x - sin(b) y + sin(a)cos(b) z, with a the Angle of Attack "
              "and b the Sideslip.")
    x, y, z, v, (ca, sa, cb, sb) = freestream(o, inp)
    g.comment("**2. Stagnation pressure coefficient.** Modified: Cp_max = 2/(gamma M^2) (p02/p_inf - 1), with p02/p_inf "
              "from the Rayleigh pitot formula above Mach 1 and the isentropic total pressure below it. Classical: "
              "Cp_max = 2. User Defined: the Cp Max Field, which may vary over the body. The Stagnation Pressure "
              "dropdown picks one (0/1 weightings, so it works on lists).")
    cp_max = o.pick(inp("Stagnation Pressure"), [cp_max_stagnation(o, inp("Mach"), inp("Gamma")), o.real(2.0),
                                                 o.real(0.0)], name="Cp Max")
    user = o.select(inp("Stagnation Pressure"), [o.real(0.0), o.real(0.0), o.real(1.0)])
    cp_max_field = o.fadd(o.fmul(inp("Cp Max Field"), user), cp_max, name="Cp Max Field Used")
    body = inp("Body")
    v_up = o.vscale(v, o.real(-1.0))
    g.comment("**3. Surface pressure field.** The unit normal is n = grad f / |grad f| of the Body's field, so any implicit "
              "works. Impact Sine = max(-n.V, 0) is sin of the local inclination on windward faces and 0 on faces turned "
              "away from the flow. Pressure Coefficient = Cp_max * (Impact Sine)^2 (Lees' modified Newtonian). There is "
              "no self-shadowing.")
    normal, sine = surface_fields(o, body, v_up)
    cp = o.fmul(o.fsq(sine), cp_max_field, name="Pressure Coefficient")

    coefficients = integrate(o, inp, body, cp, normal, (x, y, z), (ca, sa, cb, sb), with_moments) if integrate_forces else {}
    coefficients["Cp_max"] = cp_max
    return Newtonian(cp, coefficients, cp_max, v, (x, y, z), normal)


def integrate(o: Ops, inp, body, cp, normal, frame, trig, with_moments=True):
    """Force and moment coefficients of a Cp field over the body's surface.

    The scalar forms -oint Cp (n . e) dA and -oint Cp n . (e x (r - c)) dA, with n the gradient normal,
    are used rather than fluxes: nTop's vector Surface Integral ignores its tolerance and runs 3.5% low
    on a 15 degree cone's axial force, where the scalar form is within 0.4%.
    """
    x, y, z = frame
    ca, sa, cb, sb = trig
    area, length, tol = inp("Reference Area"), inp("Reference Length"), inp("Integration Tolerance")
    o.g.comment("**4. Forces.** Each component is one Surface Integral over the Body: C.e = -(1/S_ref) * integral of "
                "Cp (n.e) dA, for e = x, y, z (CA, CY, CN). CD and CL follow by rotating into the wind axes. The scalar "
                "integrand is used rather than the vector (flux) form of Surface Integral, which proved less accurate.")
    minus_cp = o.fmul(cp, o.field("-1"))
    coefficients = {}
    for label, e in (("CA", x), ("CY", y), ("CN", z)):
        force = o.surface_integral(body, o.fmul(minus_cp, o.fdot(normal, o.vconst(e))), tol)
        coefficients[label] = o.div(force, area, name=label)
    coefficients["CD"] = o.add(o.add(o.mul(o.mul(ca, cb), coefficients["CA"]), o.mul(o.mul(sb, o.real(-1.0)), coefficients["CY"])),
                               o.mul(o.mul(sa, cb), coefficients["CN"]), name="CD")
    coefficients["CL"] = o.add(o.mul(o.mul(sa, o.real(-1.0)), coefficients["CA"]), o.mul(ca, coefficients["CN"]), name="CL")
    if with_moments:
        o.g.comment("**5. Moments** about the Moment Center c: -(1/(S_ref L_ref)) * integral of Cp n.(e x (r - c)) dA, "
                    "by the scalar triple product, for e = x, y, z (Cl, Cm, Cn; Cm nose-up positive).")
        centre = inp("Moment Center")
        position = o.vfield(o.field("x", LENGTH), o.field("y", LENGTH), o.field("z", LENGTH))
        arm = o.vfsub(position, o.vfield(o.comp(centre, "x"), o.comp(centre, "y"), o.comp(centre, "z")))
        for label, e in (("Cl", x), ("Cm", y), ("Cn", z)):
            moment = o.surface_integral(body, o.fmul(minus_cp, o.fdot(normal, o.fcross(o.vconst(e), arm))), tol)
            coefficients[label] = o.div(moment, o.mul(area, length), name=label)
    return coefficients


METHOD_DOC = """
**Method.** Modified Newtonian impact (Lees): `Cp = Cp_max * max(-n.V, 0)^2`, with `n` the unit gradient of the Body
and `V` the unit freestream velocity. Surfaces facing away from the flow carry no pressure; there is no self-shadowing.
Cp_max is the stagnation value behind a normal shock (Rayleigh pitot formula), 2 for classical Newtonian, or a
User Defined scalar field (Cp Max Field), which may vary over the body.

**Forces.** Each component is one Surface Integral: `C_F.e = -(1/S_ref) * integral(Cp n.e dA)`, and moments about the
Moment Center are `-(1/(S_ref L_ref)) * integral(Cp n.(e x (r - c)) dA)`.

**Axes.** x along Aft Axis (nose to tail), z toward Up Axis, y = z x x. The freestream is
`(cos a cos b, -sin b, sin a cos b)`: positive Angle of Attack has the flow coming from below, positive Sideslip
from +y. CA, CY, CN are along x, y, z; Cm is nose-up positive. CD is along the freestream, CL normal to it in the
x-z plane.

**Accuracy.** Set by nTop's Surface Integral: about 1e-4 on smooth bodies, 0.1 to 1% across sharp edges. The
Body's field should behave like a distance near its surface.
"""
