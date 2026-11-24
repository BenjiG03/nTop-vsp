"""Local-inclination pressure laws for the mesh driver (``run_newtonian_mesh``), after hypanel's ``methods.py``.

Each face's pressure depends only on its impact angle delta, the angle between the face and the freestream
(sin delta = -n . V, positive windward), and the freestream Mach number. Supersonic only (Mach > 1).

* **Tangent wedge** (windward): the exact two-dimensional oblique shock that turns the flow through delta (the weak
  solution), Cp = 2 / (gamma M^2) (p2 / p1 - 1) with p2 / p1 = 1 + 2 gamma / (gamma + 1) (M^2 sin^2 beta - 1). Above
  the detachment angle delta_max (no attached shock), Cp is scaled with sin^2 delta from the detachment value to
  Cp_max at delta = 90 degrees, so it rises continuously to the stagnation value (hypanel's ``_detached_bridge``).
* **Tangent cone** (windward): the exact cone of half angle delta (Taylor-Maccoll, weak shock), tabulated per Mach
  number from the shock angle (hypanel's ``conical``), with the same bridge to Cp_max above the cone's detachment.
* **Prandtl-Meyer** (leeward): the isentropic expansion from the freestream through -delta, nu(M2) = nu(M) - delta,
  limited at vacuum. Faces turned more than 30 degrees away blend (smoothstep) to a separated-flow value by 60 degrees:
  freestream pressure (Cp 0) or Gabeaud's base pressure -1 / M^2.
* **Newtonian shadow** (leeward): Cp 0.

Numpy, vectorised over faces.
"""
from __future__ import annotations

import math

SEPARATION = (math.radians(30.0), math.radians(60.0))


_NUMPY = []


def _np():
    if not _NUMPY:  # once: the mesh driver's loader also puts conda's DLL folders on PATH
        from .run_newtonian_mesh import _np as np_
        _NUMPY.append(np_())
    return _NUMPY[0]


def max_deflection(mach, gamma=1.4):
    """Largest wedge angle with an attached shock (radians): the maximum over the shock angle of the theta-beta-M
    relation, tan theta = 2 cot beta (M^2 sin^2 beta - 1) / (M^2 (gamma + cos 2 beta) + 2)."""
    np = _np()
    mu = math.asin(1.0 / mach)
    beta = np.linspace(mu, 0.5 * math.pi, 20001)
    th = np.arctan(2.0 / np.tan(beta) * (mach**2 * np.sin(beta) ** 2 - 1.0) / (mach**2 * (gamma + np.cos(2 * beta)) + 2))
    return float(th.max())


def oblique_shock_cp(delta, mach, gamma=1.4):
    """Cp behind the weak oblique shock of deflection ``delta`` (array, radians, 0 <= delta <= delta_max).

    Shock angle from Thompson's closed form of the cubic in tan beta (the weak root)."""
    np = _np()
    d = np.clip(delta, 1e-10, None)
    m2 = mach * mach
    a = 1.0 + 0.5 * (gamma - 1.0) * m2
    t = np.tan(d)
    lam = np.sqrt(np.maximum((m2 - 1.0) ** 2 - 3.0 * a * (1.0 + 0.5 * (gamma + 1.0) * m2) * t * t, 0.0))
    chi = ((m2 - 1.0) ** 3 - 9.0 * a * (a + 0.25 * (gamma + 1.0) * m2 * m2) * t * t) / np.maximum(lam**3, 1e-300)
    tan_beta = (m2 - 1.0 + 2.0 * lam * np.cos((4.0 * math.pi + np.arccos(np.clip(chi, -1.0, 1.0))) / 3.0)) / (3.0 * a * t)
    beta = np.arctan(tan_beta)
    p21 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (m2 * np.sin(beta) ** 2 - 1.0)
    return 2.0 / (gamma * m2) * (p21 - 1.0)


def prandtl_meyer(m, gamma=1.4):
    np = _np()
    g = math.sqrt((gamma + 1.0) / (gamma - 1.0))
    r = np.sqrt(np.maximum(m * m - 1.0, 0.0))
    return g * np.arctan(r / g) - np.arctan(r)


def inverse_prandtl_meyer(nu, gamma=1.4):
    """Mach number with Prandtl-Meyer angle ``nu`` (array; at most nu_max), by bisection in log Mach."""
    np = _np()
    lo = np.zeros_like(nu)
    hi = np.full_like(nu, math.log(1e4))
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        low = prandtl_meyer(np.exp(mid), gamma) < nu
        lo, hi = np.where(low, mid, lo), np.where(low, hi, mid)
    return np.exp(0.5 * (lo + hi))


_TABLES = {}


def _expansion_table(mach, gamma):
    """Prandtl-Meyer Cp against turning angle, 0 to 90 degrees on 4001 points (interpolated per face)."""
    key = (round(mach, 9), gamma)
    if key not in _TABLES:
        np = _np()
        turn = np.linspace(0.0, 0.5 * math.pi, 4001)
        nu_max = 0.5 * math.pi * (math.sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0)
        nu = np.minimum(float(prandtl_meyer(np.asarray(mach), gamma)) + turn, nu_max * (1.0 - 1e-9))
        m2 = inverse_prandtl_meyer(nu, gamma)
        ratio = ((1.0 + 0.5 * (gamma - 1.0) * mach**2) / (1.0 + 0.5 * (gamma - 1.0) * m2**2)) ** (gamma / (gamma - 1.0))
        _TABLES[key] = (turn, 2.0 / (gamma * mach**2) * (ratio - 1.0))
    return _TABLES[key]


def expansion_cp(delta, mach, gamma=1.4, separated="freestream"):
    """Leeward Cp (delta <= 0): Prandtl-Meyer from the freestream, limited at vacuum, blended to the separated value
    between 30 and 60 degrees of turning."""
    np = _np()
    turn = np.maximum(-delta, 0.0)
    grid, table = _expansion_table(mach, gamma)
    cp = np.interp(turn, grid, table)
    floor = 0.0 if separated == "freestream" else -1.0 / mach**2
    s = np.clip((turn - SEPARATION[0]) / (SEPARATION[1] - SEPARATION[0]), 0.0, 1.0)
    w = s * s * (3.0 - 2.0 * s)
    return cp + (floor - cp) * w


def _normal_shock_total_pressure_ratio(mn, gamma):
    """p02 / p01 across a normal shock of upstream Mach ``mn`` (array, >= 1)."""
    np = _np()
    m2 = np.maximum(mn, 1.0) ** 2
    a = (gamma + 1.0) * m2 / (2.0 + (gamma - 1.0) * m2)
    b = (gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0))
    return a ** (gamma / (gamma - 1.0)) * b ** (1.0 / (gamma - 1.0))


def _cones_from_shocks(beta, mach, gamma, steps=400):
    """(cone half angle, surface Cp) behind conical shocks of angle ``beta`` (array): the Taylor-Maccoll equation
    integrated from the shock to the cone surface by RK4 in V_theta (hypanel's ``conical.cone_from_shock``).

    In theta the equation is singular at the Mach cone; in V_theta it is regular, and V_theta grows monotonically from
    its value behind the shock to zero at the surface, so the end point is known. Both ends get small steps:
    V_theta = V_shock (1 - S^2), S the smoothstep of t in [0, 1]. Velocities are over the limiting speed."""
    np = _np()
    sb = np.sin(beta)
    mn = np.maximum(mach * sb, 1.0)
    delta = np.arctan2(2.0 * ((mach * sb) ** 2 - 1.0) / np.tan(beta), mach**2 * (gamma + np.cos(2.0 * beta)) + 2.0)
    mn2 = np.sqrt((1.0 + 0.5 * (gamma - 1.0) * mn**2) / (gamma * mn**2 - 0.5 * (gamma - 1.0)))
    m2 = mn2 / np.maximum(np.sin(beta - delta), 1e-6)
    speed = (2.0 / ((gamma - 1.0) * m2**2) + 1.0) ** -0.5
    vt0 = -speed * np.sin(beta - delta)

    def rates(theta, vr, vt):
        a2 = 0.5 * (gamma - 1.0) * (1.0 - vr**2 - vt**2)
        dtheta = (a2 - vt**2) / (vt**2 * vr - a2 * (2.0 * vr + vt / np.tan(theta)))
        return dtheta, vt * dtheta

    def slope(theta, vr, t):
        s = t * t * (3.0 - 2.0 * t)
        dt, dv = rates(theta, vr, vt0 * (1.0 - s * s))
        k = -vt0 * 2.0 * s * 6.0 * t * (1.0 - t)
        return dt * k, dv * k

    h = 1.0 / steps
    theta, vr = beta.copy(), speed * np.cos(beta - delta)
    for i in range(steps):
        t = i * h
        k1 = slope(theta, vr, t)
        k2 = slope(theta + 0.5 * h * k1[0], vr + 0.5 * h * k1[1], t + 0.5 * h)
        k3 = slope(theta + 0.5 * h * k2[0], vr + 0.5 * h * k2[1], t + 0.5 * h)
        k4 = slope(theta + h * k3[0], vr + h * k3[1], t + h)
        theta = theta + h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        vr = vr + h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
    v2 = np.clip(vr**2, 1e-6, 1.0 - 1e-6)
    mach_c = np.sqrt(2.0 / (gamma - 1.0) * v2 / (1.0 - v2))
    p0 = lambda m: (1.0 + 0.5 * (gamma - 1.0) * m**2) ** (gamma / (gamma - 1.0))  # noqa: E731
    ratio = p0(mach) * _normal_shock_total_pressure_ratio(mn, gamma) / p0(mach_c)
    return theta, 2.0 / (gamma * mach**2) * (ratio - 1.0)


def cone_table(mach, gamma=1.4, points=320, thin_points=64):
    """Weak-branch exact cone solutions from half angle 0 to the detachment angle: (theta, Cp), theta increasing.

    Thin cones need shock angles very close to the Mach angle (geometric spacing there), the rest quadratic spacing
    (hypanel's ``conical.cone_table``). Cached per Mach number."""
    key = ("cone", round(mach, 9), gamma)
    if key not in _TABLES:
        np = _np()
        mu = math.asin(1.0 / mach)
        span = 0.5 * math.pi - 1e-3 - mu
        wide = span * np.linspace(0.0, 1.0, points + 1)[1:] ** 2
        thin = np.exp(np.linspace(math.log(1e-8), math.log(wide[0]), thin_points, endpoint=False))
        theta, cp = _cones_from_shocks(mu + np.concatenate([thin, wide]), mach, gamma)
        last = int(np.argmax(theta))  # detachment: the largest cone angle on the weak branch
        _TABLES[key] = (np.concatenate([[0.0], theta[:last + 1]]), np.concatenate([[0.0], cp[:last + 1]]))
    return _TABLES[key]


def cone_cp(delta, mach, gamma=1.4):
    """Surface Cp of the exact cone of half angle ``delta`` (array, clamped to 0..detachment); sqrt(Cp) is nearly
    linear in the angle, so it is what is interpolated."""
    np = _np()
    theta, cp = cone_table(mach, gamma)
    return np.interp(np.clip(delta, 0.0, theta[-1]), theta, np.sqrt(np.maximum(cp, 0.0))) ** 2


IMPACTS = ("wedge", "cone")


def tangent_wedge_cp(delta, mach, cp_max, gamma=1.4, expansion="prandtl_meyer", separated="freestream", limit=None,
                     impact="wedge"):
    """Face Cp from the impact angle ``delta`` (array, radians): windward the tangent wedge (``impact = wedge``, the
    oblique shock) or tangent cone (``cone``, the exact Taylor-Maccoll cone), bridged to ``cp_max`` above detachment;
    leeward Prandtl-Meyer or Newtonian shadow."""
    np = _np()
    if mach <= 1.0:
        raise ValueError("the tangent wedge and cone need a supersonic freestream (Mach > 1)")
    if impact not in IMPACTS:
        raise ValueError(f"impact = {impact}: expected one of {', '.join(IMPACTS)}")
    if impact == "cone":
        theta, _ = cone_table(mach, gamma)
        limit = float(theta[-1])
        surface = cone_cp
    else:
        if limit is None:
            key = ("limit", round(mach, 9), gamma)
            if key not in _TABLES:
                _TABLES[key] = max_deflection(mach, gamma)
            limit = _TABLES[key]
        surface = oblique_shock_cp
    wind = np.maximum(delta, 0.0)
    cp_lim = float(surface(np.asarray([limit * (1.0 - 1e-9)]), mach, gamma)[0])
    attached = surface(np.minimum(wind, limit * (1.0 - 1e-9)), mach, gamma)
    s2, s2l = np.sin(wind) ** 2, math.sin(limit) ** 2
    bridge = cp_lim + (cp_max - cp_lim) * np.clip((s2 - s2l) / max(1.0 - s2l, 1e-12), 0.0, 1.0)
    cp = np.where(wind <= limit, attached, bridge)
    cp = np.where(delta > 0.0, cp, 0.0)
    if expansion == "prandtl_meyer":
        cp = np.where(delta < 0.0, expansion_cp(delta, mach, gamma, separated), cp)
    return cp
