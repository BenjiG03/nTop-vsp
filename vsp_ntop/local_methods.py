"""Local-inclination pressure laws for the mesh driver (``run_newtonian_mesh``), after hypanel's ``methods.py``.

Each face's pressure depends only on its impact angle delta, the angle between the face and the freestream
(sin delta = -n . V, positive windward), and the freestream Mach number. Supersonic only (Mach > 1).

* **Tangent wedge** (windward): the exact two-dimensional oblique shock that turns the flow through delta (the weak
  solution), Cp = 2 / (gamma M^2) (p2 / p1 - 1) with p2 / p1 = 1 + 2 gamma / (gamma + 1) (M^2 sin^2 beta - 1). Above
  the detachment angle delta_max (no attached shock), Cp is scaled with sin^2 delta from the detachment value to
  Cp_max at delta = 90 degrees, so it rises continuously to the stagnation value (hypanel's ``_detached_bridge``).
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


def tangent_wedge_cp(delta, mach, cp_max, gamma=1.4, expansion="prandtl_meyer", separated="freestream", limit=None):
    """Face Cp from the impact angle ``delta`` (array, radians): tangent wedge windward (bridged to ``cp_max`` above
    detachment), Prandtl-Meyer or Newtonian shadow leeward."""
    np = _np()
    if mach <= 1.0:
        raise ValueError("the tangent wedge needs a supersonic freestream (Mach > 1)")
    if limit is None:
        key = ("limit", round(mach, 9), gamma)
        if key not in _TABLES:
            _TABLES[key] = max_deflection(mach, gamma)
        limit = _TABLES[key]
    wind = np.maximum(delta, 0.0)
    cp_lim = float(oblique_shock_cp(np.asarray([limit * (1.0 - 1e-9)]), mach, gamma)[0])
    attached = oblique_shock_cp(np.minimum(wind, limit * (1.0 - 1e-9)), mach, gamma)
    s2, s2l = np.sin(wind) ** 2, math.sin(limit) ** 2
    bridge = cp_lim + (cp_max - cp_lim) * np.clip((s2 - s2l) / max(1.0 - s2l, 1e-12), 0.0, 1.0)
    cp = np.where(wind <= limit, attached, bridge)
    cp = np.where(delta > 0.0, cp, 0.0)
    if expansion == "prandtl_meyer":
        cp = np.where(delta < 0.0, expansion_cp(delta, mach, gamma, separated), cp)
    return cp
