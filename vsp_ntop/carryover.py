"""Lift a fin carries onto the body, after Pitts, Nielsen and Kaattari (NACA TR 1307).

TR 1307 writes the lift of a wing-body combination as factors on the lift of the wing alone (the exposed panels
joined at the root). ``K_B(W)`` is the lift the fins put on the body. A vortex-lattice wing ending at the body surface
has none of it, and one extended through the body to the axis over-predicts it on a slender body (its inboard strips
behave as more wing).

Supersonic, TR 1307's planar model (its Figure 4) flattens the body into a strip of width d beside each fin root. The
fin's pressure reaches the strip between the Mach lines from the root's leading and trailing edges; equations 23 and
25 give it, from the conical flow of the swept leading edge, in units of the two-dimensional 4 alpha / beta:

    beta m > 1:  beta m / (pi sqrt(beta^2 m^2 - 1)) arccos((xi / beta + beta m eta) / (eta + m xi))
    beta m < 1:  2 (beta m)^(3/2) / (pi (beta m + 1)) sqrt(xi / beta - eta) / sqrt(m xi + eta)

with xi aft of the root leading edge, eta across the strip from the root and m the cotangent of the leading-edge
sweep. The report integrates it in closed form for a body that runs on behind the fin (equations 24, 26) or ends at
the trailing edge (30, 31); here it is integrated numerically up to wherever the body ends. Its integral over the
fin's planform area is K_B(W), its first moments place the load (TR 1307 equations 66 to 70).

When TR 1307's rule beta A (1 + lambda) (1 / (beta m) + 1) <= 4 says the fin's tip would reach the strip, and at
subsonic speed, slender-body theory (equation 21) gives K_B(W) instead; the load is still placed where the planar
model puts it.

The 4 alpha / beta of the planar model stands for the fin's own mean lifting pressure, so the carried lift is K_B(W)
times the fin's lift from VSPAERO, along the fin's normal. Nothing here assumes a circular body: root distance, body
width across the fin's plane and the body's length behind the root come from the OpenVSP model.
"""
from __future__ import annotations

import math

import numpy as np

SONIC_BAND = 0.02  # beta m within this of 1 blends the two leading-edge solutions
SELECTION = 4.0  # TR 1307's rule between slender-body theory and the planar model
QUADRATURE = 48  # Gauss points in each direction over the strip

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(QUADRATURE)
_U = 0.5 * (_NODES + 1.0)
_W = 0.5 * _WEIGHTS


def slender_factors(ratio):
    """(K_W(B), K_B(W)) of slender-body theory, TR 1307 equations 14 and 21, at root distance / semispan ``ratio``."""
    tau = min(max(float(ratio), 1e-9), 0.999)
    brace = ((1.0 + tau**4) * (0.5 * math.atan(0.5 * (1.0 / tau - tau)) + 0.25 * math.pi)
             - tau**2 * ((1.0 / tau - tau) + 2.0 * math.atan(tau)))
    k_wb = 2.0 / math.pi * brace / (1.0 - tau) ** 2
    return k_wb, (1.0 + tau) ** 2 - k_wb


def strip_pressure(xi, eta, beta, m):
    """Lifting pressure on the body strip in units of 4 alpha / beta, TR 1307 equations 23 and 25."""
    bm = beta * m
    tiny = 1e-30

    def supersonic(bm_edge):
        slope = bm_edge / beta
        arg = (xi / beta + bm_edge * eta) / np.maximum(eta + slope * xi, tiny)
        return bm_edge / (math.pi * math.sqrt(bm_edge**2 - 1.0)) * np.arccos(np.clip(arg, -1.0, 1.0 - 1e-12))

    def subsonic(bm_edge):
        slope = bm_edge / beta
        root = np.sqrt(np.maximum(xi / beta - eta, tiny)) / np.sqrt(np.maximum(slope * xi + eta, tiny))
        return 2.0 * bm_edge**1.5 / (math.pi * (bm_edge + 1.0)) * root

    blend = min(max((bm - (1.0 - SONIC_BAND)) / (2.0 * SONIC_BAND), 0.0), 1.0)
    value = 0.0
    if blend > 0.0:
        value = value + blend * supersonic(max(bm, 1.0 + SONIC_BAND))
    if blend < 1.0:
        value = value + (1.0 - blend) * subsonic(min(bm, 1.0 - SONIC_BAND))
    return np.where((xi >= beta * eta) & (eta >= 0.0), value, 0.0)


def planar_integrals(beta, m, width, root_chord, length):
    """(lift, xi moment, eta moment) of one fin's strip pressure, as areas in units of 4 alpha / beta.

    The strip runs across the body from eta = 0 to ``width``, between the Mach lines from the root's leading and
    trailing edges, and stops where the body does, ``length`` aft of the leading edge. Moments are about the root
    leading edge (xi aft, eta across the body)."""
    top = min(width, length / beta)
    knee = min(max((length - root_chord) / beta, 0.0), top)  # aft of here the strip stops at the body's end
    total = np.zeros(3)
    for start, stop in ((0.0, knee), (knee, top)):
        if stop <= start:
            continue
        eta = start + (stop - start) * _U
        low = beta * eta
        span = np.maximum(np.minimum(root_chord + beta * eta, length) - low, 0.0)
        xi = low[:, None] + span[:, None] * _U[None, :] ** 2  # u^2 takes out the square-root rise behind the Mach line
        weight = _W[None, :] * 2.0 * _U[None, :] * span[:, None]
        p = strip_pressure(xi, eta[:, None], beta, m) * weight
        rows = (stop - start) * _W
        total += [np.sum(rows * p.sum(axis=1)), np.sum(rows * (p * xi).sum(axis=1)),
                  np.sum(rows * eta * p.sum(axis=1))]
    return total


def body_factor(mach, root, semispan, root_chord, tip_chord, slope, planform, width, length):
    """K_B(W) of one fin, where its carried lift acts (xi aft of the root leading edge, eta across the body from the
    root) and which method set it ("planar" or "slender")."""
    if mach <= 1.0 + 1e-6:  # no Mach lines: slender-body theory, at the root's mid-chord and the body's middle
        return slender_factors(root / semispan)[1], 0.5 * root_chord, 0.5 * width, "slender"
    beta = math.sqrt(mach * mach - 1.0)
    exposed = semispan - root
    aspect = 4.0 * exposed / (root_chord + tip_chord)  # the exposed pair joined at the root
    taper = tip_chord / root_chord
    lift, m_xi, m_eta = planar_integrals(beta, slope, width, root_chord, length)
    planar = beta * aspect * (1.0 + taper) * (1.0 / (beta * slope) + 1.0) > SELECTION
    k = lift / planform if planar else slender_factors(root / semispan)[1]
    if lift <= 0.0:
        return k, 0.5 * root_chord, 0.5 * width, "planar" if planar else "slender"
    return k, m_xi / lift, m_eta / lift, "planar" if planar else "slender"
