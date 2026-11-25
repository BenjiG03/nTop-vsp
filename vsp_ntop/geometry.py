"""Shared vector, polygon and polyline helpers (numpy)."""
from __future__ import annotations

import numpy as np


def unit(v, error=ValueError, what="direction"):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise error(f"zero-length {what}")
    return v / n


def wind_axes(alpha_deg, beta_deg):
    """Unit freestream, lift and side-force directions in x aft, y right, z up."""
    a, b = np.radians(alpha_deg), np.radians(beta_deg)
    ca, sa, cb, sb = np.cos(a), np.sin(a), np.cos(b), np.sin(b)
    return np.array([ca * cb, -sb, sa * cb]), np.array([-sa, 0.0, ca]), np.array([ca * sb, cb, sa * sb])


def ls_slope(x, y):
    """Least-squares slope of y on x, or None with fewer than 2 distinct x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    dx = x - x.mean()
    sxx = dx @ dx
    if len(x) < 2 or sxx < 1e-12:
        return None
    return float(dx @ (y - y.mean()) / sxx)


def shoelace(loop):
    """Signed area (positive counter-clockwise) and centroid of a closed polygon (first point not repeated)."""
    p = np.asarray(loop, float)
    x0, y0 = np.roll(p, 1, axis=0).T
    x1, y1 = p.T
    w = x0 * y1 - x1 * y0
    a = w.sum()
    if abs(a) < 1e-300:
        return 0.5 * a, p.mean(axis=0)
    return 0.5 * a, np.array([((x0 + x1) * w).sum(), ((y0 + y1) * w).sum()]) / (3.0 * a)


def arc_length(points, closed=False):
    """Cumulative arc length along a polyline (back to the start when ``closed``)."""
    p = np.asarray(points, float)
    if closed:
        p = np.vstack([p, p[:1]])
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))])


def at_arc(points, s, closed=False):
    """Points at arc lengths ``s`` along a polyline (linear between its points)."""
    p = np.asarray(points, float)
    if closed:
        p = np.vstack([p, p[:1]])
    cum = arc_length(p)
    return np.column_stack([np.interp(s, cum, p[:, k]) for k in range(p.shape[1])])


def polygon_distance(q, loop):
    """Distance from point ``q`` to the closed polygon ``loop`` (any dimension)."""
    p1 = np.asarray(loop, float)
    p0 = np.roll(p1, 1, axis=0)
    d = p1 - p0
    l2 = (d * d).sum(axis=1)
    t = np.clip(((q - p0) * d).sum(axis=1) / np.where(l2 == 0.0, 1.0, l2), 0.0, 1.0)
    return float(np.sqrt(((q - p0 - t[:, None] * d) ** 2).sum(axis=1).min()))
