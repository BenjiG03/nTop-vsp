"""A smooth B-spline surface through an implicit body's nTop slices: the CAD surface without OpenVSP's skinning."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import fit as F
from . import slices as S

DEGREE = 3


# --- B-spline basics

def basis(knots, t, degree=DEGREE):
    """Matrix of B-spline basis functions: row i = all basis functions at t[i] (Cox-de Boor). ``knots`` is the full
    knot vector; t = the last knot belongs to the last span."""
    knots = np.asarray(knots, float)
    t = np.atleast_1d(np.asarray(t, float))
    n = len(knots) - degree - 1
    # degree 0: the span containing t (the last non-empty span takes t = end)
    last = np.max(np.nonzero(knots[:-1] < knots[-1])[0]) if np.any(knots[:-1] < knots[-1]) else len(knots) - 2
    N = np.zeros((len(t), len(knots) - 1))
    for j in range(len(knots) - 1):
        if knots[j] < knots[j + 1]:
            inside = (t >= knots[j]) & (t < knots[j + 1])
            if j == last:
                inside |= t == knots[j + 1]
            N[inside, j] = 1.0
    for d in range(1, degree + 1):
        M = np.zeros((len(t), len(knots) - 1 - d))
        for j in range(len(knots) - 1 - d):
            a = knots[j + d] - knots[j]
            b = knots[j + d + 1] - knots[j + 1]
            if a > 0:
                M[:, j] += (t - knots[j]) / a * N[:, j]
            if b > 0:
                M[:, j] += (knots[j + d + 1] - t) / b * N[:, j + 1]
        N = M
    return N[:, :n]


def ring_knots(nw, creases=(), degree=DEGREE):
    """Knots of one turn round the body, w in [0, 1), repeats included: ``nw`` even spans, a quarter's worth starting
    at each of w = 0, 1/4, 1/2, 3/4."""
    per = nw // 4
    out = []
    for q in range(4):
        a = q / 4.0
        out += [a] * (degree if q in creases else 1)
        out += [a + k / nw for k in range(1, per)]
    return np.array(out)


def unrolled(ring, extra=0, degree=DEGREE):
    """The periodic knot sequence written out over [0, 1] (plus ``extra`` more knots): the last ``degree`` knots of
    the turn before, the turn, and the first ``degree + 1 + extra`` of the turn after."""
    ring = list(ring)
    return np.array([t - 1.0 for t in ring[len(ring) - degree:]] + ring
                    + [t + 1.0 for t in ring[:degree + 1 + extra]])


def periodic_basis(ring, w, degree=DEGREE):
    """Basis matrix of the periodic cubic with knots ``ring`` (one control point per knot): the unrolled basis folded
    back onto the turn's control points."""
    n = len(ring)
    B = basis(unrolled(ring, 0, degree), np.mod(w, 1.0), degree)
    out = B[:, :n].copy()
    out[:, :degree] += B[:, n:n + degree]
    return out


def clamped_knots(u, n_ctrl, degree=DEGREE):
    """Clamped knot vector on [u[0], u[-1]] for ``n_ctrl`` control points, interior knots spread like the data
    (de Boor's averaging over the data parameters), so every span holds data."""
    u = np.asarray(u, float)
    n_int = n_ctrl - degree - 1
    if n_int <= 0:
        inner = np.array([])
    else:
        q = np.linspace(0.0, 1.0, n_int + 2)[1:-1]
        inner = np.interp(q, np.linspace(0.0, 1.0, len(u)), u)
    return np.concatenate([[u[0]] * (degree + 1), inner, [u[-1]] * (degree + 1)])


def second_difference(n, periodic=False):
    D = np.zeros((n if periodic else n - 2, n))
    for i in range(D.shape[0]):
        for k, c in enumerate((1.0, -2.0, 1.0)):
            D[i, (i + k) % n] = c
    return D


def insert_knot(knots, ctrl, t, axis, degree=DEGREE):
    """Boehm knot insertion of ``t`` once into a B-spline whose control points run along ``axis`` of ``ctrl``.
    Returns (knots, ctrl): the same curve or surface, one more control point."""
    knots = np.asarray(knots, float)
    P = np.moveaxis(np.asarray(ctrl, float), axis, 0)
    k = int(np.searchsorted(knots, t, side="right") - 1)  # knots[k] <= t < knots[k + 1]
    Q = np.empty((P.shape[0] + 1,) + P.shape[1:])
    Q[:k - degree + 1] = P[:k - degree + 1]
    for i in range(k - degree + 1, k + 1):
        a = (t - knots[i]) / (knots[i + degree] - knots[i])
        Q[i] = (1.0 - a) * P[i - 1] + a * P[i]
    Q[k + 1:] = P[k:]
    return np.insert(knots, k + 1, t), np.moveaxis(Q, 0, axis)


def evaluate(patch, u, w):
    """Points of a clamped patch at parameter arrays u (along the body) and w (round it): shape (len(u), len(w), 3)."""
    Bu = basis(patch.knots_u, u)
    Bw = basis(patch.knots_w, w)
    return np.einsum("ia,abk,jb->ijk", Bu, patch.ctrl, Bw)


# --- The body surface

@dataclass
class Patch:
    """One clamped bicubic B-spline patch: ``ctrl`` (n_u, n_w, 3) metres, full knot vectors."""
    ctrl: np.ndarray
    knots_u: np.ndarray
    knots_w: np.ndarray
    w_range: tuple


@dataclass
class End:
    x: float          # station along the slicing axis (frame z)
    pointed: bool
    point: tuple      # the tip (pointed) or the cap's centre (blunt)


@dataclass
class Loft:
    patches: list                 # Patch per face round the body
    start: End
    end: End
    frame: F.Frame
    notes: list = field(default_factory=list)
    rows: int = 0                 # data rows fitted (slices, plus the ends)
    residual: tuple = (0.0, 0.0)  # rms, max distance of the contour points from the surface (metres)
    row_residual: tuple = ()      # (station, largest distance) per data row, metres
    creases: tuple = ()           # quarter points with a sharp edge along the body (0 right, 1 bottom, 2 left, 3 top)


def _row_samples(loop, per_quarter):
    """Contour points of one slice in the surface's w order, (4 per_quarter, 2) in the slice's (u, v), split at the
    widest points (OpenVSP Model's sections, fit._section_curve)."""
    uc, vc, w, h, pts, params = F._section_curve(loop, per_quarter)
    return pts[:-1] * [w, h] + [uc, vc], (uc, vc)


def _crossings(pts, axis, value):
    """Where the closed polyline crosses the line (coordinate ``axis`` = value): (segment indices, points)."""
    a, b = pts, np.roll(pts, -1, axis=0)
    da, db = a[:, axis] - value, b[:, axis] - value
    i = np.flatnonzero((da == 0.0) | ((da < 0.0) != (db < 0.0)))
    t = np.where(da[i] == 0.0, 0.0, da[i] / np.where(da[i] == db[i], 1.0, da[i] - db[i]))
    return i, a[i] + t[:, None] * (b[i] - a[i])


def _ring_samples(loop, per_quarter, widest=(False, False)):
    """Contour points of one slice in the surface's w order, (4 per_quarter, 2) in the slice's (u, v)."""
    u0, u1, v0, v1 = S.bounds(loop)
    uc, vc = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
    pts = np.asarray(loop) if S.area(loop) < 0 else np.asarray(loop)[::-1]  # clockwise in (u, v)
    n = len(pts)
    (hi, hp), (vi, vp) = _crossings(pts, 1, vc), _crossings(pts, 0, uc)
    if not len(hi) or not len(vi):
        return _row_samples(loop, per_quarter)

    def widest_point(sign):
        i = int(np.argmax(sign * pts[:, 0]))
        return i, pts[i]

    def pick(idx, cp, k):
        return int(idx[k]), cp[k]

    right = widest_point(1.0) if widest[0] else pick(hi, hp, np.argmax(hp[:, 0]))
    left = widest_point(-1.0) if widest[1] else pick(hi, hp, np.argmin(hp[:, 0]))
    anchors = [right, pick(vi, vp, np.argmin(vp[:, 1])), left, pick(vi, vp, np.argmax(vp[:, 1]))]
    # walking forward from the right
    pos = [(a[0] - right[0]) % n for a in anchors]
    if not pos[0] <= pos[1] <= pos[2] <= pos[3]:
        return _row_samples(loop, per_quarter)
    out = []
    for q in range(4):
        (ia, pa), (ib, pb) = anchors[q], anchors[(q + 1) % 4]
        poly = np.vstack([pa, pts[(ia + 1 + np.arange((ib - ia) % n)) % n], pb])
        out.append(S.resample(poly, per_quarter + 1)[:-1])
    return np.vstack(out), (uc, vc)


CREASE_ANGLE = 30.0  # degrees: a quarter point where the contours typically turn by more is a crease


def corner_angles(rows, per_quarter):
    """Turning angle (degrees) of each contour at the four quarter points, from the chords two samples either side:
    a few degrees on a smooth section, the corner angle at a sharp edge. ``rows``: (n, 4 per_quarter, 2)."""
    m = rows.shape[1]
    out = np.zeros((rows.shape[0], 4))
    for q in range(4):
        j = q * per_quarter
        a = rows[:, j] - rows[:, (j - 2) % m]
        b = rows[:, (j + 2) % m] - rows[:, j]
        c = np.sum(a * b, axis=1) / np.maximum(np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-300)
        out[:, q] = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    return out


def loft_body(layers, frame: F.Frame, nw=48, per_quarter=24, rows_per_ctrl=2.0, smooth=1e-4, bounds=None,
              creases="auto", faces=4):
    """Fit the body's surface to its slices."""
    if nw % 4 or nw < 8:
        raise F.FitError("control points round the body must be a multiple of 4, at least 8")
    if faces % 4 or faces < 4 or nw % faces:
        raise F.FitError(f"faces round the body ({faces}) must be a multiple of 4 that divides the control points "
                         f"round it ({nw})")
    notes, used = [], []
    for l in layers:
        outer = S.outer_loops(l)
        if not outer:
            continue
        if len(outer) > 1:
            notes.append(f"z = {l.z:.4f} m: {len(outer)} separate contours; the largest is used")
        used.append(S.Layer(l.z, [outer[0]]))
    if len(used) < 6:
        raise F.FitError(f"only {len(used)} slices cut the body; use more slices")
    zs = np.array([l.z for l in used])
    spacing = (zs[-1] - zs[0]) / (len(used) - 1)
    areas = np.array([abs(S.area(l.loops[0])) for l in used])
    b0, b1 = bounds if bounds else (None, None)
    x0, pointed0, _ = F._end(zs, areas, True, spacing, b0)
    x1, pointed1, _ = F._end(zs, areas, False, spacing, b1)

    # sharp edges along the body at the four split points (chines)
    first_pass = np.stack([_row_samples(l.loops[0], per_quarter)[0] for l in used])
    angles = corner_angles(first_pass, per_quarter)
    k0, k1 = int(0.1 * len(angles)), max(int(0.9 * len(angles)), int(0.1 * len(angles)) + 1)
    typical = np.median(angles[k0:k1], axis=0)
    if creases == "auto":
        creases = tuple(q for q in range(4) if typical[q] > CREASE_ANGLE)
    creases = tuple(sorted(set(creases)))
    if creases:
        names = ("right", "bottom", "left", "top")
        notes.append("sharp edge along the body at the " + ", ".join(f"{names[q]} ({typical[q]:.0f} deg)" for q in creases))

    # data rows: (z, points (m, 2) in the slice plane), with each end added
    rows = []
    centres = []
    for l in used:
        pts, c = _ring_samples(l.loops[0], per_quarter, (0 in creases, 2 in creases))
        rows.append((l.z, pts))
        centres.append(c)
    m = len(rows[0][1])

    def end_row(at_start):
        z, pts = rows[0] if at_start else rows[-1]
        x = x0 if at_start else x1
        pointed = pointed0 if at_start else pointed1
        if pointed:
            c = centres[0] if at_start else centres[-1]
            return (x, np.tile(np.array(c), (m, 1))), pointed
        return (x, pts.copy()), pointed  # the side wall carried out to the end plane

    first, _ = end_row(True)
    last, _ = end_row(False)
    rows = [first] + rows + [last]
    z = np.array([r[0] for r in rows])
    P2 = np.stack([r[1] for r in rows])                              # (n_rows, m, 2)
    P = np.empty(P2.shape[:2] + (3,))
    P[:] = frame.point(P2[..., 0], P2[..., 1], z[:, None])

    # u: distance along the body on the centres and the equivalent radii
    cen = P.mean(axis=1)
    rad = np.sqrt(np.concatenate([[0.0 if pointed0 else areas[0]], areas, [0.0 if pointed1 else areas[-1]]]) / np.pi)
    step = np.sqrt(np.sum(np.diff(cen, axis=0) ** 2, axis=1) + np.diff(rad) ** 2)
    step = np.maximum(step, 1e-9)
    u = np.concatenate([[0.0], np.cumsum(step)])
    u /= u[-1]

    # w of the samples
    w = np.arange(m) / m

    # round the body
    ring = ring_knots(nw, creases)
    n_ring = len(ring)
    Bw = periodic_basis(ring, w)
    Dw = second_difference(n_ring, periodic=True)
    corner_ctrl = [int(np.argmax(periodic_basis(ring, [q / 4.0])[0])) for q in creases]
    Dw = Dw[[i for i in range(n_ring) if (i + 1) % n_ring not in corner_ctrl]]
    lam_w = smooth * np.trace(Bw.T @ Bw) / max(np.trace(Dw.T @ Dw), 1e-30)
    Aw = Bw.T @ Bw + lam_w * Dw.T @ Dw
    R = np.einsum("ij,rjk->rik", np.linalg.solve(Aw, Bw.T), P)         # (n_rows, n_ring, 3)

    # along the body
    n_rows = len(rows)
    nu = int(max(DEGREE + 3, min(n_rows, round(n_rows / rows_per_ctrl))))
    ku = clamped_knots(u, nu)
    Bu = basis(ku, u)
    Du = second_difference(nu)
    lam_u = smooth * np.trace(Bu.T @ Bu) / max(np.trace(Du.T @ Du), 1e-30)
    fixed = {0: R[0], nu - 1: R[-1]}
    free = [i for i in range(nu) if i not in fixed]
    rhs = R.reshape(n_rows, -1) - sum(np.outer(Bu[:, i], fixed[i].ravel()) for i in fixed)
    A = Bu[:, free]
    Dfree = Du[:, free]
    Dfix = sum(np.outer(Du[:, i], fixed[i].ravel()) for i in fixed)
    lhs = A.T @ A + lam_u * Dfree.T @ Dfree
    sol = np.linalg.solve(lhs, A.T @ rhs - lam_u * Dfree.T @ Dfix)
    C = np.empty((nu, n_ring, 3))
    for i, v in fixed.items():
        C[i] = v
    C[free] = sol.reshape(len(free), n_ring, 3)
    if pointed0:
        C[0] = P[0, 0]
    if pointed1:
        C[-1] = P[-1, 0]
    # a blunt end's boundary must lie in the end plane exactly (the cap is planar)
    for idx, pointed, x in ((0, pointed0, x0), (-1, pointed1, x1)):
        if not pointed:
            C[idx] += np.outer(x - (C[idx] - frame.origin) @ frame.ez, frame.ez)

    # how closely the surface passes through the contour points
    # (across the surface: the gap along the local normal, so sliding along the surface is not counted)
    fitted = np.einsum("ia,abk,jb->ijk", Bu, C, Bw)
    tw = np.roll(fitted, -1, axis=1) - np.roll(fitted, 1, axis=1)
    tu = np.gradient(fitted, axis=0)
    nrm = np.cross(tu, tw)
    nl = np.linalg.norm(nrm, axis=2, keepdims=True)
    gap = fitted - P
    d = np.where(nl[..., 0] > 1e-15, np.abs(np.sum(gap * nrm, axis=2)) / np.maximum(nl[..., 0], 1e-300),
                 np.linalg.norm(gap, axis=2))
    residual = (float(np.sqrt(np.mean(d ** 2))), float(d.max()))

    patches = split_faces(C, ku, ring, faces)
    start = End(x0, pointed0, tuple(P[0].mean(axis=0)))
    end = End(x1, pointed1, tuple(P[-1].mean(axis=0)))
    return Loft(patches, start, end, frame, notes, n_rows, residual, tuple(zip(z.tolist(), d.max(axis=1).tolist())),
                creases)


def split_faces(C, ku, ring, faces=4):
    """Cut the periodic surface (control net C, one control point per knot of ``ring`` round the body) into ``faces``
    clamped patches at w = k / faces: unroll it, insert each cut knot until it is triple (a cubic then passes through
    a control point there), and take the control points between the cuts."""
    # unrolled one control point further than one turn
    kw = unrolled(ring, extra=1)
    Q = np.concatenate([C, C[:, :DEGREE + 1]], axis=1)
    cuts = []
    for k in range(faces + 1):   # the knots' own values
        t = k / faces
        near = kw[np.argmin(np.abs(kw - t))]
        if abs(near - t) > 1e-9:
            raise ValueError(f"w = {t} is not a knot round the body")
        cuts.append(float(near))
    for t in cuts:
        while int(np.sum(np.abs(kw - t) < 1e-12)) < DEGREE:
            kw, Q = insert_knot(kw, Q, t, axis=1)
    patches = []
    for a, b in zip(cuts, cuts[1:]):
        s = int(np.searchsorted(kw, a - 1e-12))   # first copy of a
        e = int(np.searchsorted(kw, b - 1e-12))   # first copy of b
        ctrl = Q[:, s - 1:e].copy()
        knots = np.concatenate([[a], kw[s:e + DEGREE], [b]])
        patches.append(Patch(ctrl, ku.copy(), knots, (a, b)))
    # neighbours share their edge exactly
    for p, q in zip(patches, patches[1:] + patches[:1]):
        q.ctrl[:, 0] = p.ctrl[:, -1]
    return patches
