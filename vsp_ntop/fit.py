"""Fit OpenVSP fuselage and wing definitions to the section contours of an implicit body."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import slices as S
from .geometry import unit


class FitError(ValueError):
    pass


def cluster(t, a, b):
    """OpenVSP's clustering function: 0..1 -> 0..1, spacing scaled by ``a`` at the start and ``b`` at the end
    (1: uniform; smaller: finer there)."""
    mt = 1.0 - t
    return mt * mt * t * a + mt * t * t * (3.0 - b) + t * t * t


@dataclass(frozen=True, eq=False)
class Frame:
    """The slicing frame: contour point (u, v) on the layer at height z is origin + u ex + v ey + z ez."""

    origin: np.ndarray
    ex: np.ndarray
    ey: np.ndarray

    @property
    def ez(self):
        return np.cross(self.ex, self.ey)

    def point(self, u, v, z):
        """World point(s) of frame coordinates (scalars or arrays of one shape)."""
        u, v, z = (np.asarray(c, float)[..., None] for c in (u, v, z))
        return self.origin + u * self.ex + v * self.ey + z * self.ez


def fuselage_frame(origin):
    return Frame(np.asarray(origin, float), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))


def wing_frame(origin, span_axis):
    s = unit(span_axis, FitError, "axis")
    if abs(s[0]) > 1e-6:
        raise FitError("the span axis must be normal to x (the aft axis)")
    x = np.array([1.0, 0.0, 0.0])
    return Frame(np.asarray(origin, float), np.cross(x, s), x)


def _cycle(n, i, j, step=1):
    """Indices from ``i`` to ``j`` (inclusive) around a closed loop of ``n`` points."""
    return (i + step * np.arange((step * (j - i)) % n + 1)) % n


# --- Fuselage

@dataclass
class FuseSection:
    x: float                 # station along the slicing direction (frame z)
    yc: float = 0.0          # centre of the section's bounding box (frame u, v)
    zc: float = 0.0
    width: float = 0.0
    height: float = 0.0
    points: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))  # (u/width, v/height) about the centre,
    params: np.ndarray = field(default_factory=lambda: np.empty(0))      # closed; curve parameter 0..1 of each
    loop: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))   # the source contour (frame u, v)
    tip_angle: float = None   # a rounded point's skinning angle (degrees from the axis: 90 nose, -90 tail)
    slopes: tuple = None      # d(extent)/dx of the top, right, bottom and left

    @property
    def is_point(self):
        return len(self.points) == 0


@dataclass
class FuselageFit:
    frame: Frame
    sections: list
    notes: list


def _section_curve(loop, per_quarter):
    """Closed curve through a contour in OpenVSP's order: right-hand point, down through the bottom, round the left
    side and over the top. The four quarters meet at the widest points (right, left) and on the centreline (bottom,
    top). Returns centre, width, height, points about the centre in units of width and height, and parameters."""
    u0, u1, v0, v1 = S.bounds(loop)
    w, h = u1 - u0, v1 - v0
    if w <= 0 or h <= 0:
        raise FitError("a cross-section has no width or height")
    uc, vc = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
    pts = np.asarray(loop) if S.area(loop) < 0 else np.asarray(loop)[::-1]  # clockwise in (y, z)
    n = len(pts)
    u = pts[:, 0]

    def extreme(sign):
        """Middle of the longest run of points within 0.4% of the width of the extreme (runs may wrap)."""
        idx = np.flatnonzero(sign * u >= (sign * u).max() - 0.004 * w)
        runs = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1)
        if len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == n - 1:
            runs[0] = np.concatenate([runs.pop(), runs[0]])
        longest = max(runs, key=len)
        return int(longest[len(longest) // 2])

    def centre_crossing(start, stop):
        """Index after which the contour crosses u = uc, walking forward from ``start`` to ``stop``."""
        i = _cycle(n, start, stop)[:-1]
        a, b = u[i] - uc, u[(i + 1) % n] - uc
        hit = np.flatnonzero((a == 0.0) | ((a > 0) != (b > 0)))
        return int(i[hit[0]]) if len(hit) else start

    def cut(i):
        a, b = pts[i], pts[(i + 1) % n]
        d = b[0] - a[0]
        t = 0.0 if abs(d) < 1e-300 else np.clip((uc - a[0]) / d, 0.0, 1.0)
        return np.array([uc, a[1] + t * (b[1] - a[1])])

    right, left = extreme(1.0), extreme(-1.0)
    bottom, top = centre_crossing(right, left), centre_crossing(left, right)
    p_bottom, p_top = cut(bottom), cut(top)
    quarters = [np.vstack([pts[_cycle(n, right, bottom)], p_bottom]),
                np.vstack([p_bottom, pts[_cycle(n, (bottom + 1) % n, left)]]),
                np.vstack([pts[_cycle(n, left, top)], p_top]),
                np.vstack([p_top, pts[_cycle(n, (top + 1) % n, right)]])]
    points = np.vstack([S.resample(q, per_quarter + 1)[:-1] for q in quarters])
    points = (points - [uc, vc]) / [w, h]
    points = np.vstack([points, points[:1]])
    params = np.append(np.arange(4 * per_quarter) / (4.0 * per_quarter), 1.0)
    return uc, vc, w, h, points, params


ROUND_TIP = 0.75  # a tip whose sqrt(area) grows as (x - x_tip)^p with p below this is rounded (sphere: 0.5)


def _tip(zs, areas, order, spacing):
    """Tip position and shape of a pointed end from the slices nearest it: sqrt(area) = c |x - x_tip|^p fitted over
    the slices with less than a quarter of the largest area (3 to 8 of them), p searched from 0.35 to 1.3."""
    idx = [i for i in order[:8] if areas[i] < 0.25 * areas.max()]
    if len(idx) < 3:
        idx = order[:3]
    z, r = zs[idx], np.sqrt(areas[idx])
    p = 0.35 + 0.01 * np.arange(96)[:, None]
    y = r ** (1.0 / p)
    dz = z - z.mean()
    szz = dz @ dz
    z0 = zs[order[0]]
    lo, hi = (z0 - spacing, z0) if order[0] < order[-1] else (z0, z0 + spacing)
    if szz <= 0.0:
        return np.clip(z0, lo, hi), 1.0
    slope = (y - y.mean(axis=1, keepdims=True)) @ dz / szz
    ok = slope != 0.0
    if not ok.any():
        return np.clip(z0, lo, hi), 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        tip = z.mean() - y.mean(axis=1) / slope
        res = ((np.abs(slope[:, None] * (z - tip[:, None])) ** p - r) ** 2).sum(axis=1)  # in sqrt(area)
    k = np.argmin(np.where(ok & np.isfinite(res), res, np.inf))
    return np.clip(tip[k], lo, hi), float(p[k, 0])


def _end(zs, areas, at_start, spacing, bound=None):
    """Where the body ends beyond its first (or last) layer, whether it ends in a point, and the point's shape
    exponent p (sqrt(area) ~ distance^p: 1 for a cone or ogive, 0.5 for a sphere)."""
    order = np.arange(len(zs)) if at_start else np.arange(len(zs))[::-1]
    i0, i1 = order[0], order[min(3, len(order) - 1)]
    r0, r1 = np.sqrt(areas[i0]), np.sqrt(areas[i1])
    z0 = zs[i0]
    blunt = areas[i0] > 0.03 * areas.max() and not r1 > 1.6 * r0
    if blunt or r1 <= r0:
        sign = -1.0 if at_start else 1.0
        reach = spacing  # the next slice out is empty
        if bound is not None and 0.0 <= sign * (bound - z0) < reach:
            reach = sign * (bound - z0)  # nTop's bounding box contains the body (it may be loose)
        return z0 + sign * 0.5 * reach, False, None
    tip, p = _tip(zs, areas, order, spacing)
    return tip, True, p


def fit_fuselage(layers, frame: Frame, n_sections=24, nose_cluster=0.5, tail_cluster=0.5, per_quarter=10, bounds=None):
    """``bounds``: the body's bounding-box extent along the slicing direction (frame z), if known."""
    notes = []
    used = []
    for l in layers:
        outer = S.outer_loops(l)
        if not outer:
            continue
        if len(outer) > 1:
            notes.append(f"x = {l.z:.4f} m: {len(outer)} separate contours; the largest is used")
        used.append(S.Layer(l.z, [outer[0]]))
    if len(used) < 4:
        raise FitError(f"only {len(used)} sections cut the body; use more slices")
    n = len(used)
    zs = np.array([l.z for l in used])
    spacing = (zs[-1] - zs[0]) / (n - 1)
    areas = np.array([abs(S.area(l.loops[0])) for l in used])
    b0, b1 = bounds if bounds else (None, None)
    x_start, pointed_start, p_start = _end(zs, areas, True, spacing, b0)
    x_end, pointed_end, p_end = _end(zs, areas, False, spacing, b1)
    n_sections = max(4, int(n_sections))
    box = np.array([S.bounds(l.loops[0]) for l in used])  # u0, u1, v0, v1
    ext = np.column_stack([box[:, 3], box[:, 1], -box[:, 2], -box[:, 0]])  # top, right, bottom, left
    centre = np.column_stack([box[:, :2].mean(axis=1), box[:, 2:].mean(axis=1)])
    a, b = np.maximum(np.arange(n) - 1, 0), np.minimum(np.arange(n) + 1, n - 1)
    slopes = (ext[b] - ext[a]) / (zs[b] - zs[a])[:, None]  # central differences, one-sided at the ends

    def tip_extent(i):
        uc, vc = centre[i]
        return np.array([vc, uc, -vc, -uc])

    # nodes -1 and n are the two ends: a pointed end is a point at the tip
    node_x = {-1: x_start, n: x_end}
    node_e = {-1: tip_extent(0) if pointed_start else ext[0], n: tip_extent(n - 1) if pointed_end else ext[n - 1]}
    candidates = np.arange(int(not pointed_start), n - int(not pointed_end))
    budget = n_sections - 2
    # half the stations clustered toward the ends ...
    half = budget // 2
    targets = x_start + (x_end - x_start) * cluster(np.arange(1, half + 1) / (half + 1), nose_cluster, tail_cluster)
    chosen = {int(candidates[np.argmin(np.abs(zs[candidates] - x))]) for x in targets}
    # ... and the rest where interpolating the outline between the chosen stations misses it most
    while len(chosen) < min(budget, len(candidates)):
        nodes = sorted(chosen | {-1, n}, key=lambda j: node_x[j] if j in node_x else zs[j])
        nx = np.array([node_x[j] if j in node_x else zs[j] for j in nodes])
        ne = np.array([node_e[j] if j in node_e else ext[j] for j in nodes])
        free = np.array([j for j in candidates if j not in chosen], int)
        k = np.searchsorted(nx, zs[free], side="right") - 1
        between = (k >= 0) & (k < len(nx) - 1)
        free, k = free[between], k[between]
        between = (nx[k] < zs[free]) & (zs[free] < nx[k + 1])
        free, k = free[between], k[between]
        if not len(free):
            break
        t = ((zs[free] - nx[k]) / (nx[k + 1] - nx[k]))[:, None]
        err = np.abs(ext[free] - (ne[k] + t * (ne[k + 1] - ne[k]))).max(axis=1)
        if err.max() <= 0.0:
            break
        chosen.add(int(free[np.argmax(err)]))

    def section(i, x=None):
        uc, vc, w, h, pts, params = _section_curve(used[i].loops[0], per_quarter)
        return FuseSection(zs[i] if x is None else x, uc, vc, w, h, pts, params, used[i].loops[0],
                           slopes=tuple(slopes[i]))

    def tip(i, x, p, angle, sign):
        uc, vc = centre[i]
        s = sign * (ext[i] - node_e[-1 if sign > 0 else n]) / max(sign * (zs[i] - x), 1e-12)
        return FuseSection(x, uc, vc, tip_angle=angle if p < ROUND_TIP else None, slopes=tuple(s))

    sections = [tip(0, x_start, p_start, 90.0, 1.0) if pointed_start else section(0, x_start)]
    sections += [section(i) for i in sorted(chosen)]
    sections.append(tip(n - 1, x_end, p_end, -90.0, -1.0) if pointed_end else section(n - 1, x_end))
    clean = [sections[0]]  # stations must increase
    for s in sections[1:]:
        if s.x > clean[-1].x + 1e-9:
            clean.append(s)
    if len(clean) < 3:
        raise FitError("fewer than 3 usable sections")

    def kind(sec):
        return "blunt" if not sec.is_point else ("rounded" if sec.tip_angle is not None else "pointed")
    notes.append(f"{len(clean)} sections from x = {clean[0].x:.4f} to {clean[-1].x:.4f} m "
                 f"({kind(clean[0])} nose, {kind(clean[-1])} tail); "
                 f"{len(used)} slices {spacing * 1000:.1f} mm apart")
    return FuselageFit(frame, clean, notes)


# --- Wing

@dataclass
class WingStation:
    z: float            # along the span axis (frame z)
    le: np.ndarray      # leading edge (n, x) in the frame
    chord: float
    twist: float        # degrees, leading edge up positive
    upper: np.ndarray   # (x/c, y/c) leading edge to trailing edge
    lower: np.ndarray
    thickness: float    # t/c
    te_thickness: float  # trailing-edge base height / chord
    loop: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))


@dataclass
class WingFit:
    frame: Frame
    stations: list
    notes: list


def _farthest(pts, p):
    return int(np.argmax(((pts - p) ** 2).sum(axis=1)))


def _airfoil(loop, n_points):
    """Chord line, twist and upper/lower surfaces of a section contour (u = wing normal, v = aft)."""
    loop = np.asarray(loop)
    n = len(loop)
    if n < 12:
        raise FitError("an airfoil contour has too few points")
    spacing = S.perimeter(loop) / n
    a = _farthest(loop, S.centroid(loop))
    b = _farthest(loop, loop[a])
    a = _farthest(loop, loop[b])
    te, le = (a, b) if loop[a, 1] > loop[b, 1] else (b, a)
    te_mid, le_pt = loop[te], loop[le]

    def chord_axes():
        d = te_mid - le_pt
        chord = np.hypot(*d)
        e = d / chord
        return chord, e, np.array([e[1], -e[0]])  # chord, along it, normal with a positive u (upward) part

    up_i = low_i = te
    for _ in range(4):
        chord, e, nrm = chord_axes()
        proj, perp = (loop - le_pt) @ e, (loop - le_pt) @ nrm
        base = np.flatnonzero(proj >= proj.max() - max(0.002 * chord, 1.5 * spacing))
        up_i, low_i = base[np.argmax(perp[base])], base[np.argmin(perp[base])]
        te_mid = 0.5 * (loop[up_i] + loop[low_i])
        le_pt = loop[_farthest(loop, te_mid)]
    le = _farthest(loop, te_mid)
    chord, e, nrm = chord_axes()
    te_thick = np.linalg.norm(loop[up_i] - loop[low_i]) / chord
    if te_thick < max(0.0008, 2.5 * spacing / chord):  # a sharp trailing edge
        te_thick = 0.0
    local = (loop - le_pt) @ np.column_stack([e, nrm]) / chord
    step = 1 if local[(le + np.arange(1, max(2, n // 8))) % n, 1].sum() > 0 else -1
    upper, lower = local[_cycle(n, le, up_i, step)], local[_cycle(n, le, low_i, -step)]
    if te_thick == 0.0:
        upper[-1] = lower[-1] = (1.0, 0.0)
    fractions = 0.5 * (1.0 - np.cos(np.pi * np.arange(n_points) / (n_points - 1)))
    upper, lower = S.at_fractions(upper, fractions), S.at_fractions(lower, fractions)
    upper[0] = lower[0] = (0.0, 0.0)
    twist = np.degrees(np.arctan2(-e[0], e[1]))
    thickness = upper[:, 1].max() - lower[:, 1].min()
    return le_pt, chord, twist, upper, lower, thickness, te_thick


def _select(metric, n_max, tolerance, min_gap=0.0):
    """Indices of at most ``n_max`` rows of ``metric`` (first column the station) such that linear interpolation
    between them stays within ``tolerance`` of every row (greedy: add the worst row until done)."""
    metric = np.asarray(metric)
    x = metric[:, 0]
    keep = [0, len(metric) - 1]
    while len(keep) < n_max:
        worst, where = tolerance, None
        for a, b in zip(keep, keep[1:]):
            k = np.arange(a + 1, b)
            k = k[(x[k] - x[a] >= min_gap) & (x[b] - x[k] >= min_gap)]
            if not len(k):
                continue
            t = ((x[k] - x[a]) / (x[b] - x[a]))[:, None]
            dev = np.abs(metric[k, 1:] - (metric[a, 1:] + t * (metric[b, 1:] - metric[a, 1:]))).max(axis=1)
            if dev.max() > worst:
                worst, where = dev.max(), int(k[np.argmax(dev)])
        if where is None:
            break
        keep = sorted(keep + [where])
    return keep


def fit_wing(layers, frame: Frame, n_sections=8, airfoil_points=41, tolerance=0.01, min_span=0.04):
    notes = []
    stations = []
    for l in layers:
        outer = S.outer_loops(l)
        if not outer:
            continue
        if len(outer) > 1:
            notes.append(f"span station {l.z:.4f} m: {len(outer)} separate contours; the largest is used")
        try:
            stations.append(WingStation(l.z, *_airfoil(outer[0], airfoil_points), outer[0]))
        except FitError:
            continue
    if len(stations) < 2:
        raise FitError(f"only {len(stations)} span stations cut the wing; use more slices")
    # A rounded or faired tip: stop where the chord collapses.
    ref = max(s.chord for s in stations[: max(2, len(stations) // 4)])
    while len(stations) > 2 and stations[-1].chord < 0.5 * stations[-2].chord:
        stations.pop()
    metric = [(s.z, s.le[1] / ref, s.le[0] / ref, s.chord / ref, np.radians(s.twist) * s.chord / ref,
               s.thickness * s.chord / ref) for s in stations]
    keep = _select(metric, max(2, int(n_sections) + 1), tolerance, min_span * (stations[-1].z - stations[0].z))
    chosen = [stations[i] for i in keep]
    spacing = (stations[-1].z - stations[0].z) / (len(stations) - 1)
    notes.append(f"{len(chosen) - 1} sections over span stations {chosen[0].z:.4f} to {chosen[-1].z:.4f} m; "
                 f"{len(stations)} slices {spacing * 1000:.1f} mm apart")
    return WingFit(frame, chosen, notes)


def extend_wing(fit: WingFit, z0):
    """Extend the wing inboard to span station ``z0`` (between its root and the plane through the origin): a new root
    station whose leading edge, chord and thickness continue the first section's linear variation, with the root
    airfoil and twist."""
    a, b = fit.stations[0], fit.stations[1]
    exposed = a.z - z0
    if exposed <= 1e-9 * max(1.0, abs(a.z)):
        return 0.0  # the wing already reaches the plane
    t = (z0 - a.z) / (b.z - a.z)
    chord = a.chord + t * (b.chord - a.chord)
    if chord <= 0.0:
        raise FitError("extending the wing this far into the body makes its chord vanish")
    fit.stations.insert(0, WingStation(z0, a.le + t * (b.le - a.le), chord, a.twist, a.upper.copy(), a.lower.copy(),
                                       a.thickness + t * (b.thickness - a.thickness), a.te_thickness))
    fit.notes.append(f"extended into the body: root moved {exposed:.4f} m inboard, chord {chord:.4f} m there")
    return exposed


def wing_sections(fit: WingFit):
    """OpenVSP's section parameters between consecutive stations: span (in the section's own plane), leading
    edge sweep and dihedral in degrees, root and tip chords, tip twist."""
    out = []
    for a, b in zip(fit.stations, fit.stations[1:]):
        dy, (dz, dx) = b.z - a.z, b.le - a.le
        span = np.hypot(dy, dz)
        out.append({"span": span, "sweep": np.degrees(np.arctan2(dx, span)),
                    "dihedral": np.degrees(np.arctan2(dz, dy)), "root_chord": a.chord, "tip_chord": b.chord,
                    "twist": b.twist})
    return out
