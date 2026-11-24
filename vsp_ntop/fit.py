"""Fit OpenVSP fuselage and wing definitions to the section contours of an implicit body."""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from . import slices as S


class FitError(ValueError):
    pass


def cluster(t, a, b):
    """OpenVSP's clustering function: 0..1 -> 0..1, spacing scaled by ``a`` at the start and ``b`` at the end
    (1: uniform; smaller: finer there)."""
    mt = 1.0 - t
    return mt * mt * t * a + mt * t * t * (3.0 - b) + t * t * t


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def unit(a):
    n = math.sqrt(sum(c * c for c in a))
    if n < 1e-12:
        raise FitError("zero-length axis")
    return tuple(c / n for c in a)


@dataclass(frozen=True)
class Frame:
    """The slicing frame: contour point (u, v) on the layer at height z is origin + u ex + v ey + z ez."""

    origin: tuple
    ex: tuple
    ey: tuple

    @property
    def ez(self):
        return cross(self.ex, self.ey)

    def point(self, u, v, z):
        ez = self.ez
        return tuple(self.origin[i] + u * self.ex[i] + v * self.ey[i] + z * ez[i] for i in range(3))


def fuselage_frame(origin):
    return Frame(tuple(origin), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def wing_frame(origin, span_axis):
    s = unit(span_axis)
    if abs(s[0]) > 1e-6:
        raise FitError("the span axis must be normal to x (the aft axis)")
    x = (1.0, 0.0, 0.0)
    return Frame(tuple(origin), cross(x, s), x)


# --- Fuselage

@dataclass
class FuseSection:
    x: float                 # station along the slicing direction (frame z)
    yc: float = 0.0          # centre of the section's bounding box (frame u, v)
    zc: float = 0.0
    width: float = 0.0
    height: float = 0.0
    points: list = field(default_factory=list)   # [(u/width, v/height)] about the centre, closed (first = last);
    params: list = field(default_factory=list)   # curve parameter 0..1 of each point
    loop: list = field(default_factory=list)     # the source contour (frame u, v), for the deviation check
    tip_angle: float = None   # a rounded point's skinning angle (degrees from the axis: 90 nose, -90 tail)
    slopes: tuple = None      # d(extent)/dx of the top

    @property
    def is_point(self):
        return not self.points


@dataclass
class FuselageFit:
    frame: Frame
    sections: list
    notes: list


def _section_curve(loop, per_quarter):
    """Closed curve through a contour in OpenVSP's order: right-hand point, down through the bottom, round the left
    side and over the top. The four quarters meet at the widest points (right, left) and on the centreline (bottom,
    top)."""
    u0, u1, v0, v1 = S.bounds(loop)
    w, h = u1 - u0, v1 - v0
    if w <= 0 or h <= 0:
        raise FitError("a cross-section has no width or height")
    uc, vc = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
    pts = list(loop) if S.area(loop) < 0 else list(reversed(loop))  # clockwise in (y, z)
    n = len(pts)

    def extreme(sign):
        best = max(sign * p[0] for p in pts)
        tol = 0.004 * w
        idx = [i for i in range(n) if sign * pts[i][0] >= best - tol]
        # the tied run may wrap around the end of the list.
        runs, run = [], [idx[0]]
        for i in idx[1:]:
            if i == run[-1] + 1:
                run.append(i)
            else:
                runs.append(run)
                run = [i]
        runs.append(run)
        if len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == n - 1:
            runs[0] = runs.pop() + runs[0]
        longest = max(runs, key=len)
        return longest[len(longest) // 2]

    def centre_crossing(start, stop):
        """Index after which the contour crosses u = uc, walking forward from ``start`` to ``stop``."""
        i = start
        while i != stop:
            j = (i + 1) % n
            a, b = pts[i][0] - uc, pts[j][0] - uc
            if a == 0.0 or (a > 0) != (b > 0):
                return i
            i = j
        return start

    right, left = extreme(1.0), extreme(-1.0)
    bottom, top = centre_crossing(right, left), centre_crossing(left, right)

    def cut(i):
        j = (i + 1) % n
        a, b = pts[i], pts[j]
        d = b[0] - a[0]
        t = 0.0 if abs(d) < 1e-300 else min(1.0, max(0.0, (uc - a[0]) / d))
        return (uc, a[1] + t * (b[1] - a[1]))

    def walk(i, j):
        out = []
        k = i
        while True:
            out.append(pts[k])
            if k == j:
                return out
            k = (k + 1) % n

    p_bottom, p_top = cut(bottom), cut(top)
    quarters = [
        walk(right, bottom) + [p_bottom],
        [p_bottom] + walk((bottom + 1) % n, left),
        walk(left, top) + [p_top],
        [p_top] + walk((top + 1) % n, right),
    ]
    points, params = [], []
    for q, poly in enumerate(quarters):
        sampled = S.resample(poly, per_quarter + 1)
        for k, p in enumerate(sampled[:-1]):
            points.append(((p[0] - uc) / w, (p[1] - vc) / h))
            params.append((q + k / per_quarter) / 4.0)
    points.append(points[0])
    params.append(1.0)
    return uc, vc, w, h, points, params


ROUND_TIP = 0.75  # a tip whose sqrt(area) grows as (x - x_tip)^p with p below this is rounded (sphere: 0.5)


def _tip(layers, areas, order, spacing):
    """Tip position and shape of a pointed end from the slices nearest it: sqrt(area) = c |x - x_tip|^p fitted over
    the slices with less than a quarter of the largest area (3 to 8 of them), p searched from 0.35 to 1.3."""
    big = max(areas)
    idx = [i for i in order[:8] if areas[i] < 0.25 * big] or list(order[:3])
    if len(idx) < 3:
        idx = list(order[:3])
    z = [layers[i].z for i in idx]
    r = [math.sqrt(areas[i]) for i in idx]
    best = None
    for k in range(96):
        p = 0.35 + 0.01 * k
        y = [ri ** (1.0 / p) for ri in r]
        n = len(z)
        mz, my = sum(z) / n, sum(y) / n
        szz = sum((a - mz) ** 2 for a in z)
        if szz <= 0.0:
            continue
        slope = sum((a - mz) * (b - my) for a, b in zip(z, y)) / szz
        if slope == 0.0:
            continue
        tip = mz - my / slope
        # residual in sqrt(area), relative to its size
        res = sum((abs(slope * (a - tip)) ** p - ri) ** 2 for a, ri in zip(z, r))
        if best is None or res < best[0]:
            best = (res, tip, p)
    z0 = layers[order[0]].z
    at_start = order[0] < order[-1]
    lo, hi = (z0 - spacing, z0) if at_start else (z0, z0 + spacing)
    if best is None:
        return min(hi, max(lo, z0)), 1.0
    return min(hi, max(lo, best[1])), best[2]


def _end(layers, areas, at_start, spacing, bound=None):
    """Where the body ends beyond its first (or last) layer, whether it ends in a point, and the point's shape
    exponent p (sqrt(area) ~ distance^p: 1 for a cone or ogive, 0.5 for a sphere)."""
    order = list(range(len(layers))) if at_start else list(range(len(layers) - 1, -1, -1))
    i0, i1 = order[0], order[min(3, len(order) - 1)]
    r0, r1 = math.sqrt(areas[i0]), math.sqrt(areas[i1])
    z0 = layers[i0].z
    blunt = areas[i0] > 0.03 * max(areas) and not r1 > 1.6 * r0
    if blunt or r1 <= r0:
        sign = -1.0 if at_start else 1.0
        reach = spacing  # the next slice out is empty
        if bound is not None and 0.0 <= sign * (bound - z0) < reach:
            reach = sign * (bound - z0)  # nTop's bounding box contains the body (it may be loose)
        return z0 + sign * 0.5 * reach, False, None
    tip, p = _tip(layers, areas, order, spacing)
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
    spacing = (used[-1].z - used[0].z) / (len(used) - 1)
    areas = [abs(S.area(l.loops[0])) for l in used]
    b0, b1 = bounds if bounds else (None, None)
    x_start, pointed_start, p_start = _end(used, areas, True, spacing, b0)
    x_end, pointed_end, p_end = _end(used, areas, False, spacing, b1)
    n_sections = max(4, int(n_sections))
    n = len(used)
    zs = [l.z for l in used]
    ext = []
    for l in used:
        u0, u1, v0, v1 = S.bounds(l.loops[0])
        ext.append((v1, u1, -v0, -u0))  # top, right, bottom, left

    def slope(i):
        a, b = max(i - 1, 0), min(i + 1, n - 1)
        return tuple((ext[b][k] - ext[a][k]) / (zs[b] - zs[a]) for k in range(4)) if zs[b] > zs[a] else (0.0,) * 4

    def centre(i):
        u0, u1, v0, v1 = S.bounds(used[i].loops[0])
        return 0.5 * (u0 + u1), 0.5 * (v0 + v1)

    # End nodes: a pointed end is a point at the tip.
    def end_node(at_start):
        i = 0 if at_start else n - 1
        x, pointed, p_tip = (x_start, pointed_start, p_start) if at_start else (x_end, pointed_end, p_end)
        if not pointed:
            return x, ext[i], None
        uc, vc = centre(i)
        return x, (vc, uc, -vc, -uc), p_tip

    nodes = {-1: end_node(True), n: end_node(False)}  # virtual indices for the two ends
    blunt0, blunt1 = not pointed_start, not pointed_end
    candidates = [i for i in range(n) if not (blunt0 and i == 0) and not (blunt1 and i == n - 1)]
    budget = n_sections - 2
    # half the stations as before, clustered toward the ends ...
    chosen = set()
    for k in range(1, budget // 2 + 1):
        x = x_start + (x_end - x_start) * cluster(k / (budget // 2 + 1), nose_cluster, tail_cluster)
        chosen.add(min(candidates, key=lambda j: abs(zs[j] - x)))

    def node(j):
        return (nodes[j][0], nodes[j][1]) if j in nodes else (zs[j], ext[j])

    # ... and the rest where interpolating the outline between the chosen stations misses it most
    while len(chosen) < min(budget, len(candidates)):
        order = sorted(chosen | {-1, n}, key=lambda j: node(j)[0])
        worst, pick = 0.0, None
        for a, b in zip(order, order[1:]):
            (xa, ea), (xb, eb) = node(a), node(b)
            for j in candidates:
                if j in chosen or not xa < zs[j] < xb:
                    continue
                t = (zs[j] - xa) / (xb - xa)
                err = max(abs(ext[j][q] - (ea[q] + t * (eb[q] - ea[q]))) for q in range(4))
                if err > worst:
                    worst, pick = err, j
        if pick is None:
            break
        chosen.add(pick)
    sections = []
    x0, e0, p0 = nodes[-1]
    if pointed_start:
        uc, vc = centre(0)
        rounded = p0 < ROUND_TIP
        tip_slopes = tuple((ext[0][q] - e0[q]) / max(zs[0] - x0, 1e-12) for q in range(4))
        sections.append(FuseSection(x0, uc, vc, tip_angle=90.0 if rounded else None, slopes=tip_slopes))
    else:
        uc, vc, w, h, pts, params = _section_curve(used[0].loops[0], per_quarter)
        sections.append(FuseSection(x0, uc, vc, w, h, pts, params, used[0].loops[0], slopes=slope(0)))
    for i in sorted(chosen):
        uc, vc, w, h, pts, params = _section_curve(used[i].loops[0], per_quarter)
        sections.append(FuseSection(zs[i], uc, vc, w, h, pts, params, used[i].loops[0], slopes=slope(i)))
    x1, e1, p1 = nodes[n]
    if pointed_end:
        uc, vc = centre(n - 1)
        rounded = p1 < ROUND_TIP
        tip_slopes = tuple((e1[q] - ext[n - 1][q]) / max(x1 - zs[n - 1], 1e-12) for q in range(4))
        sections.append(FuseSection(x1, uc, vc, tip_angle=-90.0 if rounded else None, slopes=tip_slopes))
    else:
        uc, vc, w, h, pts, params = _section_curve(used[-1].loops[0], per_quarter)
        sections.append(FuseSection(x1, uc, vc, w, h, pts, params, used[-1].loops[0], slopes=slope(n - 1)))
    # stations must increase.
    clean = [sections[0]]
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
    le: tuple           # leading edge (n, x) in the frame
    chord: float
    twist: float        # degrees, leading edge up positive
    upper: list         # [(x/c, y/c)] leading edge to trailing edge
    lower: list
    thickness: float    # t/c
    te_thickness: float  # trailing-edge base height / chord
    loop: list = field(default_factory=list)


@dataclass
class WingFit:
    frame: Frame
    stations: list
    notes: list


def _farthest(pts, p):
    return max(range(len(pts)), key=lambda i: (pts[i][0] - p[0]) ** 2 + (pts[i][1] - p[1]) ** 2)


def _airfoil(loop, n_points):
    """Chord line, twist and upper/lower surfaces of a section contour (u = wing normal, v = aft)."""
    n = len(loop)
    if n < 12:
        raise FitError("an airfoil contour has too few points")
    spacing = S.perimeter(loop) / n
    a = _farthest(loop, S.centroid(loop))
    b = _farthest(loop, loop[a])
    a = _farthest(loop, loop[b])
    te, le = (a, b) if loop[a][1] > loop[b][1] else (b, a)
    te_mid, le_pt = loop[te], loop[le]
    up_i = low_i = te
    for _ in range(4):
        d = (te_mid[0] - le_pt[0], te_mid[1] - le_pt[1])
        chord = math.hypot(*d)
        e = (d[0] / chord, d[1] / chord)
        nrm = (e[1], -e[0])  # normal to the chord with a positive u (upward) part
        proj = [(p[0] - le_pt[0]) * e[0] + (p[1] - le_pt[1]) * e[1] for p in loop]
        perp = [(p[0] - le_pt[0]) * nrm[0] + (p[1] - le_pt[1]) * nrm[1] for p in loop]
        t_max = max(proj)
        eps = max(0.002 * chord, 1.5 * spacing)
        base = [i for i in range(n) if proj[i] >= t_max - eps]
        up_i, low_i = max(base, key=lambda i: perp[i]), min(base, key=lambda i: perp[i])
        te_mid = (0.5 * (loop[up_i][0] + loop[low_i][0]), 0.5 * (loop[up_i][1] + loop[low_i][1]))
        le_pt = loop[_farthest(loop, te_mid)]
    le = _farthest(loop, te_mid)
    d = (te_mid[0] - le_pt[0], te_mid[1] - le_pt[1])
    chord = math.hypot(*d)
    e = (d[0] / chord, d[1] / chord)
    nrm = (e[1], -e[0])
    te_thick = math.dist(loop[up_i], loop[low_i]) / chord
    if te_thick < max(0.0008, 2.5 * spacing / chord):  # a sharp trailing edge
        te_thick = 0.0

    def local(p):
        return (((p[0] - le_pt[0]) * e[0] + (p[1] - le_pt[1]) * e[1]) / chord,
                ((p[0] - le_pt[0]) * nrm[0] + (p[1] - le_pt[1]) * nrm[1]) / chord)

    def path(i, j, step):
        out = []
        k = i
        while True:
            out.append(loop[k])
            if k == j:
                return out
            k = (k + step) % n

    fwd_up = sum(local(loop[(le + k) % n])[1] for k in range(1, max(2, n // 8)))
    step_up = 1 if fwd_up > 0 else -1
    upper = [local(p) for p in path(le, up_i, step_up)]
    lower = [local(p) for p in path(le, low_i, -step_up)]
    if te_thick == 0.0:
        upper[-1] = lower[-1] = (1.0, 0.0)
    fractions = [0.5 * (1.0 - math.cos(math.pi * k / (n_points - 1))) for k in range(n_points)]
    upper, lower = S.at_fractions(upper, fractions), S.at_fractions(lower, fractions)
    upper[0] = lower[0] = (0.0, 0.0)
    twist = math.degrees(math.atan2(-e[0], e[1]))
    thickness = max(p[1] for p in upper) - min(p[1] for p in lower)
    return le_pt, chord, twist, upper, lower, thickness, te_thick


def _select(metric, n_max, tolerance, min_gap=0.0):
    """Indices of at most ``n_max`` stations such that linear interpolation between them stays within ``tolerance`` of
    every row of ``metric`` (greedy: add the worst station until done)."""
    keep = [0, len(metric) - 1]
    while len(keep) < n_max:
        worst, where = tolerance, None
        for a, b in zip(keep, keep[1:]):
            for k in range(a + 1, b):
                if metric[k][0] - metric[a][0] < min_gap or metric[b][0] - metric[k][0] < min_gap:
                    continue
                t = (metric[k][0] - metric[a][0]) / (metric[b][0] - metric[a][0])
                dev = max(abs(metric[k][c] - (metric[a][c] + t * (metric[b][c] - metric[a][c])))
                          for c in range(1, len(metric[k])))
                if dev > worst:
                    worst, where = dev, k
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
            le, chord, twist, upper, lower, thick, te = _airfoil(outer[0], airfoil_points)
        except FitError:
            continue
        stations.append(WingStation(l.z, le, chord, twist, upper, lower, thick, te, outer[0]))
    if len(stations) < 2:
        raise FitError(f"only {len(stations)} span stations cut the wing; use more slices")
    # A rounded or faired tip: stop where the chord collapses.
    root_chord = max(s.chord for s in stations[: max(2, len(stations) // 4)])
    while len(stations) > 2 and stations[-1].chord < 0.5 * stations[-2].chord:
        stations.pop()
    ref = root_chord
    metric = [(s.z, s.le[1] / ref, s.le[0] / ref, s.chord / ref, math.radians(s.twist) * s.chord / ref,
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
    le = tuple(p + t * (q - p) for p, q in zip(a.le, b.le))
    thick = a.thickness + t * (b.thickness - a.thickness)
    fit.stations.insert(0, WingStation(z0, le, chord, a.twist, list(a.upper), list(a.lower), thick, a.te_thickness, []))
    fit.notes.append(f"extended into the body: root moved {exposed:.4f} m inboard, chord "
                     f"{chord:.4f} m there")
    return exposed


def wing_sections(fit: WingFit):
    """OpenVSP's section parameters between consecutive stations: span (in the section's own plane), leading
    edge sweep and dihedral in degrees, root and tip chords, tip twist."""
    out = []
    for a, b in zip(fit.stations, fit.stations[1:]):
        dy, dz, dx = b.z - a.z, b.le[0] - a.le[0], b.le[1] - a.le[1]
        span = math.hypot(dy, dz)
        out.append({"span": span, "sweep": math.degrees(math.atan2(dx, span)),
                    "dihedral": math.degrees(math.atan2(dz, dy)), "root_chord": a.chord, "tip_chord": b.chord,
                    "twist": b.twist})
    return out
