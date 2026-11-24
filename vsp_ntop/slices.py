"""Section contours of an implicit body, as nTop's Slice Body block exports them (CLI, ASCII)."""
from __future__ import annotations

import math


class Layer:
    """One slicing plane: ``z`` along the slicing direction, ``loops`` closed polylines [(u, v), ...] (first
    point not repeated) with the outside on the right of a clockwise loop (holes have negative area)."""

    def __init__(self, z, loops):
        self.z = z
        self.loops = loops


def read_cli(path):
    """Layers of an ASCII CLI file, in metres, sorted by z. Empty layers are kept (no loops)."""
    unit_mm = 1.0
    layers = []
    with open(path, encoding="ascii", errors="replace") as f:
        for line in f:
            if line.startswith("$$UNITS/"):
                unit_mm = float(line.split("/", 1)[1])
            elif line.startswith("$$BINARY"):
                raise ValueError(f"{path}: binary CLI; export the slices as ASCII")
            elif line.startswith("$$LAYER/"):
                layers.append(Layer(float(line.split("/", 1)[1]) * unit_mm * 1e-3, []))
            elif line.startswith("$$POLYLINE/"):
                v = line.split("/", 1)[1].split(",")
                n = int(v[2])
                s = unit_mm * 1e-3
                pts = [(float(v[3 + 2 * i]) * s, float(v[4 + 2 * i]) * s) for i in range(n)]
                if len(pts) > 1 and pts[0] == pts[-1]:
                    pts.pop()
                if len(pts) >= 3 and layers:
                    layers[-1].loops.append(pts)
    layers.sort(key=lambda l: l.z)
    return layers


def area(loop):
    """Signed area (positive counter-clockwise in u, v)."""
    s = 0.0
    for i in range(len(loop)):
        x0, y0 = loop[i - 1]
        x1, y1 = loop[i]
        s += x0 * y1 - x1 * y0
    return 0.5 * s


def centroid(loop):
    a = cx = cy = 0.0
    for i in range(len(loop)):
        x0, y0 = loop[i - 1]
        x1, y1 = loop[i]
        w = x0 * y1 - x1 * y0
        a += w
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    if abs(a) < 1e-300:
        n = len(loop)
        return sum(p[0] for p in loop) / n, sum(p[1] for p in loop) / n
    return cx / (3.0 * a), cy / (3.0 * a)


def bounds(loop):
    us = [p[0] for p in loop]
    vs = [p[1] for p in loop]
    return min(us), max(us), min(vs), max(vs)


def perimeter(loop):
    return sum(math.dist(loop[i - 1], loop[i]) for i in range(len(loop)))


def outer_loops(layer: Layer, min_fraction=0.02):
    """The layer's outer contours (not holes), largest first, dropping slivers below ``min_fraction`` of the
    largest area."""
    sized = sorted(((abs(area(l)), l) for l in layer.loops), key=lambda t: -t[0])
    if not sized:
        return []
    keep = []
    for a, l in sized:
        if a < min_fraction * sized[0][0]:
            break
        cu, cv = centroid(l)
        if any(inside(o, cu, cv) for o in keep):  # a hole or an island inside a larger contour
            continue
        keep.append(l)
    return keep


def inside(loop, u, v):
    n = False
    for i in range(len(loop)):
        x0, y0 = loop[i - 1]
        x1, y1 = loop[i]
        if (y0 > v) != (y1 > v) and u < x0 + (v - y0) * (x1 - x0) / (y1 - y0):
            n = not n
    return n


def resample(points, n, closed=False):
    """``n`` points evenly spaced by arc length along a polyline (ends included when open)."""
    pts = list(points) + ([points[0]] if closed else [])
    cum = [0.0]
    for a, b in zip(pts, pts[1:]):
        cum.append(cum[-1] + math.dist(a, b))
    total = cum[-1]
    out = []
    k = 0
    count = n if closed else n - 1
    for i in range(n):
        s = total * i / count
        while k < len(cum) - 2 and cum[k + 1] < s:
            k += 1
        span = cum[k + 1] - cum[k]
        t = 0.0 if span <= 0 else min(1.0, max(0.0, (s - cum[k]) / span))
        out.append((pts[k][0] + t * (pts[k + 1][0] - pts[k][0]), pts[k][1] + t * (pts[k + 1][1] - pts[k][1])))
    return out


def at_fractions(points, fractions):
    """Points at the given arc-length fractions (0 to 1) of an open polyline."""
    cum = [0.0]
    for a, b in zip(points, points[1:]):
        cum.append(cum[-1] + math.dist(a, b))
    total = cum[-1]
    out = []
    k = 0
    for f in fractions:
        s = total * min(1.0, max(0.0, f))
        while k < len(cum) - 2 and cum[k + 1] < s:
            k += 1
        while k > 0 and cum[k] > s:
            k -= 1
        span = cum[k + 1] - cum[k]
        t = 0.0 if span <= 0 else (s - cum[k]) / span
        out.append((points[k][0] + t * (points[k + 1][0] - points[k][0]), points[k][1] + t * (points[k + 1][1] - points[k][1])))
    return out
