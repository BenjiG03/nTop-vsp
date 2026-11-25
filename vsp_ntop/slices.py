"""Section contours of an implicit body, as nTop's Slice Body block exports them (CLI, ASCII)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geometry import arc_length, at_arc, shoelace


@dataclass
class Layer:
    """One slicing plane: ``z`` along the slicing direction, ``loops`` closed (n, 2) contours in (u, v) (first point
    not repeated), holes with negative area."""

    z: float
    loops: list = field(default_factory=list)


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
                layers.append(Layer(float(line.split("/", 1)[1]) * unit_mm * 1e-3))
            elif line.startswith("$$POLYLINE/"):
                v = line.split("/", 1)[1].split(",")
                n = int(v[2])
                pts = np.array(v[3:3 + 2 * n], float).reshape(n, 2) * (unit_mm * 1e-3)
                if n > 1 and (pts[0] == pts[-1]).all():
                    pts = pts[:-1]
                if len(pts) >= 3 and layers:
                    layers[-1].loops.append(pts)
    layers.sort(key=lambda l: l.z)
    return layers


def area(loop):
    """Signed area (positive counter-clockwise in u, v)."""
    return shoelace(loop)[0]


def centroid(loop):
    return shoelace(loop)[1]


def bounds(loop):
    p = np.asarray(loop)
    (u0, v0), (u1, v1) = p.min(axis=0), p.max(axis=0)
    return u0, u1, v0, v1


def perimeter(loop):
    return arc_length(loop, closed=True)[-1]


def inside(loop, u, v):
    """Whether (u, v) is inside the closed polygon (even-odd rule)."""
    p = np.asarray(loop)
    x0, y0 = np.roll(p, 1, axis=0).T
    x1, y1 = p.T
    spans = (y0 > v) != (y1 > v)
    dy = np.where(spans, y1 - y0, 1.0)
    return bool(np.count_nonzero(spans & (u < x0 + (v - y0) * (x1 - x0) / dy)) % 2)


def outer_loops(layer: Layer, min_fraction=0.02):
    """The layer's outer contours (not holes), largest first, dropping slivers below ``min_fraction`` of the
    largest area."""
    sized = sorted(((abs(area(l)), l) for l in layer.loops), key=lambda t: -t[0])
    keep = []
    for a, l in sized:
        if a < min_fraction * sized[0][0]:
            break
        if not any(inside(o, *centroid(l)) for o in keep):  # not a hole or island in a larger contour
            keep.append(l)
    return keep


def resample(points, n, closed=False):
    """``n`` points evenly spaced by arc length along a polyline (ends included when open)."""
    total = arc_length(points, closed)[-1]
    return at_arc(points, total * np.arange(n) / (n if closed else n - 1), closed)


def at_fractions(points, fractions):
    """Points at the given arc-length fractions (0 to 1) of an open polyline."""
    return at_arc(points, arc_length(points)[-1] * np.clip(fractions, 0.0, 1.0))
