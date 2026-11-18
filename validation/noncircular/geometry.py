"""The hobby rocket of ``validation/rasaero`` with a noncircular body: a trapezoid section with rounded corners.

The section is a trapezoid 70 mm wide at the bottom, 30 mm at the top and 48 mm high, centred on the axis (so its
sides cross z = 0 at y = +-25 mm), with 6 mm corner radii: area 2367 mm^2 against the round tube's 2290. Along the
body it scales with the round rocket's 5:1 tangent-ogive radius (zero at the tip, full size from x = 270 mm), so the
nose and tube lengths are the round rocket's. It is top-bottom asymmetric: lift and pitching moment at zero angle of
attack, and different at +alpha and -alpha.

The same four fins (160 / 60 mm chords, 70 mm span from the body surface, 100 mm sweep, 4 mm hexagonal section) sit
at +y, +z, -y, -z with their roots on the body: 25 mm from the axis at the sides, 24 mm at the top and bottom.

The body is an nTop field written in its expression language (``set_bounding_box`` makes it an implicit body): with
d_i the distance past each edge of the trapezoid inset by the corner radius, the rounded section is
sqrt(sum max(d_i, 0)^2) + min(max d_i, 0) - r_c (exact for right-angled corners, a smooth fillet otherwise), and the
body is s(x) F(y / s(x), z / s(x)) with s(x) the ogive's radius over 27 mm.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import replace

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))
import rocket as rk  # noqa: E402

BOTTOM, TOP, HEIGHT, CORNER = 70.0, 30.0, 48.0, 6.0  # mm
R_ROUND = 27.0  # the round rocket's radius: the scale of the section along the nose
SIDE_ROOT = 0.25 * (BOTTOM + TOP)  # where the sides cross z = 0
ROOTS = {"+y": SIDE_ROOT, "+z": 0.5 * HEIGHT, "-y": SIDE_ROOT, "-z": 0.5 * HEIGHT}
TURNS = {"+y": 0, "+z": 1, "-y": 2, "-z": 3}  # quarter turns about x from the +y fin
CFD_TAIL = (1240.0, 1440.0)  # the CFD body runs on to here and closes, as hypanel's SU2 rocket


def edges():
    """The trapezoid's edges as (n_y, n_z, c): inside where n . (y, z) <= c (mm)."""
    run, rise = 0.5 * (BOTTOM - TOP), HEIGHT
    k = math.hypot(run, rise)
    ny, nz = rise / k, run / k
    c = ny * 0.5 * BOTTOM - nz * 0.5 * HEIGHT
    return [(0.0, 1.0, 0.5 * HEIGHT), (0.0, -1.0, 0.5 * HEIGHT), (ny, nz, c), (-ny, nz, c)]


def section(y, z, scale=1.0):
    """The section's field at (y, z) mm (negative inside), for the section scaled by ``scale``."""
    d = [ny * y + nz * z - scale * (c - CORNER) for ny, nz, c in edges()]
    return math.sqrt(sum(max(v, 0.0) ** 2 for v in d)) + min(max(d), 0.0) - scale * CORNER


def section_area(n=1600):
    """Section area at full scale, mm^2 (grid count)."""
    h = (BOTTOM + 2.0) / n
    return sum(h * h for i in range(n) for j in range(n)
               if section(-0.5 * BOTTOM - 1 + (i + 0.5) * h, -0.5 * HEIGHT - 1 + (j + 0.5) * h) < 0.0)


def scale(x, cfd=False):
    """s(x): the ogive radius over 27 mm, 1 along the tube; with ``cfd``, falling to 0 over the closed tail."""
    r = rk.R
    rho = (R_ROUND**2 + r.nose_length**2) / (2.0 * R_ROUND)
    xc = min(max(x, 0.0), r.nose_length)
    s = (math.sqrt(rho**2 - (r.nose_length - xc) ** 2) + R_ROUND - rho) / R_ROUND
    if cfd:
        s = min(s, (CFD_TAIL[1] - x) / (CFD_TAIL[1] - CFD_TAIL[0]))
    return s


def body_field(x, y, z, cfd=False):
    """The body's field at (x, y, z) mm, as the nTop expression (for classifying CFD wall faces)."""
    s = max(scale(x, cfd), 1e-4)
    end = CFD_TAIL[1] if cfd else rk.R.total_length
    return max(s * section(y / s, z / s), -x, x - end)


def _expression(cfd):
    """nTop expression of the body field, in millimetres."""
    r = rk.R
    rho = (R_ROUND**2 + r.nose_length**2) / (2.0 * R_ROUND)
    xc = f"min(max(X, 0), {r.nose_length!r})"
    s = f"(sqrt({rho**2!r} - ({r.nose_length!r} - {xc})^2) + {R_ROUND - rho!r}) / {R_ROUND!r}"
    if cfd:
        s = f"min({s}, ({CFD_TAIL[1]!r} - X) / {CFD_TAIL[1] - CFD_TAIL[0]!r})"
    s = f"max({s}, 0.0001)"
    terms = [f"({ny!r} * Y + {nz!r} * Z - S * {c - CORNER!r})" for ny, nz, c in edges()]
    pos = " + ".join(f"max({t}, 0)^2" for t in terms)
    big = f"max(max({terms[0]}, {terms[1]}), max({terms[2]}, {terms[3]}))"
    field = f"sqrt({pos}) + min({big}, 0) - S * {CORNER!r}"
    end = CFD_TAIL[1] if cfd else r.total_length
    expr = f"max(max({field}, -X), X - {end!r})"
    # S is the scale, X, Y, Z coordinates in mm (nTop's x, y, z carry metres).
    expr = expr.replace("S", f"({s})")
    for k, v in (("X", "(x/1mm)"), ("Y", "(y/1mm)"), ("Z", "(z/1mm)")):
        expr = expr.replace(k, v)
    return expr


def body(g, cfd=False, name="Body"):
    """The body as an nTop implicit (metres)."""
    end = CFD_TAIL[1] if cfd else rk.R.total_length
    field = g.leaf("real_field", {"expression": _expression(cfd), "units": {}})
    scaled = g.block("multiply<real_field,real_field>", "real_field", [field, g.real(1e-3, rk.LENGTH)], top=False)
    pad = 5.0
    bounds = g.block("create_bounding_box<point,point>", "bounding_box",
                     [g.point(-pad, -0.5 * BOTTOM - pad, -0.5 * HEIGHT - pad),
                      g.point(end + pad, 0.5 * BOTTOM + pad, 0.5 * HEIGHT + pad)], top=False)
    return g.block("set_bounding_box<real_field,bounding_box>[1.1.0]", "implicit", [scaled, bounds], name)


def fins(g, inner_fraction=1.0):
    """The four fins in the order +y, +z, -y, -z, each a +y fin at its root distance turned about x. The plate
    starts at ``inner_fraction`` of the root distance (1: at the body surface)."""
    out = []
    for label in ("+y", "+z", "-y", "-z"):
        r = replace(rk.R, diameter=2.0 * ROOTS[label])
        fin = rk.fin(g, r, inner_fraction)
        turn = TURNS[label]
        if turn:
            fin = g.block("rotate<spatial3d,point,vector,real>[1.1.0]", "implicit",
                          [fin, g.point(0, 0, 0), g.vector(1, 0, 0), g.real(turn * rk.FULL_TURN / 4, rk.ANGLE)],
                          f"Fin {label}", top=False)
        out.append(fin)
    return out


if __name__ == "__main__":
    print(f"section area {section_area():.0f} mm^2 (round {math.pi * R_ROUND**2:.0f}); side root {SIDE_ROOT} mm")
    print(_expression(False)[:400], "...")
