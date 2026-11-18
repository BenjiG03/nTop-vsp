"""The hobby rocket of hypanel's RASAero II comparison, as nTop implicit bodies, and helpers to run study notebooks.

The dimensions come from ``hypanel/validation/rocket/design.py`` (the same file that wrote the RASAero II design,
``hobby_rocket.CDX1``), so both tools see the same vehicle: a 5:1 tangent-ogive nose, a 54 mm tube 1.2 m long and four
hexagonal-section fins on the +y, +z, -y and -z sides. Nose tip at the origin, axis along +x (aft), z up.
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
import time
import json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PLAYGROUND = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(ROOT, "build"))
sys.path.insert(0, os.path.join(PLAYGROUND, "hypanel", "validation", "rocket"))

from design import Rocket  # noqa: E402
from ntopfile import Graph, save_host  # noqa: E402

LENGTH = {"length": 1}
AREA = {"length": 2}
ANGLE = {"angle": 1}
FULL_TURN = 2.0 * math.pi
NTOPCL = r"C:\Program Files\nTopology\nTopology\ntopcl.exe"
RUNS = os.path.join(ROOT, "runs", "rasaero")
PYTHON = os.path.join(ROOT, "dev", "venv", "Scripts", "python.exe")
VSP = os.path.join(ROOT, "OpenVSP-3.56.0-win64")
R = Rocket()


def _points(g, pts):
    return g.block("core.list<point>", "list<point>", [g.point(*p) for p in pts], top=False)


def _profile(g, pts, name):
    return g.block("profile_from_points<list<point>>[5.20.0]", "new_profile", [_points(g, pts)], name, top=False)


def _extrude(g, profile, distance_mm, direction, name):
    return g.block("extrude<new_profile,real,real,bool,vector>[5.20.0]", "implicit",
                   [profile, g.real(distance_mm * 1e-3, LENGTH), g.real(0.0, ANGLE), g.boolean(True),
                    g.vector(*direction)], name, top=False)


def boolean(g, kind, items, name, top=True):
    parts = g.block("core.list<implicit>", "list<implicit>", items, top=False)
    return g.block(f"boolean_{kind}<blend_enum,real_field,list<implicit>>[5.44.0]", "implicit",
                   [g.enum("blend_enum", 0), g.real(0.0, LENGTH), parts], name, top=top)


def ogive_points(r=R, count=160):
    radius = 0.5 * r.diameter
    rho = (radius ** 2 + r.nose_length ** 2) / (2.0 * radius)
    return [(x, 0.0, math.sqrt(rho ** 2 - (r.nose_length - x) ** 2) + radius - rho)
            for x in (r.nose_length * i / count for i in range(count + 1))]


def body(g, r=R):
    radius = 0.5 * r.diameter
    outline = ogive_points(r) + [(r.total_length, 0.0, radius), (r.total_length, 0.0, 0.0)]
    axis = g.block("axis<point,vector>", "axis", [g.point(0, 0, 0), g.vector(1, 0, 0)], "Rocket Axis", top=False)
    return g.block("revolve<new_profile,axis,real>[5.20.0]", "implicit",
                   [_profile(g, outline, "Body Profile"), axis, g.real(FULL_TURN, ANGLE)], "Body")


def fin(g, r=R, inner_fraction=0.75):
    """One fin on the +y side (hypanel's build_rocket.add_fin): a flat plate cut by two wedges (hexagonal section).
    The plate starts at ``inner_fraction`` of the body radius (inside the body; 0 reaches the axis)."""
    radius = 0.5 * r.diameter
    half = 0.5 * r.fin_thickness
    inner = inner_fraction * radius
    x_root, x_end = r.fin_location, r.total_length
    y_tip = radius + r.fin_span
    slope = r.fin_sweep / r.fin_span
    x_inner = x_root - slope * (radius - inner)
    plate = _extrude(g, _profile(g, [(x_inner, inner, 0.0), (x_root + r.fin_sweep, y_tip, 0.0), (x_end, y_tip, 0.0),
                                     (x_end, inner, 0.0)], "Fin Planform"), r.fin_thickness, (0, 0, 1), "Fin Plate")
    size = math.hypot(r.fin_sweep, r.fin_span)
    along = (r.fin_sweep / size, r.fin_span / size, 0.0)
    aft = (r.fin_span / size, -r.fin_sweep / size, 0.0)
    bevel_normal = r.fin_bevel * aft[0]
    reach = 2.0 * (x_end - x_inner)
    opening = half * reach / bevel_normal
    centre = (x_root + 0.5 * r.fin_sweep, radius + 0.5 * r.fin_span, 0.0)
    far = (centre[0] + reach * aft[0], centre[1] + reach * aft[1], 0.0)
    leading = _extrude(g, _profile(g, [centre, (far[0], far[1], opening), (far[0], far[1], -opening)], "LE Section"),
                       4.0 * size, along, "Leading Bevel")
    opening = half * reach / r.fin_bevel
    y_mid = 0.5 * (inner + y_tip)
    trailing = _extrude(g, _profile(g, [(x_end, y_mid, 0.0), (x_end - reach, y_mid, opening),
                                        (x_end - reach, y_mid, -opening)], "TE Section"),
                        4.0 * (y_tip - inner), (0, 1, 0), "Trailing Bevel")
    return boolean(g, "intersect", [plate, leading, trailing], "Fin +y", top=False)


def fins(g, r=R, inner_fraction=0.75):
    """The four fins, +y, +z, -y, -z (rotations of the +y fin about x by 90 degrees)."""
    first = fin(g, r, inner_fraction)
    out = [first]
    for i, label in ((1, "+z"), (2, "-y"), (3, "-z")):
        out.append(g.block("rotate<spatial3d,point,vector,real>[1.1.0]", "implicit",
                           [first, g.point(0, 0, 0), g.vector(1, 0, 0), g.real(i * FULL_TURN / 4, ANGLE)], f"Fin {label}",
                           top=False))
    return out


def mesh(g, body, tol_mm, name=None):
    """Mesh from Implicit Body (sharp edges kept), for the on-mesh blocks."""
    from ntopfile import NONE
    return g.block("implicit_to_mesh<implicit,real,real,bool,bool>[2.5.0]", "mesh",
                   [body, g.real(tol_mm * 1e-3, LENGTH), NONE, g.boolean(True), g.boolean(False)], name, top=False)


def env():
    e = dict(os.environ)
    e["NTOP_PYTHON_PATH"] = PYTHON
    e["VSP_PATH"] = VSP
    e["RASAERO_PATH"] = os.path.join(os.environ.get("LOCALAPPDATA", ""), "Programs", "RASAero II")
    return e


def call(g, block, given, name, top=True):
    inputs = []
    for spec in block.inputs:
        if spec.name in given:
            v = given[spec.name]
            inputs.append(v(g) if callable(v) else v)
        else:
            inputs.append(spec.default(g))
    return g.block(block.signature, block.output_type, inputs, name, top=top)


def run(g, blocks, output, output_type, path, timeout=7200, settle=6):
    """Save and run a notebook with ntopcl (trusted, the driver environment set); returns (json, log text, s).

    ntopcl may never exit after Run Command blocks inside custom blocks: once the output file is written and stable
    it is stopped."""
    save_host(g, list(blocks), path, output=output, output_type=output_type)
    res = os.path.splitext(path)[0] + ".json"
    log_path = os.path.splitext(path)[0] + ".log"
    if os.path.exists(res):
        os.remove(res)
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        p = subprocess.Popen([NTOPCL, "--trustnotebook", "-v", "2", "-o", res, path], stdout=logf,
                             stderr=subprocess.STDOUT, env=env())
        last, stable = -1, 0
        while p.poll() is None and time.time() - t0 < timeout:
            time.sleep(1)
            if os.path.exists(res):
                size = os.path.getsize(res)
                stable = stable + 1 if size == last and size > 0 else 0
                last = size
                if stable >= settle:
                    break
        if p.poll() is None:
            p.kill()
            p.wait()
    log = open(log_path, encoding="utf-8", errors="replace").read()
    out = json.load(open(res)) if os.path.exists(res) else None
    return out, log, time.time() - t0
