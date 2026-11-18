"""Inviscid SU2 on the trapezoid-section rocket, with hypanel's SU2 chain (hypanel/validation/rocket/su2).

    sh build/py.sh validation/noncircular/cfd.py mesh [tolerance edge half_angle far rate]   # nTop CFD mesh
    sh build/py.sh validation/noncircular/cfd.py run [iterations ranks]                      # the SU2 runs
    sh build/py.sh validation/noncircular/cfd.py collect                                     # results/su2.json

The mesh is built as hypanel's finned rocket (build_cfd.py rocket 3 3 20 100 0.3): the body runs on past the base
(x = 1200 mm) to a closed tail, so it has no base wake, and forces are taken from x < 1200 mm only (supersonic flow
cannot carry the tail upstream). Faces are classed as fin where they lie outside the body's own field
(``geometry.body_field``). Coefficients: Sref = pi 54^2 / 4 mm^2, moments over Sref 54 mm about x = 700 mm on the axis
(Cm positive nose up). Meshes and runs go to runs/noncircular/su2 (git-ignored).
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import geometry as geo  # noqa: E402
import rocket as rk  # noqa: E402
from ntopfile import NONE  # noqa: E402
from ntopfile import Graph  # noqa: E402

sys.path.insert(0, os.path.join(rk.PLAYGROUND, "hypanel", "validation", "rocket", "su2"))
import run_su2 as su2  # noqa: E402

CFD = os.path.join(rk.ROOT, "runs", "noncircular", "su2")
su2.CFD = CFD
MM = 1e-3
L, A = rk.LENGTH, rk.ANGLE
DOMAIN = (-40.0, 1500.0)
MACHS = (2.0, 3.0, 5.0)
ALPHAS = (-4.0, 0.0, 4.0, 8.0)
X_REF, D = 700.0, 54.0
AREF = math.pi * 0.25 * D**2
OUT = os.path.join(HERE, "results", "su2.json")


def _graded(g, field, near_mm, rate, name):
    dist = g.block("abs<real_field>", "real_field", [field], name + " |d|", top=False, units=L)
    scaled = g.block("multiply<real_field,real_field>", "real_field", [g.real(rate), dist], name + " rate", top=False,
                     units=L)
    return g.block("add<real_field,real_field>", "real_field", [g.real(near_mm * MM, L), scaled], name, top=False,
                   units=L)


def _min(g, a, b, name):
    return g.block("min<real_field,real_field>", "real_field", [a, b], name, top=False, units=L)


def mesh(tol_mm=3.0, edge_mm=3.0, half_angle=20.0, far_mm=100.0, rate=0.3, growth=1.2, fin_ratio=0.6):
    tag = f"trap_t{tol_mm:g}_h{edge_mm:g}"
    g = Graph()
    body = geo.body(g, cfd=True)
    fins = geo.fins(g, inner_fraction=0.75)
    solid = rk.boolean(g, "union", [body] + fins, "Rocket")
    r0 = 40.0
    r1 = r0 + (DOMAIN[1] - DOMAIN[0]) * math.tan(math.radians(half_angle))
    outer = g.block("cone<point,point,real,real>", "cone", [g.point(DOMAIN[0], 0, 0), g.point(DOMAIN[1], 0, 0),
                                                           g.real(r0 * MM, L), g.real(r1 * MM, L)], "Far Field")
    items = g.block("core.list<implicit>", "list<implicit>", [solid], top=False)
    domain = g.block("boolean_subtract<blend_enum,real_field,implicit,list<implicit>>[5.44.0]", "implicit",
                     [g.enum("blend_enum", 0), g.real(0.0, L), outer, items], "Fluid Domain")
    surface = g.block("implicit_to_mesh<implicit,real,real,bool,bool>[2.5.0]", "mesh",
                      [domain, g.real(tol_mm * MM, L), NONE, g.boolean(True), g.boolean(False)], "Source Surface")
    edge = _min(g, g.real(far_mm * MM, L), _graded(g, solid, edge_mm, rate, "From Body"), "Edge A")
    tip = g.block("sphere<point,real>", "sphere", [g.point(0, 0, 0), g.real(0.001 * MM, L)], "Tip", top=False)
    edge = _min(g, edge, _graded(g, tip, 0.25 * edge_mm, 0.03, "From Tip"), "Edge B")
    edge = _min(g, edge, _graded(g, rk.boolean(g, "union", fins, "Fins", top=False), fin_ratio * edge_mm, rate,
                                 "From Fins"), "Edge C")
    remesh = g.block(
        "vki_parallel_remesh_surfmesh<mesh,real_field,surface_mesh_shape_enum,real,real,real,real,real,real,integer>",
        "mesh", [surface, edge, g.enum("surface_mesh_shape_enum", 0), g.real(0.35, A), g.real(growth), g.real(0.14, A),
                 g.real(0.1 * edge_mm * MM, L), g.real(0.0, L), g.real(0.0, L), g.integer(16)], "Remeshed Surface")
    volume = g.block("vki_volume_mesh<mesh,real_field,real,list<point>>[1.2.0]", "volume_mesh",
                     [remesh, edge, g.real(growth), NONE], "Volume Mesh")
    fe = g.block("fe_volume_mesh<volume_mesh,geometric_order_enum>", "fe_mesh", [volume, g.enum("geometric_order_enum", 0)],
                 "FE Mesh")
    select = g.block("fe_boundary_by_implicit<fe_mesh,fe_boundary_enum,implicit,bool,real,text>[5.52.0]", "fe_boundary",
                     [fe, g.enum("fe_boundary_enum", 2), solid, g.boolean(False), g.real(max(1.0, tol_mm) * MM, L),
                      g.text("Wall Faces")], "Wall Faces")
    wall = g.block("cfd_boundary_on_fe_mesh<boundary,cfd_boundary_enum,text>", "cfd_boundary",
                   [select, g.enum("cfd_boundary_enum", 1), g.text("wall")], "Wall")
    walls = g.block("core.list<cfd_boundary>", "list<cfd_boundary>", [wall], top=False)
    os.makedirs(CFD, exist_ok=True)
    path = os.path.join(CFD, tag + ".msh")
    out = g.block("export_cfd_mesh<file_path,fe_mesh,unit_length_enum,list<cfd_boundary>,text>", "fe_mesh_file_data",
                  [g.file_path(path), fe, g.unit_length("mm"), walls, g.text("farfield")], "Export")
    res, log, dt = rk.run(g, [], out, "fe_mesh_file_data", os.path.join(CFD, tag + ".ntop"), settle=20)
    size = os.path.getsize(path) if os.path.exists(path) else None
    print(tag, f"{dt:.0f} s", "file", size, flush=True)
    print("\n".join(l[:240] for l in log.splitlines() if "[E]" in l or "[W]" in l)[:1500])
    return tag


def body_field(p):
    """geometry.body_field (cfd) on an array of points (mm)."""
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    r = rk.R
    rho = (geo.R_ROUND**2 + r.nose_length**2) / (2.0 * geo.R_ROUND)
    xc = np.clip(x, 0.0, r.nose_length)
    s = (np.sqrt(rho**2 - (r.nose_length - xc) ** 2) + geo.R_ROUND - rho) / geo.R_ROUND
    s = np.maximum(np.minimum(s, (geo.CFD_TAIL[1] - x) / (geo.CFD_TAIL[1] - geo.CFD_TAIL[0])), 1e-4)
    d = np.stack([ny * y + nz * z - s * (c - geo.CORNER) for ny, nz, c in geo.edges()])
    f = np.sqrt((np.maximum(d, 0.0) ** 2).sum(axis=0)) + np.minimum(d.max(axis=0), 0.0) - s * geo.CORNER
    return np.maximum(np.maximum(f, -x), x - geo.CFD_TAIL[1])


def summarize(run_dir, npz, mach):
    nodes, faces, cp = su2.read_surface(run_dir, npz, mach)
    centroid, area_vec, cp_face, fx, fy, fz = su2.forces(nodes, faces, cp)
    x, z = centroid[:, 0], centroid[:, 2]
    keep = x < 1200.0
    fin = keep & (body_field(centroid) > 0.3)
    my = ((z * fx - (x - X_REF) * fz) / D)  # nose-up positive, about (700, 0, 0) mm
    regions = {"nose": keep & ~fin & (x < 270.0), "tube": keep & ~fin & (x >= 270.0) & (x < 1040.0),
               "tube at fins": keep & ~fin & (x >= 1040.0), "fins": fin}
    out = {"CN": float(fz[keep].sum()), "CA": float(fx[keep].sum()), "CY": float(fy[keep].sum()),
           "Cm": float(my[keep].sum())}
    out["x_cp_mm"] = X_REF - out["Cm"] * D / out["CN"] if abs(out["CN"]) > 1e-9 else float("nan")
    out["regions"] = {k: {"CN": float(fz[m].sum()), "Cm": float(my[m].sum())} for k, m in regions.items()}
    edges = np.arange(0.0, 1200.0 + 1e-6, 10.0)
    out["edges"] = edges.tolist()
    out["dcn_dx"] = {k: np.round(np.histogram(x[m], bins=edges, weights=fz[m])[0] / np.diff(edges), 7).tolist()
                     for k, m in (("body", keep & ~fin), ("fins", fin))}
    return out


def run(iters=1500, ranks=14, tag="trap_t3_h3"):
    npz = su2.convert(tag)
    for m in MACHS:
        for al in ALPHAS:
            run_dir = os.path.join(CFD, f"run_{tag}_M{m:g}_a{al:g}")
            if os.path.exists(os.path.join(run_dir, "summary.json")):
                continue
            start = time.time()
            su2.run(tag, m, al, iters, ranks)
            try:
                s = summarize(run_dir, npz, m)
                json.dump(s, open(os.path.join(run_dir, "summary.json"), "w"))
                print(f"M{m:g} a{al:g}: CN {s['CN']:.4f} Cm {s['Cm']:.4f} x_cp {s['x_cp_mm']:.1f} mm "
                      f"({time.time() - start:.0f} s)", flush=True)
            except Exception as exc:  # the log keeps what went wrong
                print(f"M{m:g} a{al:g}: post failed: {exc}", flush=True)


def collect(tag="trap_t3_h3"):
    out = {}
    for path in sorted(glob.glob(os.path.join(CFD, f"run_{tag}_M*_a*", "summary.json"))):
        name = os.path.basename(os.path.dirname(path))[len(f"run_{tag}_"):]
        s = json.load(open(path))
        hist = os.path.join(os.path.dirname(path), "history.csv")
        if os.path.exists(hist):
            rows = [l.split(",") for l in open(hist).read().splitlines()[1:] if l.strip()]
            s["iterations"] = len(rows)
            s["rms_density"] = float(rows[-1][1]) if rows else float("nan")
        out[name] = s
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump({"mesh": tag, "runs": out}, open(OUT, "w"))
    for k, s in out.items():
        print(k, f"CN {s['CN']:.4f}, Cm {s['Cm']:.4f}, x_cp {s['x_cp_mm']:.1f} mm", s.get("iterations"))


if __name__ == "__main__":
    what, args = sys.argv[1], [float(v) for v in sys.argv[2:]]
    if what == "mesh":
        mesh(*args)
    elif what == "run":
        run(*(int(v) for v in args))
    elif what == "collect":
        collect()
