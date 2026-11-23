"""Check the CAD files of test_cad.py against the nTop implicit bodies.

1. Each CAD file is read with OpenCASCADE through gmsh: solids, shells, faces, total area, solid volumes. Its
   faces are meshed (about 2 mm) and the nodes kept as sample points of the CAD surface.
2. nTop mesh vertices of each body (STL from the build notebook) are sampled; the duct's end-face vertices are left
   out (the duct is open), and so are vertices inside another body (for the trimmed union).
3. One nTop notebook evaluates the three bodies' fields at all those points (Import Scalar Point Map, Evaluate
   Field, Export Point Map), so the "implicit" is nTop's own.
4. CAD -> implicit: for a CAD point, the smallest |field| over the bodies. Implicit -> CAD: for an nTop vertex,
   the distance to the nearest CAD face (OpenCASCADE closest point).

Writes runs/cad_export/check.json and prints a table.
"""
from __future__ import annotations

import json
import math
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))
import rocket as rk  # noqa: E402
import test_cad as T  # noqa: E402
from ntopfile import Graph, Ref  # noqa: E402

W = T.WORK
BODIES = ("Body", "Fin", "Duct")
MAX_VERTS = 3000
GRID = 60  # sample points per face: GRID x GRID in its parameter range, trimmed


def read_stl(path):
    data = open(path, "rb").read()
    if data[:5] == b"solid" and b"facet" in data[:400]:
        v = [tuple(map(float, l.split()[1:4])) for l in data.decode(errors="replace").splitlines()
             if l.strip().startswith("vertex")]
        return np.unique(np.array(v), axis=0)
    n = struct.unpack("<I", data[80:84])[0]
    rec = np.frombuffer(data[84:84 + 50 * n], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
    return np.unique(rec["v"].reshape(-1, 3).astype(float), axis=0)


def cad_stats(path):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.importShapes(path, highestDimOnly=False)
    gmsh.model.occ.synchronize()
    ents = {d: len(gmsh.model.getEntities(d)) for d in range(4)}
    vols = [gmsh.model.occ.getMass(3, t) for _, t in gmsh.model.getEntities(3)]
    faces = [t for _, t in gmsh.model.getEntities(2)]
    area = sum(gmsh.model.occ.getMass(2, t) for t in faces)
    # sample points: a parameter grid on every face, kept where it is inside the face's trimmed region (meshing the
    # STEP faces with gmsh returned only their corner nodes)
    chunks = []
    for t in faces:
        lo, hi = gmsh.model.getParametrizationBounds(2, t)
        n = GRID
        uv = [(lo[0] + (hi[0] - lo[0]) * (i + 0.5) / n, lo[1] + (hi[1] - lo[1]) * (j + 0.5) / n)
              for i in range(n) for j in range(n)]
        keep = [c for p in uv if gmsh.model.isInside(2, t, list(p), parametric=True) for c in p]
        if keep:
            chunks.append(np.array(gmsh.model.getValue(2, t, keep)).reshape(-1, 3))
    pts = np.concatenate(chunks) if chunks else np.zeros((0, 3))
    gmsh.finalize()
    return {"entities": ents, "volumes_mm3": vols, "area_mm2": area, "faces": len(faces)}, pts


def closest(path, pts):
    """Distance from each point to the nearest face of the CAD file."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.importShapes(path, highestDimOnly=False)
    gmsh.model.occ.synchronize()
    faces = [t for _, t in gmsh.model.getEntities(2)]
    best = np.full(len(pts), np.inf)
    flat = pts.reshape(-1).tolist()
    for t in faces:
        c, _ = gmsh.model.getClosestPoint(2, t, flat)
        d = np.linalg.norm(np.array(c).reshape(-1, 3) - pts, axis=1)
        best = np.minimum(best, d)
    gmsh.finalize()
    return best


def write_points(path, pts):
    with open(path, "w") as f:  # no header: nTop's point-map import rejects one
        for p in pts:
            f.write(f"{p[0]:.6f},{p[1]:.6f},{p[2]:.6f},0\n")


def evaluate(sets):
    """nTop fields of the three bodies at each point set: {name: csv path} -> {name: (n, 3) array of fields}."""
    g = Graph()
    bodies = T.bodies(g)
    outs = []
    one = g.leaf("real_field", {"expression": "1", "units": {}})
    for name, path in sets.items():
        pm = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map",
                     [g.file_path(path), g.unit_length("mm"), one], f"{name} Points")
        pts = Ref(pm, ("points",))
        for b in BODIES:
            vals = g.block("evaluate_field<real_field,list<point>>", "list<real>", [bodies[b], pts], top=False)
            m = g.block("point_map<list<point>,list<real>>", "real_point_map", [pts, vals], top=False)
            outs.append(g.block("export_point_map<file_path,real_point_map,unit_length_enum>", "point_map_file",
                                [g.file_path(f"{W}/eval/{name}.{b}.csv"), m, g.unit_length("mm")], f"{name} {b}"))
    out = g.block("core.group<any>", "group", outs, "Outputs")
    os.makedirs(f"{W}/eval", exist_ok=True)
    res, log, dt = rk.run(g, [], out, "group", f"{W}/eval/evaluate.ntop", timeout=7200)
    print(f"nTop evaluation: {dt:.0f} s")
    for line in log.splitlines():
        if "[E]" in line:
            print("  ", line[:300])
    vals = {}
    for name in sets:
        cols = []
        for b in BODIES:
            a = np.loadtxt(f"{W}/eval/{name}.{b}.csv", delimiter=",", ndmin=2)
            cols.append(a[:, 3] * 1000.0)  # Export Point Map writes metres whatever the unit
        vals[name] = np.stack(cols, axis=1)
    return vals


def on_cap(p):
    """Within 0.05 mm of either end plane of the duct (its end faces, left out of the open duct)."""
    d_in = np.abs((p - np.array(T.P_IN)) @ np.array(T.N_IN))
    d_out = np.abs((p - np.array(T.P_OUT)) @ np.array(T.N_OUT))
    return (d_in < 0.05) | (d_out < 0.05)


def stats(d):
    d = np.asarray(d)
    if not len(d):
        return {}
    return {"n": int(len(d)), "mean": float(d.mean()), "p99": float(np.percentile(d, 99)), "max": float(d.max())}


def main():
    rng = np.random.default_rng(1)
    report = {"cases": {}, "implicit_to_cad": {}}
    sets, cad_pts = {}, {}
    for method, fmt, caps in T.CASES:
        name = T.case_name(method, fmt, caps)
        path = f"{W}/{name}." + ("stp" if fmt == "step" else "igs")
        st, pts = cad_stats(path)
        report["cases"][name] = st
        if len(pts) > 60000:
            pts = pts[rng.choice(len(pts), 60000, replace=False)]
        cad_pts[name] = pts
        write_points(f"{W}/eval_{name}.csv", pts)
        sets[name] = f"{W}/eval_{name}.csv"
    verts = {}
    for b in BODIES:
        v = read_stl(f"{W}/{b}.ntop_mesh.stl")
        if b == "Duct":
            v = v[~on_cap(v)]
        if len(v) > MAX_VERTS:
            v = v[rng.choice(len(v), MAX_VERTS, replace=False)]
        verts[b] = v
        write_points(f"{W}/eval_mesh_{b}.csv", v)
        sets[f"mesh_{b}"] = f"{W}/eval_mesh_{b}.csv"
    f = evaluate(sets)
    print(f"{'case':28s} {'solids':>6s} {'faces':>5s} {'area m2':>9s}  CAD->implicit by body: mean / p99 / max (mm)")
    for name, pts in cad_pts.items():
        vals = f[name]
        owner = np.argmin(np.abs(vals), axis=1)
        dev = np.abs(vals)[np.arange(len(vals)), owner]
        st = report["cases"][name]
        st["cad_to_implicit"] = {b: stats(dev[owner == i]) for i, b in enumerate(BODIES)}
        line = " ".join(f"{b} {s['mean']:.3f}/{s['p99']:.3f}/{s['max']:.3f}" for b, s in st["cad_to_implicit"].items()
                        if s)
        print(f"{name:28s} {st['entities'][3]:6d} {st['faces']:5d} {st['area_mm2'] * 1e-6:9.5f}  {line}")
        if st["volumes_mm3"]:
            print(f"{'':28s} solid volumes (cm3): {[round(v * 1e-3, 2) for v in st['volumes_mm3']]}")
    print("implicit -> CAD (nTop vertices to the nearest CAD face): mean / p99 / max (mm)")
    for name in cad_pts:
        path = f"{W}/{name}." + ("igs" if "_iges" in name else "stp")
        row = {}
        for i, b in enumerate(BODIES):
            v = verts[b]
            if "trimmed" in name:  # the union: leave out vertices inside another body
                others = [j for j in range(len(BODIES)) if j != i]
                v = v[np.all(f[f"mesh_{b}"][:, others] > 0.05, axis=1)]
            row[b] = stats(closest(path, v))
        report["implicit_to_cad"][name] = row
        print(f"{name:28s} " + " ".join(f"{b} {s['mean']:.3f}/{s['p99']:.3f}/{s['max']:.3f}" for b, s in row.items()))
    # the nTop meshes themselves, against nTop's field (the size of the reference's own error)
    print("nTop mesh vertices vs nTop field (reference noise): " + " ".join(
        f"{b} {np.abs(f['mesh_' + b][:, i]).mean():.4f}/{np.abs(f['mesh_' + b][:, i]).max():.4f}"
        for i, b in enumerate(BODIES)))
    json.dump(report, open(f"{W}/check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
