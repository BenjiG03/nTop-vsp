"""Project a thin-surface (vortex lattice) VSPAERO solution onto a thick surface mesh, for visualization.

    python run_project.py SETTINGS_FILE

A vortex lattice solves for the pressure jump across a thin sheet, not for the pressure on each face of the real,
thick vehicle. For each point of the thick mesh (its nodes, its triangle centroids, or both) this driver finds the
closest point on the thin sheet, interpolates the jump there (linear over the sheet triangle, from area-weighted node
values), and gives the thick point the half of the jump that belongs to its side of the sheet:

    jump  = Cp(side the sheet normals point to) - Cp(other side)     (as VSPAERO stores it)
    Cp    = +jump / 2   where the thick surface faces the way the sheet normal points (outward normals, dot > 0)
    Cp    = -jump / 2   where it faces the other way

This is the lifting (antisymmetric) part of the pressure only: a thin solution has no thickness effect, so the two
sides of a symmetric section at zero lift both show zero. Points farther than Maximum Distance from the sheet are
left out of the point maps. The geometry of the two meshes need not match; each thick point simply takes the value
of the nearest sheet point.

Settings: ``work_dir``, ``case`` (this projection's name; the thick mesh is ``<case>.stl`` in metres),
``vspaero_results`` (base path of the vortex lattice run), ``case_number`` (VSPAERO case, from 1),
``max_distance`` (m). Writes ``<case>.ntvsp.proj.<quantity>.<points>.csv`` for quantity ``cp`` (side Cp), ``dp``
(side p - p_inf in Pa) and ``jump`` (the jump itself), and points ``nodes``, ``centroids`` and ``both``.
"""
from __future__ import annotations

import math
import os
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DriverError, Log, Settings, main_wrapper, write_point_map  # noqa: E402


# --- vectors ---------------------------------------------------------------------------------------------------------

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


# --- meshes ----------------------------------------------------------------------------------------------------------

def read_stl(path):
    """(nodes, triangles) of a binary or ASCII STL; vertices with identical coordinates are merged."""
    with open(path, "rb") as f:
        data = f.read()
    corners = []
    n = struct.unpack("<I", data[80:84])[0] if len(data) >= 84 else -1
    if len(data) == 84 + 50 * n:
        for i in range(n):
            v = struct.unpack("<12f", data[84 + 50 * i:84 + 50 * i + 48])
            corners.append([tuple(v[3 + 3 * k:6 + 3 * k]) for k in range(3)])
    else:
        tri = []
        for line in data.decode("ascii", errors="replace").splitlines():
            a = line.split()
            if a and a[0] == "vertex":
                tri.append(tuple(float(c) for c in a[1:4]))
                if len(tri) == 3:
                    corners.append(tri)
                    tri = []
    index, nodes, tris = {}, [], []
    for c in corners:
        t = []
        for p in c:
            if p not in index:
                index[p] = len(nodes)
                nodes.append(p)
            t.append(index[p])
        if len(set(t)) == 3:
            tris.append(tuple(t))
    if not tris:
        raise DriverError(f"{path}: no triangles")
    return nodes, tris


def read_panels(path):
    """The thin run's triangles, in the order of its result files: (method, nodes, triangles)."""
    if not os.path.isfile(path):
        raise DriverError(f"{path} not found; connect the Results of a VSPAERO Analysis (VSP) run")
    with open(path) as f:
        lines = f.read().split("\n")
    method, n_nodes, n_tris = lines[0].split()
    n_nodes, n_tris = int(n_nodes), int(n_tris)
    nodes = [tuple(map(float, lines[1 + i].split())) for i in range(n_nodes)]
    tris = [tuple(map(int, lines[1 + n_nodes + i].split())) for i in range(n_tris)]
    return method, nodes, tris


def read_point_map_values(path):
    with open(path) as f:
        return [float(line.split(",")[3]) for line in f if line.strip()]


def face_normal(nodes, t):
    return cross(sub(nodes[t[1]], nodes[t[0]]), sub(nodes[t[2]], nodes[t[0]]))  # length = 2 x area


# --- closest point ---------------------------------------------------------------------------------------------------

def closest_on_triangle(p, a, b, c):
    """Closest point to p on triangle abc, as barycentric weights (u, v, w) of (a, b, c)."""
    ab, ac, ap = sub(b, a), sub(c, a), sub(p, a)
    d1, d2 = dot(ab, ap), dot(ac, ap)
    if d1 <= 0 and d2 <= 0:
        return 1.0, 0.0, 0.0
    bp = sub(p, b)
    d3, d4 = dot(ab, bp), dot(ac, bp)
    if d3 >= 0 and d4 <= d3:
        return 0.0, 1.0, 0.0
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return 1.0 - v, v, 0.0
    cp = sub(p, c)
    d5, d6 = dot(ab, cp), dot(ac, cp)
    if d6 >= 0 and d5 <= d6:
        return 0.0, 0.0, 1.0
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return 1.0 - w, 0.0, w
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return 0.0, 1.0 - w, w
    denom = 1.0 / (va + vb + vc)
    v, w = vb * denom, vc * denom
    return 1.0 - v - w, v, w


class TriangleGrid:
    """Uniform grid of triangle bounding boxes for closest-point queries."""

    def __init__(self, nodes, tris):
        self.nodes, self.tris = nodes, tris
        lengths = [math.sqrt(dot(sub(nodes[t[k]], nodes[t[(k + 1) % 3]]), sub(nodes[t[k]], nodes[t[(k + 1) % 3]])))
                   for t in tris for k in range(3)]
        self.h = max(2.0 * sum(lengths) / len(lengths), 1e-9)
        self.cells = {}
        self.boxes = []
        for i, t in enumerate(tris):
            lo = [min(nodes[v][k] for v in t) for k in range(3)]
            hi = [max(nodes[v][k] for v in t) for k in range(3)]
            self.boxes.append((lo, hi))
            a, b = self.cell(lo), self.cell(hi)
            for x in range(a[0], b[0] + 1):
                for y in range(a[1], b[1] + 1):
                    for z in range(a[2], b[2] + 1):
                        self.cells.setdefault((x, y, z), []).append(i)

    def cell(self, p):
        return tuple(int(math.floor(c / self.h)) for c in p)

    def closest(self, p, max_distance):
        """(distance, triangle, (u, v, w)) of the closest sheet point within max_distance, or None."""
        base = self.cell(p)
        best = None
        seen = set()
        r = 0
        while True:
            for x in range(base[0] - r, base[0] + r + 1):
                for y in range(base[1] - r, base[1] + r + 1):
                    # Only the shell of cells at Chebyshev distance r from the base cell.
                    if r == 0 or abs(x - base[0]) == r or abs(y - base[1]) == r:
                        zs = range(base[2] - r, base[2] + r + 1)
                    else:
                        zs = (base[2] - r, base[2] + r)
                    for z in zs:
                        for i in self.cells.get((x, y, z), ()):
                            if i in seen:
                                continue
                            seen.add(i)
                            if best is not None:
                                lo, hi = self.boxes[i]
                                gap = sum(max(lo[k] - p[k], 0.0, p[k] - hi[k]) ** 2 for k in range(3))
                                if gap >= best[0] * best[0]:
                                    continue  # its bounding box is already farther than the best hit
                            t = self.tris[i]
                            a, b, c = (self.nodes[v] for v in t)
                            w = closest_on_triangle(p, a, b, c)
                            q = tuple(w[0] * a[k] + w[1] * b[k] + w[2] * c[k] for k in range(3))
                            d = math.sqrt(dot(sub(p, q), sub(p, q)))
                            if best is None or d < best[0]:
                                best = (d, i, w)
            # Every triangle not seen yet lies at least r * h away.
            if best is not None and best[0] <= r * self.h:
                return best if best[0] <= max_distance else None
            if r * self.h > max_distance:
                return None
            r += 1


# --- driver ----------------------------------------------------------------------------------------------------------

def main(settings_path):
    s = Settings(settings_path)
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    for stale in os.listdir(work):
        if stale.startswith(case + ".ntvsp.proj.") and stale.endswith(".csv"):
            os.remove(os.path.join(work, stale))
    results = s.text("vspaero_results")
    number = s.integer("case_number")
    max_distance = s.real("max_distance")
    if max_distance <= 0:
        raise DriverError("Maximum Distance must be positive")

    method, sheet_nodes, sheet_tris = read_panels(results + ".ntvsp.panels.txt")
    if method != "vlm":
        raise DriverError("the VSPAERO results are a panel-method run; projection needs a vortex lattice (thin surface) "
                          "run")
    values_path = f"{results}.ntvsp.cp.centroids.{number}.csv"
    if not os.path.isfile(values_path):
        raise DriverError(f"no case {number} in the VSPAERO results ({values_path} not found)")
    jump_tri = read_point_map_values(values_path)
    if len(jump_tri) != len(sheet_tris):
        raise DriverError(f"{values_path} has {len(jump_tri)} values for {len(sheet_tris)} sheet triangles")
    aero = Settings(results + ".ntvsp.settings.txt")
    q = 0.5 * aero.real("rho") * aero.real("vinf") ** 2

    # Node values of the jump: area-weighted averages of the triangles around each node (as the Surface Pressure
    # block's node values), so the projected field is continuous across sheet triangles.
    sheet_normals = [face_normal(sheet_nodes, t) for t in sheet_tris]
    num, den = [0.0] * len(sheet_nodes), [0.0] * len(sheet_nodes)
    for t, n, v in zip(sheet_tris, sheet_normals, jump_tri):
        area = math.sqrt(dot(n, n))
        for i in t:
            num[i] += area * v
            den[i] += area
    jump_node = [a / b if b > 0 else 0.0 for a, b in zip(num, den)]

    thick_nodes, thick_tris = read_stl(os.path.join(work, case + ".stl"))
    log(f"projection: case {number} of {results}: {len(sheet_tris)} sheet triangles onto {len(thick_tris)} thick "
        f"triangles (maximum distance {max_distance:g} m)")
    thick_normals = [face_normal(thick_nodes, t) for t in thick_tris]
    node_normal = [[0.0, 0.0, 0.0] for _ in thick_nodes]
    for t, n in zip(thick_tris, thick_normals):
        for i in t:
            for k in range(3):
                node_normal[i][k] += n[k]
    centroids = [tuple(sum(thick_nodes[i][k] for i in t) / 3.0 for k in range(3)) for t in thick_tris]
    used = sorted({i for t in thick_tris for i in t})

    grid = TriangleGrid(sheet_nodes, sheet_tris)
    samples = {"nodes": [(thick_nodes[i], node_normal[i]) for i in used],
               "centroids": list(zip(centroids, thick_normals))}
    out = {}
    for kind, pts in samples.items():
        rows = []
        for p, n in pts:
            hit = grid.closest(p, max_distance)
            if hit is None:
                continue
            _, ti, w = hit
            t = sheet_tris[ti]
            jump = w[0] * jump_node[t[0]] + w[1] * jump_node[t[1]] + w[2] * jump_node[t[2]]
            side = 1.0 if dot(n, sheet_normals[ti]) >= 0.0 else -1.0
            rows.append((p, side * 0.5 * jump, jump))
        out[kind] = rows
        log(f"{kind}: {len(rows)} of {len(pts)} points within {max_distance:g} m of the sheet")
    if not out["nodes"] and not out["centroids"]:
        raise DriverError(f"no point of the thick mesh lies within {max_distance:g} m of the thin surface")
    out["both"] = out["centroids"] + out["nodes"]
    for kind, rows in out.items():
        pts = [r[0] for r in rows] or [(0.0, 0.0, 0.0)]
        for quantity, values in (("cp", [r[1] for r in rows]), ("dp", [q * r[1] for r in rows]),
                                 ("jump", [r[2] for r in rows])):
            write_point_map(os.path.join(work, f"{case}.ntvsp.proj.{quantity}.{kind}.csv"), pts, values or [0.0])
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
