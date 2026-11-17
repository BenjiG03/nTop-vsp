"""Read the VSPAERO panel file (``.vspgeom`` v3) that OpenVSP writes for a model.

Layout: a header (version, number of mesh levels, then ``nodes faces max_nodes_per_face``), the nodes, the
faces (``n node ...``, 1-based), one tag line per face (``part tag u v ...``), one face-map line per face, the
wakes (``count`` then ``n part node ...`` per wake, wrapped at 10 nodes a line), and, when faces have more than
three nodes, a triangulation of each face (``face ntri node ...``) followed by its tag lines.

VSPAERO solves on the faces (its level 1) and stores surface results per triangle of the triangulation (level 0),
in file order; ``level0`` lists those triangles. Standard library only.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class VSPGeom:
    nodes: list            # (x, y, z)
    faces: list            # [node index, ...], 0-based
    parts: list            # part number of each face (1-based, as in the .vkey file)
    wakes: list            # (part, [node index, ...])
    triangles: list = field(default_factory=list)  # per face: [(a, b, c), ...]; empty for an all-triangle file


def _tokens(lines, pos):
    """Yield whitespace-separated tokens from ``lines`` starting at line ``pos``, with the line each came from."""
    while pos < len(lines):
        for t in lines[pos].split():
            yield t, pos
        pos += 1


def read(path) -> VSPGeom:
    with open(path) as f:
        lines = f.read().split("\n")
    if not lines[0].lower().startswith("# vspgeom v3"):
        raise ValueError(f"{path}: not a vspgeom v3 file")
    if int(lines[1].split()[0]) != 1:
        raise ValueError(f"{path}: only single-level panel files are supported")
    head = lines[2].split()
    n_nodes, n_faces = int(head[0]), int(head[1])
    pos = 3
    nodes = [tuple(map(float, lines[pos + i].split()[:3])) for i in range(n_nodes)]
    pos += n_nodes
    if int(lines[pos].split()[0]) != n_faces:
        raise ValueError(f"{path}: face count mismatch")
    pos += 1
    faces = []
    for i in range(n_faces):
        a = lines[pos + i].split()
        faces.append([int(v) - 1 for v in a[1:1 + int(a[0])]])
    pos += n_faces
    parts = [int(lines[pos + i].split()[0]) for i in range(n_faces)]
    pos += 2 * n_faces  # tags, then the face map
    tok = _tokens(lines, pos)
    n_wakes = int(next(tok)[0])
    wakes = []
    for _ in range(n_wakes):
        count = abs(int(next(tok)[0]))
        part = int(next(tok)[0])
        wakes.append((part, [int(next(tok)[0]) - 1 for _ in range(count)]))
    triangles = []
    if any(len(f) > 3 for f in faces):
        for i in range(n_faces):
            fid, line = next(tok)
            if int(fid) != i + 1:
                raise ValueError(f"{path}: triangulation of face {i + 1} not found (line {line + 1})")
            ntri = int(next(tok)[0])
            tris = []
            for _ in range(ntri):
                tris.append(tuple(int(next(tok)[0]) - 1 for _ in range(3)))
            triangles.append(tris)
    return VSPGeom(nodes, faces, parts, wakes, triangles)


def level0(g: VSPGeom):
    """VSPAERO's surface triangles in the order of its stored results: (triangle, face index) pairs."""
    if not g.triangles:
        return [(tuple(f), i) for i, f in enumerate(g.faces)]
    return [(t, i) for i, tris in enumerate(g.triangles) for t in tris]


def thin_parts(vkey_path):
    """Part numbers OpenVSP marked thin (vortex lattice) in the ``.vkey`` file."""
    out = set()
    for line in open(vkey_path, encoding="utf-8", errors="replace"):
        v = [x.strip() for x in line.split(",")]
        if len(v) >= 7 and v[0].isdigit() and v[5] in ("0", "1"):
            if v[5] == "0":
                out.add(int(v[0]))
    return out
