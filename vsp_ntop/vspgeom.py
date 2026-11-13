"""Read a VSPAERO ``.vspgeom`` (v3) panel file and give it trailing-edge wakes.

VSP writes a ``.vspgeom`` file from an imported tessellation (a MeshGeom), but it marks no trailing edges on it:
VSP's own wake finder only knows the trailing edges of its parametric wings. Without wakes VSPAERO's panel method
sheds no vorticity and the body carries almost no lift. The user lists the mesh vertices that lie on trailing edges;
this module turns the mesh edges joining those vertices into wakes and writes them into the file's wake section, in
the same form VSP uses for wing trailing edges:

    <number of wakes>
    <number of nodes> <part number> <node 1> <node 2> ...      (one line per wake, nodes 1-based)

Only single-level files (no mesh refinement levels) are handled; that is what VSP writes for a tessellation.
Pure standard library, so any Python 3 runs it.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class VSPGeom:
    lines: list  # the file's lines, without newlines
    nodes: list  # [(x, y, z)]
    faces: list  # [[node index, ...]] 0-based
    parts: list  # part number of each face
    wake_line: int  # index in ``lines`` of the wake count line
    wakes: list = field(default_factory=list)  # [(part, [node index, ...])] 0-based

    def write(self, path, wakes=None):
        wakes = self.wakes if wakes is None else wakes
        n_old = int(self.lines[self.wake_line].split()[0])
        head = self.lines[: self.wake_line]
        tail = self.lines[self.wake_line + 1 + n_old:]
        body = [str(len(wakes))] + [
            " ".join([str(len(nodes)), str(part)] + [str(i + 1) for i in nodes]) for part, nodes in wakes
        ]
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(head + body + tail))
            if not tail or tail[-1] != "":
                f.write("\n")


def read(path) -> VSPGeom:
    with open(path) as f:
        lines = f.read().split("\n")
    if not lines[0].lower().startswith("# vspgeom v3"):
        raise ValueError(f"{path}: not a vspgeom v3 file")
    levels = int(lines[1].split()[0])
    if levels != 1:
        raise ValueError(f"{path}: {levels} mesh levels; only single-level files are supported")
    n_nodes, n_faces = map(int, lines[2].split()[:2])
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
    pos += n_faces  # tags: part, tag, then u v per vertex
    pos += n_faces  # face map: "i j"
    wake_line = pos
    n_wakes = int(lines[pos].split()[0])
    wakes = []
    for i in range(n_wakes):
        a = list(map(int, lines[pos + 1 + i].split()))
        wakes.append((a[1], [v - 1 for v in a[2:2 + a[0]]]))
    return VSPGeom(lines, nodes, faces, parts, wake_line, wakes)


# --- trailing-edge vertices -----------------------------------------------------------------------------------------

def match_vertices(g: VSPGeom, points, tolerance):
    """The mesh node at each point. Every point must lie within ``tolerance`` of a node (the nearest is taken)."""
    if tolerance <= 0:
        raise ValueError("the vertex tolerance must be positive")
    cell = {}
    for i, p in enumerate(g.nodes):
        cell.setdefault(tuple(int(c // tolerance) for c in p), []).append(i)
    found, missing = [], []
    for k, p in enumerate(points, start=1):
        base = tuple(int(c // tolerance) for c in p)
        best, best_d2 = None, tolerance * tolerance
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for i in cell.get((base[0] + dx, base[1] + dy, base[2] + dz), ()):
                        q = g.nodes[i]
                        d2 = (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2
                        if d2 <= best_d2:
                            best, best_d2 = i, d2
        if best is None:
            missing.append(f"{k} at ({p[0]:.6g}, {p[1]:.6g}, {p[2]:.6g})")
        else:
            found.append(best)
    if missing:
        more = f" and {len(missing) - 5} more" if len(missing) > 5 else ""
        raise ValueError(f"{len(missing)} trailing-edge point(s) are not vertices of the mesh (no vertex within "
                         f"{tolerance:g} m): " + "; ".join(missing[:5]) + more)
    return found


def vertex_wakes(g: VSPGeom, nodes):
    """Wakes along the mesh edges whose two ends are both trailing-edge nodes. Returns [(part, [node, ...])].

    The listed vertices must define the trailing edges without ambiguity: every vertex joins one or two others, no
    face has all its vertices listed, and no vertex is a junction of three or more trailing-edge edges.
    """
    marked = set(nodes)
    edges = {}
    for fi, f in enumerate(g.faces):
        k = len(f)
        on = [v in marked for v in f]
        if all(on):
            pts = ", ".join("({:.6g}, {:.6g}, {:.6g})".format(*g.nodes[v]) for v in f)
            raise ValueError(f"every vertex of a mesh face is a trailing-edge point, so the trailing edge through "
                             f"it is ambiguous: {pts}")
        for a in range(k):
            i, j = f[a], f[(a + 1) % k]
            if on[a] and on[(a + 1) % k]:
                edges.setdefault((min(i, j), max(i, j)), g.parts[fi])
    degree = {v: 0 for v in marked}
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    lonely = [v for v, d in degree.items() if d == 0]
    if lonely:
        p = g.nodes[lonely[0]]
        raise ValueError(f"{len(lonely)} trailing-edge vertex(es) share no mesh edge with another trailing-edge vertex, "
                         f"e.g. ({p[0]:.6g}, {p[1]:.6g}, {p[2]:.6g}); list every vertex along each trailing edge")
    forks = [v for v, d in degree.items() if d > 2]
    if forks:
        p = g.nodes[forks[0]]
        raise ValueError(f"{len(forks)} trailing-edge vertex(es) join three or more trailing-edge edges, e.g. "
                         f"({p[0]:.6g}, {p[1]:.6g}, {p[2]:.6g}); the listed vertices must form simple lines")
    return chains([(i, j, part) for (i, j), part in edges.items()])


def chains(edges):
    """Join edges into simple polylines (open or closed). Returns [(part, [node, ...])]."""
    adj = {}
    part_of = {}
    for i, j, part in edges:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
        part_of[(min(i, j), max(i, j))] = part
    used = set()

    def walk(start, nxt):
        path = [start, nxt]
        used.add((min(start, nxt), max(start, nxt)))
        while len(adj[path[-1]]) == 2:
            a, b = adj[path[-1]]
            step = b if a == path[-2] else a
            key = (min(path[-1], step), max(path[-1], step))
            if key in used:
                break
            used.add(key)
            path.append(step)
        return path

    out = []
    for v in sorted(adj):  # open chains start at their ends
        if len(adj[v]) == 1 and (min(v, adj[v][0]), max(v, adj[v][0])) not in used:
            out.append(walk(v, adj[v][0]))
    for v in sorted(adj):  # closed loops
        for w in adj[v]:
            if (min(v, w), max(v, w)) not in used:
                out.append(walk(v, w))
    return [(part_of[(min(p[0], p[1]), max(p[0], p[1]))], p) for p in out]


def orient(g: VSPGeom, wakes):
    """Order each wake's nodes so they run toward +y (VSP's convention is consistent per surface; this is stable)."""
    out = []
    for part, nodes in wakes:
        a, b = g.nodes[nodes[0]], g.nodes[nodes[-1]]
        if (b[1], b[2]) < (a[1], a[2]):
            nodes = nodes[::-1]
        out.append((part, nodes))
    return out


def add_wakes(path_in, path_out, points, tolerance, keep_existing=True):
    """Write ``path_in`` to ``path_out`` with wakes along the trailing edges through ``points`` (mesh vertices).

    Returns (number of wakes, number of wake edges, list of wake polylines as node coordinates).
    """
    g = read(path_in)
    found = orient(g, vertex_wakes(g, match_vertices(g, points, tolerance))) if points else []
    wakes = (g.wakes if keep_existing else []) + found
    g.write(path_out, wakes)
    polylines = [[g.nodes[i] for i in nodes] for _, nodes in wakes]
    return len(wakes), sum(len(n) - 1 for _, n in wakes), polylines
