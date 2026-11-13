"""Read a VSPAERO ``.vspgeom`` (v3) panel file, give it trailing-edge wakes, and prepare thin surfaces.

VSP writes a ``.vspgeom`` file from an imported tessellation (a MeshGeom), but it marks no trailing edges on it:
VSP's own wake finder only knows the trailing edges of its parametric wings. Without wakes VSPAERO sheds no
vorticity and the surface carries almost no lift. The user lists the mesh vertices that lie on trailing edges;
this module turns the mesh edges joining those vertices into wakes and writes them into the file's wake section, in
the same form VSP uses for wing trailing edges:

    <number of wakes>
    <number of nodes> <part number> <node 1> <node 2> ...      (one list per wake, nodes 1-based)

Files are written with Windows line ends, as OpenVSP writes them: VSPAERO's reader of the polygon layout (see
``polygon_layout``) crashes on files with Unix line ends.

Only single-level files (no mesh refinement levels) are handled; that is what VSP writes for a tessellation.
Pure standard library, so any Python 3 runs it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

EOL = "\r\n"  # OpenVSP's line end


@dataclass
class VSPGeom:
    lines: list  # the file's lines, without line ends
    nodes: list  # [(x, y, z)]
    faces: list  # [[node index, ...]] 0-based
    parts: list  # part number of each face
    wake_line: int  # index in ``lines`` of the wake count line
    wake_end: int  # index in ``lines`` of the first line after the wake lists
    wakes: list = field(default_factory=list)  # [(part, [node index, ...])] 0-based

    def tag_lines(self):
        """The face tag lines: ``part tag u v u v ...``."""
        start = 3 + len(self.nodes) + 1 + len(self.faces)
        return self.lines[start:start + len(self.faces)]

    def write(self, path, wakes=None):
        wakes = self.wakes if wakes is None else wakes
        head = self.lines[: self.wake_line]
        tail = self.lines[self.wake_end:]
        while tail and tail[-1] == "":
            tail = tail[:-1]
        with open(path, "w", newline="") as f:
            f.write(EOL.join(head + wake_lines(wakes) + tail) + EOL)


def wake_lines(wakes):
    return [str(len(wakes))] + [
        " ".join([str(len(nodes)), str(part)] + [str(i + 1) for i in nodes]) for part, nodes in wakes
    ]


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
    pos += 1
    # Each wake list is "count part node ..."; OpenVSP wraps long lists over several lines, so read tokens.
    wakes = []
    tokens = []
    for _ in range(n_wakes):
        while len(tokens) < 2:
            tokens += lines[pos].split()
            pos += 1
        count, part = abs(int(tokens[0])), int(tokens[1])
        while len(tokens) < 2 + count:
            tokens += lines[pos].split()
            pos += 1
        wakes.append((part, [int(v) - 1 for v in tokens[2:2 + count]]))
        tokens = tokens[2 + count:]
    return VSPGeom(lines, nodes, faces, parts, wake_line, pos, wakes)


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


# --- thin (vortex lattice) surfaces ---------------------------------------------------------------------------------

def set_thick(vkey_path, thick):
    """Mark every part in the ``.vkey`` tag key file as thick (panel method) or thin (vortex lattice).

    VSPAERO chooses the method per part from the ``thick`` column of the part table
    (``# part#,geom#,surf#,gname,gid,thick,plate,copy#,geomcopy#``). OpenVSP writes nothing for a mesh in its thin
    geometry set (it builds thin surfaces only from its parametric wings), so a thin mesh goes through the thick set
    and is re-marked here. Returns the number of parts changed.
    """
    with open(vkey_path) as f:
        lines = f.read().split("\n")
    head = next((i for i, l in enumerate(lines) if l.startswith("#") and "thick" in l.split(",")), None)
    if head is None:
        raise ValueError(f"{vkey_path}: no part table with a thick column")
    col = [c.strip() for c in lines[head].lstrip("#").split(",")].index("thick")
    changed = 0
    for i in range(head + 1, len(lines)):
        cells = lines[i].split(",")
        if len(cells) <= col:
            break
        if cells[col].strip() != str(int(thick)):
            cells[col] = str(int(thick))
            lines[i] = ",".join(cells)
            changed += 1
    with open(vkey_path, "w", newline="") as f:
        f.write(EOL.join(lines))
    return changed


def free_edges(g: VSPGeom):
    """Number of mesh edges used by one face only (an open sheet has them; a closed body has none)."""
    count = {}
    for f in g.faces:
        for a in range(len(f)):
            i, j = f[a], f[(a + 1) % len(f)]
            key = (min(i, j), max(i, j))
            count[key] = count.get(key, 0) + 1
    return sum(1 for c in count.values() if c == 1)


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _quad_score(g: VSPGeom, quad):
    """Cosine of the dihedral angle between the quad's two triangles, if the quad is strictly convex; else None."""
    p = [g.nodes[i] for i in quad]
    normals = [_cross(_sub(p[(k + 1) % 4], p[k]), _sub(p[(k + 3) % 4], p[k])) for k in range(4)]
    if any(_dot(normals[0], n) <= 0.0 for n in normals[1:]):
        return None  # not convex (or degenerate)
    n1 = _cross(_sub(p[1], p[0]), _sub(p[2], p[0]))
    n2 = _cross(_sub(p[2], p[0]), _sub(p[3], p[0]))
    l1, l2 = _dot(n1, n1) ** 0.5, _dot(n2, n2) ** 0.5
    if l1 == 0.0 or l2 == 0.0:
        return None
    return _dot(n1, n2) / (l1 * l2)


def polygon_layout(path_in, path_out):
    """Rewrite a triangle panel file in VSPAERO's polygon layout, so VSPAERO solves on the triangles as given.

    VSPAERO solves on grid level 1. For a file of triangles it builds that level itself, pairing triangles (wake
    triangles too) into quads. Behind a swept or curved trailing edge the wake cells are skewed, and the pairing can
    join two wake triangles across a trailing vortex line, which removes that line's edge from level 1: VSPAERO then
    stops with "Could not find global edge for trailing wake!" (seen on thin surfaces; VSPAERO 7.2.2,
    ``VSP_AGGLOM::CreateMixedMesh_`` and ``VSP_GEOM::StoreWakeKuttaEdges``).

    OpenVSP's own wings avoid this: their files hold polygons (quads) followed by a triangulation of each polygon,
    and VSPAERO then takes the polygons as level 1 and the triangles as level 0, and builds the wake as quads. That
    layout is only recognised when at least one face has more than three nodes, so two neighbouring, nearly coplanar
    triangles are written as one quad at level 1; level 0 keeps every triangle exactly as given.

    Returns the original face index of each level-0 triangle, in the order VSPAERO stores its results.
    """
    g = read(path_in)
    if any(len(f) != 3 for f in g.faces):
        raise ValueError(f"{path_in}: the panel file is not all triangles")
    nf = len(g.faces)
    avoid = {i for _, nodes in g.wakes for i in nodes}
    edge_faces = {}
    for fi, f in enumerate(g.faces):
        for a in range(3):
            i, j = f[a], f[(a + 1) % 3]
            edge_faces.setdefault((min(i, j), max(i, j)), []).append(fi)
    for (i, j), fs in edge_faces.items():
        if len(fs) == 1:
            avoid |= {i, j}  # nodes on free edges (sheet boundaries)
    best = None
    for (i, j), fs in edge_faces.items():
        if len(fs) != 2 or g.parts[fs[0]] != g.parts[fs[1]]:
            continue
        a, b = sorted(fs)
        fa, fb = g.faces[a], g.faces[b]
        p = next(v for v in fa if v not in (i, j))
        q = next(v for v in fb if v not in (i, j))
        if {p, q, i, j} & avoid:
            continue
        k = fa.index(p)
        quad = [p, fa[(k + 1) % 3], q, fa[(k + 2) % 3]]  # keeps face a's orientation
        score = _quad_score(g, quad)
        if score is not None and (best is None or score > best[0] + 1e-12 or (abs(score - best[0]) <= 1e-12 and b - a < best[3] - best[2])):
            best = (score, quad, a, b)
            if score > 1.0 - 1e-12 and b == a + 1:
                break
    if best is None:
        raise ValueError("no two neighbouring triangles away from trailing and free edges form a convex quad")
    _, quad, a, b = best
    tags = [t.split() for t in g.tag_lines()]
    order = [fi for fi in range(nf) if fi != b]  # level-1 faces; the quad sits at face a
    fine = []
    coarse, coarse_tags, alt, alt_tags = [], [], [], []
    for n, fi in enumerate(order, start=1):
        part, tag = tags[fi][0], tags[fi][1]
        group = [a, b] if fi == a else [fi]
        fine += group
        nodes = quad if fi == a else g.faces[fi]
        coarse.append(f"{len(nodes)} " + " ".join(str(v + 1) for v in nodes))
        # u, v are not used by these blocks (OpenVSP writes none that mean anything for a mesh): zeros.
        coarse_tags.append(f"{part} {tag} " + " ".join(["0"] * (2 * len(nodes))))
        alt.append(f"{n} {len(group)} " + " ".join(str(v + 1) for t in group for v in g.faces[t]))
        alt_tags.append(f"{n} {part} {tag} " + " ".join(["0"] * (6 * len(group))))
    nc = len(order)
    lines = ["# vspgeom v3", "1", f"{len(g.nodes)} {nc} 1"] + g.lines[3:3 + len(g.nodes)] + [str(nc)]
    lines += coarse + coarse_tags + [f"{n} {n}" for n in range(1, nc + 1)]
    lines += wake_lines(g.wakes) + alt + alt_tags
    with open(path_out, "w", newline="") as f:
        f.write(EOL.join(lines) + EOL)
    return fine
