"""Planar end caps in an OpenVSP STEP file: each flat cap's patches replaced by one PLANE face, by editing the file.
"""
from __future__ import annotations

import math
import re

_ENTITY = re.compile(r"#(\d+)\s*=\s*(.*?);\s*(?=#\d+\s*=|$)", re.S)
_REF = re.compile(r"#(\d+)")


def _split(args):
    """Top-level comma-separated items of an argument string (no outer parentheses)."""
    out, depth, cur, quote = [], 0, [], False
    for ch in args:
        if ch == "'" :
            quote = not quote
        if not quote:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "," and depth == 0:
                out.append("".join(cur).strip())
                cur = []
                continue
        cur.append(ch)
    if cur:
        out.append("".join(cur).strip())
    return out


def _args(text, name=None):
    """Arguments of a simple entity ``NAME(...)``, or of the part ``name(...)`` of a complex one."""
    if name is None:
        i = text.index("(")
        body = text[i + 1:text.rindex(")")]
        return _split(body)
    m = re.search(r"\b" + name + r"\s*\(", text)
    if not m:
        return None
    i, depth = m.end(), 1
    j = i
    while depth:
        if text[j] == "(":
            depth += 1
        elif text[j] == ")":
            depth -= 1
        j += 1
    return _split(text[i:j - 1])


def _kind(text):
    t = text.strip()
    return "(complex)" if t.startswith("(") else t.split("(")[0].strip()


def _refs(text):
    return [int(r) for r in _REF.findall(text)]


def _sub(a, b):
    return tuple(p - q for p, q in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return sum(p * q for p, q in zip(a, b))


def _norm(a):
    return math.sqrt(_dot(a, a))


class Step:
    def __init__(self, text):
        self.text = text
        start = text.index("DATA;") + len("DATA;")
        end = text.index("ENDSEC;", start)
        self.head, self.tail = text[:start], text[end:]
        self.ents = {int(m.group(1)): m.group(2).strip() for m in _ENTITY.finditer(text[start:end])}
        self.next_id = max(self.ents) + 1

    def add(self, text):
        i = self.next_id
        self.next_id += 1
        self.ents[i] = text
        return i

    def point(self, i):
        a = _args(self.ents[i])
        return tuple(float(v) for v in a[1].strip("()").split(","))

    def surface_net(self, i):
        """Control points of a B-spline surface (simple or complex entity) as rows, or None."""
        text = self.ents[i]
        a = _args(text, "B_SPLINE_SURFACE_WITH_KNOTS") if _kind(text) == "B_SPLINE_SURFACE_WITH_KNOTS" else None
        if a is not None:
            net = a[3]
        else:
            b = _args(text, "B_SPLINE_SURFACE")
            if b is None:
                return None
            net = b[2]
        rows = re.findall(r"\(([^()]*)\)", net)
        return [[self.point(int(r)) for r in _REF.findall(row)] for row in rows]

    def curve_points(self, i):
        """Control points of an edge's 3D curve (through SURFACE_CURVE / SEAM_CURVE), in its parameter order."""
        text = self.ents[i]
        k = _kind(text)
        if k in ("SURFACE_CURVE", "SEAM_CURVE"):
            return self.curve_points(_refs(_args(text)[1])[0])
        a = _args(text, "B_SPLINE_CURVE_WITH_KNOTS") if k == "B_SPLINE_CURVE_WITH_KNOTS" else _args(text, "B_SPLINE_CURVE")
        if a is None:
            return []
        lst = a[2] if k == "B_SPLINE_CURVE_WITH_KNOTS" else a[1]
        return [self.point(r) for r in _refs(lst)]

    def write(self):
        body = "".join(f"#{i}={t};\n" for i, t in sorted(self.ents.items()))
        return self.head + "\n" + body + self.tail


def _plane_of(pts, tol):
    """(unit normal, offset) of the plane through ``pts`` if they all lie in one (within ``tol``), else None."""
    c = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
    far = max(pts, key=lambda p: _norm(_sub(p, c)))
    best, n = 0.0, None
    for p in pts:
        w = _cross(_sub(far, c), _sub(p, c))
        m = _norm(w)
        if m > best:
            best, n = m, tuple(x / m for x in w)
    if n is None:
        return None
    off = _dot(n, c)
    if max(abs(_dot(n, p) - off) for p in pts) > tol:
        return None
    return n, off


def _collapsed(net, tol):
    """Whether a boundary row or column of the control net is a single point."""
    rows = [net[0], net[-1], [r[0] for r in net], [r[-1] for r in net]]
    return any(max(_norm(_sub(p, row[0])) for p in row) <= tol for row in rows)


def planar_caps(text, log=print):
    """The STEP file ``text`` with every flat end cap one plane face; returns (new text, caps rebuilt)."""
    s = Step(text)
    ents = s.ents
    referenced = set()
    for t in ents.values():
        referenced.update(_refs(t))
    roots = [i for i in ents if i not in referenced]
    caps, removed_surfaces = 0, set()
    for sid in [i for i, t in ents.items() if _kind(t) in ("CLOSED_SHELL", "OPEN_SHELL")]:
        sargs = _args(ents[sid])
        faces = _refs(sargs[1])
        info = {}
        for f in faces:
            fa = _args(ents[f])
            if _kind(ents[f]) != "ADVANCED_FACE":
                continue
            bounds = _refs(fa[1])
            surf = _refs(fa[2])[0]
            net = s.surface_net(surf)
            if not net:
                continue
            pts = [p for row in net for p in row]
            ext = max(_norm(_sub(p, pts[0])) for p in pts)
            tol = 1e-7 * max(ext, 1e-12)
            pl = _plane_of(pts, tol * 10)
            if pl is None:
                continue
            oes = [oe for b in bounds for oe in _refs(_args(ents[_refs(_args(ents[b])[1])[0]])[1])]
            edges = [_refs(_args(ents[oe])[3])[0] for oe in oes]
            info[f] = {"plane": pl, "singular": _collapsed(net, tol), "oes": oes, "edges": edges, "surf": surf,
                       "ext": ext}
        # how each face of the shell runs each of its edges
        runs = {}
        for f in faces:
            if _kind(ents[f]) != "ADVANCED_FACE":
                continue
            for b in _refs(_args(ents[f])[1]):
                for oe in _refs(_args(ents[_refs(_args(ents[b])[1])[0]])[1]):
                    a = _args(ents[oe])
                    runs.setdefault(_refs(a[3])[0], []).append((f, a[4].strip() == ".T."))
        groups, seen = [], set()
        for f in [f for f in info if info[f]["singular"]]:
            if f in seen:
                continue
            group, stack = [], [f]
            seen.add(f)
            while stack:
                g = stack.pop()
                group.append(g)
                for h in info:
                    if h in seen:
                        continue
                    (ng, og), (nh, oh) = info[g]["plane"], info[h]["plane"]
                    scale = max(info[g]["ext"], info[h]["ext"])
                    coplanar = abs(abs(_dot(ng, nh)) - 1.0) < 1e-6 and abs(og - (oh if _dot(ng, nh) > 0 else -oh)) < 1e-6 * scale
                    if coplanar and set(info[g]["edges"]) & set(info[h]["edges"]):
                        seen.add(h)
                        stack.append(h)
            groups.append(group)
        new_faces = {}
        in_caps = {f for g in groups for f in g}
        for group in groups:
            use = {}
            for f in group:
                for e in info[f]["edges"]:
                    use[e] = use.get(e, 0) + 1
            rim_edges = list(dict.fromkeys(e for f in group for e in info[f]["edges"] if use[e] == 1))
            verts = {}
            for e in rim_edges:
                ec = _args(ents[e])
                verts[e] = (_refs(ec[1])[0], _refs(ec[2])[0])
            # chain the rim into closed loops by shared vertices (one, or an outer loop and holes where another part
            # passes through the cap), each as [(edge, runs from its start vertex)]
            loops, left, ok = [], list(rim_edges), True
            while left and ok:
                e0 = left.pop(0)
                chain = [(e0, True)]
                start_v, tail = verts[e0]
                while tail != start_v:
                    nxt = next((e for e in left if tail in verts[e]), None)
                    if nxt is None:
                        ok = False
                        break
                    fwd = verts[nxt][0] == tail
                    chain.append((nxt, fwd))
                    tail = verts[nxt][1] if fwd else verts[nxt][0]
                    left.remove(nxt)
                loops.append(chain)
            if not ok or not loops:
                log(f"planar caps: the rim of a cap in shell #{sid} does not close; left as it is")
                continue
            n, off = info[group[0]]["plane"]

            def travel(chain):
                """Area vector and centroid of a loop (points along its edges, in its direction of travel)."""
                poly = []
                for e, fwd in chain:
                    ec = _args(ents[e])
                    pts = s.curve_points(_refs(ec[3])[0])
                    if (ec[4].strip() == ".T.") != fwd:
                        pts = pts[::-1]
                    poly += pts[:-1] if len(pts) > 1 else pts
                c = tuple(sum(p[i] for p in poly) / len(poly) for i in range(3))
                area = (0.0, 0.0, 0.0)
                for p, q in zip(poly, poly[1:] + poly[:1]):
                    area = tuple(x + y for x, y in zip(area, _cross(_sub(p, c), _sub(q, c))))
                return area, c

            def settle(chain):
                """Turn a loop to run each edge opposite to a side face that shares it (side faces, unlike OpenVSP's
                cap patches, are oriented consistently); None when no edge of the loop has one."""
                for e, fwd in chain:
                    side = [f_fwd for g, f_fwd in runs.get(e, []) if g not in in_caps]
                    if len(side) == 1:
                        return chain if fwd != side[0] else [(x, not d) for x, d in reversed(chain)]
                return None

            settled = [settle(ch) for ch in loops]
            measured = [travel(ch) for ch in loops]
            outer = max(range(len(loops)), key=lambda k: abs(_dot(measured[k][0], n)))
            if settled[outer] is None:
                log(f"planar caps: a cap in shell #{sid} touches no side face; left as it is")
                continue
            loops[outer] = settled[outer]
            area, c = travel(loops[outer])
            if _dot(area, n) < 0:  # the outer bound runs counter-clockwise about the face normal
                n = tuple(-x for x in n)
            for k in range(len(loops)):
                if k != outer:  # holes run clockwise about it
                    loops[k] = settled[k] if settled[k] is not None else loops[k]
                    if _dot(travel(loops[k])[0], n) > 0:
                        loops[k] = [(x, not d) for x, d in reversed(loops[k])]
            ref = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
            xd = _cross(_cross(n, ref), n)
            xd = tuple(x / _norm(xd) for x in xd)
            name = _args(ents[info[group[0]]["surf"]])[0] if _kind(ents[info[group[0]]["surf"]]) != "(complex)" else "''"
            name = name[:-1] + " planar cap'" if name.endswith("'") and len(name) > 2 else "'planar cap'"
            p0 = s.add(f"CARTESIAN_POINT('',({c[0]!r},{c[1]!r},{c[2]!r}))")
            d0 = s.add(f"DIRECTION('',({n[0]!r},{n[1]!r},{n[2]!r}))")
            d1 = s.add(f"DIRECTION('',({xd[0]!r},{xd[1]!r},{xd[2]!r}))")
            ax = s.add(f"AXIS2_PLACEMENT_3D('',#{p0},#{d0},#{d1})")
            pl = s.add(f"PLANE({name},#{ax})")
            bounds = []
            for k, chain in enumerate(loops):
                oes = [s.add(f"ORIENTED_EDGE('',*,*,#{e},{'.T.' if fwd else '.F.'})") for e, fwd in chain]
                lp = s.add("EDGE_LOOP('',(" + ",".join(f"#{oe}" for oe in oes) + "))")
                kind = "FACE_OUTER_BOUND" if k == outer else "FACE_BOUND"
                bounds.append(s.add(f"{kind}('',#{lp},.T.)"))
            face = s.add("ADVANCED_FACE('',(" + ",".join(f"#{b}" for b in bounds) + f"),#{pl},.T.)")
            new_faces[group[0]] = face
            for f in group[1:]:
                new_faces[f] = None
            removed_surfaces.update(info[f]["surf"] for f in group)
            caps += 1
        if new_faces:
            kept = []
            for f in faces:
                if f in new_faces:
                    if new_faces[f] is not None:
                        kept.append(new_faces[f])
                else:
                    kept.append(f)
            sargs[1] = "(" + ",".join(f"#{f}" for f in kept) + ")"
            ents[sid] = _kind(ents[sid]) + "(" + ",".join(sargs) + ")"
    if not caps:
        return text, 0
    # rim edges' SURFACE_CURVEs
    bare = {}  # SURFACE_CURVEs left with no parameter-space curve (an edge between two caps)
    for i, t in list(ents.items()):
        if _kind(t) in ("SURFACE_CURVE", "SEAM_CURVE"):
            a = _args(t)
            pcs = [r for r in _refs(a[2]) if not (set(_refs(ents[r])) & removed_surfaces)]
            if not pcs:
                bare[i] = _refs(a[1])[0]
            elif len(pcs) != len(_refs(a[2])):
                a[2] = "(" + ",".join(f"#{r}" for r in pcs) + ")"
                ents[i] = _kind(t) + "(" + ",".join(a) + ")"
    for i, t in list(ents.items()):
        if _kind(t) == "EDGE_CURVE":
            a = _args(t)
            g = _refs(a[3])[0]
            if g in bare:
                a[3] = f"#{bare[g]}"
                ents[i] = "EDGE_CURVE(" + ",".join(a) + ")"
    # keep what the original roots still reach
    keep, stack = set(), list(roots)
    while stack:
        i = stack.pop()
        if i in keep or i not in ents:
            continue
        keep.add(i)
        stack.extend(_refs(ents[i]))
    s.ents = {i: t for i, t in ents.items() if i in keep}
    return s.write(), caps


def planar_caps_file(src, dst, log=print):
    text = open(src, encoding="utf-8", errors="replace").read()
    out, caps = planar_caps(text, log)
    with open(dst, "w", encoding="utf-8") as f:
        f.write(out)
    return caps
