"""Write a lofted body (``loft.Loft``) as a STEP file: one closed solid of B-spline faces and flat caps."""
from __future__ import annotations

import datetime
import math

MM = 1000.0  # file unit: millimetres.


def real(v):
    """A STEP REAL: always with a decimal point ("1." not "1"), exponent as E."""
    s = repr(float(v)).upper()
    if s in ("NAN", "INF", "-INF"):
        raise ValueError("non-finite coordinate")
    if "." not in s:
        s = s.replace("E", ".E") if "E" in s else s + "."
    return s


def text(s):
    """A STEP STRING: single quotes, a quote inside doubled."""
    return "'" + str(s).replace("'", "''") + "'"


def ref_list(ids):
    return "(" + ",".join(f"#{i}" for i in ids) + ")"


class StepFile:
    """The DATA section being built: ``add`` appends a record and returns its number."""

    def __init__(self):
        self.records = []

    def add(self, body):
        self.records.append(body)
        return len(self.records)

    # --- geometry
    def point(self, p):
        return self.add(f"CARTESIAN_POINT('',({','.join(real(c * MM) for c in p)}))")

    def direction(self, d):
        n = math.sqrt(sum(c * c for c in d))
        return self.add(f"DIRECTION('',({','.join(real(c / n) for c in d)}))")

    def placement(self, origin, axis, ref):
        """AXIS2_PLACEMENT_3D: a local frame (origin, z axis, x axis); planes and the model origin use it."""
        return self.add(f"AXIS2_PLACEMENT_3D('',#{self.point(origin)},#{self.direction(axis)},#{self.direction(ref)})")

    @staticmethod
    def _knots(knots):
        """Full knot vector -> (distinct knots, multiplicities), the form STEP stores."""
        distinct, mult = [], []
        for k in knots:
            if distinct and abs(k - distinct[-1]) <= 1e-14:
                mult[-1] += 1
            else:
                distinct.append(float(k))
                mult.append(1)
        return distinct, mult

    def bspline_curve(self, ctrl, knots, degree=3):
        """B_SPLINE_CURVE_WITH_KNOTS(name, degree, control points, curve form, closed?, self-intersecting?,
        multiplicities, knots, knot type)."""
        pts = [self.point(p) for p in ctrl]
        k, m = self._knots(knots)
        return self.add(f"B_SPLINE_CURVE_WITH_KNOTS('',{degree},{ref_list(pts)},.UNSPECIFIED.,.F.,.F.,"
                        f"({','.join(map(str, m))}),({','.join(real(v) for v in k)}),.UNSPECIFIED.)")

    def bspline_surface(self, ctrl, knots_u, knots_w, degree=3):
        """B_SPLINE_SURFACE_WITH_KNOTS(name, u degree, v degree, control net as rows of u, surface form, u closed?,
        v closed?, self-intersecting?, u multiplicities, v multiplicities, u knots, v knots, knot type)."""
        rows = []
        for row in ctrl:
            rows.append(ref_list(self.point(p) for p in row))
        ku, mu = self._knots(knots_u)
        kw, mw = self._knots(knots_w)
        return self.add(f"B_SPLINE_SURFACE_WITH_KNOTS('',{degree},{degree},({','.join(rows)}),.UNSPECIFIED.,"
                        f".F.,.F.,.F.,({','.join(map(str, mu))}),({','.join(map(str, mw))}),"
                        f"({','.join(real(v) for v in ku)}),({','.join(real(v) for v in kw)}),.UNSPECIFIED.)")

    # --- topology
    def vertex(self, p):
        return self.add(f"VERTEX_POINT('',#{self.point(p)})")

    def edge(self, v0, v1, curve):
        """EDGE_CURVE(name, start vertex, end vertex, curve, same sense): the edge runs along the curve from v0 to
        v1 (.T. = in the curve's own direction)."""
        return self.add(f"EDGE_CURVE('',#{v0},#{v1},#{curve},.T.)")

    def face(self, name, edges, surface, outer=True):
        """ADVANCED_FACE(name, bounds, surface, same sense). ``edges``: [(edge, forward?)] in loop order, the
        loop running anticlockwise seen from outside the solid (the surface normal must point outwards)."""
        oe = [self.add(f"ORIENTED_EDGE('',*,*,#{e},{'.T.' if fwd else '.F.'})") for e, fwd in edges]
        loop = self.add(f"EDGE_LOOP('',{ref_list(oe)})")
        bound = self.add(f"FACE_OUTER_BOUND('',#{loop},.T.)")
        return self.add(f"ADVANCED_FACE({text(name)},(#{bound}),#{surface},.T.)")

    def write(self, path, name, solids, origin_placement):
        """The product wrapper and units around the solids (one part holding them all), then the file."""
        ctx_app = self.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
        self.add(f"APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{ctx_app})")
        ctx_prod = self.add(f"PRODUCT_CONTEXT('',#{ctx_app},'mechanical')")
        product = self.add(f"PRODUCT({text(name)},{text(name)},'',(#{ctx_prod}))")
        self.add(f"PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{product}))")
        formation = self.add(f"PRODUCT_DEFINITION_FORMATION('','',#{product})")
        ctx_def = self.add(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{ctx_app},'design')")
        definition = self.add(f"PRODUCT_DEFINITION('design','',#{formation},#{ctx_def})")
        shape = self.add(f"PRODUCT_DEFINITION_SHAPE('','',#{definition})")
        # units: complex instances (several entity types at once, listed alphabetically inside one record)
        length = self.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
        angle = self.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
        solid_angle = self.add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
        tol = self.add(f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),#{length},'distance_accuracy_value',"
                       "'closure tolerance')")
        context = self.add(f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{tol}))"
                           f"GLOBAL_UNIT_ASSIGNED_CONTEXT((#{length},#{angle},#{solid_angle}))"
                           "REPRESENTATION_CONTEXT('Context3D','3D'))")
        rep = self.add(f"ADVANCED_BREP_SHAPE_REPRESENTATION({text(name)},{ref_list(list(solids) + [origin_placement])},"
                       f"#{context})")
        self.add(f"SHAPE_DEFINITION_REPRESENTATION(#{shape},#{rep})")
        stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        with open(path, "w", encoding="ascii", newline="\n") as f:
            f.write("ISO-10303-21;\nHEADER;\n")
            f.write(f"FILE_DESCRIPTION(('nTop_vsp lofted body: B-spline faces fitted to nTop slices'),'2;1');\n")
            f.write(f"FILE_NAME({text(path.replace(chr(92), '/').split('/')[-1])},'{stamp}',(''),(''),'',"
                    "'nTop_vsp step_write','');\n")
            f.write("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n")
            for i, r in enumerate(self.records, 1):
                f.write(f"#{i}={r};\n")
            f.write("ENDSEC;\nEND-ISO-10303-21;\n")


def write_loft(path, lf, name="Body"):
    """Write ``lf`` (loft.Loft) as one closed solid. Returns counts for the log."""
    return write_lofts(path, [(lf, name)], name)


def write_lofts(path, bodies, part_name="Bodies"):
    """Write several lofted bodies [(loft.Loft, name), ...] into one STEP file: one part, one closed solid each.
    Returns counts for the log."""
    sf = StepFile()
    solids, info = [], {"faces": 0, "caps": 0, "flipped": []}
    for lf, name in bodies:
        solid, faces, caps, flip = _solid(sf, lf, name)
        solids.append(solid)
        info["faces"] += faces
        info["caps"] += caps
        info["flipped"].append(flip)
    origin = sf.placement((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    sf.write(path, part_name, solids, origin)
    info["records"] = len(sf.records)
    info["solids"] = len(solids)
    return info


def _solid(sf, lf, name):
    """Add one lofted body's MANIFOLD_SOLID_BREP to ``sf``; returns (record, faces, caps, w reversed?)."""
    patches = lf.patches
    ez = lf.frame.ez
    # orientation: the side faces' normal (dS/du x dS/dw) must point out of the body.
    flip = _normal_points_in(lf)
    if flip:
        patches = _reverse_w(patches)

    # vertices: the corners of the patches at each end (one shared tip vertex at a pointed end)
    def end_vertices(row, end):
        if end.pointed:
            v = sf.vertex(patches[0].ctrl[row][0])
            return [v] * len(patches)
        return [sf.vertex(p.ctrl[row][0]) for p in patches]

    v_start = end_vertices(0, lf.start)
    v_end = end_vertices(-1, lf.end)
    # seam q: column w = start of patch q, along the body
    seams = [sf.edge(v_start[q], v_end[q], sf.bspline_curve([r[0] for r in p.ctrl], p.knots_u))
             for q, p in enumerate(patches)]

    def end_edges(row, verts, end):
        if end.pointed:
            return None
        return [sf.edge(verts[q], verts[(q + 1) % len(patches)], sf.bspline_curve(p.ctrl[row], p.knots_w))
                for q, p in enumerate(patches)]

    e_start = end_edges(0, v_start, lf.start)
    e_end = end_edges(-1, v_end, lf.end)

    faces = []
    for q, p in enumerate(patches):
        surf = sf.bspline_surface(p.ctrl, p.knots_u, p.knots_w)
        # anticlockwise in (u, w) = anticlockwise about the outward normal dS/du x dS/dw
        loop = [(seams[q], True)]
        if e_end:
            loop.append((e_end[q], True))
        loop.append((seams[(q + 1) % len(patches)], False))
        if e_start:
            loop.append((e_start[q], False))
        faces.append(sf.face(f"{name} side {q + 1}", loop, surf))
    caps = 0
    for end, edges, sign, label in ((lf.end, e_end, 1.0, "end"), (lf.start, e_start, -1.0, "start")):
        if not edges:
            continue
        # the cap's normal points out of the body (+ez at the tail, -ez at the nose).
        axis = tuple(sign * c for c in ez)
        ref = _perpendicular(axis)
        plane = sf.add(f"PLANE('',#{sf.placement(end.point, axis, ref)})")
        if sign > 0:
            loop = [(e, False) for e in reversed(edges)]
        else:
            loop = [(e, True) for e in edges]
        faces.append(sf.face(f"{name} {label} cap", loop, plane))
        caps += 1
    shell = sf.add(f"CLOSED_SHELL('',{ref_list(faces)})")
    solid = sf.add(f"MANIFOLD_SOLID_BREP({text(name)},#{shell})")
    return solid, len(faces), caps, flip


def _perpendicular(a):
    t = (1.0, 0.0, 0.0) if abs(a[0]) < 0.9 else (0.0, 1.0, 0.0)
    c = (a[1] * t[2] - a[2] * t[1], a[2] * t[0] - a[0] * t[2], a[0] * t[1] - a[1] * t[0])
    n = math.sqrt(sum(x * x for x in c))
    return tuple(x / n for x in c)


def _normal_points_in(lf):
    """True when dS/du x dS/dw points into the body, judged at the middle of every patch against the section's
    centre."""
    from . import loft as L
    import numpy as np
    score = 0.0
    for p in lf.patches:
        a, b = p.w_range
        u, w, h = 0.5, 0.5 * (a + b), 1e-4
        s = L.evaluate(p, [u - h, u, u + h], [w - h, w, w + h])
        su = s[2, 1] - s[0, 1]
        sw = s[1, 2] - s[1, 0]
        n = np.cross(su, sw)
        ring = np.concatenate([L.evaluate(q, [u], np.linspace(*q.w_range, 9))[0] for q in lf.patches])
        score += float(np.dot(n, s[1, 1] - ring.mean(axis=0)))
    return score < 0.0


def _reverse_w(patches):
    """The same surface with w running the other way round: patches in reverse order, each mirrored in w."""
    from . import loft as L
    out = []
    for p in reversed(patches):
        a, b = p.w_range
        knots = [1.0 - k for k in p.knots_w[::-1]]  # w -> 1 - w
        out.append(L.Patch(p.ctrl[:, ::-1].copy(), p.knots_u, __import__("numpy").array(knots), (1.0 - b, 1.0 - a)))
    return out
