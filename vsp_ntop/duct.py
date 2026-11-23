"""Ducts: the side wall of an implicit body between two planar ends, as an open OpenVSP Fuselage.

The body is closed (a flow passage, or a tube with a wall) and its two ends lie in planes the user gives: the
inlet plane and the outlet plane, not necessarily parallel. The duct is the body's side surface between them; the
end faces are left out, so the OpenVSP Fuselage is open at both ends (no caps).

Sections. nTop cuts the body with ``n`` planes that turn gradually from the inlet plane to the outlet plane: plane
``k`` passes through the point ``(1 - s) p_in + s p_out`` with the normal ``(1 - s) n_in + s n_out`` (normalised),
``s`` running 0 to 1 with OpenVSP's clustering function. Slicing normal to one axis instead would leave pieces
of the surface between the last full slice and a tilted end that are very short on one side of the duct and long
on the other, and OpenVSP's skinning spline overshoots there (5 cm on a 1 m tube with a 25 degree end). The end
sections are cut a little inside the body (``inset``, a fraction of the distance between the end planes, since a
cut exactly on an end face is degenerate) and moved onto their end plane along the duct's axis (inlet point to
outlet point): along the cutting normal instead, a tilted end's contour moved off the wall (1.35 mm on the 1 m
tube with a 25 degree end).

Each section is one contour of its cut: the outer one, or the hole for the inner wall of a tube. OpenVSP places a
section at its centre and turns it by XRotate, YRotate, ZRotate (rotation matrix Rx Ry Rz, found by testing
OpenVSP 3.56); the fit uses XRotate = 0, so the section's normal fixes the other two, and writes the contour in the
section's own y and z axes as a file cross-section.

Skinning. Every side of every section is left free (no angles, strengths or slews) with the top/bottom and
left/right symmetries and the left/right equalities off. OpenVSP's defaults couple the sides (each side copies the
top), which bowed a round tube with tilted ends out by 0.1 m on a 1 m radius; uncoupled, the same tube came back
within 0.2 mm. Standard library only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from . import fit as F
from . import slices as S

SIDES = ("Top", "Right", "Bottom", "Left")


def dot(a, b):
    return sum(p * q for p, q in zip(a, b))


def add(a, b, s=1.0):
    return tuple(p + s * q for p, q in zip(a, b))


def unit(a):
    n = math.sqrt(dot(a, a))
    if n == 0:
        raise F.FitError("zero-length direction")
    return tuple(c / n for c in a)


def rotation(n):
    """OpenVSP section rotation (XRotate, YRotate, ZRotate in degrees, XRotate = 0) turning the section's x axis
    onto ``n``, and its columns (the section's x, y and z axes in the world). OpenVSP's matrix is Rx Ry Rz, so with
    XRotate = 0 the x axis goes to (cos b cos g, sin g, -sin b cos g)."""
    n = unit(n)
    g = math.asin(max(-1.0, min(1.0, n[1])))
    b = math.atan2(-n[2], n[0])
    cb, sb, cg, sg = math.cos(b), math.sin(b), math.cos(g), math.sin(g)
    ex = (cb * cg, sg, -sb * cg)
    ey = (-cb * sg, cg, sb * sg)
    ez = (sb, 0.0, cb)
    return (0.0, math.degrees(b), math.degrees(g)), ex, ey, ez


def orient(n_in, n_out, p_in, p_out):
    """Both end normals pointing from the inlet towards the outlet."""
    d = add(p_out, p_in, -1.0)
    n_in = unit(n_in if dot(n_in, d) >= 0 else tuple(-c for c in n_in))
    n_out = unit(n_out if dot(n_out, d) >= 0 else tuple(-c for c in n_out))
    return n_in, n_out


def plane_fractions(n, a=1.0, b=1.0, inset=1e-3):
    """``s`` of the ``n`` cutting planes, clustered (OpenVSP's function: ``a`` scales the spacing at the inlet, ``b``
    at the outlet), the ends ``inset`` inside."""
    if n < 3:
        raise F.FitError("a duct needs at least 3 sections")
    s = [F.cluster(k / (n - 1), a, b) for k in range(n)]
    s[0], s[-1] = inset, 1.0 - inset
    return s


def section_planes(p_in, n_in, p_out, n_out, s_list):
    """(point, unit normal) of each cutting plane."""
    n_in, n_out = orient(n_in, n_out, p_in, p_out)
    out = []
    for s in s_list:
        p = add(tuple((1 - s) * c for c in p_in), p_out, s)
        n = unit(add(tuple((1 - s) * c for c in n_in), n_out, s))
        out.append((p, n))
    return out


@dataclass
class DuctSection:
    centre: tuple          # world position of the section's centre (bounding box in its own y, z)
    angles: tuple          # (XRotate, YRotate, ZRotate) degrees
    width: float
    height: float
    points: list = field(default_factory=list)  # [(y/width, z/height)] about the centre, closed
    loop3d: list = field(default_factory=list)  # the contour in the world (for checks)


@dataclass
class DuctFit:
    sections: list
    notes: list


def pick_loop(loops, wall):
    """The outer contour (largest area), or for ``wall == 'inner'`` the largest hole inside it."""
    if not loops:
        raise F.FitError("a cut has no contour (is the plane inside the body?)")
    ordered = sorted(loops, key=lambda l: abs(S.area(l)), reverse=True)
    if wall == "outer":
        return ordered[0]
    if len(ordered) < 2:
        raise F.FitError("Wall = Inner, but a cut has no hole (the body is not a tube there)")
    return ordered[1]


def fit_duct(cuts, end_planes, per_quarter=24):
    """``cuts``: per section, the chosen contour as world points [(x, y, z)] and its cutting plane (point, normal),
    inlet to outlet. ``end_planes``: ((p_in, n_in), (p_out, n_out)); the first and last contours are moved along
    their cutting normal onto these planes."""
    sections, notes = [], []
    axis = unit(add(end_planes[1][0], end_planes[0][0], -1.0))
    for k, (loop, (p, n)) in enumerate(cuts):
        if k in (0, len(cuts) - 1):
            pe, ne = end_planes[0 if k == 0 else 1]
            # along the duct's axis (inlet point to outlet point) onto the end plane, which becomes the section plane
            loop = [add(q, axis, dot(add(pe, q, -1.0), ne) / dot(axis, ne)) for q in loop]
            p, n = pe, unit(ne)
        angles, ex, ey, ez = rotation(n)
        local = [(dot(add(q, p, -1.0), ey), dot(add(q, p, -1.0), ez)) for q in loop]
        off = sum(dot(add(q, p, -1.0), ex) for q in loop) / len(loop)  # 0 up to rounding: the loop is in the plane
        uc, vc, w, h, pts, _ = F._section_curve(local, per_quarter)
        centre = add(add(add(p, ex, off), ey, uc), ez, vc)
        sections.append(DuctSection(centre, angles, w, h, pts, loop))
    return DuctFit(sections, notes)


def duct_lines(k, name, f: DuctFit, set_name, tess_w, sets):
    """Script lines making the open fuselage with one file section per cut (``set_duct_sections`` fills them)."""
    n = len(f.sections)
    gid = f"g{k}"
    L = [f"    // duct: {name}",
         f'    string {gid} = AddGeom( "FUSELAGE", "" );',
         f'    SetGeomName( {gid}, "{name}" );',
         f'    SetParmVal( {gid}, "OrderPolicy", "Design", FUSE_FREE );',
         f'    SetParmVal( {gid}, "Sym_Planar_Flag", "Sym", 0 );',
         f"    string xs{k} = GetXSecSurf( {gid}, 0 );",
         f"    while ( GetNumXSec( xs{k} ) > {n} ) {{ CutXSec( {gid}, 1 ); }}",
         f"    while ( GetNumXSec( xs{k} ) < {n} ) {{ InsertXSec( {gid}, GetNumXSec( xs{k} ) - 2, XS_ELLIPSE ); }}",
         f'    SetParmVal( {gid}, "Tess_W", "Shape", {int(tess_w)} );']
    for i in range(n):
        L.append(f"    ChangeXSecShape( xs{k}, {i}, XS_FILE_FUSE );")
    L.append(f'    SetParmVal( {gid}, "CapUMinOption", "EndCap", 0 );')
    L.append(f'    SetParmVal( {gid}, "CapUMaxOption", "EndCap", 0 );')
    L += sets(gid, set_name)
    return L


def placement(f: DuctFit):
    """Fuselage origin and length so that every section's XLocPercent lies in [0, 1]: origin at the smallest
    centre x (y and z of the first centre), length the largest extent of the centres along any axis."""
    cs = [s.centre for s in f.sections]
    x0 = min(c[0] for c in cs)
    origin = (x0, cs[0][1], cs[0][2])
    length = max(max(c[i] for c in cs) - min(c[i] for c in cs) for i in range(3))
    return origin, max(length, 1e-6)


def axial_points(f: DuctFit, tess_w, aspect):
    """Tessellation points along each section interval, as ``vspmodel.axial_points``: panels about ``aspect``
    times as long as wide; the length is the largest distance between the two contours' matching points."""
    def perim(s):
        pts = [(a * s.width, b * s.height) for a, b in s.points]
        return sum(math.dist(p, q) for p, q in zip(pts, pts[1:]))
    out = []
    for a, b in zip(f.sections, f.sections[1:]):
        width = min(perim(a), perim(b)) / max(1, tess_w - 1)
        gap = max(math.dist(a.centre, b.centre), 1e-12)
        out.append(max(2, int(math.ceil(gap / (aspect * width))) + 1))
    return out


def write_sections(geom, f: DuctFit, tess_w, aspect, continuity=2, parm=None):
    """Fill the fuselage ``geom`` (an XML ``Geom`` element made by ``duct_lines``) with the fitted sections."""
    origin, length = placement(f)
    pc = geom.find("ParmContainer")
    xf = pc.find("XForm")
    for key, v in zip(("X_Rel_Location", "Y_Rel_Location", "Z_Rel_Location"), origin):
        parm(xf, key, v)
    parm(pc.find("Design"), "Length", length)
    xsecs = geom.find("FuselageGeom").find("XSecSurf").findall("XSec")
    if len(xsecs) != len(f.sections):
        raise ValueError(f"OpenVSP model: duct has {len(xsecs)} sections, expected {len(f.sections)}")
    u_pts = axial_points(f, tess_w, aspect)
    for i, (xe, s) in enumerate(zip(xsecs, f.sections)):
        g = xe.find("ParmContainer/XSec")
        rel = add(s.centre, origin, -1.0)
        parm(g, "XLocPercent", rel[0] / length)
        parm(g, "YLocPercent", rel[1] / length)
        parm(g, "ZLocPercent", rel[2] / length)
        for key, v in zip(("XRotate", "YRotate", "ZRotate"), s.angles):
            parm(g, key, v)
        parm(g, "Spin", 0.0)
        for key in ("AllSym", "TBSym", "RLSym"):
            parm(g, key, 0)
        end = i in (0, len(f.sections) - 1)
        for side in SIDES:
            parm(g, f"Continuity{side}", 0 if end else continuity)
            for what in ("Angle", "Slew", "Strength", "Curve"):
                parm(g, f"{side}LR{what}Eq", 0)
                for lr in "LR":
                    parm(g, f"{side}{lr}{what}Set", 0)
        if i > 0:
            parm(g, "SectTess_U", u_pts[i - 1])
        curve = xe.find("XSec/XSecCurve")
        fx = curve.find("FileXSec")
        if fx is None:
            raise ValueError(f"OpenVSP model: duct section {i} is not a file section")
        fx.find("UnityFilePnts").text = "".join(f"{a:.12e}, {b:.12e}, 0.0, " for a, b in s.points)
        cg = curve.find("ParmContainer/XSecCurve")
        parm(cg, "Width", s.width)
        parm(cg, "Height", s.height)
