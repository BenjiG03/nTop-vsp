"""Ducts: the side wall of an implicit body between two planar ends, as an open OpenVSP Fuselage."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import fit as F
from . import slices as S
from .geometry import arc_length, unit as _unit

SIDES = ("Top", "Right", "Bottom", "Left")


def unit(a):
    return _unit(a, F.FitError)


def rotation(n):
    """OpenVSP section rotation (XRotate, YRotate, ZRotate in degrees, XRotate = 0) turning the section's x axis onto
    ``n``, and its columns (the section's x, y and z axes in the world)."""
    n = unit(n)
    g = np.arcsin(np.clip(n[1], -1.0, 1.0))
    b = np.arctan2(-n[2], n[0])
    cb, sb, cg, sg = np.cos(b), np.sin(b), np.cos(g), np.sin(g)
    return ((0.0, np.degrees(b), np.degrees(g)), np.array([cb * cg, sg, -sb * cg]), np.array([-cb * sg, cg, sb * sg]),
            np.array([sb, 0.0, cb]))


def orient(n_in, n_out, p_in, p_out):
    """Both end normals pointing from the inlet towards the outlet."""
    d = np.subtract(p_out, p_in)
    return tuple(unit(np.asarray(nn, float) if np.dot(nn, d) >= 0 else -np.asarray(nn, float)) for nn in (n_in, n_out))


def plane_fractions(n, a=1.0, b=1.0, inset=1e-3):
    """``s`` of the ``n`` cutting planes, clustered (OpenVSP's function: ``a`` scales the spacing at the inlet, ``b``
    at the outlet), the ends ``inset`` inside."""
    if n < 3:
        raise F.FitError("a duct needs at least 3 sections")
    s = F.cluster(np.arange(n) / (n - 1), a, b)
    s[0], s[-1] = inset, 1.0 - inset
    return s


def section_planes(p_in, n_in, p_out, n_out, s_list):
    """(point, unit normal) of each cutting plane."""
    n_in, n_out = orient(n_in, n_out, p_in, p_out)
    p_in, p_out = np.asarray(p_in, float), np.asarray(p_out, float)
    return [((1 - s) * p_in + s * p_out, unit((1 - s) * n_in + s * n_out)) for s in s_list]


@dataclass
class DuctSection:
    centre: np.ndarray     # world position of the section's centre (bounding box in its own y, z)
    angles: tuple          # (XRotate, YRotate, ZRotate) degrees
    width: float
    height: float
    points: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))  # (y/width, z/height) about the centre
    loop3d: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))  # the contour in the world (for checks)


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
    """``cuts``: per section, the chosen contour as (n, 3) world points and its cutting plane (point, normal),
    inlet to outlet. ``end_planes``: ((p_in, n_in), (p_out, n_out))."""
    sections = []
    axis = unit(np.subtract(end_planes[1][0], end_planes[0][0]))
    for k, (loop, (p, n)) in enumerate(cuts):
        loop = np.asarray(loop, float)
        if k in (0, len(cuts) - 1):  # along the duct's axis onto the end plane
            pe, ne = (np.asarray(v, float) for v in end_planes[0 if k == 0 else 1])
            loop = loop + np.outer((pe - loop) @ ne / (axis @ ne), axis)
            p, n = pe, unit(ne)
        angles, ex, ey, ez = rotation(n)
        rel = loop - p
        uc, vc, w, h, pts, _ = F._section_curve(rel @ np.column_stack([ey, ez]), per_quarter)
        centre = p + (rel @ ex).mean() * ex + uc * ey + vc * ez  # the mean along ex is 0 up to rounding
        sections.append(DuctSection(centre, angles, w, h, pts, loop))
    return DuctFit(sections, [])


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
    cs = np.array([s.centre for s in f.sections])
    origin = np.array([cs[:, 0].min(), cs[0, 1], cs[0, 2]])
    return origin, max(np.ptp(cs, axis=0).max(), 1e-6)


def axial_points(f: DuctFit, tess_w, aspect):
    """Tessellation points along each section interval, as ``vspmodel.axial_points``: panels about ``aspect``
    times as long as wide; the length is the distance between the two sections' centres."""
    per = np.array([arc_length(s.points * [s.width, s.height])[-1] for s in f.sections])
    gap = np.maximum(np.linalg.norm(np.diff([s.centre for s in f.sections], axis=0), axis=1), 1e-12)
    width = np.minimum(per[:-1], per[1:]) / max(1, tess_w - 1)
    return np.maximum(2, np.ceil(gap / (aspect * width)).astype(int) + 1).tolist()


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
        rel = s.centre - origin
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
