"""Build OpenVSP Fuselage and Wing components from fitted sections (``fit.py``).

Fuselages are built in two steps, because every change to a section in OpenVSP's API rebuilds the whole surface
(about 1.7 s per section for a 30-section fuselage): a script makes the component with the right number of
file-defined sections (``fuselage_lines``), and ``set_fuselage_sections`` then writes each section's position,
shape and size straight into the saved ``.vsp3`` XML. Wings are few sections and are built by the script alone
(``wing_lines``). Components go into the user sets "Thick" and "Thin" (OpenVSP sets 3 and 4), which the VSPAERO
analysis hands to VSPAERO as its thick (panel) and thin (vortex lattice) sets. Standard library only.
"""
from __future__ import annotations

import math
import os
import xml.etree.ElementTree as ET

from . import fit as F

THICK_SET, THIN_SET = 3, 4  # OpenVSP's SET_FIRST_USER and the next one


def _f(x):
    return repr(float(x))


def clean_name(s):
    return "".join(c for c in s if c.isalnum() or c in " _-").strip() or "Component"


def _sets(gid, set_name):
    """Put the component in OpenVSP's Thick or Thin set (VSPAERO), or neither ("none": only Wave Drag and Parasite
    Drag, which take every component)."""
    return [f"    SetSetFlag( {gid}, {THICK_SET}, {'true' if set_name == 'thick' else 'false'} );",
            f"    SetSetFlag( {gid}, {THIN_SET}, {'true' if set_name == 'thin' else 'false'} );"]


def tess_split(lengths, total, minimum=2):
    """Share ``total`` tessellation points between sections in proportion to their lengths."""
    whole = sum(lengths) or 1.0
    return [max(minimum, int(round(total * l / whole)) + 1) for l in lengths]


def fuselage_lines(k, name, f: F.FuselageFit, set_name, tess_w):
    """Script lines making the fuselage with its sections (shapes only; ``set_fuselage_sections`` fills them)."""
    secs = f.sections
    n = len(secs)
    o = f.frame.origin
    gid = f"g{k}"
    L = [f"    // fuselage: {name}",
         f'    string {gid} = AddGeom( "FUSELAGE", "" );',
         f'    SetGeomName( {gid}, "{clean_name(name)}" );',
         f'    SetParmVal( {gid}, "X_Rel_Location", "XForm", {_f(o[0] + secs[0].x)} );',
         f'    SetParmVal( {gid}, "Y_Rel_Location", "XForm", {_f(o[1])} );',
         f'    SetParmVal( {gid}, "Z_Rel_Location", "XForm", {_f(o[2])} );',
         f'    SetParmVal( {gid}, "Sym_Planar_Flag", "Sym", 0 );',
         f'    SetParmVal( {gid}, "Length", "Design", {_f(secs[-1].x - secs[0].x)} );',
         f"    string xs{k} = GetXSecSurf( {gid}, 0 );",
         f"    while ( GetNumXSec( xs{k} ) > {n} ) {{ CutXSec( {gid}, 1 ); }}",
         f"    while ( GetNumXSec( xs{k} ) < {n} ) {{ InsertXSec( {gid}, GetNumXSec( xs{k} ) - 2, XS_ELLIPSE ); }}",
         f'    SetParmVal( {gid}, "Tess_W", "Shape", {int(tess_w)} );']
    for i, s in enumerate(secs):
        L.append(f"    ChangeXSecShape( xs{k}, {i}, {'XS_POINT' if s.is_point else 'XS_FILE_FUSE'} );")
    # A blunt end (a base, or a cut-off nose) is closed with a flat cap, so the surface stays closed for the panel
    # method and for Wave Drag's area distribution; OpenVSP leaves fuselage ends open by default.
    L.append(f'    SetParmVal( {gid}, "CapUMinOption", "EndCap", {0 if secs[0].is_point else 1} );')
    L.append(f'    SetParmVal( {gid}, "CapUMaxOption", "EndCap", {0 if secs[-1].is_point else 1} );')
    L += _sets(gid, set_name)
    return L


def _parm(group, name, value):
    e = group.find(name)
    if e is None:
        raise ValueError(f"OpenVSP model: parameter {name} not found")
    e.set("Value", f"{float(value):.18e}")


def axial_points(secs, tess_w, aspect):
    """Tessellation points along each section interval, so panels are about ``aspect`` times as long (along x) as
    they are wide (around the body): the width is the smaller of the two sections' perimeters over ``tess_w - 1``.
    VSPAERO's panel method goes wrong next to a wake when a body's panels are much longer than wide."""
    def perimeter(s):
        if s.is_point:
            return 0.0
        pts = [(a * s.width, b * s.height) for a, b in s.points]
        return sum(math.dist(p, q) for p, q in zip(pts, pts[1:]))
    out = []
    for a, b in zip(secs, secs[1:]):
        per = [p for p in (perimeter(a), perimeter(b)) if p > 0]
        width = min(per) / max(1, tess_w - 1)
        out.append(max(2, int(math.ceil((b.x - a.x) / (aspect * width))) + 1))
    return out


def skin_angles(secs):
    """Skinning angle (degrees, from the x axis, positive where the side moves away from the axis going aft) of the
    top, right, bottom and left sides at every section: the body's own slopes from the nTop slices around it (the fit's
    ``slopes``), or else central differences of the sections' extents; a rounded point (``tip_angle``) meets the axis
    square.

    Every angle is set explicitly: OpenVSP put back its default fixed angles (45 degrees at the first section, -45 at
    the last) when they were left free, which flared every blunt end out by about a tenth of its radius."""
    def sides(s):
        if s.is_point:
            return (s.zc, s.yc, -s.zc, -s.yc)
        return (s.zc + 0.5 * s.height, s.yc + 0.5 * s.width, -(s.zc - 0.5 * s.height), -(s.yc - 0.5 * s.width))
    ext = [sides(s) for s in secs]
    n = len(secs)
    out = []
    for i, s in enumerate(secs):
        if s.is_point and s.tip_angle is not None:
            out.append((s.tip_angle,) * 4)
            continue
        if s.slopes is not None:  # the body's own slopes, from the nTop slices around the section
            out.append(tuple(math.degrees(math.atan(v)) for v in s.slopes))
            continue
        a, b = (i - 1 if i > 0 else i), (i + 1 if i < n - 1 else i)
        dx = secs[b].x - secs[a].x
        out.append(tuple(math.degrees(math.atan2(ext[b][k] - ext[a][k], dx)) if dx > 0 else 0.0 for k in range(4)))
    return out


def set_fuselage_sections(vsp3_in, vsp3_out, fits, tess):
    """Write the fitted sections of each fuselage into a ``.vsp3`` made by ``fuselage_lines``.

    ``fits``: {component name: FuselageFit}; ``tess``: {name: (tess_w, panel aspect ratio, nose cluster,
    tail cluster)}."""
    tree = ET.parse(vsp3_in)
    done = set()
    for geom in tree.getroot().iter("Geom"):
        pc = geom.find("ParmContainer")
        fg = geom.find("FuselageGeom")
        if pc is None or fg is None:
            continue
        name = (pc.findtext("Name") or "").strip()
        if name not in fits:
            continue
        f, (tess_w, aspect, fwd, aft) = fits[name], tess[name]
        secs = f.sections
        xsecs = fg.find("XSecSurf").findall("XSec")
        if len(xsecs) != len(secs):
            raise ValueError(f"OpenVSP model: {name} has {len(xsecs)} sections, expected {len(secs)}")
        x0, length = secs[0].x, secs[-1].x - secs[0].x
        u_pts = axial_points(secs, tess_w, aspect)
        angles = skin_angles(secs)
        for i, (xe, s) in enumerate(zip(xsecs, secs)):
            g = xe.find("ParmContainer/XSec")
            _parm(g, "XLocPercent", (s.x - x0) / length)
            _parm(g, "YLocPercent", s.yc / length)
            _parm(g, "ZLocPercent", s.zc / length)
            # each side keeps its own angle: AllSym (OpenVSP's "all sides like the top") made the sides and bottom copy
            # the top's slope; on a body with different slopes per side (validation/fuselage_sides) the boattail's
            # mean error was 3x larger with it. NTVSP_SKIN_SIDES=allsym restores it, for comparison.
            sym = 1 if os.environ.get("NTVSP_SKIN_SIDES", "independent") == "allsym" else 0
            _parm(g, "AllSym", sym)
            if not sym:
                _parm(g, "TBSym", 0)
                _parm(g, "RLSym", 0)
            for side, angle in zip(("Top", "Right", "Bottom", "Left"), angles[i]):
                for lr in ("L", "R"):
                    _parm(g, f"{side}{lr}AngleSet", 1)
                    _parm(g, f"{side}{lr}Angle", angle)
                    for what in ("Strength", "Slew", "Curve"):
                        _parm(g, f"{side}{lr}{what}Set", 0)
            if i > 0:
                _parm(g, "SectTess_U", u_pts[i - 1])
                _parm(g, "FwdCluster", fwd if i == 1 else 1.0)
                _parm(g, "AftCluster", aft if i == len(secs) - 1 else 1.0)
            if s.is_point:
                continue
            curve = xe.find("XSec/XSecCurve")
            fx = curve.find("FileXSec")
            if fx is None:
                raise ValueError(f"OpenVSP model: {name} section {i} is not a file section")
            fx.find("UnityFilePnts").text = "".join(f"{a:.12e}, {b:.12e}, 0.0, " for a, b in s.points)
            cg = curve.find("ParmContainer/XSecCurve")
            _parm(cg, "Width", s.width)
            _parm(cg, "Height", s.height)
        done.add(name)
    missing = set(fits) - done
    if missing:
        raise ValueError(f"OpenVSP model: fuselage(s) {', '.join(sorted(missing))} not found in {vsp3_in}")
    tree.write(vsp3_out, encoding="utf-8", xml_declaration=True)


def set_duct_sections(vsp3_in, vsp3_out, fits):
    """Write each duct's fitted sections (``duct.write_sections``) into a ``.vsp3`` made by ``duct.duct_lines``.
    ``fits``: {component name: (DuctFit, tess_w, panel aspect ratio)}."""
    from . import duct as D
    tree = ET.parse(vsp3_in)
    done = set()
    for geom in tree.getroot().iter("Geom"):
        pc = geom.find("ParmContainer")
        if pc is None or geom.find("FuselageGeom") is None:
            continue
        name = (pc.findtext("Name") or "").strip()
        if name in fits:
            f, tess_w, aspect = fits[name]
            D.write_sections(geom, f, tess_w, aspect, parm=_parm)
            done.add(name)
    missing = set(fits) - done
    if missing:
        raise ValueError(f"OpenVSP model: duct(s) {', '.join(sorted(missing))} not found in {vsp3_in}")
    tree.write(vsp3_out, encoding="utf-8", xml_declaration=True)


SYMMETRY_FLAGS = {"none": "0", "xz": "SYM_XZ", "xy": "SYM_XY", "xzxy": "SYM_XZ | SYM_XY"}


def wing_lines(k, name, f: F.WingFit, set_name, symmetry, tess_u, tess_w, le_cluster=0.25, te_cluster=0.25,
               root_cluster=1.0, tip_cluster=1.0, sharp_te=True):
    """Script lines for one Wing. ``symmetry``: none, xz, xy or xzxy (OpenVSP's planar symmetry: mirrored about the
    x-z plane, the x-y plane, or both). ``tess_u``: spanwise points (shared between sections by span); ``tess_w``:
    points around each airfoil."""
    st = f.stations
    secs = F.wing_sections(f)
    fr = f.frame
    s_axis = fr.ez
    root = fr.point(st[0].le[0], st[0].le[1], st[0].z)
    cant = math.degrees(math.atan2(s_axis[2], s_axis[1]))
    gid = f"g{k}"
    L = [f"    // wing: {name}",
         f'    string {gid} = AddGeom( "WING", "" );',
         f'    SetGeomName( {gid}, "{clean_name(name)}" );',
         f'    SetParmVal( {gid}, "X_Rel_Location", "XForm", {_f(root[0])} );',
         f'    SetParmVal( {gid}, "Y_Rel_Location", "XForm", {_f(root[1])} );',
         f'    SetParmVal( {gid}, "Z_Rel_Location", "XForm", {_f(root[2])} );',
         f'    SetParmVal( {gid}, "X_Rel_Rotation", "XForm", {_f(cant)} );',
         f'    SetParmVal( {gid}, "Sym_Planar_Flag", "Sym", {SYMMETRY_FLAGS[symmetry]} );',
         f'    SetParmVal( {gid}, "RotateAirfoilMatchDideralFlag", "WingGeom", 0 );',
         f'    SetParmVal( {gid}, "CorrectAirfoilthicknessFlag", "WingGeom", 0 );',
         f'    SetParmVal( {gid}, "RelativeDihedralFlag", "WingGeom", 0 );',
         f'    SetParmVal( {gid}, "RelativeTwistFlag", "WingGeom", 0 );',
         f"    for ( int i = 1; i < {len(secs)}; i++ ) {{ InsertXSec( {gid}, 1, XS_FILE_AIRFOIL ); }}",
         f"    string xs{k} = GetXSecSurf( {gid}, 0 );",
         f'    SetParmVal( {gid}, "Tess_W", "Shape", {int(tess_w)} );',
         f'    SetParmVal( {gid}, "LECluster", "WingGeom", {_f(le_cluster)} );',
         f'    SetParmVal( {gid}, "TECluster", "WingGeom", {_f(te_cluster)} );',
         f'    SetParmVal( {gid}, "Twist_Location", "XSec_0", 0.0 );',
         f'    SetParmVal( {gid}, "Twist", "XSec_0", {_f(st[0].twist)} );']
    u_pts = tess_split([s["span"] for s in secs], tess_u)
    for i, s in enumerate(secs, start=1):
        g = f'"XSec_{i}"'
        L += [f"    SetDriverGroup( {gid}, {i}, SPAN_WSECT_DRIVER, ROOTC_WSECT_DRIVER, TIPC_WSECT_DRIVER );",
              f'    SetParmVal( {gid}, "Sweep_Location", {g}, 0.0 );',
              f'    SetParmVal( {gid}, "Twist_Location", {g}, 0.0 );',
              f'    SetParmVal( {gid}, "Root_Chord", {g}, {_f(s["root_chord"])} );',
              f'    SetParmVal( {gid}, "Tip_Chord", {g}, {_f(s["tip_chord"])} );',
              f'    SetParmVal( {gid}, "Span", {g}, {_f(s["span"])} );',
              f'    SetParmVal( {gid}, "Sweep", {g}, {_f(max(-89.0, min(89.0, s["sweep"])))} );',
              f'    SetParmVal( {gid}, "Dihedral", {g}, {_f(s["dihedral"])} );',
              f'    SetParmVal( {gid}, "Twist", {g}, {_f(s["twist"])} );',
              f'    SetParmVal( {gid}, "SectTess_U", {g}, {u_pts[i - 1]} );',
              f'    SetParmVal( {gid}, "InCluster", {g}, {_f(root_cluster if i == 1 else 1.0)} );',
              f'    SetParmVal( {gid}, "OutCluster", {g}, {_f(tip_cluster if i == len(secs) else 1.0)} );',
              "    Update();"]
    for i, a in enumerate(st):
        L.append(f"    ChangeXSecShape( xs{k}, {i}, XS_FILE_AIRFOIL );")
        L.append(f"    xc = GetXSec( xs{k}, {i} );")
        L.append("    { array< vec3d > up; array< vec3d > lo;")
        close = sharp_te and a.te_thickness > 0.0
        for (xu, yu), (xl, yl) in zip(a.upper, a.lower):
            if close:  # shear both surfaces onto the chord line at the trailing edge
                yu -= xu * a.upper[-1][1]
                yl -= xl * a.lower[-1][1]
            L.append(f"      up.push_back( vec3d( {xu:.9g}, {yu:.9g}, 0.0 ) ); lo.push_back( vec3d( {xl:.9g}, {yl:.9g}, 0.0 ) );")
        L.append("      SetAirfoilPnts( xc, up, lo ); }")
    L += _sets(gid, set_name)
    return L


def build_script(component_lines, vsp3_path):
    """The first script: make every component and save the model."""
    lines = ["void main()", "{", f'    SetSetName( {THICK_SET}, "Thick" );', f'    SetSetName( {THIN_SET}, "Thin" );',
             "    string xc;"]
    lines += component_lines
    lines += ["    Update();", f'    WriteVSPFile( "{vsp3_path}", SET_ALL );', '    Print( "NTVSP_OK" );', "}"]
    return "\n".join(lines) + "\n"


def finish_script(vsp3_in, vsp3_out, stl_path, samples):
    """The second script: load the edited model, save it, export its tessellation (STL) and print a grid of
    surface points of each component: ``samples`` = [(name, n_u, n_w)], printed as ``NTVSP_PT name a b x y z``."""
    lines = ["void main()", "{", f'    ReadVSPFile( "{vsp3_in}" );', "    Update();",
             f'    WriteVSPFile( "{vsp3_out}", SET_ALL );', f'    ExportFile( "{stl_path}", SET_ALL, EXPORT_STL );']
    for name, n_u, n_w in samples:
        nm = clean_name(name)
        lines += [f'    {{ array< string > @ids = FindGeomsWithName( "{nm}" );',
                  f'      if ( ids.size() == 0 ) {{ Print( "ERR component {nm} missing after reload" ); return; }}',
                  f"      for ( int a = 0; a <= {n_u}; a++ ) for ( int b = 0; b < {n_w}; b++ ) {{",
                  f"        vec3d q = CompPnt01( ids[0], 0, a / {float(n_u)!r}, b / {float(n_w)!r} );",
                  f'        Print( "NTVSP_PT {nm} " + a + " " + b + " " + formatFloat( q.x(), "", 0, 6 ) + " " + '
                  f'formatFloat( q.y(), "", 0, 6 ) + " " + formatFloat( q.z(), "", 0, 6 ) ); }} }}']
    lines += ['    Print( "NTVSP_OK" );', "}"]
    return "\n".join(lines) + "\n"
