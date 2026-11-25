"""Analytic round tube r=1 with tilted ends through vsp_ntop.duct; deviation of the OpenVSP surface."""
import os, sys, xml.etree.ElementTree as ET
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))  # nTop_vsp
sys.path.insert(0, ROOT)
from vsp_ntop import duct as D, vspmodel as M
from vsp_ntop.common import Log, run_vspscript
VSPDIR = os.path.join(ROOT, "OpenVSP-3.56.0-win64")
WORK = os.path.join(ROOT, "runs", "duct_tube").replace("\\", "/")
R, L = 1.0, 3.0
WAIST = int(os.environ.get("WAIST", "0"))


def rot(axis, deg):
    """Rotation matrix about the y or z axis."""
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]) if axis == "y" else np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def radius(x):  # a waisted duct: r(x)
    return R * (1.0 - 0.15 * np.sin(np.pi * x / L) ** 2) if WAIST else R * np.ones_like(x)


P_IN, N_IN = np.zeros(3), rot("y", 15) @ [1.0, 0.0, 0.0]
P_OUT, N_OUT = np.array([L, 0.0, 0.0]), rot("z", 25) @ [1.0, 0.0, 0.0]


def wall(x, t):
    r = radius(x)
    return np.stack([x, r * np.cos(t), r * np.sin(t)], axis=-1)


def cut(p, n, m=720):
    """Tube wall and plane: for each angle, n.(P(x) - p) = 0 solved for x by bisection (all angles at once)."""
    t = 2 * np.pi * np.arange(m) / m
    f = lambda x: (wall(x, t) - p) @ n  # noqa: E731
    a, b = np.full(m, -2.0), np.full(m, L + 2.0)
    for _ in range(80):
        c = 0.5 * (a + b)
        left = f(a) * f(c) <= 0
        a, b = np.where(left, a, c), np.where(left, c, b)
    return wall(0.5 * (a + b), t)


def sdf_err(q):
    """(distance to the side wall, farthest distance beyond an end) of points q (n, 3)."""
    beyond = np.maximum(-q @ D.unit(N_IN), (q - P_OUT) @ D.unit(N_OUT))
    return np.abs(np.hypot(q[:, 1], q[:, 2]) - radius(q[:, 0])), beyond


def main(nsec, cont):
    os.makedirs(WORK, exist_ok=True)
    s = D.plane_fractions(nsec)
    planes = D.section_planes(P_IN, N_IN, P_OUT, N_OUT, s)
    cuts = [(cut(p, n), (p, n)) for p, n in planes]
    f = D.fit_duct(cuts, ((P_IN, N_IN), (P_OUT, N_OUT)), int(os.environ.get("PQ", "24")))
    lines = D.duct_lines(0, "Duct", f, "none", 33, M._sets)
    log = Log(os.path.join(WORK, "log.txt"))
    st1 = WORK + "/a.vsp3"
    run_vspscript(VSPDIR, M.build_script(lines, st1), WORK, "build", log, 300)
    tree = ET.parse(st1)
    geom = [g for g in tree.getroot().iter("Geom") if g.find("FuselageGeom") is not None][0]
    D.write_sections(geom, f, 33, 1.5, continuity=cont, parm=M._parm)
    tree.write(WORK + "/b.vsp3", encoding="utf-8", xml_declaration=True)
    sc = ["void main()", "{", f'    ReadVSPFile( "{WORK}/b.vsp3" );', "    Update();",
          '    string si = FindContainer("SurfaceIntersectSettings", 0);',
          '    SetParmVal(FindParm(si, "CADLenUnit", "ExportIntersect"), LEN_M);',
          f'    ExportFile("{WORK}/duct_{nsec}_{cont}.stp", SET_ALL, EXPORT_STEP_STITCH);',
          '    array< string > @ids = FindGeomsWithName( "Duct" );',
          "    for ( int a = 0; a <= 600; a++ ) for ( int b = 0; b < 48; b++ ) { vec3d q = CompPnt01( ids[0], 0, a / 600.0, b / 48.0 );",
          '      Print( "PT " + q.x() + " " + q.y() + " " + q.z() ); }', '    Print("NTVSP_OK");', "}"]
    text = run_vspscript(VSPDIR, "\n".join(sc) + "\n", WORK, "sample", log, 300)
    pts = np.array([l.split()[1:4] for l in text.splitlines() if l.strip().startswith("PT")], float)
    side, out = sdf_err(pts)
    bad = np.argsort(-side, kind="stable")[:6]
    print([(int(i) // 48, int(i) % 48, round(float(side[i]) * 1000, 2)) for i in bad])
    print(f"sections {nsec} C{cont} waist {WAIST}: wall max {side.max()*1000:.2f} mm, mean {1000*side.mean():.3f} mm; "
          f"farthest beyond an end {out.max()*1000:.2f} mm, ends: {[round(1000*float(v), 2) for v in out[[0, -1]]]}")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        n, c = a.split(":"); main(int(n), int(c))
