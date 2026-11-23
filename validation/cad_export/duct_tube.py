"""Analytic round tube r=1 with tilted ends through vsp_ntop.duct; deviation of the OpenVSP surface."""
import math, os, sys, xml.etree.ElementTree as ET
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))  # nTop_vsp
sys.path.insert(0, ROOT)
from vsp_ntop import duct as D, vspmodel as M
from vsp_ntop.common import Log, run_vspscript
VSPDIR = os.path.join(ROOT, "OpenVSP-3.56.0-win64")
WORK = os.path.join(ROOT, "runs", "duct_tube").replace("\\", "/")
R, L = 1.0, 3.0
def rot_y(v, a):
    c, s = math.cos(math.radians(a)), math.sin(math.radians(a)); return (c*v[0] + s*v[2], v[1], -s*v[0] + c*v[2])
def rot_z(v, a):
    c, s = math.cos(math.radians(a)), math.sin(math.radians(a)); return (c*v[0] - s*v[1], s*v[0] + c*v[1], v[2])
def radius(x):  # a waisted duct: r(x)
    return R * (1.0 - 0.15 * math.sin(math.pi * x / L) ** 2) if WAIST else R
WAIST = int(os.environ.get("WAIST", "0"))
P_IN, N_IN = (0.0, 0.0, 0.0), rot_y((1, 0, 0), 15)
P_OUT, N_OUT = (L, 0.0, 0.0), rot_z((1, 0, 0), 25)

def cut(p, n, m=720):
    """Tube wall ∩ plane: for each angle, solve n.(P(x)-p)=0 for x by bisection."""
    out = []
    for i in range(m):
        t = 2 * math.pi * i / m
        f = lambda x: D.dot(n, D.add((x, radius(x) * math.cos(t), radius(x) * math.sin(t)), p, -1.0))
        a, b = -2.0, L + 2.0
        for _ in range(80):
            c = 0.5 * (a + b)
            if f(a) * f(c) <= 0: b = c
            else: a = c
        x = 0.5 * (a + b)
        out.append((x, radius(x) * math.cos(t), radius(x) * math.sin(t)))
    return out

def sdf_err(q):
    x, y, z = q
    a = -D.dot(D.unit(N_IN), q); b = D.dot(D.unit(N_OUT), D.add(q, P_OUT, -1.0))
    # distance to the side wall, for points between the ends
    return abs(math.hypot(y, z) - radius(x)), max(a, b)

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
    pts = [tuple(map(float, l.split()[1:4])) for l in text.splitlines() if l.strip().startswith("PT")]
    errs = [sdf_err(q) for q in pts]
    wall = max(e[0] for e in errs); out = max(e[1] for e in errs)
    bad = sorted(range(len(errs)), key=lambda i: -errs[i][0])[:6]; print([(i // 48, i % 48, round(errs[i][0]*1000, 2)) for i in bad])
    print(f"sections {nsec} C{cont} waist {WAIST}: wall max {wall*1000:.2f} mm, mean {1000*sum(e[0] for e in errs)/len(errs):.3f} mm; "
          f"farthest beyond an end {out*1000:.2f} mm, ends: {[round(1000*sdf_err(q)[1],2) for q in (pts[0], pts[-1])]}")

if __name__ == "__main__":
    for a in sys.argv[1:]:
        n, c = a.split(":"); main(int(n), int(c))
