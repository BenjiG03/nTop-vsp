"""Lofted CAD (vsp_ntop.loft + step_write) against OpenVSP's surface, on bodies with known exact surfaces.

    sh build/py.sh validation/loft_cad/test_loft.py [build] [measure] [bodies...]

Each body is one nTop expression (mm), so the same text is the nTop implicit and, in numpy, the exact field here.
build: one nTop notebook per body (OpenVSP Fuselage + OpenVSP Model blocks: OpenVSP's surface, and the slices
as CLI). measure: loft the slices, write runs/loft_cad/<body>.stp, and score the lofted surface and OpenVSP's against
the exact field: distance |f|/|grad f| and the angle between the surface normal and grad f, by region along the body.

Bodies (x from 0, mm):
  sides   elliptic sections, rounded plan, conical side view, drooped nose, top-only boattail, flat base at 600
          (validation/fuselage_sides; reuses its notebook)
  chine   lens sections (two intersecting circles) with sharp chines along the whole body, rounded nose, flat base
  canopy  round fuselage with a canopy (ellipsoid) blended on with a 6 mm smooth union: concave sections
  blob    rounded-square sections (superellipse n = 4) on a centreline curving in y and z, round at both ends
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "validation", "rasaero"))
sys.path.insert(0, os.path.join(ROOT, "build"))
sys.path.insert(0, ROOT)
import rocket as rk  # noqa: E402
from vsp_ntop import fit as F, loft as L, slices as S, step_write as SW  # noqa: E402

WORK = os.path.join(ROOT, "runs", "loft_cad").replace("\\", "/")
VSP = os.path.join(ROOT, "OpenVSP-3.56.0-win64", "vspscript.exe")


def clamp01(e):
    return f"min(max({e}, 0), 1)"


def smin(a, b, k):
    """Polynomial smooth minimum (blend radius about k)."""
    return f"(min({a}, {b}) - max({k} - sqrt(({a} - {b})^2), 0)^2/(4*{k}))"


def sides_expr():
    t1, tb, tt = clamp01("X/200"), clamp01("X/150"), clamp01("(X-450)/150")
    a = f"max(40*sqrt(max(1-(1-{t1})^2, 0)), 0.001)"
    b = f"max(25*{tb} - 10*{tt}, 0.001)"
    zc = f"(-12*(1-{t1})^2 - 10*{tt})"
    f = f"(sqrt((Y/{a})^2 + ((Z-{zc})/{b})^2) - 1)*min({a}, {b})"
    return f"max(max({f}, -X), X - 600)"


def chine_expr():
    s = f"max(sqrt(max(1-(1-{clamp01('X/300')})^2, 0)), 0.002)"
    R, c = f"(40*{s})", f"(24*{s})"   # half-width 0.8 R = 32, half-height 0.4 R = 16 at full size
    lens = f"max(sqrt(Y^2 + (Z+{c})^2) - {R}, sqrt(Y^2 + (Z-{c})^2) - {R})"
    return f"max(max({lens}, -X), X - 500)"


def canopy_expr():
    r = f"(30*sqrt(max(1-(1-{clamp01('X/150')})^2, 0))*(1 - 0.33*{clamp01('(X-500)/100')}))"
    body = f"(sqrt(Y^2 + Z^2) - {r})"
    q = "sqrt(((X-170)/90)^2 + (Y/14)^2 + ((Z-22)/16)^2)"
    canopy = f"(({q} - 1)*14)"
    return f"max(max({smin(body, canopy, 6)}, -X), X - 600)"


def blob_expr():
    s = f"max(sqrt(max(1-(2*X/600-1)^2, 0)), 0.002)"
    a, b = f"(25*{s})", f"(18*{s})"
    yc, zc = "(8*sin(2*3.14159265358979*X/600))", "(15*sin(3.14159265358979*X/600))"
    f = f"((((Y-{yc})/{a})^4 + ((Z-{zc})/{b})^4)^0.25 - 1)*min({a}, {b})"
    return f"max(max({f}, -X), X - 600)"


BODIES = {  # name: (expression, bounding box mm, length mm)
    "sides": (sides_expr, ((-5, -45, -40), (605, 45, 30)), 600.0),
    "chine": (chine_expr, ((-5, -40, -25), (505, 40, 25)), 500.0),
    "canopy": (canopy_expr, ((-5, -35, -35), (605, 35, 45)), 600.0),
    "blob": (blob_expr, ((-5, -40, -25), (605, 40, 40)), 600.0),
}


def ntop_text(expr):
    for k, v in (("X", "(x/1mm)"), ("Y", "(y/1mm)"), ("Z", "(z/1mm)")):
        expr = expr.replace(k, v)
    return expr


def numpy_field(expr):
    """The expression as a numpy function of mm arrays (two-argument min/max as in nTop)."""
    code = expr.replace("^", "**")
    env = {"min": np.minimum, "max": np.maximum, "sqrt": lambda v: np.sqrt(np.maximum(v, 0.0)), "sin": np.sin,
           "pi": math.pi}
    return eval("lambda X, Y, Z: " + code, env)


def build(name):
    import model_blocks as mb
    from ntopfile import Graph
    expr, (lo, hi), length = BODIES[name]
    os.makedirs(WORK, exist_ok=True)
    fus, model = mb.fuselage_block(), mb.model_block()
    g = Graph()
    field = g.leaf("real_field", {"expression": ntop_text(expr()), "units": {}})
    scaled = g.block("multiply<real_field,real_field>", "real_field", [field, g.real(1e-3, rk.LENGTH)], top=False)
    box = g.block("create_bounding_box<point,point>", "bounding_box", [g.point(*lo), g.point(*hi)], top=False)
    body = g.block("set_bounding_box<real_field,bounding_box>[1.1.0]", "implicit", [scaled, box], f"{name} body")
    comp = rk.call(g, fus, {"Bodies": g.block("core.list<implicit>", "list<implicit>", [body], top=False),
                            "Names": g.block("core.list<text>", "list<text>", [g.text(name)], top=False),
                            "Slice Spacing": g.real(0.0025)}, f"{name} component")
    out = rk.call(g, model, {"Components": g.block("core.list<text>", "list<text>", [comp], top=False),
                             "Working Directory": g.file_path(WORK), "Case Name": g.text(name)}, "Model")
    res, log, dt = rk.run(g, [fus, model], out, "text", os.path.join(WORK, f"{name}.ntop"))
    print(f"{name}: build {dt:.0f} s, output {'yes' if res else 'no'}")
    for line in log.splitlines():
        if "[E]" in line:
            print("  ", line[:300])


def ntop_block(names=("chine", "canopy")):
    """End to end in nTop: the bodies in one OpenVSP Fuselage block, into Lofted CAD Export: one STEP file."""
    import model_blocks as mb
    import cad_blocks as cb
    from ntopfile import Graph
    fus, lb = mb.fuselage_block(), cb.loft_block()
    g = Graph()
    bodies = []
    for name in names:
        expr, (lo, hi), _ = BODIES[name]
        field = g.leaf("real_field", {"expression": ntop_text(expr()), "units": {}})
        scaled = g.block("multiply<real_field,real_field>", "real_field", [field, g.real(1e-3, rk.LENGTH)], top=False)
        box = g.block("create_bounding_box<point,point>", "bounding_box", [g.point(*lo), g.point(*hi)], top=False)
        bodies.append(g.block("set_bounding_box<real_field,bounding_box>[1.1.0]", "implicit", [scaled, box],
                              f"{name} body"))
    comp = rk.call(g, fus, {"Bodies": g.block("core.list<implicit>", "list<implicit>", bodies, top=False),
                            "Names": g.block("core.list<text>", "list<text>", [g.text(n) for n in names], top=False),
                            "Slice Spacing": g.real(0.0025)}, "Fuselages")
    out = rk.call(g, lb, {"Components": g.block("core.list<text>", "list<text>", [comp], top=False),
                          "Working Directory": g.file_path(WORK + "/block"), "Case Name": g.text("bodies")},
                  "Lofted CAD")
    res, log, dt = rk.run(g, [fus, lb], out, "text", os.path.join(WORK, "block", "bodies.ntop"))
    print(f"nTop block: {dt:.0f} s, output {res}")
    for l in log.splitlines():
        if "[E]" in l:
            print("  ", l[:300])


def settings_path(name):
    if name == "sides":
        return os.path.join(ROOT, "runs", "fuselage_sides", "sides.ntvsp.settings.txt")
    return os.path.join(WORK, f"{name}.ntvsp.settings.txt")


def load(name):
    comp = [l for l in open(settings_path(name)) if l.startswith("component")][0].split("=", 1)[1]
    d = {k.strip(): v.strip() for k, _, v in (kv.partition(":") for kv in comp.split(";"))}
    o = (float(d["ox"]), float(d["oy"]), float(d["oz"]))
    return S.read_cli(d["files"]), F.fuselage_frame(o), (float(d["bx0"]) - o[0], float(d["bx1"]) - o[0])


SAMPLE = """void main()
{{
    ReadVSPFile( "{model}" );
    Update();
    array< string > @ids = FindGeomsWithName( "{name}" );
    for ( int a = 0; a <= {nu}; a++ ) for ( int b = 0; b < {nw}; b++ ) {{
        vec3d q = CompPnt01( ids[0], 0, a / {nu}.0, b / {nw}.0 );
        Print( "PT " + q.x() + " " + q.y() + " " + q.z() ); }}
}}
"""


def openvsp_grid(name, nu=1200, nw=96):
    model = settings_path(name).replace(".ntvsp.settings.txt", ".vsp3").replace("\\", "/")
    script = os.path.join(WORK, f"{name}.sample.vspscript")
    open(script, "w").write(SAMPLE.format(model=model, name="Sides" if name == "sides" else name, nu=nu, nw=nw))
    out = subprocess.run([VSP, "-script", script], capture_output=True, text=True).stdout
    pts = [[float(v) for v in l.split()[1:4]] for l in out.splitlines() if l.strip().startswith("PT")]
    return np.array(pts).reshape(nu + 1, nw, 3) * 1000.0


def loft_grid(lf, nu=1200, per=24):
    u = np.linspace(0, 1, nu + 1)
    return np.concatenate([L.evaluate(p, u, np.linspace(*p.w_range, per + 1)[:-1]) for p in lf.patches], axis=1) * 1000


def chine_edges(X, Y, Z):
    """Within 1.5 mm of the chine body's sharp edges (y = +-32 s(x), z = 0), where the normal is undefined."""
    s = np.sqrt(np.maximum(1 - (1 - np.clip(X / 300, 0, 1)) ** 2, 0))
    return np.hypot(np.abs(Y) - 32 * s, Z) < 1.5


EXCLUDE = {"chine": chine_edges}


def score(grid, field, length, exclude=None):
    """Distance and normal angle against the exact field over a (nu, nw, 3) grid of surface points (mm), interior
    rows only (finite differences need neighbours), by region along the body."""
    P = grid[1:-1]
    du = grid[2:] - grid[:-2]
    dw = np.roll(grid, -1, axis=1)[1:-1] - np.roll(grid, 1, axis=1)[1:-1]
    n = np.cross(du, dw)
    nn = np.linalg.norm(n, axis=2)
    ok = nn > 1e-12
    h = 1e-3
    X, Y, Z = P[..., 0], P[..., 1], P[..., 2]
    f = field(X, Y, Z)
    g = np.stack([(field(X + h, Y, Z) - field(X - h, Y, Z)), (field(X, Y + h, Z) - field(X, Y - h, Z)),
                  (field(X, Y, Z + h) - field(X, Y, Z - h))], axis=-1) / (2 * h)
    gn = np.linalg.norm(g, axis=2)
    dist = np.abs(f) / np.maximum(gn, 1e-12)
    cosang = np.abs(np.einsum("ijk,ijk->ij", n, g)) / np.maximum(nn * gn, 1e-30)
    ang = np.degrees(np.arccos(np.clip(cosang, 0, 1)))
    out = {}
    for label, lo, hi in (("ends 2%", 0.0, 0.02), ("body", 0.02, 0.98), ("ends 2%", 0.98, 1.0)):
        m = ok & (X >= lo * length) & (X < hi * length) & (X > 0.3) & (X < length - 0.3)
        if exclude is not None:
            m &= ~exclude(X, Y, Z)
        key = "ends" if label.startswith("ends") else "body"
        prev = out.get(key, (np.array([]), np.array([])))
        out[key] = (np.concatenate([prev[0], dist[m]]), np.concatenate([prev[1], ang[m]]))
    return out


def line(tag, s):
    parts = []
    for key in ("body", "ends"):
        d, a = s[key]
        parts.append(f"{key}: dist mean {d.mean():.3f} p99 {np.percentile(d, 99):.3f} max {d.max():.3f} mm, "
                     f"normal rms {np.sqrt(np.mean(a ** 2)):.2f} p99 {np.percentile(a, 99):.2f} deg")
    return f"   {tag:8s} " + " | ".join(parts)


def measure(name, nw=48):
    expr, _, length = BODIES[name]
    field = numpy_field(expr())
    layers, frame, bounds = load(name)
    t = time.time()
    lf = L.loft_body(layers, frame, nw=nw, bounds=bounds)
    dt = time.time() - t
    stp = os.path.join(WORK, f"{name}.stp")
    info = SW.write_loft(stp, lf, name)
    print(f"{name}: loft {dt:.1f} s, {lf.rows} rows, fit residual rms {lf.residual[0] * 1e3:.3f} max "
          f"{lf.residual[1] * 1e3:.3f} mm; ends {'pointed' if lf.start.pointed else 'blunt'}/"
          f"{'pointed' if lf.end.pointed else 'blunt'} at x = {lf.start.point[0] * 1e3:.2f}, {lf.end.point[0] * 1e3:.2f} mm; "
          f"STEP {info['faces']} faces, {os.path.getsize(stp) / 1e6:.1f} MB" + (f"; {'; '.join(lf.notes)}" if lf.notes else ""))
    lg, vg = loft_grid(lf), openvsp_grid(name)
    print(line("loft", score(lg, field, length)))
    print(line("OpenVSP", score(vg, field, length)))
    if name in EXCLUDE:
        print("   away from the sharp edges (1.5 mm):")
        print(line("loft", score(lg, field, length, EXCLUDE[name])))
        print(line("OpenVSP", score(vg, field, length, EXCLUDE[name])))


if __name__ == "__main__":
    args = sys.argv[1:]
    if "block" in args:
        ntop_block()
        sys.exit()
    steps = [a for a in args if a in ("build", "measure")] or ["build", "measure"]
    names = [a for a in args if a in BODIES] or list(BODIES)
    for n in names:
        if "build" in steps and n != "sides":
            build(n)
        if "measure" in steps:
            measure(n)
