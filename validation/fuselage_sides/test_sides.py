"""OpenVSP fuselage fit of a body whose top, bottom and sides have different slopes, against its exact surface.

    sh build/py.sh validation/fuselage_sides/test_sides.py [build] [fit]

The body (mm, x from 0 to 600, flat base): elliptic sections with half-width a(x) = 40 sqrt(1 - (1 - x/200)^2)
(rounded in plan) up to x = 200, half-height b(x) = 25 x/150 (conical in side view) up to x = 150, the centre drooped
by -12 (1 - x/200)^2, and over the last 150 mm a top-only boattail (the top comes down 20 mm, the bottom stays).
So at any station the four sides have different slopes, and a fit that copies the top's skinning angle to every
side (OpenVSP's AllSym) gets the sides and bottom wrong.

build: one nTop notebook (OpenVSP Fuselage + OpenVSP Model). fit: re-runs the model driver on the same slices with
each skinning variant and measures OpenVSP's surface against the exact ellipse at each station. Runs:
runs/fuselage_sides.
"""
from __future__ import annotations

import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))
import rocket as rk  # noqa: E402
import model_blocks as mb  # noqa: E402
from ntopfile import Graph  # noqa: E402

WORK = os.path.join(rk.ROOT, "runs", "fuselage_sides").replace("\\", "/")
L = 600.0


def shape(x):
    """(half-width, half-height, centre z) at x mm."""
    t1 = min(max(x / 200.0, 0.0), 1.0)
    tb = min(max(x / 150.0, 0.0), 1.0)
    tt = min(max((x - 450.0) / 150.0, 0.0), 1.0)
    a = 40.0 * math.sqrt(max(1.0 - (1.0 - t1) ** 2, 0.0))
    b = 25.0 * tb - 10.0 * tt
    zc = -12.0 * (1.0 - t1) ** 2 - 10.0 * tt
    return a, b, zc


def expression():
    t1 = "min(max(X/200, 0), 1)"
    tb = "min(max(X/150, 0), 1)"
    tt = "min(max((X-450)/150, 0), 1)"
    a = f"max(40*sqrt(max(1-(1-{t1})^2, 0)), 0.001)"
    b = f"max(25*{tb} - 10*{tt}, 0.001)"
    zc = f"(-12*(1-{t1})^2 - 10*{tt})"
    f = f"(sqrt((Y/{a})^2 + ((Z-{zc})/{b})^2) - 1)*min({a}, {b})"
    expr = f"max(max({f}, -X), X - {L!r})"
    for k, v in (("X", "(x/1mm)"), ("Y", "(y/1mm)"), ("Z", "(z/1mm)")):
        expr = expr.replace(k, v)
    return expr


def body(g):
    field = g.leaf("real_field", {"expression": expression(), "units": {}})
    scaled = g.block("multiply<real_field,real_field>", "real_field", [field, g.real(1e-3, rk.LENGTH)], top=False)
    bounds = g.block("create_bounding_box<point,point>", "bounding_box", [g.point(-5, -45, -40), g.point(L + 5, 45, 30)],
                     top=False)
    return g.block("set_bounding_box<real_field,bounding_box>[1.1.0]", "implicit", [scaled, bounds], "Sides Body")


def build():
    os.makedirs(WORK, exist_ok=True)
    fus, model = mb.fuselage_block(), mb.model_block()
    g = Graph()
    comp = rk.call(g, fus, {"Bodies": g.block("core.list<implicit>", "list<implicit>", [body(g)], top=False),
                            "Names": g.block("core.list<text>", "list<text>", [g.text("Sides")], top=False),
                            "Slice Spacing": g.real(0.0025)}, "Sides Component")
    out = rk.call(g, model, {"Components": g.block("core.list<text>", "list<text>", [comp], top=False),
                             "Working Directory": g.file_path(WORK), "Case Name": g.text("sides")}, "Model")
    res, log, dt = rk.run(g, [fus, model], out, "text", os.path.join(WORK, "sides.ntop"))
    print(f"build: {dt:.0f} s, output {'yes' if res else 'no'}")
    for line in log.splitlines():
        if "[E]" in line:
            print("  ", line[:300])


SAMPLE = """void main()
{{
    ReadVSPFile( "{model}" );
    Update();
    array< string > @ids = FindGeomsWithName( "Sides" );
    for ( int a = 0; a <= 600; a++ ) for ( int b = 0; b < 72; b++ ) {{
        vec3d q = CompPnt01( ids[0], 0, a / 600.0, b / 72.0 );
        Print( "PT " + q.x() + " " + q.y() + " " + q.z() ); }}
}}
"""


def error(p):
    """Distance (mm) from a point to the exact surface, in its station's plane (normal to x): along the ray from the
    ellipse centre, which for these aspect ratios is within a few percent of the true distance."""
    x, y, z = p
    if x < 0.5 or x > L - 0.5:
        return None  # the tip and the base cap
    a, b, zc = shape(x)
    th = math.atan2((z - zc) / b, y / a)
    q = (a * math.cos(th), zc + b * math.sin(th))
    r_p = math.hypot(y, z - zc)
    r_q = math.hypot(q[0], q[1] - zc)
    return abs(r_p - r_q)


def fit():
    sys.path.insert(0, rk.ROOT)
    settings = os.path.join(WORK, "sides.ntvsp.settings.txt")
    env = dict(os.environ)
    py = os.path.join(rk.ROOT, "build", "py.sh")
    vsp = os.path.join(rk.ROOT, "OpenVSP-3.56.0-win64", "vspscript.exe")
    for variant in ("allsym", "independent"):  # the old default, the new one
        env["NTVSP_SKIN_SIDES"] = variant
        p = subprocess.run(["sh", py, "-m", "vsp_ntop.run_model", settings], capture_output=True, text=True, env=env,
                           cwd=rk.ROOT)
        log = open(os.path.join(WORK, "sides.ntvsp.log")).read()
        check = [l for l in log.splitlines() if "OpenVSP surface vs nTop" in l]
        model = os.path.join(WORK, "sides.vsp3").replace("\\", "/")
        script = os.path.join(WORK, "sample.vspscript")
        open(script, "w").write(SAMPLE.format(model=model))
        out = subprocess.run([vsp, "-script", script], capture_output=True, text=True).stdout
        pts = [tuple(1000.0 * float(v) for v in l.split()[1:4]) for l in out.splitlines() if l.strip().startswith("PT")]
        err = [(e, p) for p in pts for e in [error(p)] if e is not None]
        regions = {"nose x<200": (0, 200), "mid 200-450": (200, 450), "boattail 450-600": (450, 600)}
        print(f"{variant}: {check[-1].split(']')[-1].strip() if check else p.stderr[-300:]}")
        for name, (lo, hi) in regions.items():
            e = [v for v, q in err if lo <= q[0] < hi]
            worst = max(err, key=lambda t: t[0] if lo <= t[1][0] < hi else -1)
            print(f"   exact surface, {name}: mean {sum(e) / len(e):.3f} mm, max {max(e):.3f} mm at "
                  f"({worst[1][0]:.1f}, {worst[1][1]:.1f}, {worst[1][2]:.1f})")


if __name__ == "__main__":
    steps = sys.argv[1:] or ["build", "fit"]
    if "build" in steps:
        build()
    if "fit" in steps:
        fit()
