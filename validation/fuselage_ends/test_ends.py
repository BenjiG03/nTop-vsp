"""OpenVSP fuselage fits of bodies with blunt ends, against their exact radius.

    sh build/py.sh validation/fuselage_ends/test_ends.py

Bodies of revolution about x (mm): a cylinder (flat at both ends), the hobby rocket's body (pointed ogive nose, flat
base), a hemisphere-nosed cylinder, and a flat-nosed cylinder with a conical boattail to a smaller flat base. Each is
an OpenVSP Fuselage (Set None) made by the blocks in one notebook; OpenVSP's surface is then sampled in rings
(CompPnt01, as Body Lift does) and its radius compared with the exact one along x, and the ends' positions with the
body's. Runs: runs/fuselage_ends.
"""
from __future__ import annotations

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))
import rocket as rk  # noqa: E402
import model_blocks as mb  # noqa: E402
from ntopfile import Graph  # noqa: E402

sys.path.insert(0, rk.ROOT)
from vsp_ntop import run_body_lift as bl  # noqa: E402
from vsp_ntop.common import Log, run_vspscript  # noqa: E402

WORK = os.path.join(rk.ROOT, "runs", "fuselage_ends").replace("\\", "/")
R, L = 27.0, 1200.0


def revolve(g, outline, name, y=0.0):
    prof = rk._profile(g, [(x, y, z) for x, _, z in outline], name + " Profile")
    axis = g.block("axis<point,vector>", "axis", [g.point(0, y, 0), g.vector(1, 0, 0)], name + " Axis", top=False)
    return g.block("revolve<new_profile,axis,real>[5.20.0]", "implicit", [prof, axis, g.real(rk.FULL_TURN, rk.ANGLE)],
                   name)


def bodies():
    """name: (outline in the x-z half plane from x = 0 to the end, exact radius r(x))."""
    out = {}
    out["cylinder"] = ([(0, 0, 0), (0, 0, R), (L, 0, R), (L, 0, 0)], lambda x: R)
    og = rk.ogive_points()
    rho = (R**2 + 270.0**2) / (2 * R)
    out["ogive_base"] = (og + [(L, 0, R), (L, 0, 0)],
                         lambda x: math.sqrt(max(rho**2 - (270 - x) ** 2, 0)) + R - rho if x < 270 else R)
    hemi = [(R - R * math.cos(t), 0, R * math.sin(t)) for t in (0.5 * math.pi * i / 60 for i in range(61))]
    out["hemisphere"] = (hemi + [(L, 0, R), (L, 0, 0)],
                         lambda x: math.sqrt(max(R * R - (R - x) ** 2, 0)) if x < R else R)
    xb, rb = 1000.0, 15.0
    out["boattail"] = ([(0, 0, 0), (0, 0, R), (xb, 0, R), (L, 0, rb), (L, 0, 0)],
                       lambda x: R if x < xb else R + (rb - R) * (x - xb) / (L - xb))
    return out


def main():
    os.makedirs(WORK, exist_ok=True)
    fus, model = mb.fuselage_block(), mb.model_block()
    g = Graph()
    bs = bodies()
    comps = []
    for i, (name, (outline, _)) in enumerate(bs.items()):
        body = revolve(g, outline, name, 200.0 * i)  # side by side in y, so the components do not overlap
        comps.append(rk.call(g, fus, {"Bodies": g.block("core.list<implicit>", "list<implicit>", [body], top=False),
                                      "Names": g.block("core.list<text>", "list<text>", [g.text(name)], top=False),
                                      "Set": g.choice(list(mb.SETS), 2)}, name + " Component"))
    lst = g.block("core.list<text>", "list<text>", comps, "Components")
    out = rk.call(g, model, {"Components": lst, "Working Directory": g.file_path(WORK), "Case Name": g.text("ends")},
                  "Model")
    res, log, dt = rk.run(g, [fus, model], out, "text", os.path.join(WORK, "ends.ntop"))
    print(f"model: {dt:.0f} s, output {'yes' if res else 'no'}")
    for line in log.splitlines():
        if "[E]" in line or "contours:" in line or "sections from" in line:
            print("  ", line[line.find("]", 12) + 1:][:200])
    report(bs)


def report(bs):
    """Radius of OpenVSP's surface about each body's own axis (y = 200 i mm, z = 0) along x, against the exact one;
    rings in the plane of an end cap are left out (they cross the cap)."""
    lg = Log(os.path.join(WORK, "rings.log"))
    for i, (name, (_, exact)) in enumerate(bs.items()):
        text = run_vspscript(rk.VSP, bl.ring_script(f"{WORK}/ends.vsp3", [name], 400, 64), WORK, f"rings_{name}", lg, 600)
        rings = []
        for line in text.splitlines():
            p = line.split()
            if len(p) < 9 or p[0] != "NTVSP_RING":
                continue
            v = [float(t) for t in p[3:]]
            pts = [v[j:j + 3] for j in range(0, len(v) - 2, 3)]
            rings.append((sum(q[0] for q in pts) / len(pts) * 1000,
                          max(math.hypot(q[1] * 1000 - 200.0 * i, q[2] * 1000) for q in pts)))
        xs = [x for x, _ in rings]
        keep = [(x, r) for k, (x, r) in enumerate(rings)
                if not any(abs(x - rings[j][0]) < 1e-6 for j in (k - 1, k + 1) if 0 <= j < len(rings))]
        errs = [(r - exact(min(max(x, 0.0), L)), x) for x, r in keep]
        worst = max(errs, key=lambda e: abs(e[0]))
        near = lambda lo, hi: max([abs(e) for e, x in errs if lo <= x <= hi] or [0.0])  # noqa: E731
        print(f"{name:11s}: OpenVSP x {min(xs):8.2f} to {max(xs):8.2f} mm (body 0 to {L:g}); radius error largest "
              f"{worst[0]:+.2f} mm at x = {worst[1]:.1f}; x < 60: {near(-10, 60):.2f}, 60 to {L - 60:g}: "
              f"{near(60, L - 60):.2f}, x > {L - 60:g}: {near(L - 60, L + 10):.2f}")


if __name__ == "__main__":
    if sys.argv[1:] == ["report"]:
        report(bodies())
    else:
        main()
