"""Modified Newtonian on Mesh (VSP) over the full grid: the whole rocket, the body alone and the exposed fins alone.

    sh build/py.sh validation/rasaero/newtonian_sweep.py [--tolerance-mm 2]

Each runs with Stagnation Pressure User Defined and a Cp Max Field of 1, so its coefficients are per unit Cp_max:
Newtonian forces are linear in a constant Cp_max, so any Cp_max choice (Modified, Classical or a fitted value per
component and Mach number) is a scale factor afterwards. The whole rocket also runs with Modified, as a check of
that scaling against the block's own Cp_max.
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import rocket as rk  # noqa: E402
import sweep  # noqa: E402
import aero_blocks as a  # noqa: E402
from ntopfile import Graph  # noqa: E402


def build(tol_mm, fin_tol_mm):
    blk = a.newtonian_mesh_block()
    g = Graph()
    body = rk.body(g)
    whole = rk.boolean(g, "union", [body] + rk.fins(g), "Hobby Rocket")
    exposed = rk.boolean(g, "union", rk.fins(g, inner_fraction=1.0), "Exposed Fins")
    flight = g.text(f"{sweep.WORK}/flight_all")
    user = g.choice(list(a.newtonian.CP_MAX_MODES), 2)
    unit = g.leaf("real_field", {"expression": "1", "units": {}})
    outs = []
    for name, body_, tol, mode in (("newt_whole_unit", whole, tol_mm, user), ("newt_body_unit", body, tol_mm, user),
                                   ("newt_fins_unit", exposed, fin_tol_mm, user),
                                   ("newt_whole_modified", whole, tol_mm, None)):
        given = {"Body": body_, "Flight Condition": flight, "Mesh Tolerance": g.real(tol * 1e-3, rk.LENGTH),
                 "Working Directory": g.file_path(sweep.WORK), "Case Name": g.text(name), "Cp Max Field": unit}
        if mode is not None:
            given["Stagnation Pressure"] = mode
        outs.append(rk.call(g, blk, given, name))
    return g, [blk], g.block("core.group<any>", "group", outs, "Results")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tolerance-mm", type=float, default=2.0)
    p.add_argument("--fin-tolerance-mm", type=float, default=1.0)
    args = p.parse_args()
    g, blocks, out = build(args.tolerance_mm, args.fin_tolerance_mm)
    res, log, dt = rk.run(g, blocks, out, "group", os.path.join(sweep.WORK, "newtonian.ntop"))
    print(f"{dt:.0f} s, output {'yes' if res else 'no'}")
    for line in log.splitlines():
        if "[E]" in line or "Newtonian on mesh:" in line:
            print(line[:250])


if __name__ == "__main__":
    main()
