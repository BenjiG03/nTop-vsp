"""Tangent Wedge on Mesh (VSP) over the full grid: the exposed fins and the body alone.

    sh build/py.sh validation/rasaero/tangent_wedge_sweep.py

Prandtl-Meyer leeward faces, freestream pressure on separated faces (no base drag), Modified Cp_max above shock
detachment. Results: runs/rasaero/sweep/tw_fins, tw_body (``.ntvsp.coefficients.csv``).
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import rocket as rk  # noqa: E402
import sweep  # noqa: E402
import aero_blocks as a  # noqa: E402
from ntopfile import Graph  # noqa: E402


def main():
    blk = a.tangent_wedge_mesh_block()
    g = Graph()
    flight = g.text(f"{sweep.WORK}/flight_all")
    outs = []
    for name, body, tol in (("tw_fins", rk.boolean(g, "union", rk.fins(g, inner_fraction=1.0), "Exposed Fins"), 1.0),
                            ("tw_body", rk.body(g), 2.0)):
        outs.append(rk.call(g, blk, {"Body": body, "Flight Condition": flight, "Mesh Tolerance": g.real(tol * 1e-3, rk.LENGTH),
                                     "Working Directory": g.file_path(sweep.WORK), "Case Name": g.text(name)}, name))
    res, log, dt = rk.run(g, [blk], g.block("core.group<any>", "group", outs, "Results"), "group",
                          os.path.join(sweep.WORK, "tangent_wedge.ntop"))
    print(f"{dt:.0f} s, output {'yes' if res else 'no'}")
    for line in log.splitlines():
        if "[E]" in line or "on mesh:" in line:
            print(line[:250])


if __name__ == "__main__":
    main()
