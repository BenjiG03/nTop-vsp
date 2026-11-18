"""Check Modified Newtonian on Mesh (VSP) against Modified Newtonian (VSP) (Surface Integrals) on the hobby rocket, and
the User Defined Cp Max Field in both.

    sh build/py.sh validation/rasaero/test_newtonian_mesh.py [--tolerance-mm 0.5]
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import rocket as rk  # noqa: E402
import aero_blocks as a  # noqa: E402
from ntopfile import Graph  # noqa: E402

DEG = math.radians(1.0)


def lst(g, typ, items):
    return g.block(f"core.list<{typ}>", f"list<{typ}>", items, top=False)


def build(tol_mm, work):
    flight_b, si_b, mesh_b = a.flight_block(), a.newtonian_block(), a.newtonian_mesh_block()
    g = Graph()
    body = rk.body(g)
    whole = rk.boolean(g, "union", [body] + rk.fins(g), "Hobby Rocket")
    r = rk.R
    flight = rk.call(g, flight_b, {
        "Mach": lst(g, "real", [g.real(m) for m in (3.0,)]),
        "Altitude": lst(g, "real", [g.real(0.0, rk.LENGTH)]),
        "Angle of Attack": lst(g, "real", [g.real(v * DEG, rk.ANGLE) for v in (0.0, 8.0)]),
        "Sideslip": lst(g, "real", [g.real(v * DEG, rk.ANGLE) for v in (0.0, 3.0)]),
        "Reference Area": g.real(r.reference_area_mm2 * 1e-6, rk.AREA),
        "Reference Span": g.real(r.diameter * 1e-3, rk.LENGTH),
        "Reference Chord": g.real(r.diameter * 1e-3, rk.LENGTH),
        "Working Directory": g.file_path(work), "Case Name": g.text("flight")}, "Flight")
    user = lambda g: g.choice(list(a.newtonian.CP_MAX_MODES), 2)  # noqa: E731
    field = lambda g: g.leaf("real_field", {"expression": "1.5", "units": {}})  # noqa: E731
    classical = lambda g: g.choice(list(a.newtonian.CP_MAX_MODES), 1)  # noqa: E731
    body_mesh, whole_mesh = rk.mesh(g, body, tol_mm, "Body Mesh"), rk.mesh(g, whole, tol_mm, "Rocket Mesh")
    outs = []
    for name, blk, extra in (
            ("si_modified", si_b, {"Body": body}),
            ("si_user", si_b, {"Body": body, "Stagnation Pressure": user, "Cp Max Field": field}),
            ("mesh_body_modified", mesh_b, {"Mesh": body_mesh}),
            ("mesh_body_user", mesh_b, {"Mesh": body_mesh, "Stagnation Pressure": user, "Cp Max Field": field}),
            ("mesh_modified", mesh_b, {"Mesh": whole_mesh}),
            ("mesh_classical", mesh_b, {"Mesh": whole_mesh, "Stagnation Pressure": classical}),
            ("mesh_user", mesh_b, {"Mesh": whole_mesh, "Stagnation Pressure": user, "Cp Max Field": field})):
        given = {"Flight Condition": flight, "Case Name": g.text(name)}
        if blk is mesh_b:
            given["Working Directory"] = g.file_path(work)
        else:
            given["Body"] = whole
        given.update(extra)
        outs.append(rk.call(g, blk, given, name))
    group = g.block("core.group<any>", "group", outs, "Results")
    return g, [flight_b, si_b, mesh_b], group


def read(path):
    with open(path, newline="") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tolerance-mm", type=float, default=0.5)
    args = p.parse_args()
    work = os.path.join(rk.RUNS, "newtonian_mesh_test").replace("\\", "/")
    os.makedirs(work, exist_ok=True)
    g, blocks, out = build(args.tolerance_mm, work)
    res, log, dt = rk.run(g, blocks, out, "group", os.path.join(work, "test.ntop"))
    print(f"ntopcl {dt:.0f} s, output {'yes' if res else 'no'}")
    errs = [l for l in log.splitlines() if "[E]" in l]
    print("\n".join(errs[:20]))
    for line in log.splitlines():
        if "complete" in line.lower() and ("Modified" in line or "si_" in line or "mesh_" in line):
            print(line[:200])
    files = {n: (os.path.join(work, f"flight_{n}.ntvsp.coefficients.csv") if n.startswith("si") else
                 os.path.join(work, f"{n}.ntvsp.coefficients.csv")) for n in
             ("si_modified", "si_user", "mesh_body_modified", "mesh_body_user", "mesh_modified", "mesh_classical",
              "mesh_user")}
    tables = {n: read(f) for n, f in files.items() if os.path.exists(f)}
    print("tables:", list(tables))
    keys = ("CA", "CY", "CN", "Cl", "Cm", "Cn")
    for a_, b_ in (("si_modified", "mesh_body_modified"), ("si_user", "mesh_body_user"), ("mesh_modified", "mesh_user")):
        if a_ in tables and b_ in tables:
            print(f"--- {a_} vs {b_}")
            for ra, rb in zip(tables[a_], tables[b_]):
                print(f"M{ra['Mach']:.1f} a{ra['AoA']:+.0f} b{ra['Beta']:+.0f} " + "  ".join(
                    f"{k} {ra[k]:+.4f}/{rb[k]:+.4f}" for k in keys))


if __name__ == "__main__":
    main()
