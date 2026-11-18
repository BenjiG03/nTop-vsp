"""The nTop_vsp combinations and RASAero II against inviscid SU2 on the trapezoid-section rocket.

    sh build/py.sh validation/noncircular/compare.py

Sources (runs/noncircular, from build.py) enter as the Aero Database does it (run_database.Contribution); the
VSPAERO run has the fins with Exposed Strip + Carry-over (TR 1307), and the carry-over is also taken off again
(its own term, recomputed from the run's spanwise loads) to show Exposed Strip alone. SU2: results/su2.json (cfd.py).
Every coefficient on Sref = pi 54^2 / 4 mm^2, moments over Sref 54 mm about x = 0.7 m. Writes results/compare.json.
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "rasaero"))
import rocket as rk  # noqa: E402
import build as nb  # noqa: E402

sys.path.insert(0, rk.ROOT)
from vsp_ntop import run_database as db  # noqa: E402
from vsp_ntop import run_vspaero as rv  # noqa: E402
from vsp_ntop.flight import Flight  # noqa: E402

W = nb.WORK
X_REF, D = 0.7, 0.054
KEYS = ("CA", "CY", "CN", "Cl", "Cm", "Cn")
SU2 = os.path.join(HERE, "results", "su2.json")


class _S:
    def __init__(self, **kw):
        self.kw = kw

    def real(self, k, default=None):
        return self.kw.get(k, default)


def carryover_term(flight, case="vspaero", model="model"):
    """The carry-over alone, as a polar file (``<case>_carry_only``): carryover_accounting on zero rows."""
    rows = list(csv.DictReader(open(f"{W}/{case}.ntvsp.polar.csv")))
    header = list(rows[0].keys())
    keep = ("Case", "Mach", "AoA", "Beta", "Re_1e6", "Converged")
    data = [[float(r[h]) if h in keep else 0.0 for h in header] for r in rows]
    sref, bref, cref, (xr, yr, zr) = flight.ref
    s = _S(sref=sref, bref=bref, cref=cref, xcg=xr, ycg=yr, zcg=zr)
    fins = rv.carryover_fins(f"{W}/{model}.vsp3")
    rv.carryover_accounting(header, data, f"{W}/{case}.lod", fins, s, lambda *_: None)
    import shutil
    shutil.copy(f"{W}/{case}.ntvsp.settings.txt", f"{W}/{case}_carry_only.ntvsp.settings.txt")
    with open(f"{W}/{case}_carry_only.ntvsp.polar.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)
    return fins


def load(flight, kind, base, use="all"):
    con = db.Contribution({"kind": kind, "results": f"{W}/{base}", "use": use, "name": base}, lambda *_: None)
    out = []
    for c in flight.conditions:
        got = con.on(flight.ref, c)
        out.append(list(got[0]) + list(got[1]))
    return np.array(out)


def x_cp(row):
    return X_REF - row[4] * D / row[2] if abs(row[2]) > 1e-9 else float("nan")


def main():
    flight = Flight(f"{W}/flight")
    conds = flight.conditions
    fins = carryover_term(flight)
    vsp = load(flight, "vspaero", "vspaero")
    carry = load(flight, "vspaero", "vspaero_carry_only")
    import shutil
    shutil.copy(f"{W}/vspaero.ntvsp.polar_whole.csv", f"{W}/vspaero_whole.ntvsp.polar.csv")
    shutil.copy(f"{W}/vspaero.ntvsp.settings.txt", f"{W}/vspaero_whole.ntvsp.settings.txt")
    carryover_term(flight, "vspaero_plate", "plate_model")
    plate = load(flight, "vspaero", "vspaero_plate")
    plate_carry = load(flight, "vspaero", "vspaero_plate_carry_only")
    parts = {
        "VSPAERO fins, Exposed Strip + Carry-over": vsp,
        "Carry-over term": carry,
        "VSPAERO fins, Exposed Strip": vsp - carry,
        "Body Lift (slender + crossflow)": load(flight, "coefficients", "body_lift"),
        "Body Lift (section shape + crossflow)": load(flight, "coefficients", "body_lift_shape"),
        "Newtonian body (Modified)": load(flight, "coefficients", "newtonian_body"),
        "RASAero II (round body fitted)": load(flight, "coefficients", "rasaero"),
        "VSPAERO fins, Whole Wing": load(flight, "vspaero", "vspaero_whole"),
        "VSPAERO body plate + fins": plate - plate_carry,
        "VSPAERO body plate + fins + Carry-over": plate,
    }
    combos = {
        "Exposed + Carry-over + Body Lift (shape)": ("VSPAERO fins, Exposed Strip + Carry-over",
                                                     "Body Lift (section shape + crossflow)"),
        "Exposed + Body Lift (shape)": ("VSPAERO fins, Exposed Strip", "Body Lift (section shape + crossflow)"),
        "Exposed + Carry-over + Body Lift": ("VSPAERO fins, Exposed Strip + Carry-over", "Body Lift (slender + crossflow)"),
        "Exposed + Carry-over + Newtonian body": ("VSPAERO fins, Exposed Strip + Carry-over", "Newtonian body (Modified)"),
        "Exposed + Body Lift": ("VSPAERO fins, Exposed Strip", "Body Lift (slender + crossflow)"),
        "Exposed + Newtonian body": ("VSPAERO fins, Exposed Strip", "Newtonian body (Modified)"),
        "RASAero II": ("RASAero II (round body fitted)",),
        "Plate + fins (VLM)": ("VSPAERO body plate + fins",),
        "Plate + fins + Carry-over": ("VSPAERO body plate + fins + Carry-over",),
        "Whole Wing + Body Lift": ("VSPAERO fins, Whole Wing", "Body Lift (slender + crossflow)"),
        "Whole Wing + Newtonian body": ("VSPAERO fins, Whole Wing", "Newtonian body (Modified)"),
    }
    totals = {k: sum(parts[p] for p in v) for k, v in combos.items()}
    su2 = json.load(open(SU2))["runs"] if os.path.exists(SU2) else {}
    out = {"fins": fins, "conditions": [], "su2": su2}
    print("carry-over geometry: " + "; ".join(f"{f['Name']} r {f['Root'] * 1e3:.1f} mm, body {f['Body_Width'] * 1e3:.1f} "
                                              f"mm across" for f in fins))
    head = f"{'cond':>9} | {'SU2 CN':>7} {'x_cp':>6} {'Cm':>7} |" + "".join(f" {k[:24]:>24} |" for k in combos)
    print(head + "\n" + " " * 37 + "|" + " (dCN%, dx_cp mm, dCm%)  |" * len(combos))
    for i, c in enumerate(conds):
        key = f"M{c['Mach']:g}_a{c['AoA']:g}"
        row = {"Mach": c["Mach"], "AoA": c["AoA"],
               "parts": {k: dict(zip(KEYS, map(float, v[i]))) for k, v in parts.items()},
               "totals": {k: dict(zip(KEYS, map(float, v[i])), x_cp=x_cp(v[i])) for k, v in totals.items()}}
        out["conditions"].append(row)
        s = su2.get(key)
        line = f"{key:>9} | "
        if s:
            line += f"{s['CN']:7.3f} {s['x_cp_mm'] / 1e3:6.3f} {s['Cm']:7.3f} |"
            for k, v in totals.items():
                cn, cm = v[i, 2], v[i, 4]
                xc = x_cp(v[i])
                big = abs(s["CN"]) > 0.05
                dcn = f"{100 * (cn - s['CN']) / s['CN']:+4.0f}%" if big else f"{cn - s['CN']:+.3f}"
                dcm = f"{100 * (cm - s['Cm']) / s['Cm']:+4.0f}%" if abs(s["Cm"]) > 0.2 else f"{cm - s['Cm']:+.2f}"
                dx = f"{1e3 * xc - s['x_cp_mm']:+5.0f}" if big else "    -"
                line += f" {dcn:>6} {dx:>6} {dcm:>7}    |"
        else:
            line += "   (no SU2 run)        |" + "".join(f" CN {v[i, 2]:6.3f} Cm {v[i, 4]:6.3f}   |" for v in totals.values())
        print(line)
    print("\nby region, CN (x_cp m): SU2 | nTop_vsp")
    for i, c in enumerate(conds):
        s = su2.get(f"M{c['Mach']:g}_a{c['AoA']:g}")
        if not s or abs(c["AoA"]) < 1e-9:
            continue
        reg = s["regions"]
        xc = lambda cn, cm: X_REF - cm * D / cn if abs(cn) > 1e-9 else float("nan")  # noqa: E731
        body = {k: sum(reg[r][k] for r in ("nose", "tube")) for k in ("CN", "Cm")}
        p = {k: v[i] for k, v in parts.items()}
        print(f"  M{c['Mach']:g} a{c['AoA']:g}: SU2 nose {reg['nose']['CN']:.3f}, body forward of the fins "
              f"{body['CN']:.3f} ({xc(body['CN'], body['Cm']):.3f}), tube at fins {reg['tube at fins']['CN']:.3f}, "
              f"fins {reg['fins']['CN']:.3f} ({xc(reg['fins']['CN'], reg['fins']['Cm']):.3f}) | "
              f"Body Lift {p['Body Lift (slender + crossflow)'][2]:.3f} "
              f"({x_cp(p['Body Lift (slender + crossflow)']):.3f}), shape {p['Body Lift (section shape + crossflow)'][2]:.3f} "
              f"({x_cp(p['Body Lift (section shape + crossflow)']):.3f}), Newtonian body "
              f"{p['Newtonian body (Modified)'][2]:.3f} ({x_cp(p['Newtonian body (Modified)']):.3f}), carry-over "
              f"{p['Carry-over term'][2]:.3f}, exposed fins {p['VSPAERO fins, Exposed Strip'][2]:.3f} "
              f"({x_cp(p['VSPAERO fins, Exposed Strip']):.3f})")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "compare.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
