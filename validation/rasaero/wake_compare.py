"""Relaxed against fixed wake for the best VSPAERO model (thin fins through the body to the axis).

    sh build/py.sh validation/rasaero/wake_compare.py [model]

For each relaxed-wake run (3, 5 or 8 wake iterations, at sweep.RELAXED_MACHS): the largest change of CN, CY and the
moments from the fixed wake (over the largest |value| of the fixed-wake run at that Mach number), the run time and
how many cases passed the convergence check.
"""
from __future__ import annotations

import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sweep  # noqa: E402

COLS = ("CFiz", "CFiy", "CMiy", "CMiz")


def polar(case):
    path = f"{sweep.WORK}/{case}.ntvsp.polar.csv"
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return {(float(r["AoA"]), float(r["Beta"])): r for r in csv.DictReader(f)}


def seconds(case):
    text = open(f"{sweep.WORK}/{case}.ntvsp.log", encoding="utf-8").read()
    m = re.findall(r"\[\s*([0-9.]+)s\] done", text)
    return float(m[-1]) if m else float("nan")


def main():
    model = sys.argv[1] if len(sys.argv) > 1 else "thin_axis"
    print(f"{model}: relaxed vs fixed wake (largest change over the fixed run's largest |value|)")
    for wake in ("relaxed3", "relaxed5", "relaxed8"):
        for m in sweep.RELAXED_MACHS:
            fixed, relaxed = polar(sweep.vspaero_case(model, "fixed", m)), polar(sweep.vspaero_case(model, wake, m))
            if not fixed or not relaxed:
                continue
            out = []
            for c in COLS:
                scale = max(abs(float(r[c])) for r in fixed.values())
                worst = max(abs(float(relaxed[k][c]) - float(fixed[k][c])) for k in fixed)
                out.append(f"{c} {worst / scale:.2e}")
            conv = sum(int(float(r.get("Converged", 1))) for r in relaxed.values())
            print(f"  {wake} Mach {m:g}: " + ", ".join(out) + f"; {seconds(sweep.vspaero_case(model, wake, m)):.0f} s "
                  f"vs {seconds(sweep.vspaero_case(model, 'fixed', m)):.0f} s fixed; {conv}/{len(relaxed)} converged")


if __name__ == "__main__":
    main()
