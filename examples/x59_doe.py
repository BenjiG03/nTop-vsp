"""Design of experiments over the X-59's wing: one aero database with geometry as independent variables.

    python examples/x59_doe.py [--quick]

``X-59 Aero Database.ntop`` has two notebook inputs, Wing Span and Wing Inner Sweep (the starter's own parametric
variables). This script runs the notebook headlessly with nTop Automate (ntopcl) once per design point, giving the
inputs in a JSON file (``-j``). Every run rebuilds the nTop geometry, converts it to an OpenVSP model, runs VSPAERO,
Wave Drag and Parasite Drag over the same flight envelope, and its Aero Database block appends its rows, with the
two variables in the first columns, to ``runs/x59_doe.csv``. That file is the database: CL, CD, CY, Cl, Cm, Cn
against Wing Span, Wing Inner Sweep, Mach, alpha and beta.

Then it prints, per design, a few derived numbers at beta 0: lift-curve slope and pitching-moment slope at Mach 0.6
and 1.5, and zero-lift drag at Mach 1.5.

--quick runs two designs instead of five; --keep adds to an existing table, skipping designs it already has.
Needs a signed-in nTop (ntopcl), and the NTOP_PYTHON_PATH and VSP_PATH environment variables that the blocks read
(see the README's Setup); ntopcl passes this script's environment on to the notebook.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NOTEBOOK = os.path.join(HERE, "X-59 Aero Database.ntop")
TABLE = os.path.join(ROOT, "runs", "x59_doe.csv")
WORK = os.path.join(ROOT, "runs", "doe")
NTOPCL = r"C:\Program Files\nTopology\nTopology\ntopcl.exe"

# (wing span mm, inner quarter-chord sweep deg): the baseline, then one variable at a time
DESIGNS = [(9000.0, 68.3), (8000.0, 68.3), (10000.0, 68.3), (9000.0, 65.0), (9000.0, 71.0)]


def run(span_mm, sweep_deg, tag, timeout=7200):
    os.makedirs(WORK, exist_ok=True)
    inputs = {"description": "X-59 DOE point", "title": "X-59 Aero Database", "inputs": [
        {"name": "Wing Span", "type": "real", "units": "mm", "value": span_mm},
        {"name": "Wing Inner Sweep", "type": "real", "units": "deg", "value": sweep_deg}]}
    jpath = os.path.join(WORK, f"{tag}.input.json")
    json.dump(inputs, open(jpath, "w"), indent=1)
    out = os.path.join(WORK, f"{tag}.output.json")
    if os.path.exists(out):
        os.remove(out)
    log_path = os.path.join(WORK, f"{tag}.log")
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        p = subprocess.Popen([NTOPCL, "--trustnotebook", "-v", "2", "-j", jpath, "-o", out, NOTEBOOK], stdout=logf,
                             stderr=subprocess.STDOUT)
        # ntopcl can write its output and then not exit after Run Command blocks: stop it once the output settles.
        last = -1
        while p.poll() is None and time.time() - t0 < timeout:
            time.sleep(3)
            if os.path.exists(out):
                size = os.path.getsize(out)
                if size == last and size > 0:
                    break
                last = size
        if p.poll() is None:
            p.kill()
            p.wait()
    errors = [l for l in open(log_path, encoding="utf-8", errors="replace") if "[E]" in l]
    return time.time() - t0, os.path.exists(out), errors


def slope(rows, x, y):
    n = len(rows)
    mx = sum(float(r[x]) for r in rows) / n
    my = sum(float(r[y]) for r in rows) / n
    sxx = sum((float(r[x]) - mx) ** 2 for r in rows)
    return sum((float(r[x]) - mx) * (float(r[y]) - my) for r in rows) / sxx


def summary():
    rows = list(csv.DictReader(open(TABLE)))
    designs = sorted({(r["Wing_Span_m"], r["Wing_Inner_Sweep_deg"]) for r in rows}, key=lambda t: (float(t[1]), float(t[0])))
    print(f"\n{TABLE}: {len(rows)} rows, {len(designs)} designs, columns: {', '.join(rows[0].keys())}\n")
    print("span m  sweep deg | CLa/deg M0.6  Cma/deg M0.6 | CLa/deg M1.5  Cma/deg M1.5 | CD(alpha 0) M1.5")
    for span, sweep in designs:
        sel = [r for r in rows if (r["Wing_Span_m"], r["Wing_Inner_Sweep_deg"]) == (span, sweep) and abs(float(r["Beta"])) < 1e-9]
        out = []
        for m in (0.6, 1.5):
            at = [r for r in sel if abs(float(r["Mach"]) - m) < 1e-6]
            out += [slope(at, "AoA", "CL"), slope(at, "AoA", "Cm")]
        cd0 = [float(r["CD"]) for r in sel if abs(float(r["Mach"]) - 1.5) < 1e-6 and abs(float(r["AoA"])) < 1e-9]
        print(f"{float(span):6.2f}  {float(sweep):9.1f} | {out[0]:12.5f}  {out[1]:12.5f} | {out[2]:12.5f}  {out[3]:12.5f} | "
              f"{cd0[0] if cd0 else float('nan'):.5f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true", help="two designs only")
    parser.add_argument("--summary", action="store_true", help="only print the summary of an existing table")
    parser.add_argument("--keep", action="store_true", help="keep the existing table and skip the designs it has")
    args = parser.parse_args()
    missing = [v for v in ("NTOP_PYTHON_PATH", "VSP_PATH") if not os.environ.get(v)]
    if missing and not args.summary:
        raise SystemExit(f"set {' and '.join(missing)} first (see the README's Setup)")
    if not args.summary:
        done = set()
        if args.keep and os.path.exists(TABLE):
            done = {(float(r["Wing_Span_m"]), float(r["Wing_Inner_Sweep_deg"])) for r in csv.DictReader(open(TABLE))}
        elif os.path.exists(TABLE):
            os.remove(TABLE)
        for k, (span, sweep) in enumerate(DESIGNS[:2] if args.quick else DESIGNS, start=1):
            if any(abs(span / 1000 - a) < 1e-6 and abs(sweep - b) < 1e-6 for a, b in done):
                print(f"design {k}: span {span:g} mm, inner sweep {sweep:g} deg: already in the table")
                continue
            dt, ok, errors = run(span, sweep, f"design{k}")
            print(f"design {k}: span {span:g} mm, inner sweep {sweep:g} deg: {dt:.0f} s, {'ok' if ok else 'no output'}",
                  flush=True)
            for e in errors[:5]:
                print("   ", e.strip())
    summary()


if __name__ == "__main__":
    main()
