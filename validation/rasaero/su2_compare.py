"""Normal force and centre of pressure against hypanel's inviscid SU2 runs of the same rocket.

    sh build/py.sh validation/rasaero/su2_compare.py

SU2 8.2 (inviscid, tetrahedral meshes built through nTop; hypanel/validation/rocket/su2, results in
hypanel/validation/data/su2_inviscid.json) gives the normal force along the rocket, so it splits into nose (x < 270 mm),
shoulder (270-400), mid-tube (400-1040), tube between the fins (1040-1200) and fins. Its runs inside this study's range:
the finned rocket at Mach 2 alpha 2, Mach 3, 4, 5 alpha 4, Mach 5 alpha 8, and the body alone at Mach 2 alpha 2, Mach 3
alpha 4, Mach 5 alpha 8 (beta 0 throughout).

Against it: RASAero II (total, and its potential and viscous normal force separately), and the nTop_vsp pieces at the
same conditions: VSPAERO thin fins through the body to the axis (with its lift split, from VSPAERO's strip loads,
into the strips inside the body radius and the exposed ones), thin fins from the body surface, Modified Newtonian
(Modified Cp_max) on the body and on the fins.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as cmp  # noqa: E402
import rocket as rk  # noqa: E402
import sweep  # noqa: E402
from rasaero_ref import RASAero  # noqa: E402

SU2 = os.path.join(rk.PLAYGROUND, "hypanel", "validation", "data", "su2_inviscid.json")
RADIUS = 0.027
CONDS = ((2.0, 2.0), (3.0, 4.0), (4.0, 4.0), (5.0, 4.0), (5.0, 8.0))


def su2_regions(run):
    edges = np.array(run["edges"])
    w = np.diff(edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    rows = {k: np.array(v) * w for k, v in run["rows"].items()}
    body = rows["nose"] + rows["tube"]
    out = {"nose": body[mid < 270].sum(), "shoulder": body[(mid >= 270) & (mid < 400)].sum(),
           "mid-tube": body[(mid >= 400) & (mid < 1040)].sum(), "tube at fins": body[mid >= 1040].sum(),
           "fins": rows["fin"].sum()}
    out["total"] = run["cn"]
    out["x_cp"] = run["cp_mm"] / 1000.0
    return out


def lod_inner_fraction(case, alpha, beta=0.0):
    """Share of the fins' lift on strips inside the body radius (horizontal fins, from VSPAERO's .lod)."""
    path = f"{sweep.WORK}/{case}.lod"
    a = b = None
    inner = total = 0.0
    found = False
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "AoA_":
                a = float(p[1])
            elif p[0] == "Beta_":
                b = float(p[1])
                found = abs(a - alpha) < 1e-6 and abs(b - beta) < 1e-6
            elif found and p[0] == "1" and len(p) > 18:
                try:
                    y, z, area, cl = float(p[4]), float(p[5]), float(p[9]), float(p[17])
                except ValueError:
                    continue
                if abs(z) > 1e-9:  # vertical fins carry no lift at beta 0
                    continue
                total += cl * area
                if math.hypot(y, z) < RADIUS:
                    inner += cl * area
    return inner / total if total else float("nan")


def main():
    su2 = json.load(open(SU2))["runs"]
    data = cmp.Data()
    ras = RASAero()
    W = sweep.WORK
    a = np.array([c["AoA"] for c in data.conds])
    b = np.array([c["Beta"] for c in data.conds])
    cpm = cmp.cp_law("modified", data.mach)[:, None]
    parts = {
        "VSPAERO thin, fins to axis": data.source("vspaero", cmp.vsp_bases("thin_axis", "fixed")),
        "VSPAERO thin, fins from body": data.source("vspaero", cmp.vsp_bases("thin_exposed", "fixed")),
        "Newtonian body (Modified)": cpm * data.source("coefficients", f"{W}/newt_body_unit"),
        "Newtonian fins (Modified)": cpm * data.source("coefficients", f"{W}/newt_fins_unit"),
    }
    parts["VSPAERO extended, Exposed Strip"] = data.source("vspaero", cmp.vsp_bases("thin_extend_exposed", "fixed"))
    parts["VSPAERO extended, Whole Wing"] = data.source("vspaero", cmp.vsp_bases("thin_extend_whole", "fixed"))
    parts["Body Lift (slender body)"] = data.source("coefficients", f"{W}/bodylift_pot")
    parts["Body Lift (slender + crossflow)"] = data.source("coefficients", f"{W}/bodylift")
    parts["Exposed Strip + Newtonian body"] = parts["VSPAERO extended, Exposed Strip"] + parts["Newtonian body (Modified)"]
    parts["Exposed Strip + Body Lift (slender + crossflow)"] = parts["VSPAERO extended, Exposed Strip"] +         parts["Body Lift (slender + crossflow)"]
    parts["VSPAERO extended, Exposed Strip + Carry-over"] = data.source("vspaero",
                                                                        cmp.vsp_bases("thin_extend_carry", "fixed"))
    parts["Carry-over + Body Lift (slender + crossflow)"] = parts["VSPAERO extended, Exposed Strip + Carry-over"] + \
        parts["Body Lift (slender + crossflow)"]
    parts["Carry-over + Newtonian body"] = parts["VSPAERO extended, Exposed Strip + Carry-over"] + \
        parts["Newtonian body (Modified)"]
    parts["fins to axis + Newtonian body"] = parts["VSPAERO thin, fins to axis"] + parts["Newtonian body (Modified)"]
    parts["fins from body + Newtonian body"] = parts["VSPAERO thin, fins from body"] + parts["Newtonian body (Modified)"]
    out = {}
    for m, al in CONDS:
        i = int(np.where((data.mach == m) & (a == al) & (b == 0))[0][0])
        run = su2[f"rocket_t3_h3_M{m:g}_a{al:g}"]
        reg = su2_regions(run)
        t = ras.tables[ras.alphas.index(int(al))][round(m, 2)]
        row = {"SU2": reg, "RASAero": {"total": t["cn"], "potential": t["cn"] - (t["cn"] - 0) + 0, "x_cp": t["x_cp"]}}
        # the Run Test's third row: alpha, CN potential, CN viscous
        row["RASAero"] = {"total": t["cn"], "x_cp": t["x_cp"]}
        print(f"\nMach {m:g}, alpha {al:g}")
        print(f"  SU2 inviscid: CN {reg['total']:.3f}, x_cp {reg['x_cp']:.3f} m; nose {reg['nose']:.3f}, shoulder "
              f"{reg['shoulder']:.3f}, mid-tube {reg['mid-tube']:.3f}, tube at fins {reg['tube at fins']:.3f}, fins "
              f"{reg['fins']:.3f}")
        print(f"  RASAero II:   CN {t['cn']:.3f}, x_cp {t['x_cp']:.3f} m")
        row["parts"] = {}
        for name, arr in parts.items():
            cn, cmy = arr[i, 2], arr[i, 4]
            x = cmp.X_REF - cmy * cmp.D / cn if abs(cn) > 1e-9 else float("nan")
            row["parts"][name] = {"CN": float(cn), "x_cp": float(x)}
            print(f"  {name:32s} CN {cn:.3f}, x_cp {x:.3f} m")
        fr = lod_inner_fraction(sweep.vspaero_case("thin_axis", "fixed", m), al)
        cn_axis = row["parts"]["VSPAERO thin, fins to axis"]["CN"]
        row["axis_inner_fraction"] = fr
        print(f"  fins to axis: {fr * 100:.0f}% of its lift on strips inside the body radius ({fr * cn_axis:.3f}); "
              f"exposed {cn_axis * (1 - fr):.3f}")
        out[f"M{m:g}_a{al:g}"] = row
    print("\nBody alone (SU2 body_t4_h2/h4 meshes):")
    for key in ("M2_a2", "M3_a4", "M5_a8"):
        m, al = float(key[1:].split("_")[0]), float(key.split("_a")[1])
        runs = [v for k, v in su2.items() if k.startswith("body_t4") and k.endswith(key)]
        i = int(np.where((data.mach == m) & (a == al) & (b == 0))[0][0])
        nb = parts["Newtonian body (Modified)"][i]
        x = cmp.X_REF - nb[4] * cmp.D / nb[2]
        su = ", ".join(f"CN {r['cn']:.3f} x_cp {r['cp_mm'] / 1000:.3f} ({r['mesh']})" for r in runs)
        print(f"  Mach {m:g} alpha {al:g}: SU2 {su}; Newtonian body CN {nb[2]:.3f} x_cp {x:.3f}; slender-body 2 alpha "
              f"{2 * math.radians(al):.3f}")
    json.dump(out, open(os.path.join(HERE, "results", "su2_compare.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
