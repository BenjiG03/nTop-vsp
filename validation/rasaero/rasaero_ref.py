"""RASAero II 1.0.2.0 results for the hobby rocket, as force and moment coefficients at any (Mach, alpha, beta).

The Run Test files in ``rasaero_runs/`` come from ``hypanel/validation/rocket/hobby_rocket.CDX1`` (smooth, all
turbulent, sea level, power off), one file per angle of attack: 0, 1, 2, 4, 6 and 8 degrees from hypanel's comparison,
10, 12, 14, 15 and 16 run for this study with ``hypanel/validation/rocket/rasaero_driver/run_test.ps1``.

RASAero II treats the rocket as axisymmetric: it gives the normal force CN, the axial force CA and the center of
pressure at a total angle of attack. A cruciform rocket at angle of attack alpha and sideslip beta sees the total
angle alpha_t = acos(cos alpha cos beta) in the plane of the crossflow, so here (an approximation that ignores the
small effect of roll orientation on a four-fin rocket):

* the normal force CN(alpha_t) acts along the crossflow (v_y, v_z) = (-sin beta, sin alpha cos beta) of the freestream
  (cos a cos b, -sin b, sin a cos b), giving CY and CN in the model's axes (x aft, y, z up);
* CA = CA(alpha_t);
* moments about the moment point (x_ref on the axis) come from the normal force at the center of pressure:
  M = (x_cp - x_ref, 0, 0) x (CA, CY, CN) over Sref d (Cl = 0).

That is the layout of a Modified Newtonian coefficient table (moments about x aft, y, z up, all over Sref cref), so
it goes through the Aero Database like the analyses. Values are interpolated linearly in alpha_t; the Run Test files
list every 0.01 in Mach.
"""
from __future__ import annotations

import bisect
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "rasaero_runs")
INCH = 0.0254
ALPHAS = (0, 1, 2, 4, 6, 8, 10, 12, 14, 15, 16)


def read_run_test(path):
    """Run Test rows by Mach number: six numeric rows per Mach (RASAero II user's manual, Figure 108)."""
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if re.match(r"\s*\d", line):
                rows.append([float(v) for v in line.split()])
    if len(rows) % 6:
        raise ValueError(f"{path}: {len(rows)} rows do not form groups of six")
    table = {}
    for i in range(0, len(rows), 6):
        grp = rows[i:i + 6]
        e = {"alpha": grp[2][1], "cn": grp[3][2], "x_cp": grp[3][3] * INCH, "cl": grp[4][2], "cd": grp[4][3],
             "ca": grp[4][5], "cd0": grp[0][2], "cn_alpha": grp[1][1], "x_cp_low": grp[1][2] * INCH}
        if len(grp[0]) == 14:
            names = ("body_friction", "nose_wave", "base", "fin_friction", "fin_wave", "fin_interference", "fin_base",
                     "other_wave", "protuberance", "reynolds")
            e.update(dict(zip(names, grp[0][4:])))
        table[round(grp[0][0], 2)] = e
    return table


class RASAero:
    def __init__(self, runs=RUNS, tag="turbulent", alphas=ALPHAS):
        self.alphas = list(alphas)
        self.tables = [read_run_test(os.path.join(runs, f"runtest_{tag}_alpha{a:g}.txt")) for a in alphas]

    def at(self, mach, alpha_t):
        """(CN, CA, x_cp in metres from the nose tip) at a total angle of attack (degrees), interpolated."""
        m = round(mach, 2)
        a = abs(alpha_t)
        if a > self.alphas[-1] + 1e-9:
            raise ValueError(f"total angle {a:g} beyond the RASAero runs ({self.alphas[-1]:g})")
        j = min(max(bisect.bisect_right(self.alphas, a), 1), len(self.alphas) - 1)
        a0, a1 = self.alphas[j - 1], self.alphas[j]
        t = (a - a0) / (a1 - a0)
        e0, e1 = self.tables[j - 1][m], self.tables[j][m]
        cn = (1 - t) * e0["cn"] + t * e1["cn"]
        ca = (1 - t) * e0["ca"] + t * e1["ca"]
        if cn > 1e-9:  # centre of pressure: moment-weighted, so the moment is linear between the runs
            mom = (1 - t) * e0["cn"] * e0["x_cp"] + t * e1["cn"] * e1["x_cp"]
            x_cp = mom / cn
        else:
            x_cp = e0["x_cp_low"]
        return cn, ca, x_cp

    def breakdown(self, mach):
        return self.tables[0][round(mach, 2)]

    def coefficients(self, mach, alpha, beta, x_ref=0.0, d=0.054):
        """CA, CY, CN, Cl, Cm, Cn in the model's axes (x aft, y, z up), moments about (x_ref, 0, 0) over Sref d."""
        a, b = math.radians(alpha), math.radians(beta)
        vy, vz = -math.sin(b), math.sin(a) * math.cos(b)
        alpha_t = math.degrees(math.acos(max(-1.0, min(1.0, math.cos(a) * math.cos(b)))))
        cn, ca, x_cp = self.at(mach, alpha_t)
        s = math.hypot(vy, vz)
        cy, cz = (cn * vy / s, cn * vz / s) if s > 1e-12 else (0.0, 0.0)
        arm = x_cp - x_ref
        return {"CA": ca, "CY": cy, "CN": cz, "Cl": 0.0, "Cm": -arm * cz / d, "Cn": arm * cy / d,
                "alpha_t": alpha_t, "x_cp": x_cp}


def write_table(path, conditions, ref, x_ref=0.0, d=0.054):
    """A coefficient table (Modified Newtonian layout) of RASAero II at ``conditions`` [(Mach, AoA, Beta)]."""
    sref = math.pi * d * d / 4.0
    rows = []
    for m, al, be in conditions:
        c = ref.coefficients(m, al, be, x_ref, d)
        rows.append([m, al, be] + [c[k] for k in ("CA", "CY", "CN", "Cl", "Cm", "Cn")] + [sref, d, d, x_ref, 0.0, 0.0])
    with open(path, "w", encoding="utf-8") as f:
        f.write("Mach,AoA,Beta,CA,CY,CN,Cl,Cm,Cn,Sref,Bref,Cref,Xref,Yref,Zref\n")
        for r in rows:
            f.write(",".join(repr(float(v)) for v in r) + "\n")


if __name__ == "__main__":
    r = RASAero()
    for m in (1.5, 2.0, 3.0, 5.0):
        bd = r.breakdown(m)
        print(f"Mach {m}: CD0 {bd['cd0']:.3f} (base {bd.get('base', 0):.3f}, fin wave {bd.get('fin_wave', 0):.3f}, "
              f"nose wave {bd.get('nose_wave', 0):.3f}, friction {bd.get('body_friction', 0) + bd.get('fin_friction', 0):.3f})")
        for al, be in ((4, 0), (15, 0), (15, 4), (-8, -4)):
            c = r.coefficients(m, al, be)
            print(f"   a {al:+3d} b {be:+2d}: " + " ".join(f"{k} {c[k]:+.3f}" for k in ("CA", "CY", "CN", "Cm", "Cn", "x_cp")))
