"""CN and Cm against alpha on the trapezoid-section rocket: SU2 and the nTop_vsp combinations (results/cn_cm.png).

    sh build/py.sh validation/noncircular/plots.py

Reads results/compare.json (compare.py). Colours: the dataviz reference palette's first four categorical slots in
fixed order; SU2, the reference, in primary ink with markers.
"""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SERIES = (("Exposed + Carry-over + Body Lift", "#2a78d6", "-"), ("Exposed + Body Lift (shape)", "#eb6834", "--"),
          ("Exposed + Newtonian body", "#1baf7a", "-."), ("RASAero II", "#eda100", ":"))
INK, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"


def main():
    d = json.load(open(os.path.join(HERE, "results", "compare.json")))
    machs = sorted({c["Mach"] for c in d["conditions"]})
    fig, axes = plt.subplots(2, len(machs), figsize=(4.2 * len(machs), 6.6), sharex=True, facecolor=SURFACE)
    for col, m in enumerate(machs):
        rows = sorted((c for c in d["conditions"] if c["Mach"] == m), key=lambda c: c["AoA"])
        a = [c["AoA"] for c in rows]
        for row, key in enumerate(("CN", "Cm")):
            ax = axes[row, col]
            ax.set_facecolor(SURFACE)
            for name, colour, style in SERIES:
                ax.plot(a, [c["totals"][name][key] for c in rows], style, color=colour, lw=2, label=name)
            su = [(c["AoA"], d["su2"].get(f"M{m:g}_a{c['AoA']:g}")) for c in rows]
            su = [(x, s[key]) for x, s in su if s]
            if su:
                ax.plot(*zip(*su), "-o", color=INK, lw=1.2, ms=6, label="SU2 (inviscid)", zorder=5)
            ax.set_title(f"Mach {m:g}: {key}" + (" about x = 0.7 m" if key == "Cm" else ""), fontsize=10, color=INK,
                         loc="left")
            ax.grid(color=GRID, lw=0.6)
            for s in ax.spines.values():
                s.set_color(AXIS)
            ax.tick_params(colors=MUTED, labelsize=8)
            ax.axhline(0.0, color=AXIS, lw=0.8)
            if row == 1:
                ax.set_xlabel("angle of attack (deg), beta 0", color=MUTED, fontsize=9)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, fontsize=8, frameon=False, labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(HERE, "results", "cn_cm.png")
    fig.savefig(out, dpi=110, facecolor=SURFACE)
    print("written", out)


if __name__ == "__main__":
    main()
