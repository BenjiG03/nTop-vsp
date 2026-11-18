"""Plots of the best combinations against RASAero II (``results/*.png``).

    sh build/py.sh validation/rasaero/plots.py NAME [NAME ...]

NAMEs are combination names from compare.py (``results/scores.json``); with none, the three best.
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as cmp  # noqa: E402

OUT = os.path.join(HERE, "results")
SUFFIX = os.environ.get("PLOT_SUFFIX", "")


def main():
    names = sys.argv[1:]
    data = cmp.Data()
    if not names:
        names = [r["name"] for r in json.load(open(os.path.join(OUT, "scores.json")))[:3]]
    totals = {n: t for n, _, t in cmp.candidates(data) if n in names}
    missing = [n for n in names if n not in totals]
    if missing:
        raise SystemExit(f"unknown: {missing}")
    a = np.array([c["AoA"] for c in data.conds])
    b = np.array([c["Beta"] for c in data.conds])
    ref = data.rasaero
    styles = ["-", "--", "-.", ":"]

    def label(n):
        return n if len(n) < 70 else n[:67] + "..."

    # 1. CN and Cm against alpha at beta 0, four Mach numbers
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), sharex=True)
    for col, m in enumerate((1.5, 2.0, 3.0, 5.0)):
        sel = (data.mach == m) & (b == 0)
        order = np.argsort(a[sel])
        for row, (j, what) in enumerate(((2, "CN"), (4, "Cm about x = 0.7 m"))):
            ax = axes[row, col]
            ax.plot(a[sel][order], ref[sel][order, j], "k-o", ms=3, lw=2, label="RASAero II")
            for s, n in zip(styles, names):
                ax.plot(a[sel][order], totals[n][sel][order, j], s, label=label(n))
            ax.set_title(f"Mach {m:g}: {what}", fontsize=9)
            ax.grid(alpha=0.3)
            if row == 1:
                ax.set_xlabel("alpha (deg), beta 0")
    axes[0, 0].legend(fontsize=6, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"cn_cm_alpha{SUFFIX}.png"), dpi=110)

    # 2. CA against Mach at alpha 0, beta 0, and the centre of pressure at alpha 8
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
    machs = sorted(set(data.mach))
    at = lambda arr, al, j: [arr[(data.mach == m) & (a == al) & (b == 0)][0, j] for m in machs]  # noqa: E731
    ax1.plot(machs, at(ref, 0.0, 0), "k-o", lw=2, label="RASAero II CA")
    ax1.plot(machs, [x - data.base_drag[data.mach == m][0] for x, m in zip(at(ref, 0.0, 0), machs)], "k:o",
             label="RASAero II CA less base drag")
    for s, n in zip(styles, names):
        ax1.plot(machs, at(totals[n], 0.0, 0), s, label=label(n))
    ax1.set_xlabel("Mach (alpha 0)")
    ax1.set_ylabel("CA")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=6)

    def xcp(arr):
        return [cmp.X_REF - arr[(data.mach == m) & (a == 8.0) & (b == 0)][0, 4] * cmp.D /
                arr[(data.mach == m) & (a == 8.0) & (b == 0)][0, 2] for m in machs]
    ax2.plot(machs, xcp(ref), "k-o", lw=2, label="RASAero II")
    for s, n in zip(styles, names):
        ax2.plot(machs, xcp(totals[n]), s, label=label(n))
    ax2.set_xlabel("Mach (alpha 8)")
    ax2.set_ylabel("centre of pressure from the nose (m)")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"ca_xcp_mach{SUFFIX}.png"), dpi=110)
    print("written", OUT)


if __name__ == "__main__":
    main()
