"""Polar plots of aero databases with matplotlib, driven from nTop.

    python -m vsp_ntop.run_plot SETTINGS_FILE [--show]

Settings: ``results`` (``;``-separated results paths of Aero Database runs, or database CSV files), ``labels``
(``;``-separated, one per result; empty: the file names), ``x`` and ``y`` (database columns), ``along`` (the
condition that varies along each curve: AoA, Beta or Mach; every combination of the other two is one curve),
``title``, ``show`` (1: open the plot in a window), ``work_dir`` and ``case``.

Writes ``<case>.png`` in the working directory. With ``show`` the same figure opens in a window from a separate,
detached process (``--show``), so nTop does not wait for the window to be closed.
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys

from .common import DriverError, Log, Settings, conda_dll_dirs

AXIS_LABELS = {"AoA": "angle of attack (deg)", "Beta": "sideslip (deg)", "Mach": "Mach number", "Altitude": "altitude (m)",
               "L_over_D": "L/D"}
SYMBOLS = {"AoA": "α", "Beta": "β", "Mach": "M", "Altitude": "h (m)"}
FILTERS = {"beta0": "Beta", "alpha0": "AoA"}  # keep only the rows where this condition is zero
# Fixed-order categorical hues (one per result, when the results are all that differ) and, for ordered conditions,
# one blue ramp from light to dark.
CATEGORICAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
STYLES = ("-", "--", ":", "-.")
MARKERS = ("o", "s", "^", "D")


def RAMP(t):
    """Blue ramp, t from 0 (light) to 1 (dark): #86b6ef -> #0d366b."""
    a, b = (0x86, 0xb6, 0xef), (0x0d, 0x36, 0x6b)
    return "#" + "".join(f"{round(p + (q - p) * t):02x}" for p, q in zip(a, b))


def table(path):
    p = path if path.lower().endswith(".csv") else path + ".ntvsp.database.csv"
    if not os.path.isfile(p):
        raise DriverError(f"{p} not found; connect the output of Aero Database (VSP)")
    with open(p, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise DriverError(f"{p} has no rows")
    return rows


def figure(s: Settings):
    import matplotlib.pyplot as plt

    results = [r for r in s.text("results").split(";") if r.strip()]
    labels = [x for x in s.text("labels", "").split(";")]
    x, y, along = s.text("x"), s.text("y"), s.text("along", "AoA")
    others = [c for c in ("Mach", "Altitude", "Beta", "AoA") if c != along]
    curves = []
    for i, path in enumerate(results):
        rows = table(path.strip())
        for col in (x, y, along):
            if col not in rows[0]:
                raise DriverError(f"{path}: no column {col}")
        label = labels[i].strip() if i < len(labels) and labels[i].strip() else os.path.basename(path.strip())
        groups = {}
        keep = s.text("filter", "all")
        if keep in FILTERS:
            rows = [r for r in rows if abs(float(r[FILTERS[keep]])) < 1e-9]
        for r in rows:
            groups.setdefault(tuple(round(float(r.get(c, 0.0)), 9) for c in others), []).append(r)
        for key, g in groups.items():
            g.sort(key=lambda r: float(r[along]))
            curves.append((label, key, [float(r[x]) for r in g], [float(r[y]) for r in g]))
    if not curves:
        raise DriverError("no rows left to plot (check the filter)")
    # Name only what tells the curves apart.
    vary = [j for j in range(len(others)) if len({c[1][j] for c in curves}) > 1]
    names = list(dict.fromkeys(c[0] for c in curves))
    keys = sorted({c[1] for c in curves})
    if len(names) > len(STYLES) and len(keys) > 1:
        raise DriverError(f"{len(names)} results with several curves each: plot at most {len(STYLES)} at a time")
    if len(keys) == 1 and len(names) > len(CATEGORICAL):
        raise DriverError(f"{len(names)} results: plot at most {len(CATEGORICAL)} at a time")
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, key, xs, ys in curves:
        r = names.index(label)
        if len(keys) == 1:  # the results are all that differ: one hue each, in a fixed order
            color, style = CATEGORICAL[r], "-"
        else:  # the conditions are ordered: one hue, light to dark; results by line style
            t = keys.index(key) / max(1, len(keys) - 1)
            color, style = RAMP(t), STYLES[r]
        parts = ([label] if len(names) > 1 or not vary else []) + [f"{SYMBOLS[others[j]]} {key[j]:g}" for j in vary]
        ax.plot(xs, ys, style, color=color, marker=MARKERS[r % len(MARKERS)], markersize=4, linewidth=1.6,
                label=", ".join(parts))
    ax.set_xlabel(AXIS_LABELS.get(x, x))
    ax.set_ylabel(AXIS_LABELS.get(y, y))
    title = s.text("title", "").strip()
    ax.set_title(title or f"{y} against {x}" + (f", along {AXIS_LABELS.get(along, along)}" if along != x else ""))
    ax.grid(True, color="#d9d9d9", linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if len(curves) > 1:
        ax.legend(fontsize=8, ncol=2 if len(curves) > 14 else 1, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False)
    fig.tight_layout()
    return fig, len(curves)


def main(settings_path, show_only=False):
    conda_dll_dirs()
    s = Settings(settings_path)
    if show_only:
        import matplotlib.pyplot as plt
        figure(s)
        plt.show()
        return
    import matplotlib
    matplotlib.use("Agg")
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    fig, n = figure(s)
    png = os.path.join(work, case + ".png")
    fig.savefig(png, dpi=150)
    log(f"{n} curves of {s.text('y')} against {s.text('x')}: {png}")
    if s.flag("show", 1):
        show(settings_path, os.path.join(work, case + ".window.log"))
        log("opened the plot in a window")
    log.close()


def show(settings_path, err_path):
    """Start the window in its own process, detached from nTop's: nTop waits for nothing and closing it later
    does not close the window (nTop runs commands in a job, which the window breaks away from where allowed)."""
    cmd = [sys.executable, "-m", "vsp_ntop.run_plot", settings_path, "--show"]
    err = open(err_path, "w")
    kw = dict(stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=err, close_fds=True)
    if os.name != "nt":
        subprocess.Popen(cmd, start_new_session=True, **kw)
        return
    base = subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
    try:
        subprocess.Popen(cmd, creationflags=base | 0x01000000, **kw)  # CREATE_BREAKAWAY_FROM_JOB
        err.write("window process started outside nTop's job\n")
    except OSError:
        subprocess.Popen(cmd, creationflags=base, **kw)
        err.write("window process started inside nTop's job (breakaway not allowed)\n")
    err.flush()


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--show"]
    if len(args) != 1:
        sys.stderr.write("usage: python -m vsp_ntop.run_plot SETTINGS_FILE [--show]\n")
        sys.exit(2)
    try:
        main(args[0], show_only="--show" in sys.argv[1:])
    except ImportError as e:
        sys.stderr.write(f"ERROR: {e}: install matplotlib in the Python of NTOP_PYTHON_PATH (pip install matplotlib)\n")
        sys.exit(1)
    except DriverError as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
