"""Surface pressure point maps shared by the analyses, read in nTop by Surface Pressure (VSP) and Surface Pressure
Jump (VSP)."""
from __future__ import annotations

import glob
import os

from .flight import same

KINDS = ("pressure", "pressure_jump")
LINE = "%.9g,%.9g,%.9g,%.9g\n"
CHUNK = 100_000  # rows formatted at once


def path(base, case, kind="pressure"):
    return f"{base}.ntvsp.{kind}.{int(case)}.csv"


def clear(base):
    """Remove pressure files of an earlier run under ``base``; returns how many."""
    n = 0
    for kind in KINDS:
        for f in glob.glob(glob.escape(base) + f".ntvsp.{kind}.*.csv"):
            os.remove(f)
            n += 1
    return n


def cases_at(flight, mach, aoa, beta):
    """The flight condition rows (dicts) with this Mach number, angle of attack and sideslip, at any altitude."""
    return [c for c in flight.conditions if same(c["Mach"], mach) and same(c["AoA"], aoa) and same(c["Beta"], beta)]


def write(base, case, points, values, kind="pressure"):
    """Write one point map: ``points`` (n, 3) in metres, ``values`` (n) in pascals; numpy arrays or sequences."""
    target = path(base, case, kind)
    try:
        import numpy as np
        if isinstance(points, np.ndarray) or isinstance(values, np.ndarray):
            rows = np.column_stack([np.asarray(points, float), np.asarray(values, float)])
            with open(target, "w", encoding="utf-8") as f:  # one % per chunk
                for start in range(0, len(rows), CHUNK):
                    part = rows[start:start + CHUNK]
                    f.write((LINE * len(part)) % tuple(part.ravel()))
            return target
    except ImportError:  # pragma: no cover
        pass
    with open(target, "w", encoding="utf-8") as f:
        for (x, y, z), v in zip(points, values):
            f.write(f"{x:.9g},{y:.9g},{z:.9g},{v:.9g}\n")
    return target
