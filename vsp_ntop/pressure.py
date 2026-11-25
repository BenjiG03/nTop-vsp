"""Surface pressure point maps shared by the analyses, read in nTop by Surface Pressure (VSP) and Surface Pressure
Jump (VSP)."""
from __future__ import annotations

import glob
import os

from .common import write_point_map
from .flight import same

KINDS = ("pressure", "pressure_jump")


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
    """Write one point map: ``points`` (n, 3) in metres, ``values`` (n) in pascals."""
    target = path(base, case, kind)
    write_point_map(target, points, values)
    return target
