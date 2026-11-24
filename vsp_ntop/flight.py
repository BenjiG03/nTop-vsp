"""Flight conditions shared by every analysis: Mach and altitude through the US Standard Atmosphere 1976."""
from __future__ import annotations

import csv
import os

from .common import DriverError, Settings

POINT_COLUMNS = ("Point", "Mach", "Altitude", "Altitude_ft", "Temperature", "Pressure", "Density", "Speed_of_Sound",
                 "Velocity", "Dynamic_Pressure", "Dynamic_Viscosity", "Re_per_m", "Re_cref")
REFERENCE_COLUMNS = ("Sref", "Bref", "Cref", "Xref", "Yref", "Zref")
CONDITION_COLUMNS = ("Case", "Point", "Mach", "Altitude", "AoA", "Beta", "Temperature", "Pressure", "Density",
                     "Velocity", "Dynamic_Pressure", "Re_cref") + REFERENCE_COLUMNS
INVISCID_COLUMNS = ("Mach", "AoA", "Beta") + REFERENCE_COLUMNS
FEET = 0.3048
TOL = 1e-6


def same(a, b, tol=TOL):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _rows(path):
    if not os.path.isfile(path):
        raise DriverError(f"{path} not found; connect the output of Flight Condition (VSP)")
    with open(path, newline="", encoding="utf-8") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


class Flight:
    """A Flight Condition's results: ``points``, ``conditions`` (dicts of the CSV columns) and ``ref`` =
    (Sref, bref, cref, (x, y, z)) in metres."""

    def __init__(self, base):
        base = base.strip()
        if not base:
            raise DriverError("no Flight Condition connected")
        self.base = base
        s = Settings(base + ".ntvsp.settings.txt")
        self.ref = (s.real("sref"), s.real("bref"), s.real("cref"), (s.real("xref"), s.real("yref"), s.real("zref")))
        self.points = _rows(base + ".ntvsp.points.csv")
        self.conditions = _rows(base + ".ntvsp.flight.csv")

    def machs(self):
        """Distinct Mach numbers, in order of first appearance."""
        out = []
        for p in self.points:
            if not any(same(p["Mach"], m) for m in out):
                out.append(p["Mach"])
        return out

    def angles(self, key):
        out = []
        for c in self.conditions:
            if not any(same(c[key], a) for a in out):
                out.append(c[key])
        return out

    def point_at(self, mach, re_cref=None):
        """The flight point with this Mach number (and Reynolds number, when several share the Mach number)."""
        cands = [p for p in self.points if same(p["Mach"], mach)]
        if not cands:
            return None
        if re_cref is None or len(cands) == 1:
            return cands[0]
        return min(cands, key=lambda p: abs(p["Re_cref"] - re_cref))

    def describe(self):
        n = len(self.points)
        machs = ", ".join("%g" % p["Mach"] for p in self.points)
        alts = ", ".join("%g" % p["Altitude"] for p in self.points)
        return (f"{n} flight point{'s' if n != 1 else ''} (Mach {machs}; altitude {alts} m), "
                f"{len(self.angles('AoA'))} angles of attack, {len(self.angles('Beta'))} sideslips: "
                f"{len(self.conditions)} conditions")
