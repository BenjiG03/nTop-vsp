"""Aero database: the sum of any number of aerodynamic contributions on one flight condition, driven from nTop.

    python -m vsp_ntop.run_database SETTINGS_FILE
"""
from __future__ import annotations

import csv
import itertools
import os
import re

import numpy as np

from .common import DriverError, Log, Settings, main_wrapper, write_csv
from .flight import Flight, same
from .geometry import ls_slope, wind_axes

COEFFICIENTS = ("CL", "CD", "CY", "Cl", "Cm", "Cn")
TOL = 1e-6


# --- axes

def wind_to_body(cd, cy, cl, alpha, beta):
    v, lift, side = wind_axes(alpha, beta)
    return cd * v + cy * side + cl * lift


def body_to_wind(f, alpha, beta):
    v, lift, side = wind_axes(alpha, beta)
    return lift @ f, v @ f, side @ f  # CL, CD, CY


# --- interpolation

def _key(v):
    return round(v, 6)


class Grid:
    """Values on a Mach x alpha x beta tensor grid (some axes may hold one value, some points may be missing),
    interpolated multilinearly."""

    def __init__(self, name, points):
        self.name = name
        keys = np.round([p for p, _ in points], 6).reshape(-1, 3)
        seen = set()
        for k, (m, a, b) in zip(map(tuple, keys), (p for p, _ in points)):
            if k in seen:
                raise DriverError(f"{name}: two sets of values at Mach {m:g}, alpha {a:g}, beta {b:g}")
            seen.add(k)
        self.axes = [np.unique(keys[:, i]) for i in range(3)]
        idx = tuple(np.searchsorted(ax, keys[:, i]) for i, ax in enumerate(self.axes))
        self.values = np.full(tuple(len(ax) for ax in self.axes) + (len(points[0][1]) if points else 0,), np.nan)
        self.values[idx] = [v for _, v in points]

    def _bracket(self, i, x, what, cond):
        """(index, weight) of the grid values around ``x`` on axis ``i``."""
        ax = self.axes[i]
        if len(ax) == 1:
            if abs(ax[0] - x) > TOL:
                raise DriverError(f"{self.name} has {what} {ax[0]:g} only; the database asks for {what} {x:g} ({cond})")
            return [(0, 1.0)]
        if x < ax[0] - TOL or x > ax[-1] + TOL:
            raise DriverError(f"{self.name} covers {what} {ax[0]:g} to {ax[-1]:g}; the database asks for {what} {x:g} "
                              f"({cond})")
        j = int(np.clip(np.searchsorted(ax, x - TOL) - 1, 0, len(ax) - 2))  # the first interval holding x
        t = float(np.clip((x - ax[j]) / (ax[j + 1] - ax[j]), 0.0, 1.0))
        return [(j, 1.0 - t), (j + 1, t)]

    def at(self, m, a, b):
        """The values at a condition it has, else interpolated between the grid points around it."""
        cond = f"Mach {m:g}, alpha {a:g}, beta {b:g}"
        near = [np.flatnonzero(np.abs(ax - x) <= TOL) for ax, x in zip(self.axes, (m, a, b))]
        if all(len(k) for k in near) and not np.isnan(self.values[near[0][0], near[1][0], near[2][0]]).any():
            return self.values[near[0][0], near[1][0], near[2][0]].copy()
        out = 0.0
        for (i, wm), (j, wa), (k, wb) in itertools.product(*(self._bracket(n, x, w, cond) for n, (x, w) in
                                                             enumerate(((m, "Mach"), (a, "alpha"), (b, "beta"))))):
            w = wm * wa * wb
            if w == 0.0:
                continue
            vals = self.values[i, j, k]
            if np.isnan(vals).any():
                mm, aa, bb = self.axes[0][i], self.axes[1][j], self.axes[2][k]
                raise DriverError(f"{self.name} has no values at {cond} and cannot interpolate there: it lacks "
                                  f"Mach {mm:g}, alpha {aa:g}, beta {bb:g}. Interpolation needs every "
                                  "combination of the neighbouring Mach numbers and angles; or compute it at "
                                  "the database's conditions")
            out = out + w * vals
        return out


def interp_mach(name, points, m):
    x, y = np.array(sorted(points), float).reshape(-1, 2).T
    if not len(x) or m < x[0] - TOL or m > x[-1] + TOL:
        raise DriverError(f"{name} was run from Mach {x[0]:g} to {x[-1]:g} but applies at Mach {m:g}: run it "
                          "over the Mach numbers where it applies, or narrow its Mach range")
    exact = np.flatnonzero(np.abs(x - m) <= TOL)
    return float(y[exact[0]] if len(exact) else np.interp(m, x, y))


# --- contributions

def read_table(path, allow_empty=False):
    if not os.path.isfile(path):
        raise DriverError(f"{path} not found; run the analysis first and connect its output")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2 and not allow_empty:
        raise DriverError(f"{path} has no results")
    out = []
    for r in rows[1:]:
        if r:
            row = {}
            for h, v in zip(rows[0], r):
                try:
                    row[h] = float(v)
                except ValueError:
                    row[h] = v
            out.append(row)
    return out


ALT_TOL = 0.5  # metres: flight points at the same altitude


class Contribution:
    """One source of forces and moments."""

    def __init__(self, c, log):
        self.kind = c.get("kind", "")
        self.name = c.get("name", "") or self.kind
        self.mach_min = float(c.get("mach_min", 0.0))
        self.mach_max = float(c.get("mach_max", 1e9))
        self.scale = float(c.get("scale", 1.0))
        self.use = c.get("use", "all")
        self.point = None
        if self.kind == "vspaero":
            self._vspaero(c, log)
        elif self.kind in ("wave", "parasite"):
            self._drag(c, log)
        elif self.kind == "coefficients":
            self._coefficients(c, log)
        else:
            raise DriverError(f"unknown contribution kind '{self.kind}'")

    def _vspaero(self, c, log):
        base = c["results"]
        flight = Flight(Settings(base + ".ntvsp.settings.txt").text("flight"))
        self.sref, self.bref, self.cref, self.point = flight.ref
        rows = read_table(base + ".ntvsp.polar.csv")
        bad = [r for r in rows if r.get("Converged", 1) == 0]
        if bad:
            log(f"warning: {self.name}: {len(bad)} cases did not pass VSPAERO's convergence check")
        tot = c.get("terms", "inviscid") == "total"
        f = ("CFxtot", "CFytot", "CFztot", "CMxtot", "CMytot", "CMztot") if tot else \
            ("CFix", "CFiy", "CFiz", "CMix", "CMiy", "CMiz")
        # One grid per Reynolds number.
        self.grids = {}
        for re_ in sorted({r["Re_1e6"] for r in rows}):
            sub = [r for r in rows if r["Re_1e6"] == re_]
            self.grids[re_ * 1e6] = Grid(self.name, [((r["Mach"], r["AoA"], r["Beta"]), [r[k] for k in f]) for r in sub])
        log(f"{self.name}: VSPAERO {'total' if tot else 'inviscid'} coefficients, {len(rows)} cases "
            f"({len(self.grids)} Reynolds numbers) from {base}")

    def _drag(self, c, log):
        base = c["results"]
        flight = Flight(Settings(base + ".ntvsp.settings.txt").text("flight"))
        self.sref, self.bref, self.cref = flight.ref[0], 1.0, 1.0
        self.point = None  # no moment information
        if self.kind == "wave":
            rows = read_table(base + ".ntvsp.wave.csv", allow_empty=True)
            self.mach_points = {None: [(r["Mach"], r["CD_wave"]) for r in rows]}
        else:
            rows = read_table(base + ".ntvsp.parasite.csv")
            self.mach_points = {}
            for r in rows:
                self.mach_points.setdefault(r["Altitude_m"], []).append((r["Mach"], r["Total_CD"]))
        ms = [m for v in self.mach_points.values() for m, _ in v]
        log(f"{self.name}: {'Wave' if self.kind == 'wave' else 'Parasite'} Drag at "
            + (f"Mach {min(ms):g} to {max(ms):g}" if ms else "no Mach number") + f" from {base}")

    def _coefficients(self, c, log):
        base = c["results"]
        rows = read_table(base + ".ntvsp.coefficients.csv", allow_empty=True)
        self.by_altitude = rows and "Altitude" in rows[0]
        refs = {tuple(round(r[k], 9) for k in ("Sref", "Bref", "Cref", "Xref", "Yref", "Zref")) for r in rows}
        if len(refs) > 1:
            raise DriverError(f"{self.name}: its rows have different reference values or moment points")
        if refs:
            self.sref, self.bref, self.cref, *p = refs.pop()
            self.point = tuple(p)
        wind = rows and "CD" in rows[0]
        groups = {}
        for r in rows:
            if wind:  # CD, CY, CL, Cl, Cm, Cn (body-axis signs) -> x aft, y right, z up
                force = wind_to_body(r["CD"], r["CY"], r["CL"], r["AoA"], r["Beta"])
                c6 = list(force) + [-r["Cl"], r["Cm"], -r["Cn"]]
            else:
                c6 = [r[k] for k in ("CA", "CY", "CN", "Cl", "Cm", "Cn")]
            groups.setdefault(r["Altitude"] if self.by_altitude else None, []).append(((r["Mach"], r["AoA"], r["Beta"]), c6))
        self.grids = {h: Grid(self.name, pts) for h, pts in groups.items()}
        machs = [p[0][0] for pts in groups.values() for p in pts]
        log(f"{self.name}: {len(rows)} coefficient rows" + (f", Mach {min(machs):g} to {max(machs):g}" if machs else "")
            + (" by altitude" if self.by_altitude else "") + f", from {base}")

    def applies(self, m):
        return self.mach_min - TOL <= m < self.mach_max - TOL

    def _at_altitude(self, table, h):
        for key, v in table.items():
            if key is None or abs(key - h) <= ALT_TOL:
                return v
        return None

    def body(self, cond):
        m, h, a, b = cond["Mach"], cond["Altitude"], cond["AoA"], cond["Beta"]
        if self.kind == "vspaero":
            grid = self.grids[min(self.grids, key=lambda re_: abs(np.log(max(re_, 1.0)) -
                                                                  np.log(max(cond["Re_cref"], 1.0))))]
            if not grid.axes[0][0] - TOL <= m <= grid.axes[0][-1] + TOL:
                return None  # outside the Mach numbers VSPAERO ran
            v = grid.at(m, a, b)
        elif self.kind == "coefficients":
            grid = self._at_altitude(self.grids, h)
            if grid is None or not len(grid.axes[0]) or not grid.axes[0][0] - TOL <= m <= grid.axes[0][-1] + TOL:
                return None  # computed only over its own Mach range (and altitudes)
            v = grid.at(m, a, b)
        else:
            pts = self._at_altitude(self.mach_points, h)
            if pts is None:
                raise DriverError(f"{self.name} has no flight point at altitude {h:g} m: run it on the database's "
                                  "flight condition")
            if not pts or not min(x for x, _ in pts) - TOL <= m <= max(x for x, _ in pts) + TOL:
                return None  # outside the Mach numbers it ran (wave drag: subsonic)
            cd = interp_mach(self.name, pts, m)
            v = np.concatenate([wind_to_body(cd, 0.0, 0.0, a, b), np.zeros(3)])
        force, moment = np.asarray(v[:3], float), np.asarray(v[3:], float)
        vel = wind_axes(a, b)[0]
        if self.use == "drag":
            force, moment = (force @ vel) * vel, np.zeros(3)
        elif self.use == "nodrag":
            force = force - (force @ vel) * vel
        return force, moment

    def on(self, ref, cond):
        """Force and moment coefficients over the database reference ``ref`` = (S, b, c, point), about its point,
        or None where the contribution does not apply."""
        got = self.body(cond)
        if got is None:
            return None
        force, moment = got
        k, lengths = self.scale * self.sref, np.array([self.bref, self.cref, self.bref])
        f = k * force  # per dynamic pressure
        arm = np.subtract(self.point or ref[3], ref[3])
        mom = k * moment * lengths + np.cross(arm, f)
        S, B, C = ref[0], ref[1], ref[2]
        return f / S, mom / (S * np.array([B, C, B]))


def parse_contributions(path):
    """``contribution.<key> = <value>`` lines, as Aero Contribution (VSP) writes them: each ``contribution.kind``
    starts a contribution, and the following ``contribution.*`` lines are its settings."""
    out = []
    for line in open(path, encoding="utf-8"):
        key, sep, value = line.partition("=")
        key = key.strip()
        if not sep or not key.startswith("contribution."):
            continue
        field = key[len("contribution."):]
        if field == "kind":
            out.append({"kind": value.strip()})
        elif not out:
            raise DriverError(f"'{line.strip()}' comes before any contribution.kind line")
        else:
            out[-1][field] = value.strip()
    return out


def coefficients(force, moment, a, b):
    cl, cd, cy = body_to_wind(force, a, b)
    return [float(v) for v in (cl, cd, cy, -moment[0], moment[1], -moment[2])]


def column_name(name):
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_") or "part"


def centers(work, case, rows, ref, log):
    """Centre of pressure per row and aerodynamic centre per (Mach, altitude, beta) alpha sweep, from the totals
    (x_cp = x_ref - Cm cref / CN; x_ac = x_ref - dCm/dCN cref), in metres."""
    xr, yr, zr = ref[3]
    cref = ref[2]
    cp = []
    for k, r in enumerate(rows, start=1):
        ok = abs(r["CN"]) >= 1e-6
        cp.append([k, r["Mach"], r["Altitude"], r["AoA"], r["Beta"], xr - r["Cm"] * cref / r["CN"] if ok else xr, yr,
                   zr, int(ok), r["CN"], r["Cm"]])
    write_csv(os.path.join(work, case + ".ntvsp.center_of_pressure.csv"),
              ["Case", "Mach", "Altitude", "AoA", "Beta", "X", "Y", "Z", "Defined", "CN", "Cm"], cp)
    groups = {}
    for r in rows:
        groups.setdefault((r["Mach"], r["Altitude"], r["Beta"]), []).append(r)
    ac = []
    for (m, h, b), g in groups.items():
        slope = ls_slope([r["CN"] for r in g], [r["Cm"] for r in g])
        if slope is None:
            continue
        ac.append([m, h, b, xr - slope * cref, yr, zr, slope])
    write_csv(os.path.join(work, case + ".ntvsp.aerodynamic_center.csv"),
              ["Mach", "Altitude", "Beta", "X", "Y", "Z", "dCm_dCN"], ac)
    if ac:
        log("aerodynamic centre x: " + ", ".join(f"Mach {r[0]:g} ({r[1]:g} m) beta {r[2]:g}: {r[3]:.3f} m" for r in ac))


def main(settings_path):
    s = Settings(settings_path)
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    flight = Flight(s.text("flight"))
    found = parse_contributions(settings_path)
    expected = s.integer("contribution_count", len(found))
    if len(found) != expected:
        raise DriverError(f"Contributions has {expected} items but {len(found)} are contributions: connect outputs of "
                          "Aero Contribution (VSP) or the per-source contribution blocks, not the analyses themselves")
    parts = [Contribution(c, log) for c in found]
    if not parts:
        raise DriverError("no contributions connected")
    names, seen = [], {}
    for p in parts:
        n = column_name(p.name)
        seen[n] = seen.get(n, 0) + 1
        names.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    ref = flight.ref
    conds = flight.conditions
    log(f"database: {flight.describe()}; Sref {ref[0]:g} m^2, bref {ref[1]:g} m, cref {ref[2]:g} m, moments about "
        f"({ref[3][0]:g}, {ref[3][1]:g}, {ref[3][2]:g}) m")

    variables = [x.strip() for x in s.text("variables", "").split(";") if x.strip()]
    values = s.reals("values") if variables else []
    if len(variables) != len(values):
        raise DriverError(f"{len(variables)} variable names but {len(values)} values")
    header = variables + ["Mach", "Altitude", "AoA", "Beta"] + list(COEFFICIENTS) + ["L_over_D", "CA", "CN"] + \
        [f"{c}_{n}" for n in names for c in COEFFICIENTS]
    rows, table, used = [], [], {}
    for cond in conds:
        m, h, a, b = cond["Mach"], cond["Altitude"], cond["AoA"], cond["Beta"]
        tf, tm, own = np.zeros(3), np.zeros(3), []
        for p, n in zip(parts, names):
            got = p.on(ref, cond) if p.applies(m) else None
            if got is None:
                own += [0.0] * 6
                continue
            f, mo = got
            used.setdefault((m, h), []).append(n) if n not in used.get((m, h), []) else None
            tf, tm = tf + f, tm + mo
            own += coefficients(f, mo, a, b)
        tot = coefficients(tf, tm, a, b)
        rows.append(dict(zip(COEFFICIENTS, tot), Mach=m, Altitude=h, AoA=a, Beta=b, CA=float(tf[0]), CN=float(tf[2])))
        table.append(values + [m, h, a, b] + tot + [tot[0] / tot[1] if tot[1] else 0.0, float(tf[0]), float(tf[2])]
                     + own)
    for p in flight.points:
        key = next((k for k in used if same(k[0], p["Mach"]) and abs(k[1] - p["Altitude"]) <= ALT_TOL), None)
        log(f"Mach {p['Mach']:g} at {p['Altitude']:g} m: " + (", ".join(used[key]) if key else "no contribution applies"))
    write_csv(os.path.join(work, case + ".ntvsp.database.csv"), header, table)
    log(f"database: {len(table)} rows, columns {', '.join(header)}")
    centers(work, case, rows, ref, log)

    target = s.text("append_to", "").strip()
    if target:
        new = not os.path.isfile(target) or os.path.getsize(target) == 0
        if not new:
            with open(target, newline="", encoding="utf-8") as f:
                old = next(csv.reader(f), [])
            if old != header:
                raise DriverError(f"{target} has other columns ({', '.join(old)}); append to a new file")
        with open(target, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new:
                w.writerow(header)
            w.writerows(table)
        log(f"appended {len(table)} rows to {target}")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
