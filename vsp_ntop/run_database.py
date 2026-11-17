"""Aero database: the sum of any number of aerodynamic contributions, driven from nTop.

    python -m vsp_ntop.run_database SETTINGS_FILE

Each ``contribution = kind: ...; name: ...; ...`` line of the settings file is one source of forces and moments:

* ``vspaero``: a VSPAERO run (``results``), its inviscid (``terms: inviscid``, CFi*, CMi*) or total (``terms:
  total``, with VSPAERO's own skin-friction estimate) coefficients; reference and moment point from its settings.
* ``wave``: a Wave Drag run, CD_wave(Mach), a drag along the freestream through the database's moment point
  (Wave and Parasite Drag give no moments).
* ``parasite``: a Parasite Drag run, Total_CD(Mach), likewise.
* ``rows``: coefficient rows made in nTop, ``|``-separated, each ``axes, Mach, alpha, beta, c1..c6, Sref, bref,
  cref, x, y, z`` (angles in degrees, SI units). ``body`` rows give CA, CY, CN and moments about x aft, y right,
  z up (the axes of VSPAERO, Modified Newtonian and the Viscous Correction blocks); ``wind`` rows give CD, CY, CL
  and moments in aircraft body-axis signs (Cl right wing down, Cm nose up, Cn nose right). Rolling and yawing
  moments are over Sref bref, pitching over Sref cref, all about the row's moment point (x, y, z).

Every contribution also has ``mach_min`` and ``mach_max`` (it applies from Mach Min up to, not including, Mach
Max), ``scale`` (a factor) and ``use``: ``all``, ``drag`` (only the force component along the freestream) or
``nodrag`` (everything but that component).

At each database condition (Mach, alpha, beta) every contribution whose Mach range holds it is taken at that
condition if it has it, else interpolated (multilinear between the neighbouring Mach numbers, angles of attack
and sideslips, every combination of which it must have; in Mach only for wave and parasite drag),
turned into body-axis forces and moments about its own moment point, moved to the database's moment point
(M_db = M + (r - r_db) x F), and normalised by the database's reference. Their sum gives

    CL, CD, CY (wind axes; CY along the side axis, positive right)      Cl, Cm, Cn (body-axis signs)

with VSPAERO's convention for the wind axes: freestream V = (cos a cos b, -sin b, sin a cos b) in x aft, y right,
z up; lift along (-sin a, 0, cos a); side force along (cos a sin b, cos b, sin a sin b); CD = F.V.

Writes ``<case>.ntvsp.database.csv`` (the design variables, the condition, the totals, CA and CN, and each
contribution's own CL...Cn), ``<case>.ntvsp.center_of_pressure.csv`` and ``<case>.ntvsp.aerodynamic_center.csv``
(positions in metres, as the VSPAERO driver writes them). With ``append_to`` set, the rows are appended to that
file as well (a design-of-experiments table).
"""
from __future__ import annotations

import csv
import math
import os
import re

from .common import DriverError, Log, Settings, main_wrapper, write_csv

COEFFICIENTS = ("CL", "CD", "CY", "Cl", "Cm", "Cn")
TOL = 1e-6


# --- axes ----------------------------------------------------------------------------------------------------------

def axes(alpha_deg, beta_deg):
    """Unit freestream, lift and side-force directions in x aft, y right, z up."""
    a, b = math.radians(alpha_deg), math.radians(beta_deg)
    ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
    return (ca * cb, -sb, sa * cb), (-sa, 0.0, ca), (ca * sb, cb, sa * sb)


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def wind_to_body(cd, cy, cl, alpha, beta):
    v, lift, side = axes(alpha, beta)
    return tuple(cd * v[i] + cy * side[i] + cl * lift[i] for i in range(3))


def body_to_wind(f, alpha, beta):
    v, lift, side = axes(alpha, beta)
    return dot(f, lift), dot(f, v), dot(f, side)  # CL, CD, CY


# --- interpolation -------------------------------------------------------------------------------------------------

def _key(v):
    return round(v, 6)


class Grid:
    """Values on a Mach x alpha x beta tensor grid (some axes may hold one value), interpolated multilinearly."""

    def __init__(self, name, points):
        self.name = name
        self.data = {}
        for (m, a, b), vals in points:
            k = (_key(m), _key(a), _key(b))
            if k in self.data:
                raise DriverError(f"{name}: two sets of values at Mach {m:g}, alpha {a:g}, beta {b:g}")
            self.data[k] = vals
        self.axes = [sorted({k[i] for k in self.data}) for i in range(3)]
        self.order = [p for p, _ in points]

    def _bracket(self, i, x, what, cond):
        ax = self.axes[i]
        if len(ax) == 1:
            if abs(ax[0] - x) > TOL:
                raise DriverError(f"{self.name} has {what} {ax[0]:g} only; the database asks for {what} {x:g} ({cond})")
            return [(ax[0], 1.0)]
        if x < ax[0] - TOL or x > ax[-1] + TOL:
            raise DriverError(f"{self.name} covers {what} {ax[0]:g} to {ax[-1]:g}; the database asks for {what} {x:g} "
                              f"({cond})")
        for x0, x1 in zip(ax, ax[1:]):
            if x0 - TOL <= x <= x1 + TOL:
                t = min(1.0, max(0.0, (x - x0) / (x1 - x0)))
                return [(x0, 1.0 - t), (x1, t)]
        return [(ax[-1], 1.0)]

    def at(self, m, a, b):
        """The values at a condition it has, else interpolated between the grid points around it."""
        cond = f"Mach {m:g}, alpha {a:g}, beta {b:g}"
        for k, vals in self.data.items():
            if abs(k[0] - m) <= TOL and abs(k[1] - a) <= TOL and abs(k[2] - b) <= TOL:
                return list(vals)
        out = None
        for mm, wm in self._bracket(0, m, "Mach", cond):
            for aa, wa in self._bracket(1, a, "alpha", cond):
                for bb, wb in self._bracket(2, b, "beta", cond):
                    w = wm * wa * wb
                    if w == 0.0:
                        continue
                    if (mm, aa, bb) not in self.data:
                        raise DriverError(f"{self.name} has no values at {cond} and cannot interpolate there: it lacks "
                                          f"Mach {mm:g}, alpha {aa:g}, beta {bb:g}. Interpolation needs every "
                                          "combination of the neighbouring Mach numbers and angles; or compute it at "
                                          "the database's conditions")
                    vals = self.data[(mm, aa, bb)]
                    out = [w * v for v in vals] if out is None else [o + w * v for o, v in zip(out, vals)]
        return out


def interp_mach(name, points, m):
    pts = sorted(points)
    for x, y in pts:
        if abs(x - m) <= TOL:
            return y
    if not pts or m < pts[0][0] - TOL or m > pts[-1][0] + TOL:
        raise DriverError(f"{name} was run from Mach {pts[0][0]:g} to {pts[-1][0]:g} but applies at Mach {m:g}: run it "
                          "over the Mach numbers where it applies, or narrow its Mach range")
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= m <= x1:
            return y0 + (y1 - y0) * (m - x0) / (x1 - x0)
    return pts[-1][1]


# --- contributions -------------------------------------------------------------------------------------------------

def read_table(path):
    if not os.path.isfile(path):
        raise DriverError(f"{path} not found; run the analysis first and connect its output")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
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


class Contribution:
    """One source of forces and moments. ``body(m, a, b)`` gives (force, moment) coefficients in x aft, y right,
    z up, over its own reference, about its own moment point."""

    def __init__(self, c, log):
        self.kind = c.get("kind", "")
        self.name = c.get("name", "") or self.kind
        self.mach_min = float(c.get("mach_min", 0.0))
        self.mach_max = float(c.get("mach_max", 1e9))
        self.scale = float(c.get("scale", 1.0))
        self.use = c.get("use", "all")
        self.grid = None
        self.mach_points = None
        if self.kind == "vspaero":
            self._vspaero(c, log)
        elif self.kind in ("wave", "parasite"):
            self._drag(c, log)
        elif self.kind == "rows":
            self._rows(c, log)
        else:
            raise DriverError(f"unknown contribution kind '{self.kind}'")

    def _vspaero(self, c, log):
        base = c["results"]
        s = Settings(base + ".ntvsp.settings.txt")
        self.sref, self.bref, self.cref = s.real("sref"), s.real("bref"), s.real("cref")
        self.point = (s.real("xcg"), s.real("ycg"), s.real("zcg"))
        rows = read_table(base + ".ntvsp.polar.csv")
        res = sorted({r["Re_1e6"] for r in rows})
        if len(res) > 1:
            log(f"{self.name}: {len(res)} Re_cref values; the database uses Re_cref = {res[0] * 1e6:g}")
        rows = [r for r in rows if r["Re_1e6"] == res[0]]
        bad = [r for r in rows if r.get("Converged", 1) == 0]
        if bad:
            log(f"warning: {self.name}: {len(bad)} cases did not pass VSPAERO's convergence check")
        tot = c.get("terms", "inviscid") == "total"
        f = ("CFxtot", "CFytot", "CFztot", "CMxtot", "CMytot", "CMztot") if tot else \
            ("CFix", "CFiy", "CFiz", "CMix", "CMiy", "CMiz")
        self.grid = Grid(self.name, [((r["Mach"], r["AoA"], r["Beta"]), [r[k] for k in f]) for r in rows])
        log(f"{self.name}: VSPAERO {'total' if tot else 'inviscid'} coefficients, {len(rows)} cases from {base}")

    def _drag(self, c, log):
        base = c["results"]
        s = Settings(base + ".ntvsp.settings.txt")
        self.sref, self.bref, self.cref = s.real("sref"), 1.0, 1.0
        self.point = None  # no moment information: the drag acts through the database's moment point
        if self.kind == "wave":
            rows = read_table(base + ".ntvsp.wave.csv")
            self.mach_points = [(r["Mach"], r["CD_wave"]) for r in rows]
        else:
            rows = read_table(base + ".ntvsp.parasite.csv")
            self.mach_points = [(r["Mach"], r["Total_CD"]) for r in rows]
        ms = [m for m, _ in self.mach_points]
        log(f"{self.name}: {'Wave' if self.kind == 'wave' else 'Parasite'} Drag at Mach {min(ms):g} to {max(ms):g} "
            f"from {base}")

    def _rows(self, c, log):
        points, refs = [], set()
        for text in c.get("rows", "").split("|"):
            if not text.strip():
                continue
            w = [x.strip() for x in text.split(",")]
            if len(w) != 16 or w[0] not in ("body", "wind"):
                raise DriverError(f"{self.name}: a coefficient row must be 'body' or 'wind' and 15 numbers: {text}")
            v = [float(x) for x in w[1:]]
            m, a, b = v[0:3]
            c6 = v[3:9]
            if w[0] == "wind":  # CD, CY, CL, Cl, Cm, Cn (body signs) -> x aft, y right, z up
                force = wind_to_body(c6[0], c6[1], c6[2], a, b)
                c6 = list(force) + [-c6[3], c6[4], -c6[5]]
            points.append(((m, a, b), c6))
            refs.add(tuple(round(x, 9) for x in v[9:15]))
        if not points:
            raise DriverError(f"{self.name}: no coefficient rows")
        if len(refs) > 1:
            raise DriverError(f"{self.name}: its rows have different reference values or moment points; give each "
                              "contribution one reference")
        self.sref, self.bref, self.cref, *p = refs.pop()
        self.point = tuple(p)
        self.grid = Grid(self.name, points)
        log(f"{self.name}: {len(points)} coefficient rows")

    def applies(self, m):
        return self.mach_min - TOL <= m < self.mach_max - TOL

    def conditions(self):
        if self.grid is None:
            raise DriverError(f"{self.name} has no angles of attack of its own; give the database its flight "
                              "conditions, or put a VSPAERO or coefficient-row contribution first")
        return self.grid.order

    def body(self, m, a, b):
        if self.grid is not None:
            v = self.grid.at(m, a, b)
            force, moment = tuple(v[:3]), tuple(v[3:])
        else:
            cd = interp_mach(self.name, self.mach_points, m)
            force, moment = wind_to_body(cd, 0.0, 0.0, a, b), (0.0, 0.0, 0.0)
        vel = axes(a, b)[0]
        if self.use == "drag":
            d = dot(force, vel)
            force, moment = tuple(d * x for x in vel), (0.0, 0.0, 0.0)
        elif self.use == "nodrag":
            d = dot(force, vel)
            force = tuple(f - d * x for f, x in zip(force, vel))
        return force, moment

    def on(self, ref, m, a, b):
        """Force and moment coefficients over the database reference ``ref`` = (S, b, c, point), about its point."""
        force, moment = self.body(m, a, b)
        s, k = self.sref, self.scale
        f = tuple(k * x * s for x in force)  # per dynamic pressure
        mom = (k * moment[0] * s * self.bref, k * moment[1] * s * self.cref, k * moment[2] * s * self.bref)
        arm = tuple(p - q for p, q in zip(self.point or ref[3], ref[3]))
        mom = tuple(x + y for x, y in zip(mom, cross(arm, f)))
        S, B, C = ref[0], ref[1], ref[2]
        return tuple(x / S for x in f), (mom[0] / (S * B), mom[1] / (S * C), mom[2] / (S * B))


def parse_contributions(path):
    out = []
    for line in open(path, encoding="utf-8"):
        key, _, value = line.partition("=")
        if key.strip() != "contribution" or not value.strip():
            continue
        c = {}
        for item in value.split(";"):
            k, sep, v = item.partition(":")
            if sep:
                c[k.strip()] = v.strip()
        out.append(c)
    return out


def coefficients(force, moment, a, b):
    cl, cd, cy = body_to_wind(force, a, b)
    return [cl, cd, cy, -moment[0], moment[1], -moment[2]]


def column_name(name):
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_") or "part"


def centers(work, case, rows, ref, log):
    """Centre of pressure per row and aerodynamic centre per (Mach, beta) alpha sweep, from the totals, as the
    VSPAERO driver computes them (x_cp = x_ref - Cm cref / CN; x_ac = x_ref - dCm/dCN cref)."""
    xr, yr, zr = ref[3]
    cref = ref[2]
    cp = []
    for k, r in enumerate(rows, start=1):
        ok = abs(r["CN"]) >= 1e-6
        cp.append([k, r["Mach"], r["AoA"], r["Beta"], xr - r["Cm"] * cref / r["CN"] if ok else xr, yr, zr, int(ok),
                   r["CN"], r["Cm"]])
    write_csv(os.path.join(work, case + ".ntvsp.center_of_pressure.csv"),
              ["Case", "Mach", "AoA", "Beta", "X", "Y", "Z", "Defined", "CN", "Cm"], cp)
    groups = {}
    for r in rows:
        groups.setdefault((r["Mach"], r["Beta"]), []).append(r)
    ac = []
    for (m, b), g in groups.items():
        f, c = [r["CN"] for r in g], [r["Cm"] for r in g]
        n = len(f)
        mf, mc = sum(f) / n, sum(c) / n
        sff = sum((x - mf) ** 2 for x in f)
        if n < 2 or sff < 1e-12:
            continue
        slope = sum((x - mf) * (y - mc) for x, y in zip(f, c)) / sff
        ac.append([m, b, xr - slope * cref, yr, zr, slope])
    write_csv(os.path.join(work, case + ".ntvsp.aerodynamic_center.csv"),
              ["Mach", "Beta", "X", "Y", "Z", "dCm_dCN"], ac)
    if ac:
        log("aerodynamic centre x: " + ", ".join(f"Mach {r[0]:g} beta {r[1]:g}: {r[2]:.3f} m" for r in ac))


def main(settings_path):
    s = Settings(settings_path)
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    parts = [Contribution(c, log) for c in parse_contributions(settings_path)]
    if not parts:
        raise DriverError("no contributions connected")
    names, seen = [], {}
    for p in parts:
        n = column_name(p.name)
        seen[n] = seen.get(n, 0) + 1
        names.append(n if seen[n] == 1 else f"{n}_{seen[n]}")

    first = parts[0]
    if s.integer("reference_from_first", 1):
        if first.grid is None:
            raise DriverError(f"{first.name} has no reference chord and span; put a VSPAERO or coefficient-row "
                              "contribution first, or give the database its own reference")
        ref = (first.sref, first.bref, first.cref, first.point)
    else:
        ref = (s.real("sref"), s.real("bref"), s.real("cref"), (s.real("xref"), s.real("yref"), s.real("zref")))
    if s.integer("conditions_from_first", 1):
        conds = first.conditions()
    else:
        conds = [(m, a, b) for m in s.reals("mach") for a in s.reals("alpha") for b in s.reals("beta")]
    if not conds:
        raise DriverError("the database has no flight conditions")
    log(f"database: {len(conds)} conditions, Sref {ref[0]:g} m^2, bref {ref[1]:g} m, cref {ref[2]:g} m, moments "
        f"about ({ref[3][0]:g}, {ref[3][1]:g}, {ref[3][2]:g}) m")
    for m in sorted({c[0] for c in conds}):
        active = [p.name for p in parts if p.applies(m)]
        log(f"Mach {m:g}: " + (", ".join(active) if active else "no contribution applies"))

    variables = [x.strip() for x in s.text("variables", "").split(";") if x.strip()]
    values = s.reals("values") if variables else []
    if len(variables) != len(values):
        raise DriverError(f"{len(variables)} variable names but {len(values)} values")
    header = variables + ["Mach", "AoA", "Beta"] + list(COEFFICIENTS) + ["L_over_D", "CA", "CN"] + \
        [f"{c}_{n}" for n in names for c in COEFFICIENTS]
    rows, table = [], []
    for m, a, b in conds:
        tf, tm, own = [0.0] * 3, [0.0] * 3, []
        for p in parts:
            if p.applies(m):
                f, mo = p.on(ref, m, a, b)
                tf = [x + y for x, y in zip(tf, f)]
                tm = [x + y for x, y in zip(tm, mo)]
                own += coefficients(f, mo, a, b)
            else:
                own += [0.0] * 6
        tot = coefficients(tf, tm, a, b)
        r = dict(zip(COEFFICIENTS, tot), Mach=m, AoA=a, Beta=b, CA=tf[0], CN=tf[2])
        rows.append(r)
        table.append(values + [m, a, b] + tot + [tot[0] / tot[1] if tot[1] else 0.0, tf[0], tf[2]] + own)
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
