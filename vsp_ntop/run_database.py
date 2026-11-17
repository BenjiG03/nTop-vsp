"""Aero database: VSPAERO's forces and moments plus wave drag plus parasite drag, as one table, driven from nTop.

    python -m vsp_ntop.run_database SETTINGS_FILE

OpenVSP has no tool that combines VSPAERO, Wave Drag and Parasite Drag, so this driver does. Per VSPAERO case
(Mach, alpha, beta, Re_cref):

    CL = CLi + CL_viscous                     CY = CSi
    CD = CDi + CD_wave(M) + CD_parasite(M) + CD_viscous
    Cl = -CMix       Cm = CMiy + Cm_viscous   Cn = -CMiz       (about VSPAERO's moment reference point)

Moments are given in aircraft body-axis signs: Cl positive right wing down, Cm positive nose up, Cn positive nose
right; CY is positive to the right. VSPAERO's own axes are x aft, y right, z up, so its rolling and yawing moments
have the opposite signs (a fin aft gives VSPAERO a negative CMz per degree of sideslip).

CLi ... CMiz are VSPAERO's inviscid coefficients (VSPAERO's own viscous estimate, the ``o`` columns, is left out:
Parasite Drag replaces it). CD_wave and CD_parasite are interpolated linearly in Mach between the Mach numbers of
the Wave Drag and Parasite Drag runs; CD_wave is 0 below Mach 1. Optional viscous lists (one value per case or
one for all, e.g. from a viscous correction) are added as given.

Checks: the wave and parasite runs must use VSPAERO's Sref, and must cover every VSPAERO Mach number (wave drag:
every supersonic one).

Writes ``<case>.ntvsp.database.csv``. With ``append_to`` set, the rows are also appended to that file (its header
written when it is new), so repeated runs with different geometry inputs build one design-of-experiments table;
``variables`` (names) and ``values`` become its first columns.
"""
from __future__ import annotations

import csv
import os

from .common import DriverError, Log, Settings, main_wrapper, write_csv

COEFFICIENTS = ("CL", "CD", "CY", "Cl", "Cm", "Cn")
PARTS = ("L_over_D", "CL_inviscid", "CD_induced", "CD_wave", "CD_parasite", "CD_viscous", "CL_viscous", "Cm_viscous")
REL_TOL = 1e-6


def read_table(path):
    if not os.path.isfile(path):
        raise DriverError(f"{path} not found; run the analysis first and connect its output")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        raise DriverError(f"{path} has no results")
    out = []
    for r in rows[1:]:
        if not r:
            continue
        row = {}
        for h, v in zip(rows[0], r):
            try:
                row[h] = float(v)
            except ValueError:
                row[h] = v
        out.append(row)
    return out


def same(a, b):
    return abs(a - b) <= REL_TOL * max(1.0, abs(a), abs(b))


def interpolate(points, m, what):
    """Linear interpolation in Mach of [(mach, value)]; exact at the Mach numbers given."""
    pts = sorted(points)
    for x, y in pts:
        if abs(x - m) <= 1e-6:
            return y
    if not pts or m < pts[0][0] - 1e-6 or m > pts[-1][0] + 1e-6:
        raise DriverError(f"{what} was run from Mach {pts[0][0]:g} to {pts[-1][0]:g}; VSPAERO has a case at Mach "
                          f"{m:g}: run {what} over VSPAERO's Mach range")
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= m <= x1:
            return y0 + (y1 - y0) * (m - x0) / (x1 - x0)
    return pts[-1][1]


def per_case(name, values, n):
    if len(values) == 1:
        return values * n
    if len(values) != n:
        raise DriverError(f"{name}: {len(values)} values for {n} VSPAERO cases (one per case, or one for all)")
    return values


def main(settings_path):
    s = Settings(settings_path)
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))

    aero = s.text("vspaero_results")
    polar = read_table(aero + ".ntvsp.polar.csv")
    sref = Settings(aero + ".ntvsp.settings.txt").real("sref")
    n = len(polar)
    log(f"VSPAERO: {n} cases from {aero}, Sref = {sref:g} m^2")

    wave_base = s.text("wave_results", "").strip()
    wave = []
    if wave_base:
        rows = read_table(wave_base + ".ntvsp.wave.csv")
        for r in rows:
            if not same(r["Sref"], sref):
                raise DriverError(f"Wave Drag used Sref = {r['Sref']:g} m^2 but VSPAERO used {sref:g} m^2")
        wave = [(r["Mach"], r["CD_wave"]) for r in rows]
        log(f"Wave Drag: {len(wave)} Mach numbers from {min(m for m, _ in wave):g} to {max(m for m, _ in wave):g}")
    par_base = s.text("parasite_results", "").strip()
    parasite = []
    if par_base:
        rows = read_table(par_base + ".ntvsp.parasite.csv")
        for r in rows:
            if not same(r["Sref"], sref):
                raise DriverError(f"Parasite Drag used Sref = {r['Sref']:g} but VSPAERO used {sref:g}; give both the "
                                  "same reference area (Parasite Drag's is in its Model Length Unit squared)")
        parasite = [(r["Mach"], r["Total_CD"]) for r in rows]
        log(f"Parasite Drag: {len(parasite)} Mach numbers from {min(m for m, _ in parasite):g} to "
            f"{max(m for m, _ in parasite):g}")

    cd_v = per_case("Viscous CD", s.reals("viscous_cd") or [0.0], n)
    cl_v = per_case("Viscous CL", s.reals("viscous_cl") or [0.0], n)
    cm_v = per_case("Viscous Cm", s.reals("viscous_cmy") or [0.0], n)
    names = [x.strip() for x in s.text("variables", "").split(";") if x.strip()]
    values = s.reals("values") if names else []
    if len(names) != len(values):
        raise DriverError(f"{len(names)} variable names but {len(values)} values")

    header = names + ["Mach", "AoA", "Beta", "Re_cref"] + list(COEFFICIENTS) + list(PARTS)
    rows = []
    for k, row in enumerate(polar):
        m = row["Mach"]
        cdw = interpolate(wave, m, "Wave Drag") if (wave and m >= 1.0) else 0.0
        if m >= 1.0 and not wave:
            log(f"warning: case {k + 1} is supersonic (Mach {m:g}) but no Wave Drag results are connected")
        cdp = interpolate(parasite, m, "Parasite Drag") if parasite else 0.0
        cl = row["CLi"] + cl_v[k]
        cd = row["CDi"] + cdw + cdp + cd_v[k]
        cm = row["CMiy"] + cm_v[k]
        rows.append(values + [m, row["AoA"], row["Beta"], row["Re_1e6"] * 1e6,
                              cl, cd, row["CSi"], -row["CMix"], cm, -row["CMiz"],
                              cl / cd if cd else 0.0, row["CLi"], row["CDi"], cdw, cdp, cd_v[k], cl_v[k], cm_v[k]])
    write_csv(os.path.join(work, case + ".ntvsp.database.csv"), header, rows)
    log(f"database: {len(rows)} rows, columns {', '.join(header)}")
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
            w.writerows(rows)
        log(f"appended {len(rows)} rows to {target}")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
