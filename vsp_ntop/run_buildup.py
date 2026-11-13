"""Drag buildup: VSPAERO's inviscid coefficients plus wave drag plus a viscous correction, as one table.

    python run_buildup.py SETTINGS_FILE

OpenVSP has no tool that combines these (VSPAERO, Wave Drag and Parasite Drag each report on their own), so this
driver does the sums. Per VSPAERO case:

    CL  = CLi  + CL_viscous
    CD  = CDi  + CD_wave + CD_viscous
    CMy = CMiy + CMy_viscous

CLi, CDi and CMiy are VSPAERO's inviscid parts (the ``i`` columns of its polar). VSPAERO's own viscous estimate
(the ``o`` columns, CDo) is left out, since the viscous correction replaces it. CS, CMx and CMz are VSPAERO's
inviscid values with nothing added.

The three analyses must share one flight condition and one reference: every check below stops the run when they
do not.

* Wave drag: its Mach number must equal every VSPAERO case's Mach number, and its Sref VSPAERO's Sref.
* Viscous lists: one value per VSPAERO case (in case order), or a single value used for every case. They must have
  been computed at the same Mach number, angles, Sref, cref (for CMy) and moment reference point; the driver cannot
  see those, so the notebook has to wire them from the same inputs.

Writes ``<case>.ntvsp.buildup.csv``.
"""
from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DriverError, Log, Settings, main_wrapper, write_csv  # noqa: E402

COLUMNS = ("Case", "Mach", "AoA", "Beta",
           "CL_inviscid", "CL_viscous", "CL",
           "CD_inviscid", "CD_wave", "CD_viscous", "CD",
           "CMy_inviscid", "CMy_viscous", "CMy",
           "L_over_D", "CS", "CMx", "CMz")
REL_TOL = 1e-6


def read_table(path):
    if not os.path.isfile(path):
        raise DriverError(f"{path} not found; run the analysis first and connect its Results output")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        raise DriverError(f"{path} has no results")
    header = rows[0]
    return [{h: float(v) for h, v in zip(header, r)} for r in rows[1:] if r]


def same(a, b):
    return abs(a - b) <= REL_TOL * max(1.0, abs(a), abs(b))


def per_case(name, values, n):
    if len(values) == 1:
        return values * n
    if len(values) != n:
        raise DriverError(f"{name}: {len(values)} values for {n} VSPAERO cases (give one per case, in case order, "
                          f"or a single value for all)")
    return values


def main(settings_path):
    s = Settings(settings_path)
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))

    aero = s.text("vspaero_results")
    polar = read_table(aero + ".ntvsp.polar.csv")
    aero_settings = Settings(aero + ".ntvsp.settings.txt")
    sref = aero_settings.real("sref")
    n = len(polar)
    log(f"VSPAERO: {n} cases from {aero}, Sref = {sref:g} m^2")

    wave_base = s.text("wave_results", "").strip()
    cd_wave = 0.0
    if wave_base:
        wave = read_table(wave_base + ".ntvsp.wave.csv")[0]
        cd_wave = wave["CD_wave"]
        if not same(wave["Sref"], sref):
            raise DriverError(f"Wave Drag used Sref = {wave['Sref']:g} m^2 but VSPAERO used {sref:g} m^2")
        machs = sorted({row["Mach"] for row in polar})
        if any(not same(m, wave["Mach"]) for m in machs):
            raise DriverError(f"Wave Drag ran at Mach {wave['Mach']:g} but VSPAERO cases are at Mach "
                              + ", ".join(f"{m:g}" for m in machs) + "; run both at one shared Mach number")
        log(f"Wave Drag: CD_wave = {cd_wave:.6g} at Mach {wave['Mach']:g}")
    else:
        if any(row["Mach"] > 1.0 for row in polar):
            log("warning: supersonic VSPAERO cases but no Wave Drag results; CD_wave is 0")
        log("Wave Drag: none")

    cd_v = per_case("Viscous CD", s.reals("viscous_cd") or [0.0], n)
    cl_v = per_case("Viscous CL", s.reals("viscous_cl") or [0.0], n)
    cm_v = per_case("Viscous CMy", s.reals("viscous_cmy") or [0.0], n)

    rows = []
    for k, row in enumerate(polar):
        cl = row["CLi"] + cl_v[k]
        cd = row["CDi"] + cd_wave + cd_v[k]
        cm = row["CMiy"] + cm_v[k]
        rows.append([int(row["Case"]), row["Mach"], row["AoA"], row["Beta"],
                     row["CLi"], cl_v[k], cl,
                     row["CDi"], cd_wave, cd_v[k], cd,
                     row["CMiy"], cm_v[k], cm,
                     cl / cd if cd != 0.0 else 0.0, row["CSi"], row["CMix"], row["CMiz"]])
        log(f"case {k + 1}: AoA {row['AoA']:g}  CL {cl:.5f}  CD {cd:.5f} (inviscid {row['CDi']:.5f}, wave "
            f"{cd_wave:.5f}, viscous {cd_v[k]:.5f})  CMy {cm:.5f}")
    write_csv(os.path.join(work, case + ".ntvsp.buildup.csv"), COLUMNS, rows)
    log("done")


if __name__ == "__main__":
    main_wrapper(main)
