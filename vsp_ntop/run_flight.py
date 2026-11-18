"""Flight conditions from Mach number and altitude, driven from nTop.

    python -m vsp_ntop.run_flight SETTINGS_FILE

Settings: ``mach`` and ``altitude`` (``;``-separated, paired item by item into flight points; one altitude serves
every Mach number), ``alpha`` and ``beta`` (degrees), the reference ``sref``, ``bref``, ``cref`` and moment point
``xref``, ``yref``, ``zref`` (metres).

The atmosphere at each altitude (geometric) is the US Standard Atmosphere 1976 as implemented by the ``ambiance``
package: temperature, pressure, density, speed of sound and dynamic viscosity; the velocity, dynamic pressure and
Reynolds numbers follow from the Mach number. The conditions are every flight point x sideslip x angle of attack.
Files: see ``vsp_ntop.flight``.
"""
from __future__ import annotations

import os

from .common import DriverError, Log, Settings, main_wrapper, write_csv
from .flight import CONDITION_COLUMNS, FEET, INVISCID_COLUMNS, POINT_COLUMNS, same

ALTITUDE_RANGE = (-5004.0, 81020.0)  # ambiance's range, metres


def atmosphere(altitudes):
    try:
        from ambiance import Atmosphere
    except ImportError:
        raise DriverError("the ambiance package is missing from the Python of NTOP_PYTHON_PATH: pip install -e . in the "
                          "nTop_vsp folder (or pip install ambiance)")
    a = Atmosphere(list(altitudes))
    return [dict(Temperature=float(a.temperature[i]), Pressure=float(a.pressure[i]), Density=float(a.density[i]),
                 Speed_of_Sound=float(a.speed_of_sound[i]), Dynamic_Viscosity=float(a.dynamic_viscosity[i]))
            for i in range(len(altitudes))]


def main(settings_path):
    s = Settings(settings_path)
    work, case = os.path.abspath(s.text("work_dir")), s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    machs, alts = s.reals("mach"), s.reals("altitude")
    alphas, betas = s.reals("alpha"), s.reals("beta")
    if not machs:
        raise DriverError("no Mach numbers")
    if len(alts) == 1:
        alts = alts * len(machs)
    if len(alts) != len(machs):
        raise DriverError(f"{len(machs)} Mach numbers but {len(alts)} altitudes: give one altitude per Mach number (each "
                          "pair is a flight point), or one altitude for all")
    if any(m <= 0.0 for m in machs):
        raise DriverError("Mach numbers must be above 0 (the speed and Reynolds number come from them)")
    lo, hi = ALTITUDE_RANGE
    bad = [h for h in alts if not lo <= h <= hi]
    if bad:
        raise DriverError(f"altitude {bad[0]:g} m is outside the standard atmosphere ({lo:g} to {hi:g} m)")
    if not alphas or not betas:
        raise DriverError("give at least one angle of attack and one sideslip")
    sref, bref, cref = s.real("sref"), s.real("bref"), s.real("cref")
    point_ref = (s.real("xref"), s.real("yref"), s.real("zref"))
    ref = [sref, bref, cref, *point_ref]

    points = []
    for k, (m, h, atm) in enumerate(zip(machs, alts, atmosphere(alts)), start=1):
        v = m * atm["Speed_of_Sound"]
        rpm = atm["Density"] * v / atm["Dynamic_Viscosity"]
        points.append(dict(Point=k, Mach=m, Altitude=h, Altitude_ft=h / FEET, Velocity=v,
                           Dynamic_Pressure=0.5 * atm["Density"] * v * v, Re_per_m=rpm, Re_cref=rpm * cref, **atm))
    dup = [(a, b) for i, a in enumerate(points) for b in points[i + 1:]
           if same(a["Mach"], b["Mach"]) and same(a["Altitude"], b["Altitude"])]
    if dup:
        raise DriverError(f"flight point Mach {dup[0][0]['Mach']:g} at {dup[0][0]['Altitude']:g} m is given twice")
    write_csv(os.path.join(work, case + ".ntvsp.points.csv"), POINT_COLUMNS, [[p[c] for c in POINT_COLUMNS] for p in points])

    rows, inviscid = [], []
    for p in points:
        for b in betas:
            for a in alphas:
                rows.append([len(rows) + 1, p["Point"], p["Mach"], p["Altitude"], a, b] +
                            [p[c] for c in ("Temperature", "Pressure", "Density", "Velocity", "Dynamic_Pressure",
                                            "Re_cref")] + ref)
                if not any(same(r[0], p["Mach"]) and same(r[1], a) and same(r[2], b) for r in inviscid):
                    inviscid.append([p["Mach"], a, b] + ref)
    write_csv(os.path.join(work, case + ".ntvsp.flight.csv"), CONDITION_COLUMNS, rows)
    write_csv(os.path.join(work, case + ".ntvsp.inviscid.csv"), INVISCID_COLUMNS, inviscid)
    for p in points:
        log(f"point {p['Point']}: Mach {p['Mach']:g} at {p['Altitude']:g} m ({p['Altitude_ft']:.0f} ft): "
            f"T {p['Temperature']:.2f} K, p {p['Pressure']:.1f} Pa, rho {p['Density']:.5g} kg/m^3, "
            f"V {p['Velocity']:.2f} m/s, q {p['Dynamic_Pressure']:.1f} Pa, Re/m {p['Re_per_m']:.4g}")
    log(f"{len(rows)} conditions ({len(points)} points x {len(betas)} sideslips x {len(alphas)} angles of attack); "
        f"Sref {sref:g} m^2, bref {bref:g} m, cref {cref:g} m, moments about ({point_ref[0]:g}, {point_ref[1]:g}, "
        f"{point_ref[2]:g}) m")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
