"""Modified Newtonian (VSP), Modified Newtonian on Mesh (VSP) and Tangent Wedge on Mesh (VSP) against analytic
solutions, through ntopcl (see conftest.py)."""
from __future__ import annotations

import csv
import math

import pytest

from conftest import ALPHAS, MACHS, PLATE, SPHERE_R, SREF, USER_CP_MAX

GAMMA = 1.4
pytestmark = pytest.mark.e2e


# --- analytic solutions

def cp_max_rayleigh(mach, g=GAMMA):
    """Stagnation Cp behind a normal shock (Rayleigh pitot formula)."""
    m2 = mach * mach
    p02 = ((g + 1) ** 2 * m2 / (4 * g * m2 - 2 * (g - 1))) ** (g / (g - 1)) * (1 - g + 2 * g * m2) / (g + 1)
    return 2.0 / (g * m2) * (p02 - 1.0)


def oblique_shock_cp(mach, theta, g=GAMMA):
    """Cp behind the weak oblique shock turning the flow through theta: theta-beta-M solved by bisection."""
    if theta <= 0.0:
        return 0.0

    def deflection(beta):
        return math.atan(2.0 / math.tan(beta) * (mach**2 * math.sin(beta) ** 2 - 1.0)
                         / (mach**2 * (g + math.cos(2.0 * beta)) + 2.0))

    lo, hi = math.asin(1.0 / mach), 0.5 * math.pi
    for _ in range(200):  # beta at the largest deflection
        a, b = lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0
        lo, hi = (a, hi) if deflection(a) < deflection(b) else (lo, b)
    lo, hi = math.asin(1.0 / mach), 0.5 * (lo + hi)
    assert deflection(hi) > theta, "detached shock"
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if deflection(mid) < theta else (lo, mid)
    beta = 0.5 * (lo + hi)
    p21 = 1.0 + 2.0 * g / (g + 1.0) * (mach**2 * math.sin(beta) ** 2 - 1.0)
    return 2.0 / (g * mach**2) * (p21 - 1.0)


def prandtl_meyer_cp(mach, turn, g=GAMMA):
    """Cp after an isentropic expansion through ``turn`` radians from Mach ``mach``."""
    if turn <= 0.0:
        return 0.0
    k = math.sqrt((g + 1.0) / (g - 1.0))

    def nu(m):
        r = math.sqrt(m * m - 1.0)
        return k * math.atan(r / k) - math.atan(r)

    target = nu(mach) + turn
    lo, hi = mach, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if nu(mid) < target else (lo, mid)
    m2 = math.sqrt(lo * hi)
    ratio = ((1.0 + 0.5 * (g - 1.0) * mach**2) / (1.0 + 0.5 * (g - 1.0) * m2**2)) ** (g / (g - 1.0))
    return 2.0 / (g * mach**2) * (ratio - 1.0)


def wind(row):
    """(CD, CL) of a row from its body-axis CA, CN (zero sideslip)."""
    a = math.radians(row["AoA"])
    return (row["CA"] * math.cos(a) + row["CN"] * math.sin(a), row["CN"] * math.cos(a) - row["CA"] * math.sin(a))


def by_condition(rows):
    got = {(r["Mach"], r["AoA"]) for r in rows}
    want = {(m, a) for m in MACHS for a in ALPHAS}
    assert got == want, f"conditions {sorted(got)}, expected {sorted(want)}"
    return rows


SPHERE_AREA = math.pi * SPHERE_R**2 / SREF
CHORD, SPAN, THICK = PLATE
PLANFORM = CHORD * SPAN / SREF
EDGE = THICK * SPAN / SREF


# --- tests

@pytest.mark.parametrize("table, cp_max", [
    ("si_sphere_classical", lambda m: 2.0),
    ("mesh_sphere_modified", cp_max_rayleigh),
    ("mesh_sphere_user", lambda m: USER_CP_MAX),
], indirect=["table"])
def test_sphere_drag(table, cp_max):
    for r in by_condition(table):
        cd, cl = wind(r)
        expected = 0.5 * cp_max(r["Mach"]) * SPHERE_AREA
        assert cd == pytest.approx(expected, rel=5e-3), f"Mach {r['Mach']} alpha {r['AoA']}"
        assert abs(cl) < 2e-3 * expected and abs(r["CY"]) < 2e-3 * expected


@pytest.mark.parametrize("table", ["mesh_plate_classical"], indirect=True)
def test_plate_newtonian(table):
    for r in by_condition(table):
        a = math.radians(r["AoA"])
        where = f"Mach {r['Mach']} alpha {r['AoA']}"
        assert r["CN"] == pytest.approx(2.0 * math.sin(a) ** 2 * PLANFORM, rel=1e-2, abs=1e-4), where
        assert r["CA"] == pytest.approx(2.0 * math.cos(a) ** 2 * EDGE, rel=1e-2, abs=1e-4), where


@pytest.mark.parametrize("table", ["tw_plate"], indirect=True)
def test_plate_tangent_wedge(table):
    for r in by_condition(table):
        a = math.radians(r["AoA"])
        expected = (oblique_shock_cp(r["Mach"], a) - prandtl_meyer_cp(r["Mach"], a)) * PLANFORM
        assert r["CN"] == pytest.approx(expected, rel=1e-2, abs=1e-4), f"Mach {r['Mach']} alpha {r['AoA']}"


def test_cp_max_file_only_for_user_defined(results):
    files = results["_files"]
    assert "mesh_sphere_user.ntvsp.mesh_cp_max.csv" in files
    for case in ("mesh_sphere_modified", "mesh_plate_classical", "tw_plate"):
        assert f"{case}.ntvsp.coefficients.csv" in files, case
        assert f"{case}.ntvsp.mesh_cp_max.csv" not in files, case


@pytest.mark.parametrize("case", ["si_sphere_classical", "mesh_sphere_modified", "mesh_sphere_user",
                                  "mesh_plate_classical", "tw_plate", "tw_cone", "mesh_sphere_export"])
def test_aero_database_matches_table(results, case):
    """Coefficients Contribution -> Aero Database on the same flight condition gives back the table's CL and CD."""
    for key in (case, "db_" + case):
        if key not in results:
            pytest.fail(f"{key}: no table written")
    table = {(r["Mach"], r["AoA"]): wind(r) for r in by_condition(results[case])}
    rows = results["db_" + case]
    assert len(rows) == len(table)
    for r in rows:
        cd, cl = table[(r["Mach"], r["AoA"])]
        where = f"Mach {r['Mach']} alpha {r['AoA']}"
        assert r["CD"] == pytest.approx(cd, rel=1e-6, abs=1e-9), where
        assert r["CL"] == pytest.approx(cl, rel=1e-6, abs=1e-9), where


def test_aero_database_sphere_drag(results):
    rows = results.get("db_mesh_sphere_modified")
    assert rows, "no database written"
    for r in rows:
        assert r["CD"] == pytest.approx(0.5 * cp_max_rayleigh(r["Mach"]) * SPHERE_AREA, rel=5e-3)


# --- tangent cone

def taylor_maccoll_cone_cp(mach, cone_deg, g=GAMMA, step=2e-5):
    """Surface Cp of the sharp cone of half angle ``cone_deg`` at zero incidence: the Taylor-Maccoll equation
    integrated in theta (RK4) from a trial conical shock to where V_theta vanishes, the shock angle found by
    bisection."""
    target = math.radians(cone_deg)

    def surface(beta):
        mn = mach * math.sin(beta)
        delta = math.atan(2.0 / math.tan(beta) * (mn**2 - 1.0) / (mach**2 * (g + math.cos(2.0 * beta)) + 2.0))
        mn2 = math.sqrt((1.0 + 0.5 * (g - 1.0) * mn**2) / (g * mn**2 - 0.5 * (g - 1.0)))
        m2 = mn2 / math.sin(beta - delta)
        v = (2.0 / ((g - 1.0) * m2**2) + 1.0) ** -0.5
        y = (v * math.cos(beta - delta), -v * math.sin(beta - delta))  # (V_r, V_theta)

        def f(theta, y):
            vr, vt = y
            a = 0.5 * (g - 1.0) * (1.0 - vr * vr - vt * vt)
            return vt, (vt * vt * vr - a * (2.0 * vr + vt / math.tan(theta))) / (a - vt * vt)

        theta, h = beta, -step
        while True:
            k1 = f(theta, y)
            k2 = f(theta + h / 2, tuple(a + h / 2 * b for a, b in zip(y, k1)))
            k3 = f(theta + h / 2, tuple(a + h / 2 * b for a, b in zip(y, k2)))
            k4 = f(theta + h, tuple(a + h * b for a, b in zip(y, k3)))
            new = tuple(a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4))
            if new[1] >= 0.0:  # V_theta reaches zero at the surface
                t = y[1] / (y[1] - new[1])
                vr = y[0] + t * (new[0] - y[0])
                cone = theta + t * h
                break
            y, theta = new, theta + h
        m_c2 = 2.0 / (g - 1.0) * vr * vr / (1.0 - vr * vr)
        p0 = lambda m2: (1.0 + 0.5 * (g - 1.0) * m2) ** (g / (g - 1.0))  # noqa: E731
        total = ((g + 1) * mn**2 / ((g - 1) * mn**2 + 2)) ** (g / (g - 1)) * ((g + 1) / (2 * g * mn**2 - (g - 1))) ** (1 / (g - 1))
        return cone, 2.0 / (g * mach**2) * (p0(mach**2) * total / p0(m_c2) - 1.0)

    lo, hi = math.asin(1.0 / mach) + 1e-6, target + math.radians(25.0)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if surface(mid)[0] < target else (lo, mid)
    return surface(0.5 * (lo + hi))[1]


@pytest.mark.parametrize("table", ["tw_cone"], indirect=True)
def test_tangent_cone_on_a_cone(table):
    """A sharp cone at zero incidence: the tangent cone gives the exact cone pressure on every face, the base none."""
    from conftest import CONE_HALF_ANGLE, CONE_LENGTH
    base = math.pi * (CONE_LENGTH * math.tan(math.radians(CONE_HALF_ANGLE))) ** 2 / SREF
    rows = [r for r in by_condition(table) if r["AoA"] == 0.0]
    for r in rows:
        expected = taylor_maccoll_cone_cp(r["Mach"], CONE_HALF_ANGLE) * base
        assert r["CA"] == pytest.approx(expected, rel=1e-2), f"Mach {r['Mach']}"
        assert abs(r["CN"]) < 1e-3 * expected


# --- pressure exports and readers

def _flight_rows():
    from conftest import WORK
    with open(f"{WORK}/flight.ntvsp.flight.csv", newline="") as f:
        return {int(float(r["Case"])): {k: float(v) for k, v in r.items()} for r in csv.DictReader(f)}


def _point_map(path):
    with open(path) as f:
        return [tuple(float(v) for v in line.split(",")) for line in f if line.strip()]


def test_pressure_export_newtonian_sphere(results):
    """Modified Newtonian on Mesh with Export Pressures: every face centroid of every flight condition case carries
    p = p_inf + q Cp_max max(-n.V, 0)^2, n the sphere's radial normal."""
    from conftest import WORK
    flight = _flight_rows()
    assert set(flight) == set(range(1, len(MACHS) * len(ALPHAS) + 1))
    for case, row in flight.items():
        name = f"mesh_sphere_export.ntvsp.pressure.{case}.csv"
        assert name in results["_files"], name
        a = math.radians(row["AoA"])
        v = (math.cos(a), 0.0, math.sin(a))
        cp_max = cp_max_rayleigh(row["Mach"])
        worst = 0.0
        for x, y, z, p in _point_map(f"{WORK}/{name}"):
            r = math.sqrt(x * x + y * y + z * z)
            sine = max(-(x * v[0] + y * v[1] + z * v[2]) / r, 0.0)
            worst = max(worst, abs(p - (row["Pressure"] + row["Dynamic_Pressure"] * cp_max * sine * sine)))
        assert worst < 1e-5 * row["Dynamic_Pressure"], f"case {case}: off by {worst:.3g} Pa"


def test_pressure_files_only_when_exported(results):
    exported = {f.split(".ntvsp.")[0] for f in results["_files"] if ".ntvsp.pressure" in f}
    assert exported == {"mesh_sphere_export", "fake"}


def _reads(results):
    from conftest import CASES
    out = results["_json"]
    assert out, "ntopcl wrote no output"
    comps = out[0]["components"][len(CASES):]

    def plain(c):
        v = c["value"]
        return [p["val"] for p in v] if isinstance(v, list) else v["val"]
    return [[plain(c) for c in group["components"]] if group["type"] == "group" else
            [x["val"] for x in group["value"]["value"]] for group in comps]


def test_surface_pressure_reader(results):
    from conftest import READ_CASE, READ_ROWS, WORK
    rows = _point_map(f"{WORK}/{READ_CASE[0]}.ntvsp.pressure.{READ_CASE[1]}.csv")
    got = _reads(results)[0]
    n = len(READ_ROWS)
    for k, i in enumerate(READ_ROWS):
        assert got[k] == pytest.approx(rows[i][3], rel=1e-8)
        assert got[n + k] == pytest.approx(list(rows[i][:3]), abs=1e-9)


def test_surface_pressure_jump_reader(results):
    from conftest import FAKE_JUMP, REFERENCE_PRESSURE
    got = _reads(results)[1]
    n = len(FAKE_JUMP)
    for k, (x, y, z, dp) in enumerate(FAKE_JUMP):
        assert got[k] == pytest.approx(REFERENCE_PRESSURE + dp)
        assert got[n + k] == pytest.approx([x, y, z], abs=1e-12)


def test_newtonian_pressure_field(results):
    """Modified Newtonian Pressure (VSP), Classical Cp_max, at points on the sphere."""
    from conftest import FIELD_CASE, FIELD_POINTS_MM
    row = _flight_rows()[FIELD_CASE]
    a = math.radians(row["AoA"])
    got = _reads(results)[2]
    for (x, y, z), p in zip(FIELD_POINTS_MM, got):
        r = math.sqrt(x * x + y * y + z * z)
        sine = max(-(x * math.cos(a) + z * math.sin(a)) / r, 0.0)
        assert p == pytest.approx(row["Pressure"] + row["Dynamic_Pressure"] * 2.0 * sine * sine, rel=1e-6)
