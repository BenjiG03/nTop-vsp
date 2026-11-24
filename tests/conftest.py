"""End-to-end tests: build one notebook of nTop_vsp blocks on bodies with analytic answers, run it with ntopcl once
per session, and read the coefficient tables the blocks write.

    python -m pytest tests            (about 2 minutes)
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
import time

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "build"))

NTOPCL = os.environ.get("NTOPCL", r"C:\Program Files\nTopology\nTopology\ntopcl.exe")
WORK = os.path.join(ROOT, "runs", "tests").replace("\\", "/")
LENGTH = {"length": 1}
AREA = {"length": 2}
ANGLE = {"angle": 1}
DEG = math.radians(1.0)

# The flight condition every case shares.
MACHS = (3.0, 8.0)
ALPHAS = (0.0, 10.0, 30.0)
SREF = 0.01     # m^2, the plate's planform
CREF = 0.1      # m
# Bodies, small enough for Mesh from Implicit Body at MESH_TOL.
SPHERE_R = 0.05  # m
PLATE = (0.1, 0.1, 0.002)  # m: chord (x), span (y), thickness (z)
CONE_HALF_ANGLE = 15.0  # degrees.
CONE_LENGTH = 0.1  # m
MESH_TOL = 0.5e-3  # m (finer tolerances on these bodies take many GB in Mesh from Implicit Body)
COARSE_TOL = 2e-3  # m, for the pressure export (fewer faces, less to write)
USER_CP_MAX = 1.5


def _lst(g, typ, items):
    return g.block(f"core.list<{typ}>", f"list<{typ}>", items, top=False)


def _call(g, block, given, name):
    inputs = []
    for spec in block.inputs:
        v = given.get(spec.name, spec.default)
        inputs.append(v(g) if callable(v) else v)
    return g.block(block.signature, block.output_type, inputs, name)


def _mesh(g, body, name, tol=MESH_TOL):
    from ntopfile import NONE
    return g.block("implicit_to_mesh<implicit,real,real,bool,bool>[2.5.0]", "mesh",
                   [body, g.real(tol, LENGTH), NONE, g.boolean(True), g.boolean(False)], name, top=False)


def _peek(g, point_map, rows, name):
    """A group of the point map's values and points at ``rows`` (ntopcl does not write point maps to its output)."""
    from ntopfile import Ref

    def at(prop, typ, i):
        return g.block("core.list_element<list_interface,integer>", typ, [Ref(point_map, (prop,)), g.integer(i)],
                       top=False)
    items = [at("values", "real", i) for i in rows] + [at("points", "point", i) for i in rows]
    return g.block("core.group<any>", "group", items, name)


# (case name, block kind, body, Stagnation Pressure index, extra choice inputs: name -> index)
CASES = (
    ("si_sphere_classical", "si", "sphere", 1, {}),
    ("mesh_sphere_modified", "mesh", "sphere", 0, {}),
    ("mesh_sphere_export", "mesh", "sphere_coarse", 0, {"Export Pressures": 1}),
    ("mesh_sphere_user", "mesh", "sphere", 2, {}),
    ("mesh_plate_classical", "mesh", "plate", 1, {}),
    ("tw_plate", "tw", "plate", 0, {}),
    ("tw_cone", "tw", "cone", 0, {"Impact Method": 1}),
)
# Readers checked through the notebook's output (in this order)
READ_CASE = ("mesh_sphere_export", 2)
READ_ROWS = (0, 1, 2)
FAKE_JUMP = ((0.0, 0.0, 0.0, 0.0), (0.01, 0.0, 0.0, 100.0), (0.0, 0.02, 0.0, -50.0))  # x, y, z (m), Delta p (Pa)
REFERENCE_PRESSURE = 1000.0  # Pa
FIELD_CASE = 2  # Mach 3, alpha 10
FIELD_POINTS_MM = ((-50.0, 0.0, 0.0), (0.0, 0.0, -50.0), (0.0, 0.0, 50.0), (0.0, -50.0, 0.0))


def build():
    import aero_blocks as a
    import buildup_blocks as bb
    import pressure_blocks as pb
    from kit import PRESSURE
    from ntopfile import Graph

    labels = {"Export Pressures": pb.EXPORT_CHOICES, "Impact Method": a.TW_IMPACT}

    flight_b, si_b, mesh_b, tw_b = (a.flight_block(), a.newtonian_block(), a.newtonian_mesh_block(),
                                    a.tangent_wedge_mesh_block())
    contribution_b, database_b = bb.coefficients_contribution_block(), bb.database_block()
    reader_b, jump_b, field_b = pb.surface_pressure_block(), pb.surface_pressure_jump_block(), a.newtonian_pressure_block()
    g = Graph()
    bodies = {
        "sphere": g.block("sphere<point,real>", "sphere", [g.point(0, 0, 0), g.real(SPHERE_R, LENGTH)], "Sphere"),
        "plate": g.block("box<point,real,real,real>", "box", [g.point(0, 0, 0)] + [g.real(v, LENGTH) for v in PLATE],
                         "Plate"),
        "cone": g.block("cone<point,point,real,real>", "cone",
                        [g.point(0, 0, 0), g.point(CONE_LENGTH * 1e3, 0, 0), g.real(0.0, LENGTH),
                         g.real(CONE_LENGTH * math.tan(math.radians(CONE_HALF_ANGLE)), LENGTH)], "Cone"),
    }
    meshes = {k: _mesh(g, b, k.title() + " Mesh") for k, b in bodies.items()}
    bodies["sphere_coarse"] = bodies["sphere"]
    meshes["sphere_coarse"] = _mesh(g, bodies["sphere"], "Coarse Sphere Mesh", COARSE_TOL)
    work = lambda g: g.file_path(WORK)  # noqa: E731
    flight = _call(g, flight_b, {
        "Mach": _lst(g, "real", [g.real(m) for m in MACHS]),
        "Altitude": _lst(g, "real", [g.real(0.0, LENGTH)]),
        "Angle of Attack": _lst(g, "real", [g.real(v * DEG, ANGLE) for v in ALPHAS]),
        "Sideslip": _lst(g, "real", [g.real(0.0, ANGLE)]),
        "Reference Area": g.real(SREF, AREA), "Reference Span": g.real(CREF, LENGTH),
        "Reference Chord": g.real(CREF, LENGTH), "Working Directory": work, "Case Name": g.text("flight")}, "Flight")
    outs = []
    for name, kind, body, stagnation, extra in CASES:
        given = {"Body": bodies[body], "Flight Condition": flight, "Case Name": g.text(name),
                 "Stagnation Pressure": lambda g, i=stagnation: g.choice(list(a.newtonian.CP_MAX_MODES), i),
                 "Cp Max Field": lambda g: g.leaf("real_field", {"expression": repr(USER_CP_MAX), "units": {}})}
        if kind != "si":
            given.update({"Mesh": meshes[body], "Working Directory": work})
        given.update({key: (lambda g, key=key, i=i: g.choice(list(labels[key]), i)) for key, i in extra.items()})
        result = _call(g, {"si": si_b, "mesh": mesh_b, "tw": tw_b}[kind], given, name)
        # The Aero Database of that one contribution, on the same flight condition
        contribution = _call(g, contribution_b, {"Results": result}, name + " contribution")
        outs.append(_call(g, database_b, {"Flight Condition": flight, "Contributions": _lst(g, "text", [contribution]),
                                          "Working Directory": work, "Case Name": g.text("db_" + name)},
                          name + " database"))
        if name == READ_CASE[0]:
            read_from = result
    reads = [
        _peek(g, _call(g, reader_b, {"Results": read_from, "Case Number": g.integer(READ_CASE[1])}, "Read Pressure"),
              READ_ROWS, "Pressure Rows"),
        _peek(g, _call(g, jump_b, {"Results": g.text(WORK + "/fake"), "Case Number": g.integer(1),
                                   "Reference Pressure": g.real(REFERENCE_PRESSURE, PRESSURE)}, "Read Jump"),
              range(len(FAKE_JUMP)), "Jump Rows"),
        g.block("evaluate_field<real_field,list<point>>", "list<real>",
                [_call(g, field_b, {"Body": bodies["sphere"], "Flight Condition": flight,
                                    "Case Number": g.integer(FIELD_CASE),
                                    "Stagnation Pressure": lambda g: g.choice(list(a.newtonian.CP_MAX_MODES), 1)},
                       "Newtonian Pressure Field"),
                 _lst(g, "point", [g.point(*p) for p in FIELD_POINTS_MM])], "Field Values"),
    ]
    group = g.block("core.group<any>", "group", outs + reads, "Results")
    return g, [flight_b, si_b, mesh_b, tw_b, contribution_b, database_b, reader_b, jump_b, field_b], group


MEMORY_CAP_GB = float(os.environ.get("NTOPCL_MEMORY_CAP_GB", "8"))
MAX_THREADS = int(os.environ.get("NTOPCL_MAX_THREADS", max(1, (os.cpu_count() or 2) // 2)))


def _working_set_gb(pid):
    """Resident memory of a process in GB (Windows), or 0 when it cannot be read."""
    if os.name != "nt":
        return 0.0
    import ctypes
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + \
                   [(n, ctypes.c_size_t) for n in ("PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                                                   "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                                                   "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.OpenProcess.restype = wintypes.HANDLE
    handle = k32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
    if not handle:
        return 0.0
    try:
        c = Counters()
        c.cb = ctypes.sizeof(c)
        ok = ctypes.WinDLL("psapi").GetProcessMemoryInfo(wintypes.HANDLE(handle), ctypes.byref(c), c.cb)
        return max(c.WorkingSetSize, c.PagefileUsage) / 1e9 if ok else 0.0
    finally:
        k32.CloseHandle(wintypes.HANDLE(handle))


def run(g, blocks, output, path, timeout=900, settle=6):
    """Save and run the notebook."""
    from ntopfile import save_host
    save_host(g, blocks, path, output=output, output_type="group")
    res = os.path.splitext(path)[0] + ".json"
    log_path = os.path.splitext(path)[0] + ".log"
    if os.path.exists(res):
        os.remove(res)
    env = dict(os.environ)
    env.setdefault("NTOP_PYTHON_PATH", sys.executable)
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        p = subprocess.Popen([NTOPCL, "--trustnotebook", "-v", "2", "-m", str(MAX_THREADS), "-o", res, path], stdout=logf,
                             stderr=subprocess.STDOUT, env=env,
                             creationflags=getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0))
        last, stable = -1, 0
        while p.poll() is None and time.time() - t0 < timeout:
            time.sleep(1)
            used = _working_set_gb(p.pid)
            if used > MEMORY_CAP_GB:
                logf.write(f"\n[E]: test runner: ntopcl used {used:.1f} GB (cap {MEMORY_CAP_GB:g} GB); stopped\n")
                break
            with open(log_path, encoding="utf-8", errors="replace") as f:
                if "nTop exited" in f.read():
                    break
            if os.path.exists(res):
                size = os.path.getsize(res)
                stable = stable + 1 if size == last and size > 0 else 0
                last = size
                if stable >= settle:
                    break
        if p.poll() is None:
            p.kill()
            p.wait()
    return open(log_path, encoding="utf-8", errors="replace").read()


def _read(path):
    with open(path, newline="") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def _clean():
    os.makedirs(WORK, exist_ok=True)
    for f in os.listdir(WORK):
        os.remove(os.path.join(WORK, f))


@pytest.fixture(scope="session")
def results():
    """Case name -> rows of its coefficient table (dicts keyed by column), ``db_<case>`` -> rows of the Aero Database
    of that case alone, plus ``_files`` (the run folder's files) and ``_log``."""
    if not os.path.isfile(NTOPCL):
        pytest.skip(f"ntopcl not found at {NTOPCL}")
    _clean()
    with open(os.path.join(WORK, "fake.ntvsp.pressure_jump.1.csv"), "w") as f:
        f.writelines(",".join(repr(v) for v in row) + "\n" for row in FAKE_JUMP)
    g, blocks, out = build()
    log = run(g, blocks, out, os.path.join(WORK, "e2e.ntop"))
    if "Could not login" in log:
        pytest.skip("ntopcl could not log in: sign in to the nTop app first")
    tables = {}
    for name, kind, *_ in CASES:
        for key, path in ((name, (f"flight_{name}" if kind == "si" else name) + ".ntvsp.coefficients.csv"),
                          ("db_" + name, f"db_{name}.ntvsp.database.csv")):
            if os.path.exists(os.path.join(WORK, path)):
                tables[key] = _read(os.path.join(WORK, path))
    tables["_files"] = set(os.listdir(WORK))
    res = os.path.join(WORK, "e2e.json")
    tables["_json"] = json.load(open(res)) if os.path.exists(res) else None
    tables["_log"] = log
    return tables


@pytest.fixture
def table(results, request):
    """The rows of one case (``@pytest.mark.parametrize('table', [...], indirect=True)``)."""
    name = request.param
    if name not in results:
        errs = "\n".join(l for l in results["_log"].splitlines() if "[E]" in l)[:2000]
        pytest.fail(f"{name}: no coefficient table written. ntopcl errors:\n{errs}")
    return results[name]
