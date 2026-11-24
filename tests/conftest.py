"""End-to-end tests: run the prebuilt notebook e2e.ntop (nTop_vsp blocks on bodies with analytic answers) with
ntopcl once per session, and read the coefficient tables the blocks write.

    python -m pytest tests            (about 2 minutes)
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from cases import *  # noqa: E402,F401,F403
from cases import CASES, FAKE_JUMP  # noqa: E402

ROOT = os.path.dirname(HERE)
NOTEBOOK = os.path.join(HERE, "e2e.ntop")
NTOPCL = os.environ.get("NTOPCL", r"C:\Program Files\nTopology\nTopology\ntopcl.exe")
WORK = os.path.join(ROOT, "runs", "tests").replace("\\", "/")


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


def run(timeout=900, settle=6):
    """Run the notebook on WORK."""
    jpath = os.path.join(WORK, "e2e.input.json")
    with open(jpath, "w") as f:
        json.dump({"inputs": [{"name": "Working Directory", "type": "file_path", "value": WORK}]}, f)
    path = os.path.join(WORK, "e2e.ntop")
    res = os.path.splitext(path)[0] + ".json"
    log_path = os.path.splitext(path)[0] + ".log"
    if os.path.exists(res):
        os.remove(res)
    env = dict(os.environ)
    env.setdefault("NTOP_PYTHON_PATH", sys.executable)
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        cmd = [NTOPCL, "--trustnotebook", "-v", "2", "-m", str(MAX_THREADS), "-j", jpath, "-o", res, NOTEBOOK]
        p = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT, env=env,
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
    log = run()
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
