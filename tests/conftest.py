"""End-to-end tests: run the prebuilt notebook e2e.ntop (nTop_vsp blocks on bodies with analytic answers) with
ntopcl once per session, and read the coefficient tables the blocks write.

    python -m pytest tests            (about 2 minutes)
"""
from __future__ import annotations

import csv
import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from cases import *  # noqa: E402,F401,F403
from cases import CASES, FAKE_JUMP  # noqa: E402
from vsp_ntop.ntopcl import NTOPCL, run_notebook  # noqa: E402

ROOT = os.path.dirname(HERE)
NOTEBOOK = os.path.join(HERE, "e2e.ntop")
WORK = os.path.join(ROOT, "runs", "tests").replace("\\", "/")


def run():
    """Run the notebook on WORK."""
    env = dict(os.environ)
    env.setdefault("NTOP_PYTHON_PATH", sys.executable)
    inputs = [{"name": "Working Directory", "type": "file_path", "value": WORK}]
    return run_notebook(NOTEBOOK, inputs, os.path.join(WORK, "e2e"), env=env)[0]


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
    res = os.path.join(WORK, "e2e.output.json")
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
