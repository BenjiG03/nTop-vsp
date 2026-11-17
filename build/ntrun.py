"""Run generated notebooks headlessly with ntopcl and read their output back.

Python only builds notebooks and reads files here; every value comes from nTop.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ntopfile import Graph, Ref, save_host  # noqa: E402,F401

NTOPCL = r"C:\Program Files\nTopology\nTopology\ntopcl.exe"
LENGTH = {"length": 1}
AREA = {"length": 2}
ANGLE = {"angle": 1}


def run(graph: Graph, output, output_type, path, blocks=(), timeout=3600, trust=False, verbose=2):
    """Save ``graph`` as a host notebook at ``path``, run it and return (json or None, log, seconds)."""
    save_host(graph, list(blocks), path, output=output, output_type=output_type)
    res = os.path.splitext(path)[0] + ".json"
    if os.path.exists(res):
        os.remove(res)
    cmd = [NTOPCL] + (["--trustnotebook"] if trust else []) + ["-v", str(verbose), "-o", res, path]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    dt = time.time() - t0
    log = p.stdout + p.stderr
    out = json.load(open(res)) if os.path.exists(res) else None
    return out, log, dt


def errors(log):
    return [line for line in log.splitlines() if "[E]" in line or "rror" in line]


def values(out):
    """Plain numbers from an ntopcl output: a real, a vector, or a list of them."""
    def conv(v):
        if isinstance(v, dict):
            if "val" in v:
                return v["val"]
            if v.get("value") is not None:
                return conv(v["value"])
            if "components" in v:
                return [conv(c) for c in v["components"]]
        if isinstance(v, list):
            return [conv(c) for c in v]
        return v
    return conv(out[0] if isinstance(out, list) and len(out) == 1 else out)


def run_notebook(path, timeout=3600):
    res = os.path.splitext(path)[0] + ".json"
    if os.path.exists(res):
        os.remove(res)
    t0 = time.time()
    log_path = os.path.splitext(path)[0] + ".log"
    # ntopcl writes its output file and then may never exit after notebooks that use Run Command inside custom
    # blocks, so wait for the output file to settle and stop it.
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        p = subprocess.Popen([NTOPCL, "--trustnotebook", "-v", "2", "-o", res, path], stdout=logf,
                             stderr=subprocess.STDOUT)
        last = -1
        while p.poll() is None and time.time() - t0 < timeout:
            time.sleep(2)
            if os.path.exists(res):
                size = os.path.getsize(res)
                if size == last and size > 0:
                    break
                last = size
        if p.poll() is None:
            p.kill()
            p.wait()
    log = open(log_path, encoding="utf-8", errors="replace").read()
    out = json.load(open(res)) if os.path.exists(res) else None
    return out, log, time.time() - t0
