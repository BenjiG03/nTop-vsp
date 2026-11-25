"""Run an nTop notebook headlessly with nTop Automate (ntopcl), within memory and thread caps."""
from __future__ import annotations

import json
import os
import subprocess
import time

NTOPCL = os.environ.get("NTOPCL", r"C:\Program Files\nTopology\nTopology\ntopcl.exe")
MEMORY_CAP_GB = float(os.environ.get("NTOPCL_MEMORY_CAP_GB", "8"))
MAX_THREADS = int(os.environ.get("NTOPCL_MAX_THREADS", max(1, (os.cpu_count() or 2) // 2)))


def working_set_gb(pid):
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


def run_notebook(notebook, inputs, base, meta=None, env=None, timeout=900, settle=6):
    """Run ``notebook`` with its inputs set from ``inputs`` (ntopcl -j entries), writing ``base.input.json``,
    ``base.output.json`` and ``base.log``. ntopcl runs at low priority on ``MAX_THREADS`` threads; it is stopped
    above ``MEMORY_CAP_GB``, at ``timeout`` seconds, or once its output stops growing (it may not exit after Run
    Command blocks). Returns (log text, seconds)."""
    jpath, out, log_path = (f"{base}.{e}" for e in ("input.json", "output.json", "log"))
    with open(jpath, "w") as f:
        json.dump({**(meta or {}), "inputs": inputs}, f, indent=1)
    if os.path.exists(out):
        os.remove(out)
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8", errors="replace") as logf:
        cmd = [NTOPCL, "--trustnotebook", "-v", "2", "-m", str(MAX_THREADS), "-j", jpath, "-o", out, notebook]
        p = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT, env=env,
                             creationflags=getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0))
        last, stable = -1, 0
        while p.poll() is None and time.time() - t0 < timeout:
            time.sleep(1)
            used = working_set_gb(p.pid)
            if used > MEMORY_CAP_GB:
                logf.write(f"\n[E]: ntopcl used {used:.1f} GB (cap {MEMORY_CAP_GB:g} GB); stopped\n")
                break
            with open(log_path, encoding="utf-8", errors="replace") as f:
                if "nTop exited" in f.read():
                    break
            if os.path.exists(out):
                size = os.path.getsize(out)
                stable = stable + 1 if size == last and size > 0 else 0
                last = size
                if stable >= settle:
                    break
        if p.poll() is None:
            p.kill()
            p.wait()
    return open(log_path, encoding="utf-8", errors="replace").read(), time.time() - t0
