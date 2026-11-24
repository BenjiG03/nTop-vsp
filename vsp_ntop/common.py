"""Shared helpers for the nTop -> OpenVSP drivers."""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import time

LENGTH_UNITS = {  # VSP's LEN_* enum: name, metres per unit
    0: ("mm", 0.001), 1: ("cm", 0.01), 2: ("m", 1.0), 3: ("in", 0.0254), 4: ("ft", 0.3048), 5: ("yd", 0.9144),
    6: ("unitless", 1.0),
}


class DriverError(RuntimeError):
    pass


def read_settings(path) -> dict:
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


class Settings:
    """Typed access to the settings file, with clear errors naming the missing key."""

    def __init__(self, path):
        self.path = path
        self.raw = read_settings(path)

    def text(self, key, default=None):
        if key in self.raw:
            return self.raw[key]
        if default is None:
            raise DriverError(f"setting '{key}' missing from {self.path}")
        return default

    def real(self, key, default=None):
        v = self.text(key, None if default is None else str(default))
        try:
            return float(v)
        except ValueError:
            raise DriverError(f"setting '{key}' = '{v}' is not a number")

    def integer(self, key, default=None):
        return int(round(self.real(key, default)))

    def flag(self, key, default=None):
        v = self.text(key, None if default is None else str(default)).lower()
        return v in ("1", "on", "yes", "true")

    def reals(self, key, sep=";"):
        v = self.text(key, "")
        return [float(x) for x in v.split(sep) if x.strip()]


class Log:
    """Writes progress to a log file and to stdout (nTop shows stdout in its log panel)."""

    def __init__(self, path):
        self.path = path
        self.f = open(path, "w", encoding="utf-8")
        self.t0 = time.time()

    def __call__(self, msg):
        line = f"[{time.time() - self.t0:7.1f}s] {msg}"
        self.f.write(line + "\n")
        self.f.flush()
        print(line, flush=True)

    def close(self):
        self.f.close()


def linspace(start, end, n):
    """VSP's sweep points: n values from start to end inclusive (start alone when n is 1)."""
    n = max(1, int(n))
    if n == 1:
        return [start]
    return [start + i * (end - start) / (n - 1) for i in range(n)]


def fmt(x):
    """A number as VSP writes it: shortest round-trip form."""
    if float(x).is_integer() and abs(x) < 1e15:
        return str(int(x))
    return repr(float(x))


def openvsp_dir(s: "Settings"):
    """OpenVSP's folder: the settings' ``vsp_dir`` (the blocks fill it from the VSP_PATH environment variable), or
    VSP_PATH itself. A path to one of OpenVSP's programs is taken as its folder."""
    path = s.raw.get("vsp_dir", "").strip() or os.environ.get("VSP_PATH", "").strip()
    if not path:
        raise DriverError("OpenVSP not found: set the VSP_PATH environment variable to the OpenVSP folder "
                          "(the one holding vspscript.exe and vspaero.exe)")
    path = os.path.abspath(path.strip('"'))
    if os.path.isfile(path):
        path = os.path.dirname(path)
    if not os.path.isfile(os.path.join(path, "vspscript.exe")):
        raise DriverError(f"VSP_PATH = {path}: no vspscript.exe there (point it at the OpenVSP folder)")
    return path


def vsp_path(p):
    return os.path.abspath(p).replace("\\", "/")


def run(cmd, log: Log, cwd, stdout_path, timeout=None):
    """Run a program, sending its output to ``stdout_path``. Raises DriverError on failure."""
    log("run: " + " ".join(f'"{c}"' if " " in c else c for c in cmd))
    with open(stdout_path, "w", encoding="utf-8", errors="replace") as out:
        try:
            p = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=subprocess.STDOUT, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise DriverError(f"{os.path.basename(cmd[0])} timed out after {timeout} s; see {stdout_path}")
    if p.returncode != 0:
        raise DriverError(f"{os.path.basename(cmd[0])} exited with code {p.returncode}; see {stdout_path}")
    return p.returncode


def run_vspscript(vsp_dir, script_text, work, stem, log: Log, timeout=None):
    """Write ``script_text`` to ``<stem>.vspscript`` in ``work`` and run it with vspscript.exe."""
    exe = os.path.join(vsp_dir, "vspscript.exe")
    if not os.path.isfile(exe):
        raise DriverError(f"vspscript.exe not found in OpenVSP directory '{vsp_dir}'")
    script = os.path.join(work, stem + ".vspscript")
    with open(script, "w", encoding="utf-8", newline="\n") as f:
        f.write(script_text)
    out = os.path.join(work, stem + ".vspscript.log")
    run([exe, "-script", script], log, work, out, timeout)
    text = open(out, encoding="utf-8", errors="replace").read()
    if "NTVSP_OK" not in text:
        errors = [l.strip() for l in text.splitlines() if "ERR" in l or "Error" in l]
        raise DriverError("OpenVSP script failed: " + ("; ".join(errors[:5]) or "no completion marker") + f" (see {out})")
    return text


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if header:
            w.writerow(header)
        for r in rows:
            w.writerow(r)


def write_point_map(path, points, values):
    """CSV for nTop's Import Scalar Point Map: x, y, z, value per row, no header."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        for p, v in zip(points, values):
            f.write(f"{p[0]:.9g},{p[1]:.9g},{p[2]:.9g},{v:.9g}\n")


def read_vsp_results_csv(path):
    """VSP's WriteResultsCSVFile output: blocks of ``Results_Name,...`` then ``name,v1,v2,...`` rows."""
    results = {}
    current = None
    if not os.path.isfile(path):
        raise DriverError(f"OpenVSP wrote no results file {path}")
    with open(path, encoding="utf-8", errors="replace") as f:
        for row in csv.reader(f):
            if not row:
                continue
            key = row[0].strip()
            if key == "Results_Name":
                current = results.setdefault(row[1].strip() if len(row) > 1 else "", {})
                continue
            if current is not None and key:
                current.setdefault(key, []).append([c.strip() for c in row[1:]])
    return results


def conda_dll_dirs():
    """Put a conda Python's DLL folders on PATH. nTop starts the Python of NTOP_PYTHON_PATH without activating its
    environment, so numpy's BLAS/LAPACK DLLs are not found at its first matrix product and the process dies with no
    traceback."""
    if os.name != "nt":
        return
    dirs = [os.path.join(sys.prefix, *p) for p in (("Library", "mingw-w64", "bin"), ("Library", "usr", "bin"),
                                                   ("Library", "bin"), ("Scripts",))]
    base = getattr(sys, "base_prefix", sys.prefix)
    if base != sys.prefix:  # a venv made from a conda Python
        dirs += [os.path.join(base, *p) for p in (("Library", "mingw-w64", "bin"), ("Library", "bin"))]
    dirs = [d for d in dirs if os.path.isdir(d)]
    if not dirs:
        return
    os.environ["PATH"] = os.pathsep.join(dirs + [os.environ.get("PATH", "")])
    for d in dirs:
        try:
            os.add_dll_directory(d)
        except (AttributeError, OSError):
            pass


def main_wrapper(fn):
    """Run ``fn(settings_path)`` after ``conda_dll_dirs``; print a one-line error and exit 1 on failure (nTop shows
    stderr)."""
    if len(sys.argv) != 2:
        sys.stderr.write(f"usage: python -m vsp_ntop.{os.path.splitext(os.path.basename(sys.argv[0]))[0]} SETTINGS_FILE\n")
        sys.exit(2)
    conda_dll_dirs()
    try:
        fn(sys.argv[1])
    except DriverError as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
