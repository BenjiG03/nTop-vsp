"""Surface pressure per case from VSPAERO's binary ``.adb`` (aerothermal database) file.

The layout was worked out against VSPAERO 7.2.2 (OpenVSP 3.56) outputs; the Tecplot file VSPAERO writes for
its last case is used as a check (see ``check_against_tecplot``):

* int32 words. Word 5: node count, word 6: triangle count (VSPAERO's, which is the ``.vspgeom`` triangles in
  order, followed by extra triangles it adds along wakes).
* For each solved (Mach, alpha, beta) condition, one block of ``3 * triangles`` float32 words: per triangle
  (Cp, unsteady Cp, a third quantity). The unsteady Cp column is exactly zero in steady runs, which is how the
  blocks are found. Reynolds-number variants share the same block (the potential solution does not depend on it).

The conditions are listed, in order, in ``<case>.adb.cases`` ("Mach AoA Beta Case: n ...").
"""
from __future__ import annotations

import math
from array import array


def read_cases(path):
    """[(mach, alpha, beta)] from ``<case>.adb.cases``."""
    out = []
    with open(path) as f:
        for line in f:
            a = line.split()
            if len(a) >= 3:
                try:
                    out.append(tuple(float(v) for v in a[:3]))
                except ValueError:
                    pass
    return out


def read_cp(path, n_cases):
    """List (one per case) of per-triangle Cp lists, in VSPAERO's triangle order."""
    with open(path, "rb") as f:
        raw = f.read()
    n = len(raw) // 4
    ints = array("i")
    ints.frombytes(raw[: 4 * n])
    floats = array("f")
    floats.frombytes(raw[: 4 * n])
    n_tri = ints[6]
    width = 3 * n_tri
    blocks = []
    k = 0
    while k + width <= n and len(blocks) < n_cases:
        if ints[k + 1] == 0 and ints[k + 4] == 0 and ints[k + 7] == 0 and floats[k] != 0.0:
            mid = ints[k + 1:k + width:3]
            if not any(mid):
                cp = floats[k:k + width:3]
                if all(math.isfinite(v) and abs(v) < 1e4 for v in cp):
                    blocks.append(list(cp))
                    k += width
                    continue
        k += 1
    if len(blocks) != n_cases:
        raise ValueError(f"{path}: found {len(blocks)} pressure blocks for {n_cases} cases")
    return blocks, n_tri


def read_tecplot_cp(path):
    """Cell-centred Cp of VSPAERO's Tecplot file (its last case)."""
    with open(path) as f:
        lines = f.read().split("\n")
    head = lines[1]
    n = int(head.split("N=")[1].split(",")[0])
    e = int(head.split("E=")[1].split(",")[0])
    vals = []
    i = 4
    while len(vals) < 3 * n + e:
        vals += [float(v) for v in lines[i].split()]
        i += 1
    return vals[3 * n:3 * n + e]


def check_against_tecplot(blocks, tec_path, tol=1e-4):
    tec = read_tecplot_cp(tec_path)
    last = blocks[-1]
    err = max(abs(a - b) for a, b in zip(last, tec))
    if err > tol:
        raise ValueError(f"surface pressures read from the .adb disagree with the Tecplot file by {err:.3g}")
    return err
