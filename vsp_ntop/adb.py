"""Surface pressure per case from VSPAERO's binary ``.adb`` (aerothermal database) file."""
from __future__ import annotations

import numpy as np


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


def read_cp(path, cases, n_surface):
    """List (one per case of ``cases``) of per-loop Cp arrays, in VSPAERO's level-0 loop order."""
    with open(path, "rb") as f:
        raw = f.read()
    n = len(raw) // 4
    ints = np.frombuffer(raw, "<i4", n)
    floats = np.frombuffer(raw, "<f4", n).astype(float)
    n_tri, n_nodes = int(ints[6]), int(ints[5])
    width = 3 * n_tri
    blocks = []
    pos = 0
    for mach, alpha, beta in cases:
        a, b = np.radians(alpha), np.radians(beta)
        k = np.arange(pos, max(pos, n - 5))
        hit = np.flatnonzero((np.abs(floats[k] - mach) <= 1e-5 * max(1.0, abs(mach)))
                             & (np.abs(floats[k + 1] - a) < 2e-6) & (np.abs(floats[k + 2] - b) < 2e-6))
        if not len(hit):
            raise ValueError(f"{path}: no results header for Mach {mach}, alpha {alpha}, beta {beta}")
        k = np.arange(k[hit[0]] + 5, max(k[hit[0]] + 5, n - width))
        cand = k[(ints[k + 1] == 0) & (ints[k + 4] == 0) & (ints[k + width] >= 0) & (ints[k + width] <= n_nodes)]
        found = None
        for k in cand:
            if ints[k + 1:k + width:3].any():
                continue
            cp = floats[k:k + width:3]
            if (np.isfinite(cp).all() and (np.abs(cp) < 1e4).all()
                    and np.count_nonzero(cp[:n_surface]) * 2 > n_surface):
                found = cp
                break
        if found is None:
            raise ValueError(f"{path}: no surface pressures after the header of Mach {mach}, alpha {alpha}, beta {beta}")
        blocks.append(found)
        pos = k + width
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
        vals += lines[i].split()
        i += 1
    return np.array(vals[3 * n:3 * n + e], float)


def check_against_tecplot(blocks, tec_path, tol=1e-4):
    tec = read_tecplot_cp(tec_path)
    last = blocks[-1]
    m = min(len(last), len(tec))
    err = np.abs(last[:m] - tec[:m]).max()
    if err > tol:
        raise ValueError(f"surface pressures read from the .adb disagree with the Tecplot file by {err:.3g}")
    return err
