"""Surface pressure per case from VSPAERO's binary ``.adb`` (aerothermal database) file.

Layout (VSPAERO 7.2.2 / OpenVSP 3.56, ``VSP_SOLVER::WriteOutAerothermalDatabase*``), in 4-byte words:

* A header and the geometry. Word 6 of the file is the number of level-0 loops: the triangles of the ``.vspgeom``
  file in order (for the polygon layout, the triangles of its triangulation section), followed by the triangles
  VSPAERO builds along the wakes.
* Per solved (Mach, alpha, beta) condition: Mach, alpha and beta (float32, angles in radians), Cp clip limits,
  solver-level arrays (doubles), then one (Cp, unsteady Cp, circulation) float32 triple per level-0 loop, then the
  number of trailing vortices (int32) and their data. On thick (panel) surfaces Cp is the surface pressure
  coefficient; on thin (vortex lattice) surfaces it is the pressure jump across the sheet: Cp on the side the
  triangle normals point to, minus Cp on the other side (negative where the sheet lifts toward its normals).

The conditions are listed, in order, in ``<case>.adb.cases`` ("Mach AoA Beta Case: n ..."). Each condition's header
is found by its Mach, alpha and beta; the triple block is the run of ``3 * loops`` words after it whose middle
column (unsteady Cp, zero in steady runs) is all zero and which is followed by the trailing-vortex count. The
Tecplot file VSPAERO writes for its last case checks the result (``check_against_tecplot``).
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


def _close(a, b, tol):
    return abs(a - b) <= tol * max(1.0, abs(b))


def read_cp(path, cases, n_surface):
    """List (one per case of ``cases``) of per-loop Cp lists, in VSPAERO's level-0 loop order.

    ``n_surface``: the number of surface triangles (they come first; the wake loops after them carry no pressure).
    Most of them have a non-zero Cp, which tells the pressure block from runs of zeros elsewhere in the file.
    """
    with open(path, "rb") as f:
        raw = f.read()
    n = len(raw) // 4
    ints = array("i")
    ints.frombytes(raw[: 4 * n])
    floats = array("f")
    floats.frombytes(raw[: 4 * n])
    n_tri, n_nodes = ints[6], ints[5]
    width = 3 * n_tri
    blocks = []
    pos = 0
    for mach, alpha, beta in cases:
        a, b = math.radians(alpha), math.radians(beta)
        head = next((k for k in range(pos, n - 5)
                     if _close(floats[k], mach, 1e-5) and abs(floats[k + 1] - a) < 2e-6 and abs(floats[k + 2] - b) < 2e-6),
                    None)
        if head is None:
            raise ValueError(f"{path}: no results header for Mach {mach}, alpha {alpha}, beta {beta}")
        found = None
        k = head + 5
        while k + width < n:
            if ints[k + 1] == 0 and ints[k + 4] == 0 and 0 <= ints[k + width] <= n_nodes:
                if not any(ints[k + 1:k + width:3]):
                    cp = floats[k:k + width:3]
                    if (all(math.isfinite(v) and abs(v) < 1e4 for v in cp)
                            and sum(1 for v in cp[:n_surface] if v != 0.0) * 2 > n_surface):
                        found = list(cp)
                        break
            k += 1
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
