"""Modified Newtonian forces and moments of an implicit body, integrated over nTop's mesh of it.

    python -m vsp_ntop.run_newtonian_mesh faces STL_FILE FACES_CSV      (step 1, from the block)
    python -m vsp_ntop.run_newtonian_mesh SETTINGS_FILE                  (step 2)

The Modified Newtonian on Mesh (VSP) block meshes the implicit body in nTop (Mesh from Implicit Body) and exports it
as an STL in millimetres. Step 1 writes the faces' centroids as a point map (``x, y, z, area``, mm and mm^2) that
nTop imports; nTop evaluates the body's gradient and the Cp Max Field at those centroids and writes them, one row per
face, as point maps ``<case>.ntvsp.mesh_gradient.csv`` and ``<case>.ntvsp.mesh_cp_max.csv``. Step 2 does the sums at every
(Mach, AoA, Beta) of the flight condition within Mach Min to Mach Max:

    Cp = Cp_max max(-n . V, 0)^2,  n = grad f / |grad f| at the face centroid (the implicit's own normal)
    C_F = -(1 / Sref) sum Cp S_f,   C_M = -(1 / (Sref cref)) sum Cp (r_f - c) x S_f

with S_f the face's vector area (outward, from the mesh) and r_f its centroid: the pressure is the implicit body's,
the integral the mesh's. Cp_max: Modified (behind a normal shock above Mach 1, Rayleigh pitot; isentropic below),
Classical (2) or User Defined (the Cp Max Field at each face).

Writes ``<case>.ntvsp.coefficients.csv`` with the Modified Newtonian (VSP) block's layout (Mach, AoA, Beta, CA, CY,
CN, Cl, Cm, Cn about x = Aft Axis, z = Up Axis made normal to it, y = z x x, every moment over Sref cref, and the
reference), for Coefficients Contribution (VSP). Needs numpy.
"""
from __future__ import annotations

import math
import os
import struct
import sys

from .common import DriverError, Log, Settings, main_wrapper, write_csv
from .flight import Flight, _rows

MM = 0.001
CP_MAX_MODES = ("modified", "classical", "user")
HEADER = ("Mach", "AoA", "Beta", "CA", "CY", "CN", "Cl", "Cm", "Cn", "Sref", "Bref", "Cref", "Xref", "Yref", "Zref")


def _np():
    from .run_plot import conda_dll_dirs
    conda_dll_dirs()  # a conda Python started by nTop needs its DLL folders for numpy
    try:
        import numpy
    except ImportError:  # pragma: no cover
        raise DriverError("Modified Newtonian on Mesh needs numpy in the Python named by NTOP_PYTHON_PATH")
    return numpy


def read_stl(path):
    """Triangles of an ASCII or binary STL as an (n, 3, 3) array, in the file's units."""
    np = _np()
    if not os.path.isfile(path):
        raise DriverError(f"mesh file {path} not found")
    with open(path, "rb") as f:
        data = f.read()
    if len(data) >= 84:
        n = struct.unpack_from("<I", data, 80)[0]
        if 84 + 50 * n == len(data):
            rec = np.frombuffer(data, dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]),
                                count=n, offset=84)
            return rec["v"].astype(float)
    verts = [[float(t) for t in line.split()[1:4]] for line in data.decode("ascii", "ignore").splitlines()
             if line.strip().startswith("vertex")]
    if not verts or len(verts) % 3:
        raise DriverError(f"{path} is not an STL file with triangles")
    return np.asarray(verts, dtype=float).reshape(-1, 3, 3)


def faces(stl_path):
    """(centroids, vector areas) of the mesh's non-degenerate triangles, in the file's units."""
    np = _np()
    tri = read_stl(stl_path)
    vec = 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(vec, axis=1)
    keep = area > 1e-12 * max(area.max(), 1e-300)
    return tri[keep].mean(axis=1), vec[keep]


def write_faces(stl_path, out_csv):
    """Step 1: the faces' centroids and areas as a point map for nTop (x, y, z, area; no header)."""
    np = _np()
    c, s = faces(stl_path)
    area = np.linalg.norm(s, axis=1)
    np.savetxt(out_csv, np.column_stack([c, area]), delimiter=",", fmt="%.9g")
    return len(c)


def cp_max_stagnation(mach, gamma=1.4):
    """Stagnation pressure coefficient: behind a normal shock above Mach 1 (Rayleigh pitot), isentropic below."""
    if mach <= 0:
        return 1.0
    m2 = mach * mach
    if mach > 1.0:
        p02 = ((gamma + 1) ** 2 * m2 / (4 * gamma * m2 - 2 * (gamma - 1))) ** (gamma / (gamma - 1)) * \
            (1 - gamma + 2 * gamma * m2) / (gamma + 1)
    else:
        p02 = (1 + 0.5 * (gamma - 1) * m2) ** (gamma / (gamma - 1))
    return 2.0 / (gamma * m2) * (p02 - 1.0)


def frame(aft, up):
    np = _np()
    x = np.asarray(aft, float)
    x = x / np.linalg.norm(x)
    z = np.asarray(up, float)
    z = z - x * (z @ x)
    if np.linalg.norm(z) < 1e-9:
        raise DriverError("Up Axis is parallel to Aft Axis")
    z = z / np.linalg.norm(z)
    return x, np.cross(z, x), z


def coefficients(centroids, vector_areas, normals, cp_max_field, conditions, mode, axes, ref, gamma=1.4):
    """Rows of HEADER for each condition (dict with Mach, AoA, Beta). Lengths in metres, areas in m^2."""
    np = _np()
    x, y, z = axes
    sref, _, cref, centre = ref
    arm = centroids - np.asarray(centre, float)
    moment_area = np.cross(arm, vector_areas)  # (r - c) x S per face
    rows = []
    for c in conditions:
        a, b = math.radians(c["AoA"]), math.radians(c["Beta"])
        v = x * math.cos(a) * math.cos(b) - y * math.sin(b) + z * math.sin(a) * math.cos(b)
        sine = np.maximum(-(normals @ v), 0.0)
        if mode == "user":
            cpm = cp_max_field
        else:
            cpm = cp_max_stagnation(c["Mach"], gamma) if mode == "modified" else 2.0
        cp = cpm * sine * sine
        force = -(cp[:, None] * vector_areas).sum(axis=0) / sref
        moment = -(cp[:, None] * moment_area).sum(axis=0) / (sref * cref)
        rows.append([c["Mach"], c["AoA"], c["Beta"], force @ x, force @ y, force @ z, moment @ x, moment @ y,
                     moment @ z, sref, cref, cref, centre[0], centre[1], centre[2]])
    return rows


def read_values(base, centroids_mm):
    """nTop's values at the face centroids: (Cp Max Field, gradient (n, 3)), from its point-map exports
    (x, y, z in mm, then the value or vector), checked to list the faces in order."""
    np = _np()
    out = []
    for suffix, width in (("cp_max", 4), ("gradient", 6)):
        path = f"{base}.ntvsp.mesh_{suffix}.csv"
        if not os.path.isfile(path):
            raise DriverError(f"{path} not found (nTop's evaluation of the body at the face centroids)")
        v = np.loadtxt(path, delimiter=",", ndmin=2)
        if v.shape != (len(centroids_mm), width):
            raise DriverError(f"{path}: {v.shape[0]} rows of {v.shape[1]} values, expected {len(centroids_mm)} rows of "
                              f"{width} (one per face)")
        gap = np.abs(v[:, :3] - centroids_mm).max()
        if gap > 1e-3 * max(1.0, np.abs(centroids_mm).max()):
            raise DriverError(f"{path} does not list the faces in the mesh's order (points {gap:.3g} mm apart)")
        out.append(v[:, 3] if width == 4 else v[:, 3:6])
    return out


def unit_normals(gradients, vector_areas, log):
    """The implicit's unit normals; where its gradient vanishes or is not finite, the face's own normal."""
    np = _np()
    g = np.where(np.isfinite(gradients), gradients, 0.0)
    norm = np.linalg.norm(g, axis=1)
    bad = norm < 1e-9
    face = vector_areas / np.linalg.norm(vector_areas, axis=1)[:, None]
    n = np.where(bad[:, None], face, g / np.where(bad, 1.0, norm)[:, None])
    if bad.any():
        log(f"note: {int(bad.sum())} faces had no usable gradient; their mesh normals are used")
    flipped = (n * face).sum(axis=1) < 0.0
    if flipped.mean() > 0.5:
        raise DriverError("the body's gradient points against the mesh's outward normals on most faces; is the body "
                          "an implicit that is negative inside?")
    return n


def main(settings_path):
    np = _np()
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    base = os.path.join(work, case)
    flight = s.text("flight").strip()
    fl = Flight(flight)
    rows = _rows(flight + ".ntvsp.inviscid.csv")
    lo, hi = s.real("mach_min"), s.real("mach_max")
    conds = [r for r in rows if lo <= r["Mach"] < hi]
    mode = s.text("stagnation").strip().lower()
    if mode not in CP_MAX_MODES:
        raise DriverError(f"stagnation = {mode}: expected one of {', '.join(CP_MAX_MODES)}")
    c_mm, s_mm2 = faces(base + ".ntvsp.mesh.stl")
    cp_max_field, gradient = read_values(base, c_mm)
    centroids, vector_areas = c_mm * MM, s_mm2 * MM * MM
    normals = unit_normals(gradient, vector_areas, log)
    closure = np.linalg.norm(vector_areas.sum(axis=0)) / np.linalg.norm(vector_areas, axis=1).sum()
    log(f"Modified Newtonian on mesh: {len(c_mm)} faces, wetted area {np.linalg.norm(vector_areas, axis=1).sum():.6g} "
        f"m^2, closure {closure:.1e}; {len(conds)} of {len(rows)} conditions in Mach [{lo:g}, {hi:g}); Cp_max {mode}")
    if closure > 1e-3:
        log("note: the mesh is not closed (its vector areas do not sum to zero)")
    axes = frame([s.real(f"aft_{c}") for c in "xyz"], [s.real(f"up_{c}") for c in "xyz"])
    out = coefficients(centroids, vector_areas, normals, cp_max_field, conds, mode, axes, fl.ref)
    write_csv(base + ".ntvsp.coefficients.csv", HEADER, out)
    log("done")
    log.close()


def _faces_main():
    if len(sys.argv) != 4:
        sys.stderr.write("usage: python -m vsp_ntop.run_newtonian_mesh faces STL_FILE FACES_CSV\n")
        sys.exit(2)
    try:
        n = write_faces(sys.argv[2], sys.argv[3])
        print(f"{n} faces")
    except DriverError as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "faces":
        _faces_main()
    else:
        main_wrapper(main)
