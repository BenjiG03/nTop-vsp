"""Modified Newtonian or tangent-wedge forces and moments of an implicit body, integrated over nTop's mesh of it.

    python -m vsp_ntop.run_newtonian_mesh faces STL_FILE FACES_CSV      (step 1, from the block)
    python -m vsp_ntop.run_newtonian_mesh SETTINGS_FILE                  (step 2)

The block exports its Mesh input as an STL in millimetres. Step 1 writes the faces' centroids as a point map
(``x, y, z, area``, mm and mm^2). nTop evaluates the Body's gradient there (``<case>.ntvsp.mesh_gradient.csv``) and,
only when Cp_max is User Defined, the Cp Max Field (``<case>.ntvsp.mesh_cp_max.csv``). Step 2 does the sums at every
(Mach, AoA, Beta) of the flight condition within Mach Min to Mach Max:

    Cp = Cp_max max(-n . V, 0)^2,  n = grad f / |grad f| at the face centroid (the implicit Body's own normal)
    C_F = -(1 / Sref) sum Cp S_f,   C_M = -(1 / (Sref cref)) sum Cp (r_f - c) x S_f

with S_f the face's vector area (outward, from the Mesh) and r_f its centroid: the pressure is the implicit's, the
integral the mesh's. Cp_max: Modified (behind a normal shock above Mach 1, Rayleigh pitot; isentropic below),
Classical (2) or User Defined (the Cp Max Field at each face).

With ``method = tangent_wedge`` (the Tangent Wedge on Mesh (VSP) block) the face pressure is the tangent wedge's or,
with ``impact = cone``, the tangent cone's, with Prandtl-Meyer or Newtonian-shadow leeward faces
(``local_methods.py``); Cp_max then sets the pressure above shock detachment.

With ``export_pressure = 1`` it also writes the static pressure p = p_inf + q Cp at the face centroids for every flight
condition case it computed (``pressure.py``: ``<case>.ntvsp.pressure.<Case>.csv``, metres and pascals).

Writes ``<case>.ntvsp.coefficients.csv`` with the Modified Newtonian (VSP) block's layout (Mach, AoA, Beta, CA, CY,
CN, Cl, Cm, Cn about x = Aft Axis, z = Up Axis made normal to it, y = z x x, every moment over Sref cref, and the
reference), for Coefficients Contribution (VSP). Needs numpy.
"""
from __future__ import annotations

import math
import os
import struct
import sys

from . import pressure
from .common import DriverError, Log, Settings, conda_dll_dirs, main_wrapper, write_csv
from .flight import Flight, _rows

MM = 0.001
ALIGNED_DEG = 25.0  # a face is aligned when the Body's normal at its centroid is within this of its own normal
ALIGNED_COS = math.cos(math.radians(ALIGNED_DEG))
MIN_ALIGNED_SHARE = 0.5  # below this share of the mesh area aligned, the Mesh is not a mesh of the Body
CP_MAX_MODES = ("modified", "classical", "user")
HEADER = ("Mach", "AoA", "Beta", "CA", "CY", "CN", "Cl", "Cm", "Cn", "Sref", "Bref", "Cref", "Xref", "Yref", "Zref")


def _np():
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


def coefficients(centroids, vector_areas, normals, cp_max_field, conditions, mode, axes, ref, gamma=1.4,
                 method="newtonian", expansion="prandtl_meyer", separated="freestream", impact="wedge",
                 on_cp=None):
    """Rows of HEADER for each condition (dict with Mach, AoA, Beta). Lengths in metres, areas in m^2. ``method``:
    newtonian (Cp_max sin^2 windward, 0 leeward) or tangent_wedge (``local_methods``: windward the tangent wedge or,
    with ``impact = cone``, the tangent cone; ``expansion`` and ``separated`` for the leeward faces; Cp_max sets the
    pressure above shock detachment). ``on_cp(condition, cp)``, when given, gets every condition's face Cp."""
    np = _np()
    x, y, z = axes
    sref, _, cref, centre = ref
    arm = centroids - np.asarray(centre, float)
    moment_area = np.cross(arm, vector_areas)  # (r - c) x S per face
    rows = []
    for c in conditions:
        a, b = math.radians(c["AoA"]), math.radians(c["Beta"])
        v = x * math.cos(a) * math.cos(b) - y * math.sin(b) + z * math.sin(a) * math.cos(b)
        if mode == "user":
            cpm = cp_max_field
        else:
            cpm = cp_max_stagnation(c["Mach"], gamma) if mode == "modified" else 2.0
        if method == "tangent_wedge":
            from . import local_methods as lm
            delta = np.arcsin(np.clip(-(normals @ v), -1.0, 1.0))
            cp = lm.tangent_wedge_cp(delta, c["Mach"], cpm, gamma, expansion, separated, impact=impact)
        else:
            sine = np.maximum(-(normals @ v), 0.0)
            cp = cpm * sine * sine
        if on_cp is not None:
            on_cp(c, cp)
        force = -(cp[:, None] * vector_areas).sum(axis=0) / sref
        moment = -(cp[:, None] * moment_area).sum(axis=0) / (sref * cref)
        rows.append([c["Mach"], c["AoA"], c["Beta"], force @ x, force @ y, force @ z, moment @ x, moment @ y,
                     moment @ z, sref, cref, cref, centre[0], centre[1], centre[2]])
    return rows


def read_values(base, centroids_mm, user_cp_max):
    """nTop's values at the face centroids: (Cp Max Field or None, gradient (n, 3)), from its point-map exports
    (x, y, z in mm, then the value or vector), checked to list the faces in order. The Cp Max Field is exported (and
    read) only when Stagnation Pressure is User Defined."""
    np = _np()
    out = []
    for suffix, width in (("cp_max", 4), ("gradient", 6)):
        path = f"{base}.ntvsp.mesh_{suffix}.csv"
        if suffix == "cp_max" and not user_cp_max:
            out.append(None)
            continue
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
    cos = (n * face).sum(axis=1)
    flipped = cos < 0.0
    if flipped.mean() > 0.5:
        raise DriverError("the body's gradient points against the mesh's outward normals on most faces; is the body "
                          "an implicit that is negative inside?")
    area = np.linalg.norm(vector_areas, axis=1)
    alignment = float((cos * area).sum() / area.sum())
    aligned = float(area[cos >= ALIGNED_COS].sum() / area.sum())
    log(f"Body normals against Mesh faces: area-weighted mean cosine {alignment:.6f}; {100 * (1 - aligned):.1f}% of "
        f"the mesh area has normals more than {ALIGNED_DEG:g} degrees off the Body's")
    if aligned < MIN_ALIGNED_SHARE:
        raise DriverError(
            f"the Body's normals at the face centroids do not match the Mesh's faces ({100 * aligned:.0f}% of the mesh "
            f"area within {ALIGNED_DEG:g} degrees, mean cosine {alignment:.3f}): the Mesh is not a mesh of the Body "
            "where it is, for example an imported mesh in other length units, a moved or scaled body, or face "
            "centroids that nTop read or wrote in other units than millimetres")
    if aligned < 0.99:
        log("note: faces whose normals disagree with the Body's are often parts of the mesh inside the Body (OpenVSP's "
            "wing extensions into the fuselage, overlapping components); their pressure is still summed")
    return n


def face_normals(vector_areas, centroids):
    """Unit outward normals of the faces. The mesh should be closed with outward faces; one whose faces point inward
    on the whole (negative enclosed volume) is turned round."""
    np = _np()
    n = vector_areas / np.linalg.norm(vector_areas, axis=1)[:, None]
    volume = (centroids * vector_areas).sum() / 3.0
    return (n, False) if volume >= 0.0 else (-n, True)


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
    method = s.text("method", "newtonian").strip().lower()
    impact = s.text("impact", "wedge").strip().lower()
    export = s.integer("export_pressure", 0) != 0
    expansion = s.text("expansion", "prandtl_meyer").strip().lower()
    separated = s.text("separated", "freestream").strip().lower()
    if method not in ("newtonian", "tangent_wedge"):
        raise DriverError(f"method = {method}: expected newtonian or tangent_wedge")
    if impact not in ("wedge", "cone"):
        raise DriverError(f"impact = {impact}: expected wedge or cone")
    if expansion not in ("prandtl_meyer", "shadow") or separated not in ("freestream", "gabeaud"):
        raise DriverError(f"expansion = {expansion}, separated = {separated}: expected prandtl_meyer or shadow, and "
                          "freestream or gabeaud")
    if method == "tangent_wedge":
        sub = sorted({r["Mach"] for r in conds if r["Mach"] <= 1.0})
        if sub:
            raise DriverError(f"the tangent wedge needs a supersonic freestream: the flight condition has Mach "
                              f"{', '.join(f'{m:g}' for m in sub)} in [Mach Min, Mach Max); set Mach Min above 1")
    stale = base + ".ntvsp.mesh_cp_max.csv"
    if mode != "user" and os.path.isfile(stale):  # left by an earlier User Defined run of this case
        os.remove(stale)
    c_mm, s_mm2 = faces(base + ".ntvsp.mesh.stl")
    cp_max_field, gradient = read_values(base, c_mm, mode == "user")
    centroids, vector_areas = c_mm * MM, s_mm2 * MM * MM
    _, flipped = face_normals(vector_areas, centroids)
    if flipped:
        vector_areas = -vector_areas
        log("note: the mesh's faces point inward (negative enclosed volume); they are turned round")
    normals = unit_normals(gradient, vector_areas, log)
    closure = np.linalg.norm(vector_areas.sum(axis=0)) / np.linalg.norm(vector_areas, axis=1).sum()
    law = ("Modified Newtonian" if method == "newtonian" else
           f"tangent {impact} + {'Prandtl-Meyer' if expansion == 'prandtl_meyer' else 'Newtonian shadow'}"
           + (f" (separated faces: {separated})" if expansion == "prandtl_meyer" else ""))
    log(f"{law} on mesh: {len(c_mm)} faces, wetted area {np.linalg.norm(vector_areas, axis=1).sum():.6g} "
        f"m^2, closure {closure:.1e}; {len(conds)} of {len(rows)} conditions in Mach [{lo:g}, {hi:g}); Cp_max {mode}")
    if closure > 1e-3:
        log("note: the mesh is not closed (its vector areas do not sum to zero)")
    axes = frame([s.real(f"aft_{c}") for c in "xyz"], [s.real(f"up_{c}") for c in "xyz"])
    removed = pressure.clear(base)
    on_cp = None
    if export:
        written = []

        def on_cp(c, cp):  # every flight condition (any altitude) at this Mach, alpha, beta: p = p_inf + q Cp
            for row in pressure.cases_at(fl, c["Mach"], c["AoA"], c["Beta"]):
                pressure.write(base, row["Case"], centroids, row["Pressure"] + row["Dynamic_Pressure"] * cp)
                written.append(int(row["Case"]))
    out = coefficients(centroids, vector_areas, normals, cp_max_field, conds, mode, axes, fl.ref, method=method,
                       expansion=expansion, separated=separated, impact=impact, on_cp=on_cp)
    if export:
        log(f"pressures (Pa) at the {len(centroids)} face centroids of flight condition cases "
            f"{', '.join(str(n) for n in sorted(written))}: {os.path.basename(base)}.ntvsp.pressure.<Case>.csv")
    elif removed:
        log(f"removed {removed} pressure files of an earlier run (Export Pressures is off)")
    write_csv(base + ".ntvsp.coefficients.csv", HEADER, out)
    log("done")
    log.close()


def _faces_main():
    if len(sys.argv) != 4:
        sys.stderr.write("usage: python -m vsp_ntop.run_newtonian_mesh faces STL_FILE FACES_CSV\n")
        sys.exit(2)
    conda_dll_dirs()
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
