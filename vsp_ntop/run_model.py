"""Build an OpenVSP model from nTop implicit bodies, driven from nTop.

    python -m vsp_ntop.run_model SETTINGS_FILE

Each component block (OpenVSP Fuselage, OpenVSP Wing) slices its implicit body in nTop (Slice Body) and exports
the section contours as an ASCII CLI file; its settings arrive here as one ``component = key: value; ...`` line.
This driver

1. reads each component's contours and fits OpenVSP sections to them (``fit.py``): a Fuselage's closed
   cross-sections along x, or a Wing's airfoils, chords, twist, sweep and dihedral along its span axis;
2. builds the components in OpenVSP (``vspmodel.py``), each in the Thick or Thin set, with its tessellation and
   clustering, and saves ``<case>.vsp3``;
3. checks the model against the contours: a grid of points on each OpenVSP surface is compared with the nTop
   section contours, and the mean and largest distances go to ``<case>.ntvsp.model.csv``;
4. exports OpenVSP's tessellation as ``<case>.stl`` (metres) to show in nTop.

Lengths are metres. Standard library only.
"""
from __future__ import annotations

import math
import os


from . import fit as F
from . import slices as S
from . import vspmodel as M
from .common import openvsp_dir, DriverError, Log, Settings, main_wrapper, run_vspscript, vsp_path, write_csv

X, Y, Z = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)


def dot(a, b):
    return sum(p * q for p, q in zip(a, b))


PER_BODY = ("file", "ox", "oy", "oz", "sx", "sy", "sz")


def parse_components(path):
    """``component = key: value; key: value; ...`` lines, in order. A component block given several bodies writes
    one line for all of them: ``names`` and the per-body values (``files``, ``ox`` ... ``sz``) are ``|``-separated
    lists, expanded here into one component per body. One name for several bodies numbers them (Name_1, Name_2...)."""
    out = []
    for line in open(path, encoding="utf-8"):
        key, _, value = line.partition("=")
        if key.strip() != "component" or not value.strip():
            continue
        c = {}
        for item in value.split(";"):
            k, sep, v = item.partition(":")
            if sep:
                c[k.strip()] = v.strip()
        if "files" in c:
            files = [x.strip() for x in c.pop("files").split("|") if x.strip()]
            names = [x.strip() for x in c.pop("names", "").split("|") if x.strip()] or [c.get("kind", "Component").capitalize()]
            if len(names) == 1 and len(files) > 1:
                names = [f"{names[0]}_{i + 1}" for i in range(len(files))]
            if len(names) != len(files):
                raise DriverError(f"{len(names)} names for {len(files)} bodies: give one name per body, or one for all")
            lists = {k: [x.strip() for x in c.pop(k).split("|")] for k in PER_BODY[1:] if k in c}
            for k, v in lists.items():
                if len(v) != len(files):
                    raise DriverError(f"component {names[0]}: {len(v)} values of {k} for {len(files)} bodies")
            for i, (n, f) in enumerate(zip(names, files)):
                out.append(dict(c, name=n, file=f, **{k: v[i] for k, v in lists.items()}))
            continue
        for need in ("kind", "name", "file"):
            if need not in c:
                raise DriverError(f"component line without '{need}': {line.strip()}")
        out.append(c)
    return out


def num(c, key, default=None):
    if key not in c:
        if default is None:
            raise DriverError(f"component {c.get('name')}: '{key}' missing")
        return default
    try:
        return float(c[key])
    except ValueError:
        raise DriverError(f"component {c.get('name')}: {key} = '{c[key]}' is not a number")


def slicer_basis(d, ex, ey):
    """The in-plane axes nTop's Slice Body writes contours in: the frame's own, except for slicing along a world
    axis, where it uses world axes (u, v) with u x v the slicing direction: the cyclic ones for +x, +y, +z and
    (z, y), (x, z), (y, x) for -x, -y, -z (found by testing nTop 6.3.3; the four fins of a cruciform rocket)."""
    neg = lambda a: tuple(-c for c in a)  # noqa: E731
    for axis, basis in ((X, (Y, Z)), (Y, (Z, X)), (Z, (X, Y)), (neg(X), (Z, Y)), (neg(Y), (X, Z)), (neg(Z), (Y, X))):
        if dot(d, axis) > 1.0 - 1e-9:
            return basis
    return ex, ey


def contours(c, frame: F.Frame):
    """The component's layers, in the coordinates of ``frame`` (whose origin and slicing direction the block
    used)."""
    path = c["file"]
    if not os.path.isfile(path):
        raise DriverError(f"component {c['name']}: section file {path} not found")
    layers = S.read_cli(path)
    bu, bv = slicer_basis(frame.ez, frame.ex, frame.ey)
    m = [[dot(bu, frame.ex), dot(bv, frame.ex)], [dot(bu, frame.ey), dot(bv, frame.ey)]]
    if not (abs(m[0][0] - 1) < 1e-9 and abs(m[1][1] - 1) < 1e-9):
        for l in layers:
            l.loops = [[(m[0][0] * u + m[0][1] * v, m[1][0] * u + m[1][1] * v) for u, v in loop] for loop in l.loops]
    return [l for l in layers if l.loops]


def frame_of(c):
    o = (num(c, "ox"), num(c, "oy"), num(c, "oz"))
    if c["kind"] == "fuselage":
        return F.fuselage_frame(o)
    return F.wing_frame(o, (num(c, "sx"), num(c, "sy"), num(c, "sz")))


def _dist_loop(loop, u, v):
    best = 1e300
    for i in range(len(loop)):
        ax, ay = loop[i - 1]
        bx, by = loop[i]
        dx, dy = bx - ax, by - ay
        l2 = dx * dx + dy * dy
        t = 0.0 if l2 == 0 else max(0.0, min(1.0, ((u - ax) * dx + (v - ay) * dy) / l2))
        best = min(best, (u - ax - t * dx) ** 2 + (v - ay - t * dy) ** 2)
    return math.sqrt(best)


def deviation(points, layers, frame: F.Frame, skip_ends):
    """Mean and largest distance from OpenVSP surface points to the nTop contours (interpolated between the two
    nearest slices), over points between the first and last slice. ``skip_ends``: fraction of the grid's u range
    left out at each end (a wing's root and tip caps)."""
    loops = []
    for l in layers:
        loop = S.outer_loops(l)[0]
        loops.append((l.z, S.resample(loop, min(len(loop), 240), closed=True)))
    if len(loops) < 2:
        return float("nan"), float("nan"), 0
    zs = [z for z, _ in loops]
    n_u = max(a for a, *_ in points) or 1
    out = []
    for a, b, x, y, z in points:
        if a < skip_ends * n_u - 1e-9 or a > (1 - skip_ends) * n_u + 1e-9:
            continue
        rel = (x - frame.origin[0], y - frame.origin[1], z - frame.origin[2])
        u, v, h = dot(rel, frame.ex), dot(rel, frame.ey), dot(rel, frame.ez)
        if h < zs[0] or h > zs[-1]:
            continue
        j = max(0, min(len(zs) - 2, next((i for i in range(len(zs) - 1) if zs[i + 1] >= h), len(zs) - 2)))
        t = (h - zs[j]) / (zs[j + 1] - zs[j]) if zs[j + 1] > zs[j] else 0.0
        out.append((1 - t) * _dist_loop(loops[j][1], u, v) + t * _dist_loop(loops[j + 1][1], u, v))
    if not out:
        return float("nan"), float("nan"), 0
    return sum(out) / len(out), max(out), len(out)


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = openvsp_dir(s)
    timeout = s.real("timeout_s", 7200)
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    comps = parse_components(settings_path)
    if not comps:
        raise DriverError("no components: connect OpenVSP Fuselage or OpenVSP Wing blocks to Components")
    seen = {}
    for c in comps:  # names must differ in OpenVSP: number repeats
        n = M.clean_name(c["name"])
        seen[n] = seen.get(n, 0) + 1
        if seen[n] > 1:
            log(f"note: a second component named {n} is called {n}_{seen[n]}")
            c["name"] = f"{n}_{seen[n]}"

    lines, fuse_fits, fuse_tess, checks, rows, wings = [], {}, {}, [], [], []
    for k, c in enumerate(comps):
        name, kind = M.clean_name(c["name"]), c["kind"]
        set_name = c.get("set", "thick")
        if set_name not in ("thick", "thin", "none"):
            raise DriverError(f"{name}: set '{set_name}': expected thick, thin or none")
        frame = frame_of(c)
        layers = contours(c, frame)
        try:
            if kind == "fuselage":
                if abs(dot(frame.ez, X) - 1.0) > 1e-9:
                    raise DriverError(f"{name}: a fuselage is sliced along +x")
                f = F.fit_fuselage(layers, frame, int(num(c, "sections", 24)), num(c, "nose_cluster", 0.5),
                                   num(c, "tail_cluster", 0.5))
                lines += M.fuselage_lines(k, name, f, set_name, int(num(c, "tess_w", 33)))
                fuse_fits[name] = f
                fuse_tess[name] = (int(num(c, "tess_w", 33)), num(c, "panel_aspect", 1.5),
                                   num(c, "fwd_tess_cluster", 1.0), num(c, "aft_tess_cluster", 1.0))
                n_sec = len(f.sections)
                checks.append((name, layers, frame, 0.0, 6 * n_sec, 32))
                size = f.sections[-1].x - f.sections[0].x
            elif kind == "wing":
                symmetry = c.get("symmetry", "xz" if int(num(c, "symmetric", 0)) else "none")
                if symmetry not in M.SYMMETRY_FLAGS:
                    raise DriverError(f"{name}: symmetry '{symmetry}': expected none, xz, xy or xzxy")
                if symmetry != "none":
                    o_s = dot(frame.origin, frame.ez)
                    layers = [l for l in layers if o_s + l.z >= -1e-6]
                f = F.fit_wing(layers, frame, int(num(c, "sections", 8)), int(num(c, "airfoil_points", 41)),
                               num(c, "section_tolerance", 0.01), num(c, "min_section_span", 0.04))
                z0 = -dot(frame.origin, frame.ez)  # the plane through the origin normal to the span axis
                root_span = f.stations[0].z - z0  # from that plane to where the wing body starts
                if int(num(c, "extend", 0)):
                    F.extend_wing(f, z0)
                accounting = c.get("accounting", "whole")
                if accounting not in ("whole", "exposed"):
                    raise DriverError(f"{name}: accounting '{accounting}': expected whole or exposed")
                st0 = f.stations[0]
                wings.append({"Name": name, "Accounting": accounting, "Extended": int(num(c, "extend", 0)),
                              "Symmetry": symmetry, "Root_Span": root_span,
                              **{f"S{a}": v for a, v in zip("xyz", frame.ez)},
                              **{f"N{a}": v for a, v in zip("xyz", frame.ex)},
                              **{f"Q{a}": v for a, v in zip("xyz", frame.point(st0.le[0], st0.le[1], st0.z))}})
                lines += M.wing_lines(k, name, f, set_name, symmetry, int(num(c, "tess_u", 20)),
                                      int(num(c, "tess_w", 33)), num(c, "le_cluster", 0.25), num(c, "te_cluster", 0.25),
                                      num(c, "root_cluster", 1.0), num(c, "tip_cluster", 1.0),
                                      bool(int(num(c, "sharp_te", 1))))
                n_sec = len(f.stations) - 1
                checks.append((name, layers, frame, 1.0 / (n_sec + 2), 8 * (n_sec + 2), 32))
                size = max(st.chord for st in f.stations)
            else:
                raise DriverError(f"{name}: unknown kind '{kind}'")
        except F.FitError as e:
            raise DriverError(f"{name}: {e}")
        for note in f.notes:
            log(f"{name}: {note}")
        rows.append({"Name": name, "Type": kind.capitalize(), "Set": set_name.capitalize(),
                     "Sections": n_sec, "Slices": len(layers), "Size": size})

    log(f"OpenVSP: build {len(comps)} components")
    stage1 = os.path.join(work, case + ".build.vsp3")
    run_vspscript(vsp_dir, M.build_script(lines, vsp_path(stage1)), work, case + ".build", log, timeout)
    stage2 = os.path.join(work, case + ".sections.vsp3")
    try:
        M.set_fuselage_sections(stage1, stage2, fuse_fits, fuse_tess)
    except (ValueError, OSError) as e:
        raise DriverError(str(e))
    model = os.path.join(work, case + ".vsp3")
    stl = os.path.join(work, case + ".stl")
    for p in (model, stl):
        if os.path.exists(p):
            os.remove(p)
    text = run_vspscript(vsp_dir, M.finish_script(vsp_path(stage2), vsp_path(model), vsp_path(stl),
                                                  [(n, nu, nw) for n, _, _, _, nu, nw in checks]),
                         work, case + ".finish", log, timeout)
    if not os.path.isfile(model):
        raise DriverError("OpenVSP wrote no model file")
    pts = {}
    for line in text.splitlines():
        v = line.split()
        if len(v) == 7 and v[0] == "NTVSP_PT":
            pts.setdefault(v[1], []).append((int(v[2]), int(v[3]), float(v[4]), float(v[5]), float(v[6])))
    for row, (name, layers, frame, skip, _, _) in zip(rows, checks):
        mean, worst, n = deviation(pts.get(name, []), layers, frame, skip)
        row.update({"Mean_Deviation": mean, "Max_Deviation": worst, "Max_Deviation_Fraction": worst / row["Size"],
                    "Check_Points": n})
        log(f"{name}: OpenVSP surface vs nTop contours: mean {mean * 1000:.2f} mm, largest {worst * 1000:.2f} mm "
            f"({100 * worst / row['Size']:.2f}% of {'length' if row['Type'] == 'Fuselage' else 'root chord'}) over {n} points")
    cols = ["Name", "Type", "Set", "Sections", "Slices", "Size", "Mean_Deviation", "Max_Deviation",
            "Max_Deviation_Fraction", "Check_Points"]
    write_csv(os.path.join(work, case + ".ntvsp.model.csv"), cols, [[r[c] for c in cols] for r in rows])
    wing_cols = ["Name", "Accounting", "Extended", "Symmetry", "Root_Span", "Sx", "Sy", "Sz", "Nx", "Ny", "Nz", "Qx",
                 "Qy", "Qz"]
    write_csv(os.path.join(work, case + ".ntvsp.wings.csv"), wing_cols, [[w[c] for c in wing_cols] for w in wings])
    for p in (stage1, stage2):
        os.remove(p)
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
