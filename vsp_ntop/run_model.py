"""Build an OpenVSP model from nTop implicit bodies, driven from nTop.

    python -m vsp_ntop.run_model SETTINGS_FILE
"""
from __future__ import annotations

import math
import os


from . import duct as D
from . import fit as F
from . import slices as S
from . import vspmodel as M
from .common import openvsp_dir, DriverError, Log, Settings, main_wrapper, run_vspscript, vsp_path, write_csv

X, Y, Z = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)


def dot(a, b):
    return sum(p * q for p, q in zip(a, b))


PER_BODY = ("file", "ox", "oy", "oz", "sx", "sy", "sz", "bx0", "bx1")

CARRY_COLS = ["Name", "Copy", "Body", "Sx", "Sy", "Sz", "Nx", "Ny", "Nz", "Qx", "Qy", "Qz", "Root", "Semispan",
              "Root_Chord", "Tip_Chord", "LE_Slope", "Planform", "Body_Width", "Body_Length"]
FLIPS = {"none": [(1, 1)], "xz": [(1, 1), (-1, 1)], "xy": [(1, 1), (1, -1)],
         "xzxy": [(1, 1), (-1, 1), (1, -1), (-1, -1)]}


def carryover_fins(name, frame, z0, stations, symmetry, fuse_fits):
    """The geometry TR 1307's carry-over needs for each copy of a wing (``carryover.py``), from its exposed stations
    (before any extension into the body) and the fuselage it sits on: root and tip distances from the plane through
    the origin normal to the span axis, root and tip chords, the cotangent of the leading-edge sweep, the planform
    area, the body's width across the span direction at the root and its length aft of the root leading edge."""
    if not fuse_fits:
        raise DriverError(f"{name}: Force Accounting with carry-over needs the body it sits on as an OpenVSP Fuselage "
                          "in the model (in any Set, None included)")
    root, tip = stations[0], stations[-1]
    r, semispan = root.z - z0, tip.z - z0
    x_le = frame.origin[0] + root.le[1]
    x_tip = frame.origin[0] + tip.le[1]
    slope = (semispan - r) / max(x_tip - x_le, 1e-3 * (semispan - r))
    planform = sum(0.5 * (a.chord + b.chord) * (b.z - a.z) for a, b in zip(stations, stations[1:]))
    rows = []
    for copy, (fy, fz) in enumerate(FLIPS[symmetry], start=1):
        m = lambda v: (v[0], fy * v[1], fz * v[2])  # noqa: E731
        s, n = m(frame.ez), m(frame.ex)
        q = m(frame.point(root.le[0], root.le[1], root.z))
        mid = m(frame.point(root.le[0], root.le[1] + 0.5 * root.chord, root.z))
        best = None
        for body, ff in fuse_fits.items():
            secs = [sec for sec in ff.sections if sec.loop]
            pts = [[ff.frame.point(u, v, sec.x) for u, v in sec.loop] for sec in secs]
            xs = [ff.frame.origin[0] + sec.x for sec in secs]
            near = min(range(len(secs)), key=lambda i: abs(xs[i] - mid[0]))
            gap = min(math.dist(p, mid) for p in pts[near])
            if best is None or gap < best[0]:
                over = [i for i, x in enumerate(xs) if x_le <= x <= x_le + root.chord] or [near]
                # across the body in the fin's plane
                band = [dot(p, s) for i in over for p in pts[i] if abs(dot(p, n) - dot(q, n)) <= 0.2 * r]
                low = min(band or [dot(p, s) for i in over for p in pts[i]])
                end = max(ff.frame.origin[0] + sec.x for sec in ff.sections)
                best = (gap, body, low, end)
        _, body, low, end = best
        rows.append({"Name": name, "Copy": copy, "Body": body, **{f"S{a}": v for a, v in zip("xyz", s)},
                     **{f"N{a}": v for a, v in zip("xyz", n)}, **{f"Q{a}": v for a, v in zip("xyz", q)},
                     "Root": r, "Semispan": semispan, "Root_Chord": root.chord, "Tip_Chord": tip.chord,
                     "LE_Slope": slope, "Planform": planform, "Body_Width": r - low, "Body_Length": end - x_le})
    return rows


def parse_components(path):
    """``component = key: value; key: value; ...`` lines, in order."""
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
        if c.get("kind") == "duct":  # one body, per-section lists
            for need in ("name", "cuts"):
                if need not in c:
                    raise DriverError(f"duct line without '{need}': {line.strip()}")
            out.append(c)
            continue
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


def _vec(c, key):
    v = [x.strip() for x in c[key].split(",")]
    if len(v) != 3:
        raise DriverError(f"component {c.get('name')}: {key} = '{c[key]}' is not x, y, z")
    return tuple(float(x) for x in v)


def duct_cuts(c):
    """The chosen contour of each cut of a duct, in world coordinates, with its cutting plane. Each cut is a CLI
    file sliced in its own frame (origin p, x axis u, y axis v, normal u x v); the layer nearest the frame origin is
    the cut."""
    files = [x.strip() for x in c["cuts"].split("|") if x.strip()]
    per = {k: [float(x) for x in c[k].split("|")] for k in ("px", "py", "pz", "ux", "uy", "uz", "vx", "vy", "vz")}
    for k, v in per.items():
        if len(v) != len(files):
            raise DriverError(f"duct {c['name']}: {len(v)} values of {k} for {len(files)} cuts")
    wall = c.get("wall", "outer")
    cuts = []
    for i, path in enumerate(files):
        if not os.path.isfile(path):
            raise DriverError(f"duct {c['name']}: section file {path} not found")
        o = (per["px"][i], per["py"][i], per["pz"][i])
        ex = F.unit((per["ux"][i], per["uy"][i], per["uz"][i]))
        ey = F.unit((per["vx"][i], per["vy"][i], per["vz"][i]))
        n = F.cross(ex, ey)
        bu, bv = slicer_basis(n, ex, ey)
        layers = [l for l in S.read_cli(path) if l.loops]
        if not layers:
            raise DriverError(f"duct {c['name']}: cut {i + 1} of {len(files)} is empty (the cutting plane misses the "
                              "body: do the inlet and outlet planes lie on its ends?)")
        layer = min(layers, key=lambda l: abs(l.z))
        try:
            loop = D.pick_loop(layer.loops, wall)
        except F.FitError as e:
            raise DriverError(f"duct {c['name']}, cut {i + 1}: {e}")
        pts = [tuple(o[j] + u * bu[j] + v * bv[j] + layer.z * n[j] for j in range(3)) for u, v in loop]
        cuts.append((pts, (o, n)))
    return cuts


def duct_deviation(points, f):
    """Distance from OpenVSP's surface at each section (grid rows at the section parameters) to that section's nTop
    contour, in 3D."""
    n_u = max((a for a, *_ in points), default=0)
    per = n_u // max(1, len(f.sections) - 1) if n_u else 1
    out = []
    for a, b, x, y, z in points:
        if a % per:
            continue
        loop = f.sections[a // per].loop3d
        loop = loop[::max(1, len(loop) // 720)]  # at most ~720 points
        q = (x, y, z)
        best = 1e300
        for i in range(len(loop)):
            p0, p1 = loop[i - 1], loop[i]
            d = tuple(b1 - b0 for b0, b1 in zip(p0, p1))
            l2 = sum(c * c for c in d)
            t = 0.0 if l2 == 0 else max(0.0, min(1.0, sum((qq - pp) * dd for qq, pp, dd in zip(q, p0, d)) / l2))
            best = min(best, math.dist(q, tuple(pp + t * dd for pp, dd in zip(p0, d))))
        out.append(best)
    if not out:
        return float("nan"), float("nan"), 0
    return sum(out) / len(out), max(out), len(out)


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
    nearest slices), over points between the first and last slice."""
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
        if h <= zs[0] + 1e-7 or h >= zs[-1] - 1e-7:  # outside the slices
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
    for c in comps:  # names must differ in OpenVSP
        n = M.clean_name(c["name"])
        seen[n] = seen.get(n, 0) + 1
        if seen[n] > 1:
            log(f"note: a second component named {n} is called {n}_{seen[n]}")
            c["name"] = f"{n}_{seen[n]}"

    lines, fuse_fits, fuse_tess, checks, rows, wings, carry = [], {}, {}, [], [], [], []
    duct_fits, duct_checks = {}, []
    for k, c in enumerate(comps):
        name, kind = M.clean_name(c["name"]), c["kind"]
        set_name = c.get("set", "thick")
        if set_name not in ("thick", "thin", "none"):
            raise DriverError(f"{name}: set '{set_name}': expected thick, thin or none")
        if kind == "duct":
            try:
                ends = ((_vec(c, "inlet_p"), _vec(c, "inlet_n")), (_vec(c, "outlet_p"), _vec(c, "outlet_n")))
                cuts = duct_cuts(c)
                f = D.fit_duct(cuts, ends, int(num(c, "per_quarter", 24)))
            except F.FitError as e:
                raise DriverError(f"{name}: {e}")
            tess_w = int(num(c, "tess_w", 33))
            lines += D.duct_lines(k, name, f, set_name, tess_w, M._sets)
            duct_fits[name] = (f, tess_w, num(c, "panel_aspect", 1.5))
            duct_checks.append(name)
            length = sum(math.dist(a.centre, b.centre) for a, b in zip(f.sections, f.sections[1:]))
            rows.append({"Name": name, "Type": "Duct", "Set": set_name.capitalize(), "Sections": len(f.sections),
                         "Slices": len(cuts), "Size": length})
            log(f"{name}: duct of {len(f.sections)} sections ({c.get('wall', 'outer')} wall), centre line "
                f"{length * 1000:.1f} mm, open at both ends")
            continue
        frame = frame_of(c)
        layers = contours(c, frame)
        try:
            if kind == "fuselage":
                if abs(dot(frame.ez, X) - 1.0) > 1e-9:
                    raise DriverError(f"{name}: a fuselage is sliced along +x")
                bounds = None
                if "bx0" in c and "bx1" in c:  # the bounding box's x range
                    bounds = (num(c, "bx0") - frame.origin[0], num(c, "bx1") - frame.origin[0])
                f = F.fit_fuselage(layers, frame, int(num(c, "sections", 24)), num(c, "nose_cluster", 0.5),
                                   num(c, "tail_cluster", 0.5), bounds=bounds)
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
                exposed_stations = list(f.stations)
                extension = num(c, "extend", 0.0)
                if extension < 0.0:
                    raise DriverError(f"{name}: Extension {extension:g} m is negative")
                if extension > 0.0:
                    if extension >= root_span - 1e-9:
                        raise DriverError(
                            f"{name}: an Extension of {extension * 1000:.2f} mm reaches the plane through the origin "
                            f"normal to the span axis (the root is {root_span * 1000:.2f} mm from it): stop short of it "
                            "(wings that meet or cross there make VSPAERO fail)")
                    F.extend_wing(f, f.stations[0].z - extension)
                accounting = c.get("accounting", "whole")
                if accounting not in ("whole", "exposed", "carryover"):
                    raise DriverError(f"{name}: accounting '{accounting}': expected whole, exposed or carryover")
                if accounting == "carryover":
                    carry.append((name, frame, z0, exposed_stations, symmetry))
                st0 = f.stations[0]
                wings.append({"Name": name, "Accounting": accounting, "Extended": num(c, "extend", 0.0),
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
        if duct_fits:
            M.set_duct_sections(stage2, stage2, duct_fits)
    except (ValueError, OSError) as e:
        raise DriverError(str(e))
    model = os.path.join(work, case + ".vsp3")
    stl = os.path.join(work, case + ".stl")
    for p in (model, stl):
        if os.path.exists(p):
            os.remove(p)
    text = run_vspscript(vsp_dir, M.finish_script(vsp_path(stage2), vsp_path(model), vsp_path(stl),
                                                  [(n, nu, nw) for n, _, _, _, nu, nw in checks]
                                                  + [(n, (len(duct_fits[n][0].sections) - 1) * 4, 48)
                                                     for n in duct_checks]),
                         work, case + ".finish", log, timeout)
    if not os.path.isfile(model):
        raise DriverError("OpenVSP wrote no model file")
    pts = {}
    for line in text.splitlines():
        v = line.split()
        if len(v) == 7 and v[0] == "NTVSP_PT":
            pts.setdefault(v[1], []).append((int(v[2]), int(v[3]), float(v[4]), float(v[5]), float(v[6])))
    by_name = {r["Name"]: r for r in rows}
    for name in duct_checks:
        mean, worst, n = duct_deviation(pts.get(name, []), duct_fits[name][0])
        by_name[name].update({"Mean_Deviation": mean, "Max_Deviation": worst, "Check_Points": n,
                              "Max_Deviation_Fraction": worst / by_name[name]["Size"]})
        log(f"{name}: OpenVSP section outlines vs nTop cuts: mean {mean * 1000:.3f} mm, largest {worst * 1000:.3f} mm "
            f"over {n} points")
    for row, (name, layers, frame, skip, _, _) in zip([r for r in rows if r["Type"] != "Duct"], checks):
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
    fins = [row for item in carry for row in carryover_fins(*item, fuse_fits)]
    for r in fins:
        log(f"{r['Name']} copy {r['Copy']} (TR 1307 carry-over): root {r['Root'] * 1000:.1f} mm and tip "
            f"{r['Semispan'] * 1000:.1f} mm from the origin plane, chords {r['Root_Chord'] * 1000:.1f} / "
            f"{r['Tip_Chord'] * 1000:.1f} mm, body {r['Body_Width'] * 1000:.1f} mm across, running "
            f"{r['Body_Length'] * 1000:.1f} mm aft of the root leading edge ({r['Body']})")
    write_csv(os.path.join(work, case + ".ntvsp.carryover.csv"), CARRY_COLS, [[r[c] for c in CARRY_COLS] for r in fins])
    for p in (stage1, stage2):
        os.remove(p)
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
