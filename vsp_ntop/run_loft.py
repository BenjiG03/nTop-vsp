"""Lofted CAD: STEP solids fitted straight to nTop's slices of each body, without OpenVSP's skinning.

    python -m vsp_ntop.run_loft SETTINGS_FILE
"""
from __future__ import annotations

import os

from . import loft as L
from . import step_write as W
from .common import DriverError, Log, Settings, main_wrapper, write_csv
from .run_model import contours, frame_of, num, parse_components

EDGE_NAMES = ("right", "bottom", "left", "top")


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    nw = int(s.real("points_around", 48))
    if nw % 4 or nw < 8:
        raise DriverError(f"Points Around = {nw}: give a multiple of 4, at least 8")
    faces = int(s.real("faces_around", 4))
    if faces % 4 or faces < 4 or nw % faces:
        raise DriverError(f"Faces Around = {faces}: give a multiple of 4 that divides Points Around ({nw})")
    per_ctrl = s.real("rows_per_ctrl", 2.0)
    if per_ctrl < 1.0:
        raise DriverError("Slices per Control Point must be at least 1")
    edges = s.text("sharp_edges", "auto")
    if edges not in ("auto", "none"):
        raise DriverError(f"sharp_edges '{edges}': expected auto or none")
    comps = [c for c in parse_components(settings_path)]
    if not comps:
        raise DriverError("no components: connect OpenVSP Fuselage blocks to Components")
    bodies, rows = [], []
    for c in comps:
        name = c.get("name", "Body")
        if c.get("kind") != "fuselage":
            raise DriverError(f"{name}: Lofted CAD takes OpenVSP Fuselage components only (got {c.get('kind')})")
        frame = frame_of(c)
        layers = contours(c, frame)
        bounds = None
        if "bx0" in c and "bx1" in c:
            bounds = (num(c, "bx0") - frame.origin[0], num(c, "bx1") - frame.origin[0])
        try:
            lf = L.loft_body(layers, frame, nw=nw, rows_per_ctrl=per_ctrl, bounds=bounds,
                             creases="auto" if edges == "auto" else (), faces=faces)
        except L.F.FitError as e:
            raise DriverError(f"{name}: {e}")
        for n in lf.notes:
            log(f"{name}: {n}")
        ends = "/".join("pointed" if e.pointed else "blunt" for e in (lf.start, lf.end))
        around = nw + 2 * len(lf.creases)  # a sharp edge's triple knot adds two control points
        log(f"{name}: surface through {lf.rows - 2} slices, {lf.patches[0].ctrl.shape[0]} x {around} control points "
            f"(along x around); contour points off the surface: rms "
            f"{lf.residual[0] * 1e3:.3f} mm, largest {lf.residual[1] * 1e3:.3f} mm; ends {ends}")
        bodies.append((lf, name))
        rows.append({"Name": name, "Slices": lf.rows - 2, "Control_Rows": lf.patches[0].ctrl.shape[0],
                     "Points_Around": around,
                     "Residual_RMS_mm": lf.residual[0] * 1e3, "Residual_Max_mm": lf.residual[1] * 1e3,
                     "Start": "pointed" if lf.start.pointed else "blunt", "End": "pointed" if lf.end.pointed else "blunt",
                     "Sharp_Edges": " ".join(EDGE_NAMES[q] for q in lf.creases) or "none"})
    out = os.path.join(work, case + ".stp")
    info = W.write_lofts(out, bodies, case)
    log(f"written {out}: {info['solids']} solid(s), {info['faces']} faces ({info['caps']} flat caps), millimetres")
    write_csv(os.path.join(work, case + ".ntvsp.loft.csv"), list(rows[0]), [list(r.values()) for r in rows])
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
