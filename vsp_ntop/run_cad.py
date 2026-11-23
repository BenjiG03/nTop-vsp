"""Export an OpenVSP model as CAD (STEP or IGES), driven from nTop.

    python -m vsp_ntop.run_cad SETTINGS_FILE

OpenVSP 3.56 has three ways to write its surfaces to CAD; ``method`` picks one:

* ``trimmed`` (OpenVSP's Surface Intersection): the components are intersected with each other, the parts inside
  other components are trimmed away, and the result is written as trimmed surfaces. With ``representation = brep``
  (STEP only) closed surfaces become one BREP solid (``MANIFOLD_SOLID_BREP``); ``shell`` writes a shell model. Open
  surfaces (a duct's wall) stay surfaces.
* ``stitched`` (File > Export, "split and stitch"): each component's surface patches are split along their feature
  lines, degenerate patches are dropped, and the patches are stitched back into one watertight BREP per closed
  component. Components are neither intersected nor trimmed against each other.
* ``surfaces`` (File > Export, STEP or IGES): each component's untrimmed surface patches, as in OpenVSP's plain
  STEP/IGES export (surfaces only, no topology).

Units: OpenVSP's CAD length unit does not scale anything; it declares the unit of the model's numbers (a 10-long
pod written with LEN_M reads as 10 m, with the default LEN_FT as 10 ft). OpenVSP Model writes metres, so the driver
always declares metres and the CAD program converts on import. The plain export takes the unit from OpenVSP's
STEP/IGES settings, the stitched and trimmed ones from its Surface Intersection settings (found by testing: the
stitched export ignores the STEP settings' unit), so the driver sets both. ``units`` may name another unit only for
a .vsp3 built in it.

Writes ``<case>.stp`` or ``<case>.igs`` and ``<case>.ntvsp.cad.csv`` (one row: method, format, units, the file's
entity counts read back from it: solids, shells, faces, B-spline surfaces). With ``planar_caps`` (STEP), every flat
end cap's patches are replaced by one plane face (``step_caps.py``; OpenVSP's own file is kept as
``<case>.openvsp.stp``). Standard library only.
"""
from __future__ import annotations

import os
import re

from . import step_caps
from .common import (openvsp_dir, DriverError, Log, Settings, main_wrapper, run_vspscript, vsp_path, write_csv)

UNITS = {"mm": "LEN_MM", "cm": "LEN_CM", "m": "LEN_M", "in": "LEN_IN", "ft": "LEN_FT"}
SETS = {"all": "SET_ALL", "thick": "3", "thin": "4"}


def script(s: Settings, model, out):
    method, fmt = s.text("method", "trimmed"), s.text("format", "step")
    unit = UNITS[s.text("units", "m")]
    tol = s.real("tolerance", 1e-6)
    cubic = 1 if s.flag("to_cubic", False) else 0
    cubic_tol = s.real("to_cubic_tolerance", 1e-6)
    merge = 1 if s.flag("merge_points", False) else 0
    lete = 1 if s.flag("merge_le_te", False) else 0
    label = 1 if s.flag("label_names", True) else 0
    geom_set = SETS[s.text("geom_set", "all")]
    L = ["void main()", "{", f'    ReadVSPFile( "{vsp_path(model)}" );',
         '    if ( FindGeoms().size() == 0 ) { Print( "ERR the model has no components" ); return; }', "    Update();",
         '    string veh = FindContainer( "Vehicle", 0 );',
         '    string si = FindContainer( "SurfaceIntersectSettings", 0 );']
    for grp in ("STEPSettings", "IGESSettings"):
        L += [f'    SetParmVal( FindParm( veh, "LenUnit", "{grp}" ), {unit} );',
              f'    SetParmVal( FindParm( veh, "ToCubic", "{grp}" ), {cubic} );',
              f'    SetParmVal( FindParm( veh, "ToCubicTol", "{grp}" ), {cubic_tol!r} );',
              f'    SetParmVal( FindParm( veh, "MergeLETE", "{grp}" ), {lete} );',
              f'    SetParmVal( FindParm( veh, "LabelName", "{grp}" ), {label} );']
    L += [f'    SetParmVal( FindParm( veh, "Tolerance", "STEPSettings" ), {tol!r} );',
          f'    SetParmVal( FindParm( veh, "MergePoints", "STEPSettings" ), {merge} );',
          f'    SetParmVal( FindParm( si, "CADLenUnit", "ExportIntersect" ), {unit} );',
          f'    SetParmVal( FindParm( si, "STEPTol", "ExportIntersect" ), {tol!r} );',
          f'    SetParmVal( FindParm( si, "STEPMergePoints", "ExportIntersect" ), {merge} );',
          f'    SetParmVal( FindParm( si, "CADLabelName", "ExportIntersect" ), {label} );']
    if method == "trimmed":
        rep = "STEP_BREP" if s.text("representation", "brep") == "brep" else "STEP_SHELL"
        step, iges = (1, 0) if fmt == "step" else (0, 1)
        L += ['    SetAnalysisInputDefaults( "SurfaceIntersection" );',
              f'    SetIntAnalysisInput( "SurfaceIntersection", "SelectedSetIndex", array<int> = {{{geom_set}}} );',
              f'    SetIntAnalysisInput( "SurfaceIntersection", "CADLenUnit", array<int> = {{{unit}}} );',
              f'    SetDoubleAnalysisInput( "SurfaceIntersection", "STEPTol", array<double> = {{{tol!r}}} );',
              f'    SetIntAnalysisInput( "SurfaceIntersection", "STEPRepresentation", array<int> = {{{rep}}} );',
              f'    SetIntAnalysisInput( "SurfaceIntersection", "CADLabelName", array<int> = {{{label}}} );',
              f'    SetIntAnalysisInput( "SurfaceIntersection", "STEPFileFlag", array<int> = {{{step}}} );',
              f'    SetIntAnalysisInput( "SurfaceIntersection", "IGESFileFlag", array<int> = {{{iges}}} );',
              f'    SetStringAnalysisInput( "SurfaceIntersection", "STEPFileName", array<string> = {{"{vsp_path(out)}"}} );',
              f'    SetStringAnalysisInput( "SurfaceIntersection", "IGESFileName", array<string> = {{"{vsp_path(out)}"}} );',
              '    ExecAnalysis( "SurfaceIntersection" );']
    else:
        kind = {("stitched", "step"): "EXPORT_STEP_STITCH", ("stitched", "iges"): "EXPORT_IGES_STITCH",
                ("surfaces", "step"): "EXPORT_STEP", ("surfaces", "iges"): "EXPORT_IGES"}[(method, fmt)]
        L.append(f'    ExportFile( "{vsp_path(out)}", {geom_set}, {kind} );')
    L += ["    while ( GetNumTotalErrors() > 0 ) { ErrorObj e = PopLastError(); Print( \"VSPERR \" + e.GetErrorString() ); }",
          '    Print( "NTVSP_OK" );', "}"]
    return "\n".join(L) + "\n"


def step_counts(path):
    """Entity counts of a STEP file (its text)."""
    text = open(path, encoding="utf-8", errors="replace").read()
    def n(name):
        return len(re.findall(r"=\s*" + name + r"\s*\(", text))
    return {"Solids": n("MANIFOLD_SOLID_BREP"), "Closed_Shells": n("CLOSED_SHELL"), "Open_Shells": n("OPEN_SHELL"),
            "Faces": n("ADVANCED_FACE"),
            # simple or complex (rational) entities
            "Surfaces": len(re.findall(r"B_SPLINE_SURFACE_WITH_KNOTS\s*\(", text))}


def iges_counts(path):
    """Entity counts of an IGES file (directory entries: 128 B-spline surface, 144 trimmed surface, 186 solid,
    514 shell, 510 face)."""
    counts = {}
    for line in open(path, encoding="ascii", errors="replace"):
        if len(line) >= 73 and line[72] == "D":
            try:
                t = int(line[:8])
            except ValueError:
                continue
            counts[t] = counts.get(t, 0) + 1
    half = lambda t: counts.get(t, 0) // 2  # noqa: E731  (two directory lines per entity)
    return {"Solids": half(186), "Closed_Shells": half(514), "Open_Shells": 0, "Faces": half(510) or half(144),
            "Surfaces": half(128)}


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = openvsp_dir(s)
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    model = s.text("model")
    if not model.lower().endswith(".vsp3"):
        model += ".vsp3"
    if not os.path.isfile(model):
        raise DriverError(f"OpenVSP model {model} not found")
    method, fmt, units = s.text("method", "trimmed"), s.text("format", "step"), s.text("units", "m")
    if method not in ("trimmed", "stitched", "surfaces"):
        raise DriverError(f"method '{method}': expected trimmed, stitched or surfaces")
    if fmt not in ("step", "iges"):
        raise DriverError(f"format '{fmt}': expected step or iges")
    if units not in UNITS:
        raise DriverError(f"units '{units}': expected one of {', '.join(UNITS)}")
    out = os.path.join(work, case + (".stp" if fmt == "step" else ".igs"))
    if os.path.exists(out):
        os.remove(out)
    log(f"CAD export: {method} {fmt.upper()} in {units} from {model}")
    text = run_vspscript(vsp_dir, script(s, model, out), work, case + ".cad", log, s.real("timeout_s", 7200))
    for line in text.splitlines():
        if "VSPERR" in line:
            log("OpenVSP: " + line.strip())
    if not os.path.isfile(out) or os.path.getsize(out) == 0:
        raise DriverError(f"OpenVSP wrote no {fmt.upper()} file (see {case}.cad output in the working directory; the "
                          "trimmed method needs components in the chosen set)")
    if s.flag("planar_caps", False):
        if fmt != "step":
            raise DriverError("Planar End Caps works on STEP files only: choose Format STEP")
        raw = os.path.join(work, case + ".openvsp.stp")
        os.replace(out, raw)
        caps = step_caps.planar_caps_file(raw, out, log)
        log(f"planar end caps: {caps} flat cap(s) rebuilt as one plane face each (OpenVSP's own file kept as "
            f"{os.path.basename(raw)})")
    counts = step_counts(out) if fmt == "step" else iges_counts(out)
    log("written " + out + ": " + ", ".join(f"{k.replace('_', ' ').lower()} {v}" for k, v in counts.items()))
    cols = ["File", "Method", "Format", "Units"] + list(counts)
    write_csv(os.path.join(work, case + ".ntvsp.cad.csv"), cols,
              [[out.replace("\\", "/"), method, fmt, units] + list(counts.values())])
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
