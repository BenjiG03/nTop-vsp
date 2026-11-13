"""OpenVSP Wave Drag (supersonic area rule) of a tessellated surface, driven from nTop.

    python run_wave.py SETTINGS_FILE

OpenVSP imports the surface mesh (``<case>.stl``, metres) as a MeshGeom and runs its Wave Drag analysis: the
body is cut by Mach planes at several roll angles (rotation sections), and the wave drag of each equivalent area
distribution is averaged (von Karman slender-body theory, as in OpenVSP's Wave Drag screen). Inflow and outflow
faces are not used (out of scope).

Writes ``<case>.ntvsp.wave.csv`` (CD_wave and the run's settings) and ``<case>.ntvsp.wave_areas.csv`` (for each
rotation section: roll angle, x/L station and slice area per slice).
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from common import DriverError, Log, Settings, main_wrapper, read_vsp_results_csv, run_vspscript, vsp_path, write_csv  # noqa: E402


def script(s: Settings, stl, results_out):
    return f"""void main()
{{
    string gid = ImportFile("{vsp_path(stl)}", IMPORT_STL, "");
    if (gid == "") {{ Print("ERR could not import the surface mesh"); return; }}
    Update();
    SetVSP3FileName("{vsp_path(os.path.splitext(stl)[0] + '.vsp3')}");
    string wc = FindContainer("WaveDragSettings", 0);
    SetParmVal(wc, "RefFlag", "WaveDrag", 0);
    SetParmVal(wc, "Sref", "WaveDrag", {s.real("sref")!r});
    SetIntAnalysisInput("WaveDrag", "Set", array<int> = {{0}});
    SetIntAnalysisInput("WaveDrag", "NumSlices", array<int> = {{{s.integer("num_slices")}}});
    SetIntAnalysisInput("WaveDrag", "NumRotSects", array<int> = {{{s.integer("num_rot_sects")}}});
    SetIntAnalysisInput("WaveDrag", "SymmFlag", array<int> = {{{1 if s.flag("symmetry") else 0}}});
    SetDoubleAnalysisInput("WaveDrag", "Mach", array<double> = {{{s.real("mach")!r}}});
    string rid = ExecAnalysis("WaveDrag");
    if (rid == "") {{ Print("ERR Wave Drag returned no results"); return; }}
    WriteResultsCSVFile(rid, "{vsp_path(results_out)}");
    Print("NTVSP_OK");
}}
"""


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = os.path.abspath(s.text("vsp_dir"))
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    stl = os.path.join(work, case + ".stl")
    if not os.path.isfile(stl):
        raise DriverError(f"surface mesh {stl} not found")
    if s.real("mach") < 1.0:
        raise DriverError("Wave Drag needs a supersonic Mach number (1 or more)")
    results_out = os.path.join(work, case + ".ntvsp.wave_results.csv")
    if os.path.exists(results_out):
        os.remove(results_out)
    log(f"Wave Drag at Mach {s.real('mach')}: {s.integer('num_slices')} slices, {s.integer('num_rot_sects')} rotation sections")
    run_vspscript(vsp_dir, script(s, stl, results_out), work, case + ".wave", log, s.real("timeout_s", 7200))
    res = read_vsp_results_csv(results_out).get("WaveDrag")
    if not res:
        raise DriverError("no WaveDrag results in OpenVSP's output")
    first = lambda k: res[k][0][0]  # noqa: E731
    write_csv(os.path.join(work, case + ".ntvsp.wave.csv"),
              ["CD_wave", "Mach", "Num_Slices", "Num_Rotation_Sections", "Sref"],
              [[first("CDWave"), first("Mach"), first("Num_Slices"), first("Num_Cone_Sections"), s.real("sref")]])
    thetas = res["Theta"][0]
    x_norm = res["X_Norm"][0]
    areas = res["Slice_Area"]
    rows = []
    for t, (theta, a_row) in enumerate(zip(thetas, areas), start=1):
        for x, a in zip(x_norm, a_row):
            rows.append([t, theta, x, a, res["Start_X"][0][t - 1], res["End_X"][0][t - 1], res["MaxArea"][0][t - 1]])
    write_csv(os.path.join(work, case + ".ntvsp.wave_areas.csv"),
              ["Section", "Theta", "X_Norm", "Area", "Start_X", "End_X", "Max_Area"], rows)
    # One file per rotation section, so nTop reads one area distribution with Import Table alone.
    for t in range(1, len(thetas) + 1):
        write_csv(os.path.join(work, f"{case}.ntvsp.wave_areas.{t}.csv"), ["X_Norm", "Area"],
                  [[r[2], r[3]] for r in rows if r[0] == t])
    log(f"CD_wave = {first('CDWave')}")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
