"""OpenVSP Wave Drag (supersonic area rule) of an OpenVSP model over a sweep of Mach numbers, driven from nTop.

    python -m vsp_ntop.run_wave SETTINGS_FILE
"""
from __future__ import annotations

import os


from .common import (openvsp_dir, DriverError, Log, Settings, main_wrapper, read_vsp_results_csv, run_vspscript,
                     vsp_path, write_csv)
from .flight import Flight

SETS = {"all": 0, "thick": 3, "thin": 4}


def script(s: Settings, model, machs, results_out):
    lines = ["void main()", "{", f'    ReadVSPFile( "{vsp_path(model)}" );',
             '    if ( FindGeoms().size() == 0 ) { Print( "ERR the model has no components" ); return; }', "    Update();",
             '    string wc = FindContainer( "WaveDragSettings", 0 );',
             '    SetParmVal( wc, "RefFlag", "WaveDrag", 0 );',
             f'    SetParmVal( wc, "Sref", "WaveDrag", {s.real("sref")!r} );',
             f'    SetIntAnalysisInput( "WaveDrag", "Set", array<int> = {{{SETS[s.text("geom_set", "all")]}}} );',
             f'    SetIntAnalysisInput( "WaveDrag", "NumSlices", array<int> = {{{s.integer("num_slices")}}} );',
             f'    SetIntAnalysisInput( "WaveDrag", "NumRotSects", array<int> = {{{s.integer("num_rot_sects")}}} );',
             f'    SetIntAnalysisInput( "WaveDrag", "SymmFlag", array<int> = {{{1 if s.flag("symmetry") else 0}}} );']
    for k, m in enumerate(machs, start=1):
        lines += [f'    SetDoubleAnalysisInput( "WaveDrag", "Mach", array<double> = {{{m!r}}} );',
                  '    { string rid = ExecAnalysis( "WaveDrag" );',
                  '      if ( rid == "" ) { Print( "ERR Wave Drag returned no results" ); return; }',
                  f'      WriteResultsCSVFile( rid, "{vsp_path(results_out)}.{k}.csv" ); }}']
    lines += ['    Print( "NTVSP_OK" );', "}"]
    return "\n".join(lines) + "\n"


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
    flight = Flight(s.text("flight"))
    s.raw["sref"] = repr(flight.ref[0])
    machs = [m for m in flight.machs() if m > 1.0]
    for stale in os.listdir(work):
        if stale.startswith(case + ".ntvsp.wave"):
            os.remove(os.path.join(work, stale))
    results_out = os.path.join(work, case + ".ntvsp.wave_results")
    header = ["Mach", "CD_wave", "Num_Slices", "Num_Rotation_Sections", "Sref"]
    if not machs:
        write_csv(os.path.join(work, case + ".ntvsp.wave.csv"), header, [])
        log("the flight condition has no supersonic Mach number: no wave drag")
        log.close()
        return
    log(f"Wave Drag at {len(machs)} Mach numbers from {machs[0]} to {machs[-1]}: {s.integer('num_slices')} slices, "
        f"{s.integer('num_rot_sects')} rotation sections")
    run_vspscript(vsp_dir, script(s, model, machs, results_out), work, case + ".wave", log, s.real("timeout_s", 7200))
    rows = []
    for k, m in enumerate(machs, start=1):
        path = f"{results_out}.{k}.csv"
        res = read_vsp_results_csv(path).get("WaveDrag")
        if not res:
            raise DriverError(f"no WaveDrag results for Mach {m}")
        first = lambda key: res[key][0][0]  # noqa: E731
        rows.append([first("Mach"), first("CDWave"), first("Num_Slices"), first("Num_Cone_Sections"), s.real("sref")])
        thetas, x_norm, areas = res["Theta"][0], res["X_Norm"][0], res["Slice_Area"]
        for t, a_row in enumerate(areas, start=1):
            write_csv(os.path.join(work, f"{case}.ntvsp.wave_areas.{k}.{t}.csv"), ["Theta", "X_Norm", "Area"],
                      [[thetas[t - 1], x, a] for x, a in zip(x_norm, a_row)])
        os.remove(path)
        log(f"Mach {m}: CD_wave = {first('CDWave')}")
    write_csv(os.path.join(work, case + ".ntvsp.wave.csv"), header, rows)
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
