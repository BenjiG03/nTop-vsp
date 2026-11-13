"""OpenVSP Parasite Drag (component drag buildup) of an OpenVSP model, driven from nTop.

    python run_parasite.py SETTINGS_FILE

OpenVSP's Parasite Drag tool needs OpenVSP's own components: it takes each component's wetted area from CompGeom
and its reference length, thickness ratio or fineness ratio from DegenGeom, and it leaves out mesh components
(an imported tessellation is reported as "not included in Parasite Drag calculation"). So this driver works on a
``.vsp3`` model, with the flow condition and the per-component settings from nTop.

Writes ``<case>.ntvsp.parasite.csv`` (one row of totals and flight condition), ``<case>.ntvsp.parasite_components.csv``
(one row per component) and OpenVSP's own buildup table ``<case>.ParasiteDragBuildUp.csv``.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from common import (LENGTH_UNITS, DriverError, Log, Settings, main_wrapper, read_vsp_results_csv,  # noqa: E402
                    run_vspscript, vsp_path, write_csv)

COMPONENT_FIELDS = ("name", "ff_body", "ff_wing", "ff_user", "percent_laminar", "q", "roughness_m", "te_tw", "taw_tw",
                    "incorporated_gen")

SUMMARY = [("Total_CD", "Total_CD_Total"), ("Geom_CD", "Geom_CD_Total"), ("Total_f", "Total_f_Total"),
           ("Geom_f", "Geom_f_Total"), ("Re_per_length", "Re_L"), ("Mach", "FC_Mach"), ("Vinf", "FC_Vinf"),
           ("Altitude", "FC_Alt"), ("Temperature", "FC_Temp"), ("Pressure", "FC_Pres"), ("Density", "FC_Rho"),
           ("Sref", "FC_Sref"), ("Num_Components", "Num_Comp")]
COMPONENTS = [("Label", "Comp_Label"), ("ID", "Comp_ID"), ("Surface", "Comp_SurfNum"), ("Swet", "Comp_Swet"),
              ("Lref", "Comp_Lref"), ("t_c_or_fineness", "Comp_FineRat"), ("FF_Equation", "Comp_FFEqnName"),
              ("FF", "Comp_FFOut"), ("Re", "Comp_Re"), ("Percent_Laminar", "Comp_PercLam"), ("Cf", "Comp_Cf"),
              ("Q", "Comp_Q"), ("f", "Comp_f"), ("CD", "Comp_CD"), ("Percent_Total_CD", "Comp_PercTotalCD"),
              ("Roughness", "Comp_Roughness"), ("Te_Tw", "Comp_TeTwRatio"), ("Taw_Tw", "Comp_TawTwRatio")]


def parse_components(s: Settings):
    """Lines ``component = name | ff_body | ff_wing | ff_user | %lam | Q | roughness_m | Te/Tw | Taw/Tw | gen``."""
    out = []
    for line in open(s.path, encoding="utf-8"):
        if not line.strip().startswith("component"):
            continue
        key, value = line.split("=", 1)
        if key.strip() != "component":
            continue
        parts = [p.strip() for p in value.split("|")]
        if not parts[0]:
            continue
        if len(parts) != len(COMPONENT_FIELDS):
            raise DriverError(f"component line has {len(parts)} fields, expected {len(COMPONENT_FIELDS)}: {line.strip()}")
        out.append(dict(zip(COMPONENT_FIELDS, parts)))
    return out


def script(s: Settings, components, csv_out, results_out, model_out):
    unit = s.integer("length_unit")
    to_model = 1.0 / LENGTH_UNITS[unit][1]  # metres -> model length units
    # Renaming the model keeps OpenVSP's side files (CompGeom tables) in the working directory, not beside the model.
    lines = [f'    ReadVSPFile("{vsp_path(s.text("model"))}");', f'    SetVSP3FileName("{vsp_path(model_out)}");',
             "    Update();",
             '    if (FindGeoms().size() == 0) { Print("ERR the model has no components"); return; }']
    for c in components:
        name = c["name"].replace('"', "")
        lines += [
            f'    {{ array<string> @ids = FindGeomsWithName("{name}");',
            f'      if (ids.size() == 0) {{ Print("ERR no component named {name} in the model"); return; }}',
            "      for (int i = 0; i < int(ids.size()); i++) {",
        ]
        for parm, value in (("FFBodyEqnType", c["ff_body"]), ("FFWingEqnType", c["ff_wing"]), ("FFUser", c["ff_user"]),
                            ("PercLam", c["percent_laminar"]), ("Q", c["q"]),
                            ("Roughness", float(c["roughness_m"]) * to_model), ("TeTwRatio", c["te_tw"]),
                            ("TawTwRatio", c["taw_tw"]), ("IncorporatedGen", c["incorporated_gen"])):
            lines.append(f'        SetParmVal(ids[i], "{parm}", "ParasiteDragProps", {float(value)!r});')
        lines += ["      }", "    }"]
    sref = s.real("sref") * to_model ** 2
    ints = {"GeomSet": s.integer("geom_set"), "RefFlag": 0, "FreestreamPropChoice": 0, "LengthUnit": unit,
            "LamCfEqnChoice": s.integer("lam_cf"), "TurbCfEqnChoice": s.integer("turb_cf"),
            "VelocityUnit": s.integer("speed_unit"), "AltLengthUnit": 1, "TempUnit": 0, "PresUnit": 3,
            "ExportSubCompFlag": s.integer("export_subcomp"), "RecomputeGeom": 1}
    doubles = {"Sref": sref, "Vinf": s.real("speed"), "Altitude": s.real("altitude"), "DeltaTemp": s.real("delta_temp"),
               "SpecificHeatRatio": s.real("gamma")}
    for k, v in ints.items():
        lines.append(f'    SetIntAnalysisInput("ParasiteDrag", "{k}", array<int> = {{{int(v)}}});')
    for k, v in doubles.items():
        lines.append(f'    SetDoubleAnalysisInput("ParasiteDrag", "{k}", array<double> = {{{float(v)!r}}});')
    lines += [
        f'    SetStringAnalysisInput("ParasiteDrag", "FileName", array<string> = {{"{vsp_path(csv_out)}"}});',
        '    string rid = ExecAnalysis("ParasiteDrag");',
        '    if (rid == "") { Print("ERR Parasite Drag returned no results"); return; }',
        f'    WriteResultsCSVFile(rid, "{vsp_path(results_out)}");',
        '    Print("NTVSP_OK");',
    ]
    return "void main()\n{\n" + "\n".join(lines) + "\n}\n"


def main(settings_path):
    s = Settings(settings_path)
    work = os.path.abspath(s.text("work_dir"))
    case = s.text("case")
    vsp_dir = os.path.abspath(s.text("vsp_dir"))
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(work, case + ".ntvsp.log"))
    model = s.text("model")
    if not os.path.isfile(model):
        raise DriverError(f"OpenVSP model {model} not found")
    components = parse_components(s)
    log(f"Parasite Drag of {model}: {len(components)} component settings")
    csv_out = os.path.join(work, case + ".ParasiteDragBuildUp.csv")
    results_out = os.path.join(work, case + ".ntvsp.parasite_results.csv")
    for p in (csv_out, results_out):
        if os.path.exists(p):
            os.remove(p)
    text = run_vspscript(vsp_dir, script(s, components, csv_out, results_out, os.path.join(work, case + ".vsp3")), work, case + ".parasite", log,
                         s.real("timeout_s", 7200))
    for line in text.splitlines():
        if "not included in Parasite Drag" in line:
            log("note: " + line.strip())
    res = read_vsp_results_csv(results_out).get("Parasite_Drag")
    if not res:
        raise DriverError("no Parasite_Drag results in OpenVSP's output")
    write_csv(os.path.join(work, case + ".ntvsp.parasite.csv"), [k for k, _ in SUMMARY],
              [[res.get(src, [""])[0][0] if res.get(src) and res[src][0] else "" for _, src in SUMMARY]])
    n = int(float(res["Num_Comp"][0][0]))
    rows = []
    for i in range(n):
        rows.append([res[src][0][i] if src in res and i < len(res[src][0]) else "" for _, src in COMPONENTS])
    write_csv(os.path.join(work, case + ".ntvsp.parasite_components.csv"), [k for k, _ in COMPONENTS], rows)
    log(f"Total CD = {res['Total_CD_Total'][0][0]} over {n} components")
    log("done")
    log.close()


if __name__ == "__main__":
    main_wrapper(main)
