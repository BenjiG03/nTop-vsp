"""Write the nTop custom blocks: nTop implicit bodies to an OpenVSP model, and OpenVSP's VSPAERO, Wave Drag and
Parasite Drag on it.

    sh build/py.sh build/build_blocks.py [--out DIR]

Each file in ``blocks/`` is one custom block; bring it into nTop with File > Import (or drag it into a
notebook). The blocks call the Python drivers in ``vsp_ntop/`` with Run Command, so nTop asks once to trust
the notebook's Run Command blocks.
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import aero_blocks as a  # noqa: E402
import buildup_blocks as b  # noqa: E402
import drag_blocks as d  # noqa: E402
import model_blocks as m  # noqa: E402
import vspaero_blocks as v  # noqa: E402
from ntopfile import save_block  # noqa: E402


def blocks():
    """(file stem, block): the model blocks, the flight condition, each analysis with its settings and readers, then
    the buildup."""
    return [
        ("OpenVSP Fuselage (VSP)", m.fuselage_block()),
        ("OpenVSP Wing (VSP)", m.wing_block()),
        ("OpenVSP Model (VSP)", m.model_block()),
        ("OpenVSP Model Check (VSP)", m.model_check_block()),
        ("OpenVSP Model Mesh (VSP)", m.model_mesh_block()),
        ("Flight Condition (VSP)", a.flight_block()),
        ("VSPAERO Solver Settings (VSP)", v.solver_block()),
        ("VSPAERO Analysis (VSP)", v.analysis_block()),
        ("VSPAERO Coefficient (VSP)", v.coefficient_block()),
        ("VSPAERO Surface Pressure (VSP)", v.surface_block()),
        ("Wave Drag (VSP)", d.wave_block()),
        ("Wave Drag Result (VSP)", d.wave_result_block()),
        ("Wave Drag Area Distribution (VSP)", d.wave_areas_block()),
        ("Parasite Drag Component (VSP)", d.component_block()),
        ("Parasite Drag (VSP)", d.parasite_block()),
        ("Parasite Drag Result (VSP)", d.parasite_result_block()),
        ("Parasite Drag Component Values (VSP)", d.parasite_components_block()),
        ("Modified Newtonian (VSP)", a.newtonian_block()),
        ("Modified Newtonian Pressure (VSP)", a.newtonian_pressure_block()),
        ("Viscous Correction (VSP)", a.viscous_block()),
        ("Aero Coefficients (VSP)", a.coefficients_block()),
        ("VSPAERO Contribution (VSP)", b.vspaero_contribution_block()),
        ("Wave Drag Contribution (VSP)", b.drag_contribution_block("wave")),
        ("Parasite Drag Contribution (VSP)", b.drag_contribution_block("parasite")),
        ("Coefficients Contribution (VSP)", b.coefficients_contribution_block()),
        ("Aero Database (VSP)", b.database_block()),
        ("Aero Database Coefficient (VSP)", b.database_coefficient_block()),
        ("Aero Database Table (VSP)", b.database_table_block()),
        ("Center of Pressure (VSP)", v.center_block()),
        ("Polar Plot (VSP)", b.plot_block()),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=os.path.join(os.path.dirname(HERE), "blocks"))
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    current = set()
    for stem, block in blocks():
        print(save_block(block, os.path.join(args.out, stem + ".ntop")))
        current.add(stem + ".ntop")
    for old in sorted(set(os.listdir(args.out)) - current):
        if old.endswith(".ntop"):
            os.remove(os.path.join(args.out, old))
            print("removed", old)


if __name__ == "__main__":
    main()
