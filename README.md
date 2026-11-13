# nTop_vsp

nTop custom blocks for three OpenVSP 3.56 analyses: **VSPAERO** (panel method), **Parasite Drag** and **Wave
Drag**, plus separate blocks that bring the results back into nTop, including surface pressure contours as
Scalar Point Maps. The test case is the HL-20 lifting body (`hl20.vsp3`, `h20_vsp_mesh.stl`, `h20_ntop_mesh.stl`).

## How it works

Each analysis block writes its inputs to a settings file (`key = value`, SI units), exports the geometry, and
calls a Python driver with nTop's **Run Command**. The driver runs OpenVSP (`vspscript.exe`) and VSPAERO
(`vspaero.exe`) and writes plain CSV files. The block's output is the path of those results (text); reader
blocks load them with Import Table and Import Scalar Point Map. A failed run stops the notebook with the
driver's error (Validate block), and every run leaves `<Case Name>.ntvsp.log` in the Working Directory.

```
nTop inputs ──► settings.txt + <case>.stl ──► vsp_ntop/run_*.py ──► OpenVSP / VSPAERO ──► *.ntvsp.*.csv ──► reader blocks
```

* `blocks/` — the custom blocks, one `.ntop` file each (File > Import, or drag into a notebook).
* `examples/` — HL-20 notebooks with floating comments: `HL-20 VSPAERO.ntop`, `HL-20 Parasite Drag.ntop`,
  `HL-20 Wave Drag.ntop`. They embed the blocks they use.
* `vsp_ntop/` — the drivers (standard-library Python 3.8+, no packages): `run_vspaero.py`, `run_parasite.py`,
  `run_wave.py`, `vspgeom.py` (VSPAERO panel file and trailing-edge wakes), `adb.py` (VSPAERO's binary results).
* `build/` — writes the blocks and examples (`build_blocks.py`, `build_examples.py`; `ntopfile.py` is the
  `.ntop` writer from `ntop_newtonian`). `dev/inventory.txt` lists OpenVSP's analysis inputs and parameter
  containers, as dumped by `dev/inventory.vspscript`.

nTop asks once to trust Run Command blocks. The defaults of every block point at this folder (OpenVSP in
`OpenVSP-3.56.0-win64`, drivers in `vsp_ntop`, runs in `runs`) and at `C:/Users/Benji/anaconda3/python.exe`;
change the **OpenVSP Directory**, **Python**, **Driver Directory** and **Working Directory** inputs if those move.
OpenVSP is not in the repository: unpack the Windows 64-bit OpenVSP 3.56.0 build from openvsp.org into
`OpenVSP-3.56.0-win64/`, or point **OpenVSP Directory** at an existing install.

## Blocks

| Block | Output | Purpose |
|---|---|---|
| VSPAERO Flow Condition (VSP) | text | alpha, beta, Mach, Re_cref sweeps (start, end, points); Vinf, Rho; reference velocity |
| VSPAERO Solver Settings (VSP) | text | threads, wake model, far field, ground effect, freezing, implicit wake, convergence factors, preconditioners, quad tree, restart file |
| VSPAERO Wake Settings (VSP) | text | which mesh edges shed wakes (see below) |
| VSPAERO Analysis (VSP) | results path | surface mesh, Sref, bref, cref, moment point, the three settings, run folder |
| VSPAERO Coefficient (VSP) | list of reals | any column of VSPAERO's polar, one value per case |
| VSPAERO Surface Pressure (VSP) | Scalar Point Map | Cp or p - p_inf of one case at panel centroids, nodes or both — the contour plot |
| VSPAERO Wake Lines (VSP) | Scalar Point Map | the trailing-edge nodes wakes were shed from (value = wake number) |
| Parasite Drag Component (VSP) | text | one component's form factor equations, % laminar, Q, roughness, temperature ratios |
| Parasite Drag (VSP) | results path | OpenVSP model, units, Sref, Cf equations, speed, altitude, dT, gamma, components |
| Parasite Drag Result (VSP) | real | total CD, f, Re/L and the flight condition |
| Parasite Drag Component Values (VSP) | list of reals | one column of the buildup table, per component |
| Wave Drag (VSP) | results path | surface mesh, Mach, slices, rotation sections, Sref, symmetry |
| Wave Drag Result (VSP) | real | CD_wave |
| Wave Drag Area Distribution (VSP) | list of reals | x/L stations or slice areas of one rotation section |

All geometry and lengths are sent to OpenVSP in metres, so reference areas and lengths are nTop quantities with
units and the panels match them.

## Coverage of OpenVSP's options

Every option of OpenVSP 3.56's VSPAERO, Parasite Drag and Wave Drag analyses was listed from OpenVSP itself
(`dev/inventory.txt`: parameter containers and analysis inputs) and is either a block input, out of scope as you
listed, or not applicable to a tessellation (with the reason). Plot-axis and plot-selection parameters of the
GUI's result windows (`m_*XMin`, `ThetaIndex`, `AreaPlotType`...) are display state, not analysis options.

**VSPAERO**

| OpenVSP option | Here |
|---|---|
| Alpha, Beta, Mach, Re_cref start / end / points | Flow Condition |
| Vinf, Rho, Vref (manual flag), Machref | Flow Condition |
| Sref, bref, cref; Xcg, Ycg, Zcg | Analysis (Reference Area, Span, Chord, Moment Reference Point) |
| NCPU (`-omp`) | Solver: Processors |
| Fixed wake, wake iterations, wake nodes | Solver |
| Far field distance (toggle and value) | Solver |
| Ground effect (toggle and height, `-groundheight`) | Solver |
| Freeze multipole / wake at iteration, freeze root vortices, implicit wake and its start iteration, wake relaxation | Solver |
| Forward / adjoint GMRES and nonlinear convergence factors, core size factor, multipole far-away factor | Solver |
| Update matrix preconditioner, wake node matrix preconditioner, quad tree buffer levels | Solver |
| `-save`, `-restart` | Solver: Restart File |
| Find body wakes | Wake Settings |
| Thick geometry set | the Surface Mesh |
| X-Z symmetry; 2D FEM and Tecplot files; Sref/cref from the model (RefFlag, MAC, S curve); CG calculation (mass sets, slices); stall model (Clo2D, CLMax2D); propellers, rotor blade mode, actuator disks, hover ramp, RPM, noise (PSU-WOPWOP); stability type and its time stepping (`-stab`, `-pstab`, `-qstab`, `-rstab`, `-unsteady`, time step, averaging); Cp slicer; control surfaces | out of scope (your list). Noise and the unsteady time-step settings only apply to rotors and unsteady stability runs |
| Thin geometry set (VLM surfaces) | **flagged**: built from OpenVSP's parametric wings (DegenGeom camber surfaces); a tessellation has none |
| NRef (mesh refinement levels) | not applicable: for an imported mesh OpenVSP writes identical copies of the mesh (checked) |
| Cull fraction, merge degenerate tri regions, clean intersections | not applicable: thin-surface (DegenGeom) processing |
| Modes / geometry sets | not applicable: the mesh is the geometry |
| `-fem` (FEM deformation file), `-interrogate` with survey points | **flagged**, not implemented (see below) |
| `-fs`, `-geom` | not needed: the flow condition is in the input file; geometry-only runs have no use here |

**Parasite Drag**: geometry set, model length unit, Sref, laminar and turbulent Cf equations, speed and its unit
(ft/s, m/s, mph, km/hr, KEAS, KTAS, Mach), altitude, temperature offset, specific heat ratio and sub-component
export are inputs of Parasite Drag (VSP); each component's body and wing form factor equations, manual form
factor, % laminar, Q, roughness, Te/Tw, Taw/Tw and incorporated generation are inputs of Parasite Drag
Component (VSP). Out of scope: atmospheres other than US Standard 1976 (and with them the manual pressure,
density, temperature, Re/L and viscosity inputs), excrescences, Sref from a wing. Pressure and temperature
units only change how OpenVSP labels its own report; the CSVs here use pascals and kelvin.

**Wave Drag**: Mach, number of slices, number of rotation sections, Sref and symmetry are inputs. Out of scope:
inflow and outflow faces. The geometry set is the Surface Mesh.

## Flagged (implementation cliffs, not built)

1. **Parasite Drag cannot take a tessellation.** OpenVSP excludes mesh components from the buildup (it reports
   "Geom ... not included in Parasite Drag calculation" for an imported STL), because it needs each component's
   reference length and thickness or fineness ratio from its parametric definition. Parasite Drag (VSP) therefore
   takes an OpenVSP model (`.vsp3`). Using an nTop mesh would mean computing Swet, Lref and t/c or l/d per
   component in nTop and re-implementing OpenVSP's equations outside OpenVSP.
2. **Wakes on a tessellation.** OpenVSP marks no trailing edges on an imported mesh (its body-wake finder found
   none on either HL-20 mesh), and without wakes VSPAERO's panel method gives almost no lift. Choosing the
   trailing edges is the user's job: **VSPAERO Wake Settings (VSP)** takes a list of the mesh vertices on the
   trailing edges, and every mesh edge joining two listed vertices becomes a wake edge (connected edges form one
   wake). The run stops with an error when a point is not a mesh vertex (within Vertex Tolerance), a listed vertex
   has no listed neighbour, three wake edges meet at a vertex, or all the vertices of one face are listed. See the
   blunt trailing-edge result below: the mesh must provide a good edge for the wake.
3. **VLM (thin surfaces)** needs OpenVSP's parametric wings; not possible from a mesh.
4. **`-fem` and survey points (`-interrogate`)** were not implemented: `-fem` needs a structural deformation file
   whose format was not available to test, and survey points need an output reader for VSPAERO's survey file.
   Both look moderate; say if you want them.

## HL-20 test case

**Wave Drag**, Mach 2, 20 slices, 10 rotation sections, Sref 100: CD_wave = 0.2176 on the nTop mesh,
0.2204 on OpenVSP's own STL of the model (1.3% apart). About 8 s in nTop.

**Parasite Drag** on `hl20.vsp3` (metres, Sref 100 m², Mach 0.5 at 6096 m, wings at 10% laminar and Q = 1.1):
CD0 = 0.01371 over five components. This is OpenVSP's own tool, so the numbers are OpenVSP's.

**VSPAERO** at the conditions of the run that came with the model (`hl20.vspaero`: Mach 0, alpha 1, 5.5, 10,
Sref 100, cref = bref = 1, moments about the origin). OpenVSP's own run reproduces `hl20.polar` exactly.

| Panels and wakes | CL at 1 / 5.5 / 10 deg | CMy at 1 / 5.5 / 10 deg |
|---|---|---|
| OpenVSP's own model (`hl20.polar`) | -0.108 / -0.088 / -0.065 | 1.762 / 2.364 / 2.904 |
| OpenVSP's STL through this pipeline, wakes on OpenVSP's trailing-edge edges | -0.123 / -0.097 / -0.067 | 1.736 / 2.278 / 2.741 |
| nTop mesh, wakes on one corner of each blunt base | -0.260 / -0.210 / -0.162 | 4.676 / 4.576 / 4.626 |
| nTop mesh, wakes on the other corner | +0.100 / +0.113 / +0.129 | -2.076 / -1.371 / -0.729 |

* The pipeline reproduces OpenVSP when it gets the same surface and the same wake edges: lift within 0.015 and
  pitching moment within 2-6%. Induced drag does not agree (CDi -0.007 against 0.075 from OpenVSP's
  quadrilateral panels); this was not resolved and VSPAERO's surface-integrated CDi on triangles should not be
  trusted for this body.
* **Blunt trailing edges decide the answer.** The HL-20's fins and tail end in a flat base about 0.37 m thick.
  OpenVSP sheds the wake from the middle of that base; the nTop mesh has no edge there, so the wake has to leave
  the upper or the lower corner, and the two choices bracket OpenVSP's result by a factor of two or more in CL.
  For VSPAERO the mesh needs either sharp trailing edges or a line of mesh edges along the middle of each blunt
  base (then list that line's vertices). The driver uses the mesh as it is and does not insert one.
* `examples/HL-20 VSPAERO.ntop` lists 53 vertices of the nTop mesh along the fins' and tail's trailing edges, on
  the lower-corner side of each base. The notebook gives CL = 0.116 / 0.142 /
  0.167 and CMy = -2.24 / -1.76 / -1.28 (186 s in nTop, 20,358 panels). `examples/hl20_cp_alpha5.5_preview.png`
  is a plot of its Cp point map at 5.5 degrees. A few panels reach Cp = -21 (minimum); clip the colour range in
  the point map's HUD.

## Rebuilding and running headlessly

```bash
sh build/py.sh build/build_blocks.py
sh build/py.sh build/build_examples.py --run
```

`build/py.sh` runs the `nTop_pkgs` environment (it has `implicit_io`). `--run` runs the examples with `ntopcl
--trustnotebook` and prints their outputs. ntopcl writes the outputs and then may not exit after notebooks with
Run Command inside custom blocks; the runner stops it once the output file is written. A signed-in nTop is
needed.

## Known limits

* Runs with the same Working Directory and Case Name overwrite each other. nTop's list processing over an
  analysis block (for example a list of meshes) needs a list of Case Names too.
* The surface pressures are read from VSPAERO's binary `.adb` by a layout worked out for VSPAERO 7.2.2; every run
  checks the last case against VSPAERO's Tecplot output and stops if they differ.
* Reynolds-number variants of a case share one surface solution (VSPAERO's potential solution does not depend
  on Re_cref); only the skin-friction coefficients differ.
