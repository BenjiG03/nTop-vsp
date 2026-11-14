# nTop_vsp

nTop custom blocks for three OpenVSP 3.56 analyses: **VSPAERO** (panel method or vortex lattice), **Parasite
Drag** and **Wave Drag**, plus separate blocks that bring the results back into nTop, including surface pressure
contours as Scalar Point Maps, and a drag buildup that combines them. The test case is the HL-20 lifting body (`hl20.vsp3`, `h20_vsp_mesh.stl`, `h20_ntop_mesh.stl`).

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
* `examples/` — HL-20 notebooks with floating comments: `HL-20 VSPAERO.ntop`, `HL-20 VSPAERO VLM.ntop`,
  `HL-20 Parasite Drag.ntop`, `HL-20 Wave Drag.ntop`, `HL-20 Drag Buildup.ntop`. They embed the blocks they use.
* `vsp_ntop/` — the drivers (standard-library Python 3.8+, no packages): `run_vspaero.py`, `run_parasite.py`,
  `run_wave.py`, `run_buildup.py`, `run_project.py` (thin solution onto a thick mesh), `vspgeom.py` (VSPAERO
  panel file, trailing-edge wakes, thin surfaces), `adb.py` (VSPAERO's binary results).
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
| VSPAERO Analysis (VSP) | results path | surface mesh, method (panel or vortex lattice), Sref, bref, cref, moment point, the three settings, run folder |
| VSPAERO Coefficient (VSP) | list of reals | any column of VSPAERO's polar, one value per case |
| VSPAERO Surface Pressure (VSP) | Scalar Point Map | Cp or p - p_inf of one case at panel centroids, nodes or both — the contour plot (vortex lattice: the jump across the sheet) |
| VSPAERO Wake Lines (VSP) | Scalar Point Map | the trailing-edge nodes wakes were shed from (value = wake number) |
| VSPAERO Thin Surface Projection (VSP) | Scalar Point Map | a vortex lattice solution shown on a thick mesh: half the pressure jump, signed by side, at the thick mesh's points |
| VSPAERO Mesh Faces (VSP) | Scalar Point Map | a panel run's mesh as surface points weighted by area (m^2), 1 or 4 per face, for ntop_newtonian's Viscous Correction on Mesh blocks |
| Parasite Drag Component (VSP) | text | one component's form factor equations, % laminar, Q, roughness, temperature ratios |
| Parasite Drag (VSP) | results path | OpenVSP model, units, Sref, Cf equations, speed, altitude, dT, gamma, components |
| Parasite Drag Result (VSP) | real | total CD, f, Re/L and the flight condition |
| Parasite Drag Component Values (VSP) | list of reals | one column of the buildup table, per component |
| Wave Drag (VSP) | results path | surface mesh, Mach, slices, rotation sections, Sref, symmetry |
| Wave Drag Result (VSP) | real | CD_wave |
| Wave Drag Area Distribution (VSP) | list of reals | x/L stations or slice areas of one rotation section |
| Drag Buildup (VSP) | table | VSPAERO's inviscid coefficients + CD_wave + viscous CD, CL, CMy per case; also written to `<Case Name>.ntvsp.buildup.csv` |

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
| Find body wakes (experimental mode) | Wake Settings; marks wakes behind the blunt bases of OpenVSP's own bodies, and found none on either HL-20 mesh |
| Thick geometry set | the Surface Mesh, Method = Panel |
| Thin geometry set (VLM surfaces) | the Surface Mesh, Method = Vortex lattice (see "Vortex lattice" below) |
| X-Z symmetry; 2D FEM and Tecplot files; Sref/cref from the model (RefFlag, MAC, S curve); CG calculation (mass sets, slices); stall model (Clo2D, CLMax2D); propellers, rotor blade mode, actuator disks, hover ramp, RPM, noise (PSU-WOPWOP); stability type and its time stepping (`-stab`, `-pstab`, `-qstab`, `-rstab`, `-unsteady`, time step, averaging); Cp slicer; control surfaces | out of scope (your list). Noise and the unsteady time-step settings only apply to rotors and unsteady stability runs |
| NRef (mesh refinement levels) | not applicable: for an imported mesh OpenVSP writes identical copies of the mesh (checked) |
| Cull fraction, merge degenerate tri regions, clean intersections | not applicable: OpenVSP's processing of the thin surfaces it builds from its own wings (DegenGeom); a thin mesh is used as given |
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
   none on either HL-20 mesh), and without wakes VSPAERO gives almost no lift. Choosing the
   trailing edges is the user's job: **VSPAERO Wake Settings (VSP)** takes a list of the mesh vertices on the
   trailing edges, and every mesh edge joining two listed vertices becomes a wake edge (connected edges form one
   wake). The run stops with an error when a point is not a mesh vertex (within Vertex Tolerance), a listed vertex
   has no listed neighbour, three wake edges meet at a vertex, or all the vertices of one face are listed. See the
   blunt trailing-edge result below: the mesh must provide a good edge for the wake.
3. **`-fem` and survey points (`-interrogate`)** were not implemented: `-fem` needs a structural deformation file
   whose format was not available to test, and survey points need an output reader for VSPAERO's survey file.
   Both look moderate; say if you want them.

## Vortex lattice (thin surfaces)

**VSPAERO Analysis (VSP)** has a **Method** input. Panel solves a closed, thick surface; Vortex lattice solves a
thin sheet (a camber surface or planform: one open mesh) and returns the pressure jump across it. VSPAERO picks the
method per component from the `thick` column of OpenVSP's `.vkey` file. OpenVSP's own thin geometry set builds thin
surfaces only from its parametric wings and writes nothing for an imported mesh (checked), so the driver passes the
mesh through the thick set and marks it thin in the `.vkey` file. Wakes are listed as for the panel method: the
sheet's trailing-edge vertices.

**A VSPAERO bug, and the way around it.** VSPAERO solves on its grid level 1. For a panel file of triangles it
builds that level itself by pairing triangles, wake triangles included, into quads. Behind a swept or curved
trailing edge the wake cells are skewed, and a wake triangle can be paired across a trailing vortex line; that line
then loses an edge and VSPAERO stops with "Could not find global edge for trailing wake!" (found in VSPAERO's
source, `VSP_AGGLOM::CreateMixedMesh_` and `VSP_GEOM::StoreWakeKuttaEdges`). It failed on the thin HL-20 for both
swept wing trailing edges whatever the wake settings, and on OpenVSP's own default wing once its quads were split
into triangles; only the straight body trailing edge worked. OpenVSP's own wing files avoid it with a "polygon"
layout: polygons as level 1, plus a triangulation of each polygon as level 0, and the wake built as quads. The
driver writes thin meshes in that layout: level 0 is every triangle exactly as given, and level 1 is the same
triangles, except that two neighbouring, coplanar triangles away from the edges form one quad (VSPAERO recognises
the layout only when a face has more than three nodes). VSPAERO then solves on the mesh's triangles. VSPAERO also
crashes on this layout with Unix line ends, so every panel file is written with Windows line ends, as OpenVSP
writes them. Checked on OpenVSP's default wing (Mach 0.3, 10 degrees): rewritten in this layout with its own
quads it gives OpenVSP's result exactly (CL = 0.3570); with every quad but one split into triangles, as the driver
does for a mesh, it gives 0.3654 (+2.3%, the finer solution grid). The panel method is unchanged: its example gives
identical results.

**Pressure jump.** For a thin surface VSPAERO stores, per triangle, Cp on the side the triangle's normal points to
minus Cp on the other side (found by integrating it on OpenVSP's own wing and on the HL-20 sheet: negative where
the sheet lifts toward its normals). VSPAERO Surface Pressure (VSP) shows that jump.

**VSPAERO Thin Surface Projection (VSP)** shows a vortex lattice solution on the vehicle's thick mesh. For every
node and centroid of the thick mesh the driver finds the closest point of the sheet, interpolates the jump there
(linear over the sheet triangle, from area-weighted node values), and gives the point Cp = +jump/2 where its outward
normal points the way the sheet's normals do and -jump/2 where it points the other way: the lifting part of the
pressure. A thin solution has no thickness effect, so this is not the thick body's pressure field. Points farther
than Maximum Distance from the sheet are left out (for parts the sheet does not model). The two meshes need not
match. Output quantities: Cp, p - p_inf, or the jump itself.

**HL-20.** `h20_vlm_mesh_from ntop.stl` is a planform sheet at z = 0 (2,569 triangles, 20.5 m x 8.5 m, 116 m^2),
not the same model as the thick mesh: it has no fins or tail. `examples/HL-20 VSPAERO VLM.ntop` runs it at the
panel example's conditions (Mach 0, alpha 1, 5.5, 10, Sref 100, cref = bref = 1) with wakes from the two wing
trailing edges and the body's base edge (33 vertices), and projects case 2 onto the thick nTop mesh with Maximum
Distance 3 m (the body lies within about 2.5 m of the sheet, the fins and tail reach 8 m). CL = 0.016 / 0.108 /
0.208, CMy = -0.085 / -0.821 / -1.759; 54 s in nTop, of which the projection takes 16 s.
`examples/hl20_vlm_projection_preview.png` plots the jump and the projected upper and lower Cp at 5.5 degrees.

* Slender-wing theory (CL_alpha = pi AR / 2, AR = 0.63, on the planform area, scaled to Sref 100) gives CL = 0.20
  at 10 degrees; the vortex lattice gives 0.208.
* **Wake choice at the wing-body junction.** The wing trailing edges start at |y| = 0.94 m, inboard of the body's
  side edge (which runs aft to |y| = 1.46 m), so the wakes of their two inboard vertices trail across the body part
  of the sheet, in its plane: the jump shows streaks there. Starting the wing wakes outboard of the body instead
  (|y| >= 1.53 m, 29 vertices) gives CL = 0.012 / 0.065 / 0.120 and CMy = -0.025 / -0.123 / -0.295 at Mach 0.3
  (against 0.017 / 0.102 / 0.207 and -0.088 / -0.684 / -1.692 with all 33). As for the blunt bases of the thick
  mesh, where the wakes leave decides much of the answer, and the planform needs a trailing edge that is the aft
  edge of the sheet for a clean result.
* **Surface pressures are not VSPAERO's forces.** VSPAERO computes forces from its vortex edges, not by
  integrating these pressures: integrated over the surface, the stored Cp gives 50 to 60% of CLi on the thick
  HL-20, 67% on the thin one and 85% on OpenVSP's default wing. Use the pressure maps to look at the flow and as
  edge conditions (as the viscous correction does); take coefficients from VSPAERO Coefficient (VSP).

## Drag buildup

OpenVSP has no tool that adds these analyses together: VSPAERO, Wave Drag and Parasite Drag each report on their
own, and VSPAERO's polar carries only its own rough viscous estimate (the `o` columns, CDo, driven by Re_cref).
**Drag Buildup (VSP)** is therefore a small Python driver (`vsp_ntop/run_buildup.py`). Per VSPAERO case:

    CL  = CLi  + CL_viscous
    CD  = CDi  + CD_wave + CD_viscous
    CMy = CMiy + CMy_viscous

CLi, CDi and CMiy are VSPAERO's inviscid parts; its CDo is left out because the viscous correction replaces it.
CS, CMx and CMz are VSPAERO's inviscid values. The table goes to `<Case Name>.ntvsp.buildup.csv` and is the
block's output.

**Shared conditions.** In the notebook the flight condition and reference (Mach, alpha start/end/points, Sref,
cref, moment reference point, freestream pressure and temperature) are notebook variables wired into every
analysis (`examples/HL-20 Drag Buildup.ntop`). The driver also checks what it can see: Wave Drag's Mach number
must equal every VSPAERO case's and its Sref VSPAERO's, and each viscous list must have one value per case (or
one value for all). It cannot see the viscous correction's inputs.

**Viscous correction (ntop_newtonian).** The Viscous Correction Coefficient (nTop) blocks from `ntop_newtonian`
fit VSPAERO's output:

* VSPAERO Surface Pressure gives each case's Cp as a point map; Field from Point Map turns it into the Cp field
  the viscous block takes. The body is an implicit of the same geometry (nTop's own, or Implicit Body from Mesh).
* **On-mesh blocks.** Viscous Correction on Mesh and Viscous Correction Subsonic on Mesh sum friction over
  surface points with their areas instead of integrating over the body: seconds instead of minutes per call.
  **VSPAERO Mesh Faces (VSP)** supplies the points: on panel runs the driver writes the analysed mesh's faces as
  `<Case Name>.ntvsp.faces.{1,4}.csv` (face centroids, or each face split in four at its edge midpoints, with
  the area each point stands for). On the HL-20 the six friction coefficients take about 3 s against about 140 s
  for the integrals, CD within 0.4%; Cm, a small difference of large opposite moments, needs the 4-point rule.
* Axes and coefficients agree: x aft, z up, y = z x x in both; CD and CL in wind axes; Cm nose-up about the
  moment center, so with Reference Length = cref its Cm is VSPAERO's CMy.
* **Two blocks, by Mach number.** Viscous Correction Coefficient (nTop) takes the edge state from the Rayleigh
  pitot pressure (a normal shock), so it holds above Mach 1. **Viscous Correction Subsonic Coefficient (nTop)**,
  added for this buildup, takes it from the freestream's isentropic total pressure for 0 < Mach < 1, and is
  otherwise the same method (within 0.5% of hypanel's formulas on a flat plate at Mach 0.3 and 0.8; see
  `ntop_newtonian/validation/REPORT.md`). Both normalize by q from the Mach number, so Mach must be above 0: for
  VSPAERO's incompressible Mach 0 runs, give the viscous block the flight Mach number.
* It has one reference length for all moments; VSPAERO uses bref for CMx and CMz, so only CMy gets a viscous part.
* Its method comes from hypersonic practice: friction fades out on surfaces inclined 30 to 60 degrees leeward,
  and the running length is measured along the freestream from the most upstream point.

**No supersonic buildup on the HL-20.** Wave drag exists only above Mach 1, but VSPAERO's panel method gave no
usable supersonic solution on the HL-20: on the nTop mesh, Mach 1.2 gave NaN and Mach 2 gave zero inviscid forces
with a diverging residual; on OpenVSP's own HL-20 panels (`hl20.vspgeom`), Mach 1.2 gave CL = -617 and Mach 2 did
not finish in 500 s. Mach 0.8 works (CLi = 0.103 at 5 degrees). The VSPAERO driver stops with an error on NaN
results instead of returning an empty polar.

**Subsonic buildup on the HL-20** (`examples/HL-20 Drag Buildup.ntop`): Mach 0.5 at 6096 m (US Standard
Atmosphere 1976: 46,563 Pa, 248.5 K), alpha 1, 5.5, 10, Sref 100, cref 1, moments about the origin; VSPAERO's panel
method on the nTop mesh (wakes as in the VSPAERO example), no wave drag, and the subsonic viscous correction
(on-mesh block, fully turbulent, adiabatic wall) summed over the analysed mesh's faces (4 points per face), with
normals from nTop's implicit HL-20 (`h20.implicit`) and each case's Cp at the mesh nodes.

| alpha | CL inviscid | CL viscous | CL | CD inviscid | CD viscous | CD | CMy inviscid | CMy viscous | CMy |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.1042 | 0.0001 | 0.1043 | 0.0289 | 0.0086 | 0.0376 | -2.0716 | 0.0037 | -2.0679 |
| 5.5 | 0.1341 | -0.0002 | 0.1339 | 0.0355 | 0.0086 | 0.0442 | -1.5690 | 0.0002 | -1.5688 |
| 10 | 0.1639 | -0.0006 | 0.1634 | 0.0410 | 0.0087 | 0.0497 | -1.0893 | -0.0032 | -1.0926 |

* 216 s in nTop, 184 s of it VSPAERO; the three viscous calls take seconds. With the integral block (Viscous
  Correction Subsonic Coefficients) the same notebook took 999 s, each call being six surface integrals over the
  HL-20 at 24 to 36 s each. At 5.5 degrees the on-mesh block gives friction CD 0.008649 against the integral's
  0.008686 (-0.4%) and Cm within 1% (4 points per face; 1 point per face is 0.4% on CD but 10% on Cm).
* The friction CD (0.0086) is below OpenVSP's own Parasite Drag CD0 for the same flight condition (0.0137, on the
  OpenVSP model with form factors and interference), as flat-plate friction without form factors should be.
* The inviscid CD comes from VSPAERO's panels on triangles, which disagreed with OpenVSP's own panels at Mach 0 (see
  the HL-20 test case below); treat it with care.


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
* The surface pressures are read from VSPAERO's binary `.adb` (layout from VSPAERO's source,
  `WriteOutAerothermalDatabaseSolution`; each case is found by its Mach, alpha and beta); every run checks the last
  case against VSPAERO's Tecplot output and stops if they differ.
* Reynolds-number variants of a case share one surface solution (VSPAERO's potential solution does not depend
  on Re_cref); only the skin-friction coefficients differ.
