# nTop_vsp

nTop custom blocks that turn **nTop implicit bodies into an OpenVSP model** and run three OpenVSP 3.56 analyses on
it: **VSPAERO**, **Wave Drag** and **Parasite Drag**. A last block combines them into a **six-degree-of-freedom aero
database** (CL, CD, CY, Cl, Cm, Cn against Mach, alpha and beta). Because the geometry stays parametric in nTop, its
inputs can be independent variables of the database too: `examples/x59_doe.py` runs a design of experiments over
the wing span and sweep of an X-59-style aircraft.

Test geometry: the X-59 hackathon starter in `../x59-hackathon-starter` (MIT License, Copyright (c) 2026
matthewmuellernTop; notice in `examples/X-59_STARTER_LICENSE.txt`), a parametric nTop model whose fuselage, wing,
canards and tails are separate implicit bodies. `build/x59_example.py` reads it from there.

## How it works

```
implicit body ──► OpenVSP Fuselage / Wing (VSP) ──► component line ─┐
   (nTop Slice Body → section contours, CLI)                         ├─► OpenVSP Model (VSP) ──► <case>.vsp3
implicit body ──► ...                                               ─┘     (vsp_ntop/run_model.py: fit + build)
                                                                                      │
             ┌────────────────────────────┬───────────────────────────┬───────────────┘
   VSPAERO Analysis (VSP)          Wave Drag (VSP)            Parasite Drag (VSP)
   (Thin set: vortex lattice;      (area rule, Mach sweep)    (component buildup, Mach sweep)
    Thick set: panel method)
             └────────────────────────────┴───────────────────────────┴──► Aero Database (VSP) ──► table / DOE CSV
```

Each analysis block writes its inputs to a settings file (`key = value`, SI units) and runs a driver of the
`vsp_ntop` Python package (`python -m vsp_ntop.<module> <settings file>`) with nTop's **Run Command**; the driver runs OpenVSP (`vspscript.exe`) and VSPAERO and writes plain CSV
files, which reader blocks load with Import Table and Import Scalar Point Map. A failed run stops the notebook with
the driver's error, and every run leaves `<Case Name>.ntvsp.log` in the Working Directory.

* `blocks/` — the custom blocks, one `.ntop` file each (File > Import, or drag into a notebook).
* `examples/X-59 Aero Database.ntop` — the X-59 starter's geometry with the whole chain appended, floating comments
  throughout; `examples/x59_doe.py` — the design of experiments over it.
* `vsp_ntop/` — the drivers, an installable package (`pyproject.toml`; standard-library Python 3.8+, no
  dependencies): `run_model.py` with `slices.py` (contour reader),
  `fit.py` (section fitting) and `vspmodel.py` (OpenVSP script and `.vsp3` writing); `run_vspaero.py`,
  `run_wave.py`, `run_parasite.py`, `run_database.py`; `vspgeom.py` and `adb.py` read VSPAERO's panel file and
  binary results.
* `build/` — writes the blocks and the example (`build_blocks.py`, `x59_example.py`; `x59.py` derives notebooks from
  the starter; `ntopfile.py` is the `.ntop` writer from `ntop_newtonian`).
* `validation/wave_drag_mach.py` — the wave drag study below (results in `validation/results/`).

## Setup

1. Install OpenVSP 3.56: unpack the Windows 64-bit build (for example into `OpenVSP-3.56.0-win64/` here; it is not
   in the repository).
2. Install the drivers into a Python 3.8 or newer: `pip install -e .` in this folder (or `pip install .`).
3. Set two environment variables (Windows: System Properties > Environment Variables), then restart nTop so it
   sees them:
   * `NTOP_PYTHON_PATH`: that Python's executable, e.g. `C:/Users/me/anaconda3/python.exe`;
   * `VSP_PATH`: the OpenVSP folder, the one holding `vspscript.exe` and `vspaero.exe`.

The blocks read both with nTop's **Get Environment Variable** block, so they have no Python, driver or OpenVSP path
inputs; a missing variable stops the notebook with nTop's "Environment variable ... is not set". The only path
input left is each analysis block's **Working Directory** (default: `runs/` here). nTop asks once to trust Run
Command blocks.

## Blocks

| Block | Output | Purpose |
|---|---|---|
| OpenVSP Fuselage (VSP) | text | an implicit body as an OpenVSP Fuselage: slicing, sections, clustering, tessellation, set |
| OpenVSP Wing (VSP) | text | an implicit lifting surface as an OpenVSP Wing: span axis, symmetry, sections, airfoils, tessellation, clustering, set |
| OpenVSP Model (VSP) | model path | builds `<case>.vsp3` from a list of component lines |
| OpenVSP Model Check (VSP) | list of reals | per component: largest and mean distance from the OpenVSP surface to the nTop contours, sections, slices |
| OpenVSP Model Mesh (VSP) | mesh | OpenVSP's tessellation of the model, to view against the bodies |
| VSPAERO Flow Condition (VSP) | text | alpha, beta, Mach, Re_cref sweeps; Vinf, Rho; reference velocity |
| VSPAERO Solver Settings (VSP) | text | threads, wake model and nodes, far field, ground effect, freezing, implicit wake, convergence factors, preconditioners, quad tree, restart file |
| VSPAERO Analysis (VSP) | results path | model, Sref, bref, cref, moment point, body wakes, the two settings |
| VSPAERO Coefficient (VSP) | list of reals | any column of VSPAERO's polar, one value per case |
| VSPAERO Center of Pressure (VSP) | list of points | each case's center of pressure, or each alpha sweep's aerodynamic center (neutral point), as absolute positions with length units |
| VSPAERO Surface Pressure (VSP) | Scalar Point Map | Cp or p - p_inf of one case at panel centroids, nodes or both (thin surfaces: the jump across the sheet) |
| VSPAERO Mesh Faces (VSP) | Scalar Point Map | the Thick-set panels as points weighted by area, for ntop_newtonian's Viscous Correction on Mesh blocks |
| Wave Drag (VSP) | results path | model, Mach sweep, slices, rotation sections, Sref, symmetry, geometry set |
| Wave Drag Result (VSP) | list of reals | CD_wave (or Mach, counts, Sref), one value per Mach number |
| Wave Drag Area Distribution (VSP) | list of reals | x/L or slice areas of one rotation section at one Mach number |
| Parasite Drag Component (VSP) | text | one component's form factor equations, % laminar, Q, roughness, temperature ratios |
| Parasite Drag (VSP) | results path | model, units, Sref, Cf equations, speed sweep and unit, altitude, dT, gamma, components |
| Parasite Drag Result (VSP) | list of reals | total CD, f, Re/L and the flight condition, one value per speed |
| Parasite Drag Component Values (VSP) | list of reals | one column of the buildup table at one speed, per component |
| Aero Database (VSP) | table | CL, CD, CY, Cl, Cm, Cn per VSPAERO case, drag split into induced, wave, parasite (and optional viscous); design variables; appends to a DOE file |

Every length is sent to OpenVSP in metres, so reference areas and lengths are nTop quantities with units.

## From implicit to OpenVSP

![OpenVSP model of the X-59 starter](examples/x59_openvsp_model.png)

**Slicing in nTop.** A component block cuts its body with nTop's **Slice Body** (planes normal to x for a fuselage,
normal to the span axis for a wing) and exports the contours with **Export Slices to CLI** into a temporary folder.
The cuts start just outside the body's bounding box; **Slice Spacing** (a fraction of the body's length, or of its
extent along the span axis) and **Contour Tolerance** (a fraction of its bounding-box diagonal) set the resolution.
On the X-59 the six bodies slice in about a minute. Meshing the bodies instead was not workable: the fuselage alone
became a 1.4 GB mesh after 6 minutes at 6 mm tolerance.

**Fuselage fit** (`fit.fit_fuselage`). Stations are placed with OpenVSP's clustering function between the nose and
tail (**Sections**, **Nose/Tail Clustering**); each takes the largest closed contour at its slice. A pointed end
(section area tending to zero) becomes an OpenVSP point section placed by extrapolating the square root of the area;
a blunt end keeps its last contour and gets a flat cap. Each section is the contour itself (an OpenVSP file
cross-section, normalised by its width and height, positioned at its bounding-box centre), with free skinning
between sections. Building 30 such sections through OpenVSP's API takes 51 s, because every change rebuilds the
whole surface; the driver therefore makes the fuselage with the right number of sections by script and writes each
section's points, size and position straight into the `.vsp3` XML (2 s).

**Wing fit** (`fit.fit_wing`). At every slice the chord line runs from the leading edge (the point farthest from
the trailing edge) to the trailing edge (the mid-point of the base for a blunt one); that gives the leading-edge
position, chord, twist, thickness and the airfoil's upper and lower surfaces (cosine-spaced **Airfoil Points**).
Stations become OpenVSP wing sections where leading edge, chord, twist or thickness stop varying linearly (**Section
Tolerance**, 1% of the root chord), at most **Sections**, none narrower than 4% of the span, so measurement noise of
a few millimetres does not create sliver sections. Each section's span, sweep and dihedral come from its end
stations; the surface is rotated about x for a canted or vertical span axis. **Symmetric** mirrors it about the x-z
plane and uses only the +span side of the body. **Trailing Edge: Sharp** shears a blunt base onto the chord line so
that VSPAERO's wake leaves one edge. OpenVSP marks every wing trailing edge as a wake: no wake input.

**Tessellation.** Wings: **Tess U** spanwise points (shared between sections by span), **Tess W** around the
airfoil, OpenVSP's leading- and trailing-edge clustering, root and tip clustering. Fuselages: **Tess W** around each
section and a **Panel Aspect Ratio** (panel length along the body over its width around it, default 1.5) that sets
the points along the body section by section (see below), plus nose and tail clustering.

**Fit check.** The driver samples each OpenVSP surface on a grid and measures its distance to the nTop contours
(interpolated between slices): `<case>.ntvsp.model.csv`, **OpenVSP Model Check (VSP)**. X-59 baseline:

| Component | Sections | Mean | Largest | Largest / size |
|---|---|---|---|---|
| Fuselage (fuselage + nacelle, 30.3 m) | 24 | 10.2 mm | 435 mm | 1.4% (the nacelle inlet's step) |
| Wing (13.2 m root chord) | 2 | 10.0 mm | 65 mm | 0.5% |
| Canard | 1 | 2.8 mm | 13 mm | 0.6% |
| Horizontal tail | 1 | 0.2 mm | 12 mm | 0.2% |
| T-tail | 1 | 0.1 mm | 4 mm | 0.2% |
| Fin | 1 | 0.6 mm | 16 mm | 0.3% |

The whole X-59 model (slicing, fitting, building and checking) takes about 70 s in nTop.

**Slicer conventions found by testing nTop 6.3.3**: Slice Body cuts only from the frame's origin in the +z direction
of the frame; contours come in the frame's x and y relative to its origin, except when slicing along exactly +x, +y
or +z, where they come in the world's cyclic axes (y-z, z-x, x-y). The driver handles both; the blocks place the
frame origin outside the body.

## VSPAERO defaults, and why

OpenVSP puts components in VSPAERO's **Thick** set (panel method on the closed surface) or **Thin** set (vortex
lattice on the camber surface); each component block has a **Set** input with a third choice, **None**, which leaves
the component out of VSPAERO but keeps it in Wave Drag and Parasite Drag. Defaults: **Wing: Thin, Fuselage: None.**
They come from these tests on the X-59 (Mach 0.3 to 2.0, alpha 0 to 8 degrees, beta 0 and 4 degrees) and on
OpenVSP's own default fuselage and wing:

| Configuration | Subsonic | Supersonic |
|---|---|---|
| Thick fuselage, thin wing (OpenVSP's own model) | CL 0.035 at 2 deg, fine | CDi 0.075 at Mach 1.6, wrong |
| Thick fuselage, thick nacelle and wings (X-59) | VSPAERO refuses the mesh: "edge has 3 tris" where the overlapping bodies intersect | — |
| Thick fuselage, thin lifting surfaces (X-59) | wing alone fine; full aircraft CDi 0.014 at CL 0.09, wrong | wrong |
| Thin fuselage (camber plates), thin lifting surfaces | NaN | NaN |
| **Thin lifting surfaces, fuselage in no set** | CL_alpha 1.7/rad, span efficiency 0.9 | fine to Mach 2 |

* VSPAERO's panel method (Thick set) gave wrong drag above Mach 1 in every test, on OpenVSP's own models as well,
  so it is a subsonic option. Its body also needs care: union touching bodies in nTop into one component instead of
  letting OpenVSP intersect overlapping ones, close blunt ends (the driver adds flat caps), and keep body panels
  about as long as they are wide. Panels four times longer than wide (the first tessellation tried) gave CL -0.007
  and CDi -0.024 for the fuselage with a T-tail; at aspect ratio 1.5 the same pair gives sensible values. Hence
  **Panel Aspect Ratio**.
* Thin surfaces that cross each other near a tip make VSPAERO return NaN above Mach 1: the X-59's fin rises 55 mm
  through its T-tail. Cut such a surface in nTop (the example cuts the fin 65 mm below the T-tail).
* **Wake Nodes** default 16 (OpenVSP's GUI writes 64): on the X-59 at Mach 1.6 the coefficients agree within 0.2% and
  the run is 20 times faster. A fixed wake changed CL by under 0.1% there and is about four times faster than a
  relaxed one; the example uses it.
* VSPAERO's own axes are x aft, y right, z up, so its rolling and yawing moments have the opposite sign to aircraft
  body axes; the Aero Database reports body-axis signs (Cl right wing down, Cm nose up, Cn nose right).
* Leaving the fuselage out of VSPAERO loses its lift and moment contributions; give it the Thick set for subsonic
  work where its panel mesh behaves.
* **Convergence check.** After every run the driver checks each case: VSPAERO's linear solve (GMRES) of the last
  wake iteration must reach both of its targets (residual reduction and max residual; VSPAERO otherwise stops it
  at 1500 iterations), every force and moment coefficient must be finite, and with a relaxed wake CLtot, CDtot and
  CMytot must change by less than the **Wake Change Tolerance** (10%) over the last wake iteration. By default an
  unconverged case stops the notebook with the reason (**Convergence Check** in Solver Settings can make it a
  warning). `<case>.ntvsp.convergence.csv` lists every case, and the polar has a `Converged` column. Tried on the
  53 earlier runs of this project: every run of a converted model whose numbers were sound passes (the X-59 baseline: 60
  of 60); the failures it reports are the NaN cases above (T-tail crossing, thin fuselage) and relaxed wakes that
  were still moving 15 to 290% per iteration (thick fuselage with fins, the old mesh-route vortex lattice). Sound
  relaxed-wake runs moved under 6%. VSPAERO's wake span efficiency `Ew` is infinite in every supersonic case (it
  computes no wake forces there) and is not checked.

## Wave drag: the Mach range it holds over

OpenVSP's Wave Drag is the supersonic area rule: Mach-plane cuts at several roll angles, each turned into an
equivalent body of revolution whose von Kármán wave drag is averaged. That is zero-lift linear theory.
`validation/wave_drag_mach.py` measures its error in two parts (results in `validation/results/wave_drag_mach.json`
and `.png`).

![Wave drag study](validation/results/wave_drag_mach.png)

**1. OpenVSP against linear theory, which it implements.**

* Sears-Haack body, fineness 10 (linear theory: D/q = 9π A_max² / 2L², the same at every Mach): OpenVSP is within
  1.4% from Mach 1.05 to 2.5, then over-predicts as the Mach cone opens: +5% at Mach 3, +18% at 4, +29 to 30% at
  5 and 6. More slices drift upward too (+4% at 100 slices, Mach 2), because the area distribution of a tessellated
  surface gets noisy: 20 to 40 slices are best.
* Rectangular wings of 4% biconvex sections at Mach 2 (Ackeret: CD = (16/3)(t/c)²/√(M²−1); a fine roll-angle
  integration of the area rule reproduces it): OpenVSP's result depends on the number of **rotation sections**,
  because a wing's area-rule integrand peaks sharply at the roll angles where the cuts lie along the span. Aspect
  ratio 3: 3.7x Ackeret at 10 sections, 1.6x at 30, 1.08x at 120. Aspect ratio 12: 14x at 10, 1.6x at 120. On the
  X-59 (aspect ratio 1.3), 30 and 120 sections agree within 0.3%: slender configurations converge by 30, the
  default; straight or high-aspect-ratio wings need far more.
* A step in the area distribution (a blunt base, an inlet face) makes the result grow with the number of slices
  instead of converging: on the X-59, 0.0018, 0.0021, 0.0026 and 0.0028 at 20, 30, 60 and 100 slices (Mach 1.2).
  OpenVSP's inflow and outflow faces, out of scope here, are its answer to that.

**2. Linear theory against exact inviscid flow.** Biconvex airfoils against shock-expansion theory and cones
against the Taylor-Maccoll solution:

* Wings: Ackeret is within 2% of exact at t/c ≤ 4% from Mach 1.4 to 6, and within 10% up to Mach × t/c ≈ 0.6 (10%
  thick at Mach 6: -9%). Below about Mach 1.2 (t/c 2%), 1.4 (4%) or 1.6 (10%) the leading-edge shock detaches and
  linear theory does not apply.
* Bodies: slender-body theory under-predicts cone pressure drag by 10% at Mach × half-angle ≈ 0.2 (5° at Mach 2.5,
  7.5° at Mach 1.6, 10° at Mach 1.4) and by 25 to 40% at Mach × half-angle 0.4 to 0.5.

**So:** OpenVSP's Wave Drag is a good estimate from about Mach 1.2 (above shock detachment; linear theory is
singular at Mach 1) up to where Mach × nose or forebody angle reaches about 0.2 and Mach × wing thickness ratio about
0.5 to 0.6; for a slender supersonic aircraft (forebody half-angles of 5 to 7°, 3 to 4% wings) that is roughly Mach
1.2 to 2.5. Above about Mach 2.5 to 3 its own Mach-plane cutting adds a further over-prediction on bodies, and the
answer is also zero-lift only (drag due to lift comes from VSPAERO). Within that range, use enough rotation sections
for the wing aspect ratio and check the slice count when the configuration has blunt bases or inlets.

## Aero database and design of experiments

**Aero Database (VSP)** takes the three analyses' results and, per VSPAERO case, writes CL = CLi, CD = CDi +
CD_wave(M) + CD_parasite(M), CY = CSi, Cl = -CMix, Cm = CMiy, Cn = -CMiz (plus optional viscous lists, e.g. from
`ntop_newtonian`'s viscous correction). Wave and parasite drag are interpolated linearly in Mach between the Mach
numbers they ran at; wave drag is 0 below Mach 1. VSPAERO's own viscous estimate (CDo) is left out, Parasite Drag
replacing it. The driver checks that the three used the same Sref and that the wave and parasite runs cover
VSPAERO's Mach numbers. **Variable Names** and **Variable Values** become the first columns, and **Append To**
collects the rows of repeated runs into one design-of-experiments table.

**`examples/X-59 Aero Database.ntop`**: the starter's geometry (unchanged), six component blocks (fuselage assembly
in no VSPAERO set; wing, canards, horizontal tail, T-tail and the fin cut below the T-tail, all thin), the model, and
a shared envelope: VSPAERO at Mach 0.6 to 1.8 (5 values) x alpha -2 to 8 degrees (6) x beta 0 and 4 degrees (60
cases, fixed wake), Wave Drag at Mach 1.2, 1.5, 1.8 (30 slices, 30 rotation sections), Parasite Drag at the five
Mach numbers at 16,764 m (55,000 ft), Sref 61.5 m², bref 9.04 m, cref 8.72 m (the fitted baseline wing), moments
about the starter's CG station (x = 18.27 m). **Wing Span** and **Wing Inner Sweep** are notebook inputs, recorded
in every row. One run takes about 23 minutes on 16 threads, nearly all of it VSPAERO's subsonic cases.

Baseline results at beta 0 (derivatives per degree):

| Mach | CL_alpha | Cm_alpha | CD at alpha 0 = induced + wave + parasite | CY_beta | Cl_beta | Cn_beta | best L/D |
|---|---|---|---|---|---|---|---|
| 0.6 | 0.0294 | -0.0091 | 0.0117 = 0.0013 + 0 + 0.0104 | -0.0030 | -0.0020 | +0.0030 | 9.2 |
| 0.9 | 0.0306 | -0.0096 | 0.0112 = 0.0014 + 0 + 0.0098 | -0.0030 | -0.0020 | +0.0030 | 9.6 |
| 1.2 | 0.0331 | -0.0129 | 0.0279 = 0.0044 + 0.0142 + 0.0094 | -0.0039 | -0.0023 | +0.0039 | 4.4 |
| 1.5 | 0.0322 | -0.0135 | 0.0240 = 0.0040 + 0.0109 + 0.0091 | -0.0038 | -0.0022 | +0.0038 | 4.6 |
| 1.8 | 0.0308 | -0.0133 | 0.0216 = 0.0038 + 0.0089 + 0.0089 | -0.0039 | -0.0021 | +0.0038 | 4.8 |

**VSPAERO Center of Pressure (VSP)** gives the aerodynamic center (neutral point) of each alpha sweep, x_ac =
x_ref - (dCMy/dCFz) cref: 20.95 m at Mach 0.6 and 0.9, moving aft to 21.62, 21.88 and 21.98 m at Mach 1.2, 1.5 and
1.8. The starter's own AeroSandbox estimate is 20.90 m; with the CG at 18.27 m both give a static margin of about
30% of the mean aerodynamic chord subsonic. The block also gives each case's center of pressure, x_cp = x_ref -
CMy cref / CFz, which runs away as the normal force goes to zero.

Sideslip derivatives are from beta 0 and 4 degrees at alpha 2. The aircraft is statically stable in pitch, yaw and
roll throughout, more so in pitch above Mach 1 as the aerodynamic centre moves aft. The parasite drag (0.009 to
0.010) is OpenVSP's component buildup on the fitted components (wing t/c 3.0%, tails 4.0%, fuselage fineness 20,
as the starter's NACA 0003 and 0004 sections and fuselage give). The fuselage's lift is not in these numbers (see
above), and the wave drag carries the uncertainty of the nacelle's inlet and base steps.

**`examples/x59_doe.py`** runs the notebook headlessly once per design with `ntopcl --trustnotebook -j inputs.json`
(the JSON sets Wing Span and Wing Inner Sweep) and collects every run's rows in `runs/x59_doe.csv`: one table with
Wing_Span_m, Wing_Inner_Sweep_deg, Mach, AoA, Beta, Re_cref, CL, CD, CY, Cl, Cm, Cn and the drag parts. Designs: the
baseline (9 m, 68.3°), span 8 and 10 m, inner sweep 65 and 71°. `--keep` adds designs to an existing table and
`--summary` prints the slopes per design.

Result (`examples/x59_doe_results.csv`, 300 rows: 5 designs x 60 cases; about 23 minutes per design on 16 threads,
nearly all of it VSPAERO). Derivatives per degree at beta 0, all on the same Sref 61.5 m²:

| Wing span | Inner sweep | CL_alpha, Mach 0.6 | Cm_alpha, Mach 0.6 | CL_alpha, Mach 1.5 | Cm_alpha, Mach 1.5 | CD at alpha 0, Mach 1.5 |
|---|---|---|---|---|---|---|
| 8 m | 68.3° | 0.0255 | -0.0067 | 0.0276 | -0.0102 | 0.0229 |
| 9 m | 68.3° (baseline) | 0.0294 | -0.0091 | 0.0322 | -0.0135 | 0.0240 |
| 10 m | 68.3° | 0.0348 | -0.0132 | 0.0377 | -0.0183 | 0.0256 |
| 9 m | 65° | 0.0306 | -0.0068 | 0.0342 | -0.0116 | 0.0245 |
| 9 m | 71° | 0.0281 | -0.0110 | 0.0303 | -0.0150 | 0.0235 |

More span adds wing area behind the (fixed) CG: more lift per degree, more pitch stability and more zero-lift drag
(wetted area and volume). More inner-panel sweep moves area aft and makes the wing more slender: less lift per
degree, more stability, less wave drag. These are the trends to expect, from a database whose geometry columns came
straight from nTop's parametric model.

## Coverage of OpenVSP's options

Every option of OpenVSP 3.56's VSPAERO, Parasite Drag and Wave Drag analyses was listed from OpenVSP itself
(`dev/inventory.txt`) and is a block input, out of scope as you listed, or not applicable (with the reason).

**VSPAERO**

| OpenVSP option | Here |
|---|---|
| Alpha, Beta, Mach, Re_cref start / end / points; Vinf, Rho, Vref (manual flag), Machref | Flow Condition |
| Sref, bref, cref; Xcg, Ycg, Zcg | Analysis |
| NCPU, fixed wake, wake iterations, wake nodes, far field, ground effect, freezing, implicit wake, wake relaxation, convergence factors, core size, multipole, preconditioners, quad tree, `-save`/`-restart` | Solver Settings |
| Thick and thin geometry sets | each component's Set (Thick, Thin or None) |
| Find body wakes | Analysis: Body Wakes |
| Wing trailing-edge wakes | automatic: OpenVSP marks every wing trailing edge |
| X-Z symmetry; 2D FEM and Tecplot files; Sref/cref from the model; CG calculation; stall model; propellers, rotors, actuator disks, noise; stability types and unsteady time stepping; Cp slicer; control surfaces | out of scope (your list) |
| NRef (refinement levels), cull fraction, degenerate-region merging, intersection cleaning | OpenVSP's defaults |
| `-fem`, `-interrogate` with survey points | not implemented (say if you want them) |

**Parasite Drag**: geometry set, model length unit, Sref, laminar and turbulent Cf equations, speed sweep and unit,
altitude, temperature offset, specific heat ratio and sub-component export are inputs; each component's form factor
equations, manual form factor, % laminar, Q, roughness, Te/Tw, Taw/Tw and incorporated generation are inputs of
Parasite Drag Component (VSP). Out of scope: atmospheres other than US Standard 1976, excrescences, Sref from a wing.

**Wave Drag**: Mach sweep, number of slices, number of rotation sections, Sref, symmetry and geometry set are
inputs. Out of scope: inflow and outflow faces.

## Rebuilding and running headlessly

```bash
sh build/py.sh build/build_blocks.py
sh build/py.sh build/x59_example.py --run
sh build/py.sh examples/x59_doe.py
sh build/py.sh validation/wave_drag_mach.py
```

`build/py.sh` runs the `nTop_pkgs` environment (it has `implicit_io`, which writes the `.ntop` files); the drivers
run in the Python of `NTOP_PYTHON_PATH`, and the example, DOE and study need `NTOP_PYTHON_PATH` and `VSP_PATH`
set. ntopcl may write its output and then not exit after notebooks with Run Command inside custom blocks; the
runners stop it once the output file is written. A signed-in nTop is needed.

## Known limits

* Runs with the same Working Directory and Case Name overwrite each other (the DOE runs one design at a time).
* A component block's body must have one closed contour per cut (the largest is used when there are several, with
  a note in the log); a fuselage is sliced along +x and a wing's span axis must be normal to x.
* The fuselage fit follows the body between sections with OpenVSP's own skinning, so abrupt steps (an inlet face)
  are smoothed: 0.43 m at the X-59's nacelle inlet.
* Surface pressures are read from VSPAERO's binary `.adb`; every run checks the last case against VSPAERO's Tecplot
  output. Integrated over the surface, these pressures do not give VSPAERO's forces (it computes those from its
  vortex edges): use them for viewing and as viscous-correction edge conditions.
* The HL-20 files in this folder (`h20*`, `hl20.*`) belong to the earlier mesh-based workflow, which these blocks
  replace.
