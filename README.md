# nTop_vsp

nTop custom blocks that turn **nTop implicit bodies into an OpenVSP model** and run aerodynamic analyses on one
shared **Flight Condition** (Mach numbers and altitudes through the US Standard Atmosphere 1976, angles of attack and
sideslips, and the reference area, span, chord and moment point): OpenVSP 3.56's **VSPAERO**, **Wave Drag** and
**Parasite Drag**, and, computed entirely in nTop on the implicit body itself, **Modified Newtonian** pressure and a
**Viscous Correction** (skin friction). **Aero Database (VSP)** sums any mix of them, each over its own Mach range,
into a **six-degree-of-freedom aero database** (CL, CD, CY, Cl, Cm, Cn against Mach, altitude, alpha and beta);
**Center of Pressure (VSP)** gives its centers of pressure and aerodynamic centers as positions, and **Polar Plot
(VSP)** plots it with matplotlib. Because the geometry stays parametric in nTop, its
inputs can be independent variables of the database too: `examples/x59_doe.py` runs a design of experiments over
the wing span and sweep of an X-59-style aircraft.

Test geometry: the X-59 hackathon starter in `../x59-hackathon-starter` (MIT License, Copyright (c) 2026
matthewmuellernTop; notice in `examples/X-59_STARTER_LICENSE.txt`), a parametric nTop model whose fuselage, wing,
canards and tails are separate implicit bodies. `build/x59_example.py` reads it from there.

## How it works

```
implicit bodies ─► OpenVSP Fuselage / Wing (VSP) ─► component line ─┐
  (lists; nTop Slice Body → section contours, CLI)                  ├─► OpenVSP Model (VSP) ─► <case>.vsp3
implicit bodies ─► ...                                             ─┘    (run_model.py: fit + build)
                                                                                    │
Flight Condition (VSP): Mach + altitude (US Std Atmosphere 1976), alpha, beta,     │
                        Sref, bref, cref, moment point ─── into every analysis ────┤
                                                                                    │
   ┌────────────────────┬───────────────────┬─────────────────────┬───────────────────┴───┐
VSPAERO Analysis     Wave Drag          Parasite Drag      Modified Newtonian    Viscous Correction
(vortex lattice /    (area rule)        (component         (in nTop, on the      (in nTop, on the
 panel method)                           buildup)           implicit body)        implicit body)
   │                    │                   │                     │                     │
VSPAERO              Wave Drag          Parasite Drag      Coefficients          Coefficients
Contribution         Contribution       Contribution       Contribution          Contribution
   └────────────────────┴───────────────────┴──────────┬──────────┴─────────────────────┘
                                                        ▼
                         Aero Database (VSP) on the Flight Condition ─► CSV / DOE table, readers
                                                        ▼
                         Center of Pressure (VSP) ─► points       Polar Plot (VSP) ─► PNG + window
```

The OpenVSP analyses, the Flight Condition and the database write their inputs to a settings file (`key = value`,
SI units) and run a driver of the `vsp_ntop` Python package (`python -m vsp_ntop.<module> <settings file>`) with
nTop's **Run Command**; the drivers run OpenVSP (`vspscript.exe`) and VSPAERO and write plain CSV files, which reader
blocks load with Import Table and Import Scalar Point Map. A failed run stops the notebook with the driver's error,
and every run leaves `<Case Name>.ntvsp.log` in the Working Directory. Modified Newtonian and the Viscous Correction
run no program: they are nTop block graphs evaluated once per flight condition, which write their coefficient table
with Export Text.

* `blocks/` — the custom blocks, one `.ntop` file each (File > Import, or drag into a notebook).
* `examples/X-59 Aero Database.ntop` — the X-59 starter's geometry with VSPAERO, Wave Drag and Parasite Drag on one
  flight condition, floating comments throughout; `examples/x59_doe.py` — the design of experiments over it.
* `examples/X-59 Mixed Analyses.ntop` — the same geometry, VSPAERO and Wave Drag below Mach 2.5, Modified Newtonian
  from Mach 2.5 and the Viscous Correction in place of Parasite Drag, on three flight points at three altitudes.
* `vsp_ntop/` — the drivers, an installable package (`pyproject.toml`; Python 3.8+; `ambiance` for the standard
  atmosphere, matplotlib for the plots): `run_flight.py` and `flight.py` (the flight condition), `run_model.py` with
  `slices.py` (contour reader),
  `fit.py` (section fitting) and `vspmodel.py` (OpenVSP script and `.vsp3` writing); `run_vspaero.py`,
  `run_wave.py`, `run_parasite.py`, `run_database.py` (the buildup), `run_plot.py`; `vspgeom.py` and `adb.py` read VSPAERO's panel file and
  binary results.
* `build/` — writes the blocks and the examples (`build_blocks.py` with `model_blocks.py`, `aero_blocks.py` (flight
  condition, Newtonian, viscous, coefficients), `vspaero_blocks.py`, `drag_blocks.py`, `buildup_blocks.py`;
  `newtonian.py`, `viscous.py` and `ops.py` are the Modified Newtonian and skin-friction block graphs, brought over
  from the former `ntop_newtonian` project; `x59_example.py`, `mixed_example.py`; `x59.py` derives notebooks from the
  starter; `ntopfile.py` writes `.ntop` files).
* `validation/wave_drag_mach.py` — the wave drag study below (results in `validation/results/`);
  `validation/newtonian/` — the Modified Newtonian and viscous-correction validation report from `ntop_newtonian`;
  `validation/rasaero/` — which combination of these analyses best reproduces RASAero II on a finned rocket (Mach
  1.5 to 5, alpha -15 to 15, beta -4 to 4 degrees): `validation/rasaero/REPORT.md`.

## Setup

1. Install OpenVSP 3.56: unpack the Windows 64-bit build (for example into `OpenVSP-3.56.0-win64/` here; it is not
   in the repository).
2. Install the drivers into a Python 3.8 or newer: `pip install -e .` in this folder (or `pip install .`); this
   also installs `ambiance` (the US Standard Atmosphere 1976 of Flight Condition) and matplotlib (Polar Plot). A conda
   environment's Python works without activating it.
3. Set two environment variables (Windows: System Properties > Environment Variables), then restart nTop so it
   sees them:
   * `NTOP_PYTHON_PATH`: that Python's executable, e.g. `C:/Users/me/anaconda3/python.exe`;
   * `VSP_PATH`: the OpenVSP folder, the one holding `vspscript.exe` and `vspaero.exe`;
   * `RASAERO_PATH` (only for RASAero II (VSP)): the RASAero II folder, the one holding `RASAero II.exe`. RASAero II
     (Rogers Aeroscience, free download from rasaero.com, Windows) is not in this repository: its site offers it free
     but grants no right to redistribute it, so install it yourself.

The blocks read both with nTop's **Get Environment Variable** block, so they have no Python, driver or OpenVSP path
inputs; a missing variable stops the notebook with nTop's "Environment variable ... is not set". The only path
input left is each driver block's **Working Directory**, empty by default: then the files go to a folder `ntop_vsp`
in the system temporary folder (%TEMP%). Modified Newtonian, Viscous Correction and Aero Coefficients write next to
their Flight Condition's files. nTop asks once to trust Run Command
blocks.

## Blocks

| Block | Output | Purpose |
|---|---|---|
| OpenVSP Fuselage (VSP) | text | a list of implicit bodies as OpenVSP Fuselages: names, slicing, sections, clustering, tessellation, set |
| OpenVSP Wing (VSP) | text | a list of implicit lifting surfaces as OpenVSP Wings: names, span axes (one per body), symmetry (none, x-z, x-y, both), extension into the body, force accounting (whole wing or exposed strips), sections, airfoils, tessellation, clustering, set |
| OpenVSP Model (VSP) | model path | builds `<case>.vsp3` from a list of component lines |
| OpenVSP Model Check (VSP) | list of reals | per component: largest and mean distance from the OpenVSP surface to the nTop contours, sections, slices |
| OpenVSP Model Mesh (VSP) | mesh | OpenVSP's tessellation of the model, to view against the bodies |
| Flight Condition (VSP) | results path | Mach numbers paired with altitudes (US Standard Atmosphere 1976), angles of attack, sideslips; reference area, span, chord, moment point |
| VSPAERO Solver Settings (VSP) | text | threads, wake model and nodes, far field, ground effect, freezing, implicit wake, convergence factors, preconditioners, quad tree, restart file, convergence check |
| VSPAERO Analysis (VSP) | results path | model, flight condition, Mach range, body wakes, solver settings |
| VSPAERO Coefficient (VSP) | list of reals | any column of VSPAERO's polar, one value per case |
| VSPAERO Surface Pressure (VSP) | Scalar Point Map | Cp or p - p_inf of one case at panel centroids, nodes or both (thin surfaces: the jump across the sheet) |
| Wave Drag (VSP) | results path | model, flight condition (its supersonic Mach numbers), slices, rotation sections, symmetry, geometry set |
| Wave Drag Result (VSP) | list of reals | CD_wave (or Mach, counts, Sref), one value per Mach number |
| Wave Drag Area Distribution (VSP) | list of reals | x/L or slice areas of one rotation section at one Mach number |
| Parasite Drag Component (VSP) | text | one component's form factor equations, % laminar, Q, roughness, temperature ratios |
| Parasite Drag (VSP) | results path | model (metres), flight condition (one run per flight point), Cf equations, components |
| Parasite Drag Result (VSP) | list of reals | total CD, f, Re/L and the flight condition, one value per flight point |
| Parasite Drag Component Values (VSP) | list of reals | one column of the buildup table at one flight point, per component |
| Body Lift (VSP) | results path | the fuselages' normal force from their OpenVSP area distribution: slender-body nose lift, optional viscous crossflow |
| RASAero II (VSP) | results path | RASAero II run on a rocket's OpenVSP model (design fitted to its surfaces) at every condition, as a coefficient table |
| Modified Newtonian (VSP) | results path | an implicit body's Modified Newtonian forces and moments at every condition of a flight condition in a Mach range, computed in nTop |
| Modified Newtonian on Mesh (VSP) | results path | the same Modified Newtonian pressure on the implicit body, with the force and moment integrals summed over nTop's mesh of it (much faster) |
| Modified Newtonian Pressure (VSP) | scalar field | the Modified Newtonian Cp of one condition, to view on the body |
| Viscous Correction (VSP) | results path | an implicit body's skin-friction forces and moments at every condition (standard-atmosphere edge state), computed in nTop |
| Aero Coefficients (VSP) | results path | six lists of coefficients from anywhere, one value per flight condition, as a coefficient table |
| Aero Contribution (VSP) | text | the common contribution block: source, results, name, terms, use, Mach range, scale, all with defaults; the blocks below wrap it |
| VSPAERO Contribution (VSP) | text | a VSPAERO run for the database: inviscid or total coefficients, which forces to use, Mach range, scale |
| Wave Drag Contribution (VSP) | text | a Wave Drag run for the database: CD_wave(Mach) along the freestream, Mach range (default from 1), scale |
| Parasite Drag Contribution (VSP) | text | a Parasite Drag run for the database: CD(Mach, altitude) along the freestream, Mach range, scale |
| Coefficients Contribution (VSP) | text | a Modified Newtonian, Viscous Correction, Body Lift or Aero Coefficients table for the database |
| RASAero Contribution (VSP) | text | a RASAero II run for the database |
| Aero Database (VSP) | results path | sums a list of contributions at every condition of a flight condition, on its reference; design variables; appends to a DOE file |
| Aero Database Coefficient (VSP) | list of reals | one column of the database (CL, CD, CY, Cl, Cm, Cn, L/D, CA, CN, Mach, Altitude, AoA, Beta) |
| Aero Database Table (VSP) | table | the whole database, each contribution's share included |
| Center of Pressure (VSP) | list of points | from an Aero Database: each condition's center of pressure, or each alpha sweep's aerodynamic center (neutral point), as positions with length units |
| Polar Plot (VSP) | image path | one database quantity against another for a list of databases (matplotlib); opens in a window by default |

Every length is sent to OpenVSP in metres, so reference areas and lengths are nTop quantities with units.

## From implicit to OpenVSP

![OpenVSP model of the X-59 starter](examples/x59_openvsp_model.png)

**Lists of bodies.** Each component block takes a list of bodies (**Bodies**), a list of **Names** (one per body, or
one name that the driver numbers) and, for wings, a list of **Span Axes** (one per body, in order; nTop stops with
"Input list length mismatch" otherwise); every other input applies to all of them. nTop slices them one by one
(list processing), and the driver makes one OpenVSP component of each. The X-59 example has one Wing block for the
wing, canards and both horizontal tails, and another for the fin, which is not mirrored.

**Slicing in nTop.** A component block cuts each body with nTop's **Slice Body** (planes normal to x for a fuselage,
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
stations; the surface is rotated about x for a canted or vertical span axis. **Symmetry** is OpenVSP's planar
symmetry: about the x-z plane (left and right), the x-y plane (top and bottom: a dorsal fin mirrored into a ventral
one) or both (cruciform tails); with any of them only the part of the body on the +span side of the plane through
the origin normal to its span axis is used. **Trailing Edge: Sharp** shears a blunt base onto the chord line so
that VSPAERO's wake leaves one edge. OpenVSP marks every wing trailing edge as a wake: no wake input.

**Wings in a body.** **Extend Into Body** continues each wing inboard to the plane through the origin normal to its
span axis (for fins on a body along the x axis: to the axis), with the first section's leading- and trailing-edge
lines, thickness and root airfoil. VSPAERO then solves panels through the body, which removes the free root edge (and
its tip-like loss of lift) of a wing that ends at the body surface. **Force Accounting** says which of those loads
count: Whole Wing (every strip: the inboard part's lift stands for the body's carry-over, and over-predicts it on a
slender body), or Exposed Strip (only outboard of the wing body's own root: give the wing trimmed at the fuselage
surface). For Exposed Strip, OpenVSP Model records each wing's root (`<case>.ntvsp.wings.csv`) and VSPAERO Analysis
takes the inboard strips' forces and moments off every case from VSPAERO's spanwise loads (`.lod`; the whole-wing
polar is kept as `.ntvsp.polar_whole.csv`). On a cruciform rocket (`validation/rasaero/REPORT.md`) the exposed strips'
lift is within -6% to +11% of SU2's fin lift, where fins ending at the body give 59 to 76% of it; Whole Wing puts 30 to
38% of the lift on the inboard strips, 2 to 6 times SU2's carry-over. Opposite wings that meet at the origin plane
become one continuous wing in VSPAERO, and four fins meeting on the axis are two crossing thin surfaces: VSPAERO
returned no finite result at Mach 1.5 for them (Mach 2 to 5 converged).

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
or +z, where they come in the world's cyclic axes (y-z, z-x, x-y), and along exactly -x, -y or -z, where they come in
(z, y), (x, z) and (y, x) (again world axes whose cross product is the slicing direction; found on the four fins of a
cruciform rocket, whose -y fin was mis-fitted before). The driver handles all of them; the blocks place the frame
origin outside the body.

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
* **Parallel cases.** VSPAERO solves its cases one after another, each on all its threads (**Processors**). nTop runs
  independent blocks at the same time (two Run Commands in one notebook started together in a test), so to spread a
  large sweep over cases, split it into two or more VSPAERO Analysis blocks with adjacent Mach ranges (**Mach Min**,
  **Mach Max**), separate Case Names and the threads shared out, and give each its own VSPAERO Contribution: the
  database takes each where it has values. For several designs, run separate `ntopcl` processes.
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

## Flight condition

**Flight Condition (VSP)** is the one place flight conditions and the reference are defined; every analysis and the
Aero Database take its output and use what they need.

* **Mach** and **Altitude** are lists paired item by item into flight points (one altitude serves every Mach
  number). The atmosphere at each altitude (geometric, -16,400 to 265,800 ft) is the US Standard Atmosphere 1976 as
  implemented by the `ambiance` package: temperature, pressure, density, speed of sound and viscosity; the velocity,
  dynamic pressure and Reynolds number (per metre and on cref) follow from the Mach number. Enter the altitude in feet
  (nTop takes any length unit).
* **Angle of Attack** and **Sideslip** are lists; every flight point is run at each. The conditions are flight point x
  sideslip x angle of attack (angle of attack fastest).
* **Reference Area**, **Reference Span**, **Reference Chord** and **Moment Reference Point** travel with the
  conditions: every coefficient of every analysis and of the database is over them and about that point.

What each analysis takes from it: VSPAERO its distinct Mach numbers, angles and sideslips (all combinations, in one
run), the Reynolds numbers of its flight points and the reference; Wave Drag its supersonic Mach numbers and Sref;
Parasite Drag each flight point (one OpenVSP run per Mach number and altitude) and Sref; Modified Newtonian the
distinct (Mach, alpha, beta) combinations and the reference; the Viscous Correction every condition with its
pressure and temperature. VSPAERO, Modified Newtonian and the Viscous Correction also have a **Mach Min** and **Mach
Max**, so each runs only where it is meant to (VSPAERO below Mach 3, Newtonian from Mach 3, say).

## Body lift and RASAero II

**Body Lift (VSP)** gives a fuselage's normal force from the OpenVSP model alone: the fuselages' surface is sampled in
rings at planes normal to x (OpenVSP's fuselage sections are such planes), giving the section area S(x).
Slender-body (Munk) potential lift puts dCN/dx = (dS/dx / Sref) sin(2 alpha) cos(alpha / 2) where the area grows (the
nose); optional Allen and Perkins viscous crossflow adds eta Cdn (2 r / Sref) sin^2(alpha) along the body, with
Jorgensen's (NASA TR R-474) Cdn against crossflow Mach number and eta against fineness. Both act along the crossflow at
the total angle of attack, so sideslip is included; no axial force. It is meant to go with Exposed Strip wings. On the
cruciform rocket, Exposed Strip fins with Body Lift (slender body and crossflow) put the centre of pressure within 21
mm of inviscid SU2's at Mach 2 to 5 (CN 6 to 12% low: the carry-over onto the body between the fins is in no block).
OpenVSP's tessellation-based Planar Slice was not usable for S(x): 80% off near a nose.

**RASAero II (VSP)** runs RASAero II (Rogers Aeroscience; `RASAERO_PATH`) on a rocket's OpenVSP model, with no rocket
description of its own: `vsp_ntop.run_rasaero` fits a RASAero II design to the model's surfaces (nose: the closest
of Conical, Tangent Ogive and Von Karman Ogive, to where the radius reaches the body tube's; body tube; boattail where
the radius falls by 5%; one set of alike trapezoidal fins, mirrored copies counted, with hexagonal sections whose
thickness and bevels come from the root section), starts RASAero II, opens the design and runs its Run Test through
the window with Microsoft Active Accessibility (`vsp_ntop/rasaero_ui`, the scripts of hypanel's RASAero comparison),
at every whole degree of total angle of attack up to the largest needed (its angle box drops decimal points) and
each altitude (as the launch site altitude), and closes it. RASAero II is axisymmetric: at (alpha, beta) its normal
force is taken at the total angle of attack along the crossflow, with its centre of pressure on the axis, on the
flight condition's reference. RASAero II opens a window while it runs, so it needs a desktop session. On the rocket
of `validation/rasaero` the fitted design gives CA within 1%, CN within 2% and the moments within 1.5% of hand-written
design runs (50 s). Its own methods are not published (its site describes a "Rogers Modified Barrowman" method
calibrated against wind-tunnel data), so they are not reimplemented here; Body Lift is the public-method counterpart
for the body.

## Modified Newtonian and the viscous correction, in nTop

These two analyses work on the implicit body itself, entirely with nTop blocks (no Run Command, no Python): the
method of the former `ntop_newtonian` project, now built here (`build/newtonian.py`, `build/viscous.py`) and
evaluated once per flight condition through nTop's list processing.

* **Modified Newtonian (VSP)**: Cp = Cp_max max(-n.V, 0)^2 with n the unit gradient of the body's field, so any
  implicit works; each force and moment component is one Surface Integral. Cp_max is the stagnation pressure behind
  a normal shock above Mach 1 (Rayleigh pitot) and the isentropic one below it (both 1.839 at Mach 1), 2
  (Classical), or **User Defined**: the **Cp Max Field** input, a scalar field, so Cp_max may differ over the body
  (one value on the fins, another on the nose); a constant field is one number everywhere. No self-shadowing: correct
  for convex bodies, a part in another's shadow is over-loaded.
* **Modified Newtonian on Mesh (VSP)**: the same pressure, computed on the implicit (its gradient normal, the same
  Cp_max choices including the Cp Max Field), with the integrals done over nTop's mesh of the body (Mesh from Implicit
  Body at the **Mesh Tolerance**, sharp edges kept) instead of Surface Integrals: nTop meshes the body and exports it,
  the driver `vsp_ntop.run_newtonian_mesh` lists the faces' centroids, nTop evaluates the body's gradient and the Cp
  Max Field there and exports them as point maps (Evaluate Field Gradient, Point Map, Export Point Map), and the driver
  sums -Cp S_f and -Cp (r_f - c) x S_f over the faces (S_f the face's outward vector area) at every condition (numpy).
  It writes the same coefficient table. See the run time and accuracy below.
* **Viscous Correction (VSP)**: skin friction at the inviscid edge state, hypanel's `viscous.correct`: edge pressure
  from the Modified Newtonian Cp of each condition (or the freestream, Cp = 0, for flat-plate friction), edge total
  pressure behind a normal shock above Mach 1 and isentropic below it, total temperature conserved; Blasius with
  Eckert's reference temperature (laminar) and White's van Driest II (turbulent), sudden transition at a running-length
  Reynolds number, adiabatic or fixed wall; shear along the surface flow, faded out between 30 and 60 degrees of
  leeward inclination. The freestream pressure and temperature of each condition come from the flight condition.
* **Modified Newtonian Pressure (VSP)** shows one condition's Cp on the body.

Axes: x along **Aft Axis** (the model's +x), z toward **Up Axis**, y = z x x: CA, CY, CN along them, moments about
them (the axes of VSPAERO too). Each writes a coefficient table, `<Flight Condition>_<Case Name>.ntvsp.coefficients.csv`,
with Export Text, for **Coefficients Contribution (VSP)**. Accuracy is nTop's Surface Integral (about 1e-4 on smooth
bodies, 0.1 to 1% across sharp edges, **Integration Tolerance**); validation against exact results and wind-tunnel
data, beside hypanel's implementation, is in `validation/newtonian/REPORT.md`: sphere CD within 8e-4, 15 degree cone
CA and CN within 0.002 and 0.010 from 0 to 60 degrees, the Apollo capsule's pressures at hypanel's RMS, skin friction
within 0.5% of hypanel's on a plate and a cone.

On the flight condition, a cone checked against the former single-condition blocks gives the same CA, CN and Cm at
alpha 0 and within 0.15% at alpha 4 degrees (the Surface Integral's tolerance), and finite, flat-plate-like friction
at Mach 0.6. Not carried over: the on-mesh viscous blocks (sums over VSPAERO's panel centroids instead of Surface
Integrals). Evaluating a list of fields at a list of points makes a list of lists, which nTop does not allow, so they
cannot run over a flight condition.

**Run time.** Each condition costs six Surface Integrals for each of Newtonian and friction, about 1 to 3 s each on
simple bodies, more on a full aircraft and much more on thin parts (fins); use Mach Min and Max, and few angles, where
Newtonian applies. Modified Newtonian on Mesh costs one mesh, one evaluation of the gradient at the face centroids
and a sum per condition. On the hobby rocket of `validation/rasaero` (1.2 m, four 4 mm fins), four Newtonian
runs (the whole rocket at 611,000 faces, the body, the fins, and the rocket again) over 330 conditions took 99 s in
all; Surface Integrals on the body alone took 750 s for 4 conditions. On the body at Mach 3 (alpha 0 and 8, beta 0 and
3 degrees) a 4 mm Mesh Tolerance (118,000 faces) agrees with the Surface Integral block within 0.5% in CN, CY and the
moments and to four digits in CA, and gives exact zeros where the Surface Integral's noise gives 7e-4; a User Defined
Cp Max Field of 1.5 scales every coefficient by exactly 1.5 / Cp_max(Modified).

## Aero buildup: mixing analyses in one database

Every source of forces and moments enters the database the same way, as a **contribution**, collected with other
contributions in a list and connected to **Aero Database (VSP)**, which also takes the Flight Condition (its
conditions and reference). Every contribution is made by one block, **Aero Contribution (VSP)**: Source, Results,
Name, Terms, Use, Mach Min, Mach Max and Scale, all with defaults. It checks that Results is connected and writes one
`contribution.<key> = <value>` line per setting (names and paths may hold any character but a line break). The
per-source blocks (VSPAERO, Wave Drag, Parasite Drag and Coefficients Contribution) are wrappers: custom blocks that
call it with their source fixed, their own defaults and only the inputs that matter for that source. The database
checks that every item of its Contributions list is one.

| Analysis | Contribution block | What the database takes |
|---|---|---|
| VSPAERO Analysis | VSPAERO Contribution | forces and moments of every case, inviscid (default) or with VSPAERO's own friction; at each condition the run's Reynolds number nearest the condition's |
| Wave Drag | Wave Drag Contribution | CD_wave against Mach (from Mach 1 by default), along the freestream |
| Parasite Drag | Parasite Drag Contribution | Total_CD against Mach at the condition's altitude, along the freestream |
| Modified Newtonian, Viscous Correction, Aero Coefficients | Coefficients Contribution | their coefficient tables (Newtonian's by Mach, alpha and beta; friction's and Aero Coefficients' by altitude too) |

Each contribution applies over the Mach numbers its analysis ran (Wave Drag: none subsonic; VSPAERO, Newtonian and
friction: their own Mach Min to Max) and within the contribution's own **Mach Min** and **Mach Max** (from Min up to,
not including, Max). It has a **Scale** and, where it has more than drag, a **Use** choice: all forces and moments,
drag only (the force component along the freestream), or everything but drag. So, for example:

* VSPAERO below Mach 3 and Modified Newtonian from Mach 3: give VSPAERO Analysis Mach Max 3 and Modified Newtonian
  Mach Min 3 (the analyses then only compute there, and the contributions follow).
* The Viscous Correction instead of Parasite Drag: connect its contribution (Use: Drag only, or All to keep its small
  lift and moments) and leave Parasite Drag out of the list.
* Lift and moments from one source and drag from others: Use: All but drag on the first.

**What the database does** (`vsp_ntop/run_database.py`). At each condition of the flight condition, every
contribution that applies is taken at that condition if it has it, else interpolated (multilinear between the
neighbouring Mach numbers, angles and sideslips; Wave and Parasite Drag linearly in Mach). It is turned into
body-axis force and moment coefficients (x aft, y right, z up), made dimensional with its own reference, moved to the
flight condition's moment point (M + (r - r_ref) x F) and normalised by its reference. Wave and Parasite Drag carry
no moments and act through the moment point. The sum gives CL, CD and CY in wind axes (VSPAERO's: freestream (cos a
cos b, -sin b, sin a cos b), CD along it, CL along (-sin a, 0, cos a), CY along (cos a sin b, cos b, sin a sin b)) and
Cl, Cm, Cn in aircraft body-axis signs (right wing down, nose up, nose right). `<case>.ntvsp.database.csv` has the
design variables, Mach, Altitude (m), AoA, Beta, the totals, L/D, CA, CN, and each contribution's own six
coefficients in columns named after it (CD_Wave_Drag, Cm_Newtonian, ...).

**Center of Pressure (VSP)** takes the Aero Database: each condition's center of pressure, x_cp = x_ref - Cm cref /
CN, or each alpha sweep's aerodynamic center (neutral point), x_ac = x_ref - (dCm/dCN) cref at each Mach number,
altitude and sideslip, from the totals, as points with length units on the line through the moment point along x.

**Checks.** The driver reproduced the earlier VSPAERO + Wave + Parasite database of the X-59 to 1e-12 in every
coefficient, and VSPAERO's coefficients re-expressed about another moment point and reference, in body and in wind
axes, come back unchanged (to 1e-16). A database mixing VSPAERO and Wave Drag below Mach 1.6, Modified Newtonian of a
cone and its friction (drag only) from Mach 1.6, and wind-axis coefficients on another reference and moment point
matched a hand calculation to 5e-13 at every condition, including Mach 2.5, interpolated between Newtonian at Mach 2
and 3.

**Polar Plot (VSP)** plots one database column against another (X and Y dropdowns: CL, CD, CY, Cl, Cm, Cn, L/D,
CA, CN, AoA, Beta, Mach, Altitude) for a list of databases, one curve per combination of the conditions that do not
vary along the curve (**Curves Along**: angle of attack, sideslip or Mach), optionally only at zero sideslip or zero
angle of attack. Curves of one result are shaded light to dark by condition and results differ by line style (by
colour when each has one curve). It saves `<Case Name>.png` and by default opens the plot in matplotlib's window,
which zooms, pans and saves; nTop does not wait for it. The window is a separate process: it stays open while nTop
is, and a headless `ntopcl` run closes it when it exits.

## Aero database and design of experiments

**`examples/X-59 Aero Database.ntop`**: the starter's geometry (unchanged), three component blocks (the fuselage
assembly in no VSPAERO set; one Wing block for the wing, canards, horizontal tail and T-tail; one for the fin, cut
below the T-tail; all thin), the model, and one Flight Condition: Mach 0.6, 0.9, 1.2, 1.5 and 1.8 at 55,000 ft x alpha
-2 to 8 degrees (6) x beta 0 and 4 degrees (60 conditions), Sref 61.5 m², bref 9.04 m, cref 8.72 m (the fitted
baseline wing), moments about the starter's CG station (x = 18.27 m). VSPAERO runs every condition (fixed wake),
Wave Drag the three supersonic Mach numbers (30 slices, 30 rotation sections), Parasite Drag the five flight
points. **Wing Span** and **Wing Inner Sweep** are notebook inputs, recorded
in every row. One run takes about 23 minutes on 16 threads, nearly all of it VSPAERO's subsonic cases.

Baseline results at beta 0 (derivatives per degree):

| Mach | CL_alpha | Cm_alpha | CD at alpha 0 = induced + wave + parasite | CY_beta | Cl_beta | Cn_beta | best L/D |
|---|---|---|---|---|---|---|---|
| 0.6 | 0.0294 | -0.0091 | 0.0117 = 0.0013 + 0 + 0.0104 | -0.0030 | -0.0020 | +0.0030 | 9.2 |
| 0.9 | 0.0306 | -0.0096 | 0.0112 = 0.0014 + 0 + 0.0098 | -0.0030 | -0.0020 | +0.0030 | 9.6 |
| 1.2 | 0.0331 | -0.0129 | 0.0279 = 0.0044 + 0.0142 + 0.0094 | -0.0039 | -0.0023 | +0.0039 | 4.4 |
| 1.5 | 0.0322 | -0.0135 | 0.0240 = 0.0040 + 0.0109 + 0.0091 | -0.0038 | -0.0022 | +0.0038 | 4.6 |
| 1.8 | 0.0308 | -0.0133 | 0.0216 = 0.0038 + 0.0089 + 0.0089 | -0.0039 | -0.0021 | +0.0038 | 4.8 |

**Center of Pressure (VSP)** gives the aerodynamic center (neutral point) of each alpha sweep, x_ac =
x_ref - (dCm/dCN) cref, in metres. From the VSPAERO run alone: 20.95 m at Mach 0.6 and 0.9, moving aft to 21.62,
21.88 and 21.98 m at Mach 1.2, 1.5 and 1.8. From the database's totals (the drag's normal-force component included):
20.94, 20.98, 21.58, 21.84 and 21.95 m. The starter's own AeroSandbox estimate is 20.90 m; with the CG at 18.27 m both give a static margin of about
30% of the mean aerodynamic chord subsonic. The block also gives each case's center of pressure, x_cp = x_ref -
CMy cref / CFz, which runs away as the normal force goes to zero.

Sideslip derivatives are from beta 0 and 4 degrees at alpha 2. The aircraft is statically stable in pitch, yaw and
roll throughout, more so in pitch above Mach 1 as the aerodynamic centre moves aft. The parasite drag (0.009 to
0.010) is OpenVSP's component buildup on the fitted components (wing t/c 3.0%, tails 4.0%, fuselage fineness 20,
as the starter's NACA 0003 and 0004 sections and fuselage give). The fuselage's lift is not in these numbers (see
above), and the wave drag carries the uncertainty of the nacelle's inlet and base steps.

**`examples/x59_doe.py`** runs the notebook headlessly once per design with `ntopcl --trustnotebook -j inputs.json`
(the JSON sets Wing Span and Wing Inner Sweep) and collects every run's rows in `runs/x59_doe.csv`: one table with
Wing_Span_m, Wing_Inner_Sweep_deg, Mach, Altitude, AoA, Beta, CL, CD, CY, Cl, Cm, Cn, L/D, CA, CN and each
contribution's share (CL_VSPAERO ... Cn_Parasite_Drag). `examples/x59_doe_results.csv` was made with the first
database layout (Re_cref and the drag parts as columns; the same totals). Designs: the
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
| Alpha, Beta, Mach; Re_cref, Vinf, Rho (from the standard atmosphere at each flight point); Sref, bref, cref; Xcg, Ycg, Zcg | Flight Condition |
| Vref, Machref (manual reference velocity) | not exposed: Vref = Vinf |
| NCPU, fixed wake, wake iterations, wake nodes, far field, ground effect, freezing, implicit wake, wake relaxation, convergence factors, core size, multipole, preconditioners, quad tree, `-save`/`-restart` | Solver Settings |
| Thick and thin geometry sets | each component's Set (Thick, Thin or None) |
| Find body wakes | Analysis: Body Wakes |
| Wing trailing-edge wakes | automatic: OpenVSP marks every wing trailing edge |
| X-Z symmetry; 2D FEM and Tecplot files; Sref/cref from the model; CG calculation; stall model; propellers, rotors, actuator disks, noise; stability types and unsteady time stepping; Cp slicer; control surfaces | out of scope (your list) |
| NRef (refinement levels), cull fraction, degenerate-region merging, intersection cleaning | OpenVSP's defaults |
| `-fem`, `-interrogate` with survey points | not implemented (say if you want them) |

**Parasite Drag**: geometry set, laminar and turbulent Cf equations and sub-component export are inputs; the model's
length unit is metres, fixed (OpenVSP Model writes every nTop length in metres, so it is not a choice); the default
turbulent Cf is Roughness Schlichting Avg Compressible, because OpenVSP 3.56's "Schlichting Compressible" applies no
Mach correction (it equals the incompressible 0.455/(log10 Re)^2.58: 1.7x too much friction at Mach 3); Mach, altitude
and Sref come from the Flight Condition (one run per flight point, Mach as the speed unit);
each component's form factor equations, manual form factor, % laminar, Q, roughness, Te/Tw, Taw/Tw and incorporated
generation are inputs of Parasite Drag Component (VSP). Not exposed: atmospheres other than US Standard 1976 (no
temperature offset; gamma 1.4), excrescences, Sref from a wing.

**Wave Drag**: number of slices, number of rotation sections, symmetry and geometry set are inputs; the Mach numbers
(the supersonic ones) and Sref come from the Flight Condition. Out of scope: inflow and outflow faces.

## Rebuilding and running headlessly

```bash
sh build/py.sh build/build_blocks.py
sh build/py.sh build/x59_example.py --run
sh build/py.sh build/mixed_example.py --run
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
  vortex edges): use them for viewing.
* nTop pairs lists item by item: a Wing block's Span Axes must have one vector per body, and Aero Coefficients' lists
  one value per flight condition ("Input list length mismatch" otherwise).
* Modified Newtonian and the Viscous Correction cost Surface Integrals per condition (see above); the on-mesh friction
  sums of the former `ntop_newtonian` cannot run over a flight condition and were not carried over.
