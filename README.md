# nTop_vsp

nTop custom blocks that turn **nTop implicit bodies into an OpenVSP model** and run aerodynamic analyses on one
shared **Flight Condition** (Mach numbers and altitudes through the US Standard Atmosphere 1976, angles of attack and
sideslips, and the reference area, span, chord and moment point): OpenVSP 3.56's **VSPAERO**, **Wave Drag** and
**Parasite Drag**, and, computed entirely in nTop on the implicit body itself, **Modified Newtonian** pressure and a
**Viscous Correction** (skin friction). **Aero Database (VSP)** sums any mix of them, each over its own Mach range,
into a **six-degree-of-freedom aero database** (CL, CD, CY, Cl, Cm, Cn against Mach, altitude, alpha and beta);
**Center of Pressure (VSP)** gives its centers of pressure and aerodynamic centers as positions, and **Polar Plot
(VSP)** plots it with matplotlib. Because the geometry stays parametric in nTop, its inputs can be independent
variables of the database too: `examples/x59_doe.py` runs a design of experiments over the wing span and sweep of an
X-59-style aircraft.

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
  `slices.py` (contour reader), `fit.py` (section fitting) and `vspmodel.py` (OpenVSP script and `.vsp3` writing);
  `run_vspaero.py`, `run_wave.py`, `run_parasite.py`, `run_database.py` (the buildup), `run_plot.py`;
  `run_newtonian_mesh.py` with `local_methods.py` (the on-mesh blocks); `pressure.py` (pressure point maps);
  `vspgeom.py` and `adb.py` read VSPAERO's panel file and binary results.
* `build/` — writes the blocks and the examples (`build_blocks.py` with `model_blocks.py`, `aero_blocks.py` (flight
  condition, Newtonian, viscous, coefficients), `vspaero_blocks.py`, `drag_blocks.py`, `buildup_blocks.py`;
  `newtonian.py`, `viscous.py` and `ops.py` are the Modified Newtonian and skin-friction block graphs, brought over
  from the former `ntop_newtonian` project; `x59_example.py`, `mixed_example.py`; `x59.py` derives notebooks from the
  starter; `ntopfile.py` writes `.ntop` files).
* `validation/wave_drag_mach.py` — the wave drag study below (results in `validation/results/`);
  `validation/newtonian/` — the Modified Newtonian and viscous-correction validation report from `ntop_newtonian`;
  `validation/rasaero/` — which combination of these analyses best reproduces RASAero II on a finned rocket (Mach
  1.5 to 5, alpha -15 to 15, beta -4 to 4 degrees): `validation/rasaero/REPORT.md`;
  `validation/noncircular/` — the same rocket with a trapezoid-section body, against inviscid SU2: where methods built
  for round bodies fail (`validation/noncircular/REPORT.md`).

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
their Flight Condition's files. nTop asks once to trust Run Command blocks.

## Blocks

| Block | Output | Purpose |
|---|---|---|
| OpenVSP Fuselage (VSP) | text | a list of implicit bodies as OpenVSP Fuselages: names, slicing, sections, clustering, tessellation, set |
| OpenVSP Wing (VSP) | text | a list of implicit lifting surfaces as OpenVSP Wings: names, span axes (one per body), symmetry (none, x-z, x-y, both), extension into the body, force accounting (whole wing or exposed strips), sections, airfoils, tessellation, clustering, set |
| OpenVSP Duct (VSP) | text | the side wall of an implicit body with planar ends in two given planes (not necessarily parallel), as an open OpenVSP Fuselage: outer or inner wall, sections, clustering, tessellation, set |
| OpenVSP Model (VSP) | model path | builds `<case>.vsp3` from a list of component lines |
| OpenVSP Model Check (VSP) | list of reals | per component: largest and mean distance from the OpenVSP surface to the nTop contours, sections, slices |
| OpenVSP Model Mesh (VSP) | mesh | OpenVSP's tessellation of the model, to view against the bodies |
| OpenVSP CAD Export (VSP) | file path | the model as STEP or IGES by OpenVSP's own CAD export: trimmed and intersected (one BREP solid), stitched (one BREP per component) or untrimmed surfaces; optional planar end caps |
| Lofted CAD Export (VSP) | file path | bodies (from OpenVSP Fuselage blocks) as STEP solids fitted straight to their nTop slices, without OpenVSP's skinning: smoother and closer than OpenVSP's export on non-round bodies; sharp chines kept |
| Flight Condition (VSP) | results path | Mach numbers paired with altitudes (US Standard Atmosphere 1976), angles of attack, sideslips; reference area, span, chord, moment point |
| VSPAERO Solver Settings (VSP) | text | threads, wake model and nodes, far field, ground effect, freezing, implicit wake, convergence factors, preconditioners, quad tree, restart file, convergence check |
| VSPAERO Analysis (VSP) | results path | model, flight condition, Mach range, body wakes, solver settings, pressure export |
| VSPAERO Coefficient (VSP) | list of reals | any column of VSPAERO's polar, one value per case |
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
| Tangent Wedge on Mesh (VSP) | results path | local-inclination pressure on the same mesh pipeline: exact oblique shock (tangent wedge) or exact cone (tangent cone) windward, Prandtl-Meyer or Newtonian shadow leeward; supersonic |
| Modified Newtonian Pressure (VSP) | scalar field | the Modified Newtonian static pressure (Pa) of one flight condition case, to view on the body |
| Surface Pressure (VSP) | Scalar Point Map | static pressure (Pa) of one flight condition case of VSPAERO, Modified Newtonian on Mesh or Tangent Wedge on Mesh run with Export Pressures |
| Surface Pressure Jump (VSP) | Scalar Point Map | the pressure jump across VSPAERO's thin surfaces for one case, plus a required Reference Pressure (Pa) |
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

**Fuselage fit** (`fit.fit_fuselage`). Of the **Sections**, half are placed with OpenVSP's clustering function between
the nose and tail (**Nose/Tail Clustering**); the rest go, one at a time, to the slice where interpolating the
outline (top, bottom and side extents) between the stations already chosen misses it most, so curved noses,
shoulders and boattails get the stations they need (a plain cylinder needs few). Each takes the largest closed
contour at its slice: the contour itself as an OpenVSP file cross-section, normalised by its width and height and
positioned at its bounding-box centre.

Ends:
- **Pointed** (section area tending to zero): an OpenVSP point section. Its position comes from fitting
  √area ∝ |x − x_tip|^p to the slices nearest the tip (p = 1 for a cone or ogive, 0.5 for a sphere), held within one
  slice spacing of the outermost slice. With p below 0.75 the tip is **rounded** and the surface meets the axis at 90°.
- **Blunt**: the outermost slice's contour, carried out to midway between that slice and the next (empty) one, or
  nTop's bounding box if nearer (it can be a few millimetres loose), so within half a slice spacing of the true end,
  with a flat cap.

Skinning angles are set explicitly at every section and side, from the body's own slopes in the nTop slices around it.
Left free, OpenVSP puts back its default fixed angles (45° at the first section, 0 between, −45° at the last), which
flare every blunt end out by about a tenth of its radius and bend pointed noses. Each side keeps its own angle
(OpenVSP's AllSym, "all sides like the top", is off). On a body whose sides have different slopes
(`validation/fuselage_sides`: elliptic sections, a drooped nose, a top-only boattail) that cuts the error against the
exact surface from 0.210 to 0.066 mm mean (0.80 to 0.37 mm largest) on the boattail, 0.025 to 0.009 mm along the
middle and 0.129 to 0.086 mm on the nose.

On test bodies of revolution (`validation/fuselage_ends`, 6 mm slices), the OpenVSP radius against the exact one,
before and after the end placement and skinning angles above:

| Body | Before | After |
|---|---|---|
| cylinder, flat at both ends | +2.4 mm at both ends | exact |
| ogive nose, flat base (the hobby rocket) | +2.4 mm (nose and base) | 0.09 mm |
| hemisphere nose, flat base | −5.6 mm (a conical tip) | 0.5 mm past the first 1.5 mm (tip 0.8 mm aft) |
| flat nose, conical boattail to a smaller base | +2.4 mm | 0.17 mm |

The model check leaves out end-cap points, which would read as up to a radius of error against the section outlines.
Building 30 sections through OpenVSP's API takes 51 s, because every change rebuilds the whole surface; the driver
instead writes each section's points, size and position straight into the `.vsp3` XML (2 s).

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

**Wings in a body.** **Extension** (a length, 0 for none) continues each wing that far inboard of its root, along its
span axis, with the first section's leading- and trailing-edge lines, thickness and root airfoil. It must stop short of
the plane through the origin normal to the span axis (for fins on a body along the x axis, the axis); the driver
refuses one that reaches it. Four fins touching on the axis are two crossing thin surfaces: VSPAERO then returned no
finite result (at Mach 1.5 on a round body, at every Mach number on a trapezoid body), and its spanwise loads lost a
strip of the merged pair. Stopping 1 to 2% of the root's distance short (27.6 mm for roots 27.7 mm off the axis) keeps
opposite wings' root vortices nearly cancelling: at Mach 3 the whole-fin lift was 1.026 with 0.05% left, 1.013 with
0.2% and 0.929 with 2%. VSPAERO then solves panels through the body, which removes the free root edge (and
its tip-like loss of lift) of a wing that ends at the body surface. **Force Accounting** says which of those loads
count: Whole Wing (every strip: the inboard part's lift stands for the body's carry-over, and over-predicts it on a
slender body), or Exposed Strip (only outboard of the wing body's own root: give the wing trimmed at the fuselage
surface). For Exposed Strip, OpenVSP Model records each wing's root (`<case>.ntvsp.wings.csv`) and VSPAERO Analysis
takes the inboard strips' forces and moments off every case from VSPAERO's spanwise loads (`.lod`; the whole-wing
polar is kept as `.ntvsp.polar_whole.csv`). On a cruciform rocket (`validation/rasaero/REPORT.md`) the exposed strips'
lift is 10% under to 5% over SU2's fin lift (Mach 2 to 5), where fins ending at the body give 59 to 76% of it; Whole
Wing puts 30 to 38% of the lift on the inboard strips, 2 to 6 times SU2's carry-over.

**Exposed Strip + Carry-over (TR 1307)** is Exposed Strip plus the lift each wing carries onto the body. It is
implemented in `vsp_ntop/carryover.py`, after W. C. Pitts, J. N. Nielsen and G. E. Kaattari, *Lift and Center of
Pressure of Wing-Body-Tail Combinations at Subsonic, Transonic, and Supersonic Speeds*, NACA Report 1307 (1957).

TR 1307 writes a wing-body combination's lift as factors on the lift $L_W$ of the wing alone (the exposed panels
joined at the root). $K_B(W)$ is the share the wing puts on the body:

$$L_{B(W)} = K_{B(W)}\,L_W$$

The driver takes $L_W$ to be each wing's own normal force $F_n$ from its exposed VSPAERO strips. It adds

$$\Delta\mathbf{F} = K_{B(W)}\,F_n\,\hat{\mathbf{n}}$$

along the wing's normal $\hat{\mathbf{n}}$, at the centre of the carried load. That force is added to every force and
moment column of the polar (totals and the inviscid part).

*Supersonic: the planar model* (TR 1307's Figure 4). The body is flattened into a strip of width $d$ beside the wing
root. The wing's pressure reaches the strip between the Mach lines from the root's leading and trailing edges. With
$\xi$ aft of the root leading edge and $\eta$ across the strip from the root, TR 1307's equations 23 and 25 give the
lifting pressure there, in units of the two-dimensional value $4\alpha/\beta$:

$$
\hat p(\xi,\eta) =
\begin{cases}
\dfrac{\beta m}{\pi\sqrt{\beta^2m^2-1}}\,\cos^{-1}\!\left(\dfrac{\xi/\beta+\beta m\,\eta}{\eta+m\,\xi}\right), & \beta m>1 \text{ (supersonic leading edge)}\\[3ex]
\dfrac{2(\beta m)^{3/2}}{\pi(\beta m+1)}\sqrt{\dfrac{\xi/\beta-\eta}{m\,\xi+\eta}}, & \beta m<1 \text{ (subsonic leading edge)}
\end{cases}
$$

The two are blended linearly over $|\beta m-1|<0.02$. The pressure is integrated numerically (48 × 48 Gauss points)
over

$$0\le\eta\le d,\qquad \beta\eta\le\xi\le\min(c_r+\beta\eta,\ \ell)$$

That gives

$$
K_{B(W)}=\frac{1}{S_W}\iint\hat p\,d\xi\,d\eta,\qquad
\bar\xi=\frac{\iint\hat p\,\xi}{\iint\hat p},\qquad
\bar\eta=\frac{\iint\hat p\,\eta}{\iint\hat p}
$$

The load acts at the root leading edge $+\,\bar\xi$ aft and $\bar\eta$ across the body. Integrating to $\ell$ covers
both of TR 1307's closed forms: a body that runs on behind the wing (equations 24, 26) and one that ends at its
trailing edge (equations 30, 31).

*Slender-body theory*, TR 1307's equation 21, is used instead subsonic, and where TR 1307's selection rule says the
wing's tip would reach the strip:

$$\beta A(1+\lambda)\left(\frac{1}{\beta m}+1\right)\le 4$$

It gives

$$K_{B(W)}=(1+\tau)^2-K_{W(B)}$$

$$K_{W(B)}=\frac{2}{\pi(1-\tau)^2}\left[(1+\tau^4)\left(\tfrac12\tan^{-1}\!\tfrac12\!\left(\tfrac1\tau-\tau\right)+\tfrac\pi4\right)-\tau^2\left(\tfrac1\tau-\tau+2\tan^{-1}\tau\right)\right]$$

That is equation 14. The load is still placed where the planar model puts it (subsonic: at the root mid-chord and the
body's middle).

| Symbol | Meaning | Where it comes from |
|---|---|---|
| $M$, $\beta=\sqrt{M^2-1}$ | freestream Mach number, Prandtl-Glauert factor | flight condition |
| $\alpha$ | the wing's incidence (it cancels: only $F_n$ is used) | |
| $r$ | root's distance from the plane through the origin normal to the span axis (the body axis) | first exposed station |
| $s$ | tip's distance from that plane (the semispan) | last station |
| $\tau=r/s$ | body-to-semispan ratio | |
| $c_r$, $c_t$, $\lambda=c_t/c_r$ | root and tip chords, taper ratio | first and last stations |
| $m=\cot\Lambda_{LE}$ | cotangent of the leading-edge sweep, $(s-r)/(x_{LE,tip}-x_{LE,root})$ | stations |
| $S_W$ | one wing's exposed planform area | stations (trapezoidal rule) |
| $A=4(s-r)/(c_r+c_t)$ | aspect ratio of the exposed pair joined at the root | |
| $d$ | body width across the wing's plane at the root: root distance plus the body's extent on the far side, from fuselage contour points within $0.2r$ of the plane over the root chord | fuselage sections |
| $\ell$ | body length from the root leading edge to the body's end | fuselage sections |
| $\xi$, $\eta$ | distance aft of the root leading edge, distance across the body from the root | |
| $F_n$ | the wing's normal force from its exposed VSPAERO strips (along $\hat{\mathbf n}$) | `.lod` |

These are measured per wing and mirror copy (`<case>.ntvsp.carryover.csv`). The body is the OpenVSP Fuselage with a
section nearest the root's mid-chord point, so it must be in the model (Set None is enough).

*Assumptions and limits:*
- Linear supersonic theory: thin wings at small incidence, no shocks, no viscous effects.
- The body is flat (TR 1307's planar model). Only its width $d$ and length $\ell$ enter: not its section shape, or how
  its crossflow changes the wing's upwash. Nothing assumes a round body, but a noncircular body's own effect on the
  interference is not modelled (`validation/noncircular`).
- The strip pressure is that of the root leading edge's conical flow. Tip effects are left out when the selection
  rule keeps the planar model, and the wing is taken as trapezoidal: only its root and tip chords and leading-edge
  sweep enter.
- The $4\alpha/\beta$ of linear theory stands for the wing's own mean lifting pressure, $F_n/S_W$ from VSPAERO. That
  keeps the wing and the carried load consistent, but it carries everything in $F_n$, including lift from camber or
  deflection, at the incidence factor. TR 1307 gives control deflection its own, smaller factor $k_{B(W)}$, which is
  not modelled.
- Each wing carries onto the body independently. Fin-fin interference is only what VSPAERO itself solves.
- Only $K_{B(W)}$ is applied. TR 1307's $K_{W(B)}$, the wing's own lift raised by the body's upwash, is not, because
  it assumes the wing-alone lift already includes tip and root losses.

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

## Ducts

**OpenVSP Duct (VSP)** takes one closed implicit body whose two ends are planar and lie in two given nTop Planes (an
inlet and an outlet plane, which need not be parallel), and makes its side wall an OpenVSP Fuselage **open at both
ends** (no caps, no inner/outer subtraction): only the surface. **Wall** picks the body's outside (a solid flow
passage) or, for a tube, its hole (the inside wall); for both walls of a tube use two Duct blocks.

**Cuts on turning planes.** nTop cuts the body with **Sections** planes that turn gradually from the inlet plane to
the outlet plane: plane s (0 to 1, OpenVSP's clustering function) passes through (1 - s) p_in + s p_out, the
centre of the body's bounding box projected onto each end plane, with the normal (1 - s) n_in + s n_out. Each cut
is one OpenVSP file cross-section in its own plane (XRotate 0; YRotate and ZRotate from the normal: OpenVSP turns a
section by Rx Ry Rz about its centre). The end cuts are made 0.2% of the way inside the body (a cut on an end face
is degenerate) and moved along the duct's axis onto their end plane, so the open ends lie in the given planes.
Slicing normal to one axis instead leaves pieces next to a tilted end that are very short on one side and long on
the other, where OpenVSP's skinning spline overshot (5 cm on a 1 m tube with a 25 degree end).

**Skinning.** Every side of every section is left free (no angles, slews or strengths) with OpenVSP's top/bottom and
left/right symmetries and its left/right equalities switched off. OpenVSP's defaults tie the sides together; on a
round tube with tilted ends they bowed the surface out by 0.1 m on a 1 m radius, while uncoupled it came back within
0.2 mm, and the open surface's area within 1 ppm of 2 pi R L. A waisted duct (15% contraction, tilted ends) with
12 sections: 0.8 mm worst on 1 m (`validation/cad_export/duct_tube.py`).

**nTop's in-plane axes.** For a cut plane that is not normal to a world axis, nTop's Slice Body writes the contour
in its own axes, v = n x y (normalised) and u = v x n, whatever frame it is given (found on the test duct: the
frame's axes were ignored). The block builds its frames in exactly those axes. A duct running along y has no such
axes: not supported.

## From OpenVSP to CAD

**OpenVSP CAD Export (VSP)** writes the OpenVSP model with OpenVSP 3.56's own CAD export (`vsp_ntop/run_cad.py`):

| Method | OpenVSP | Result |
|---|---|---|
| Trimmed (intersected, one solid) | Surface Intersection | components intersected, parts inside other components trimmed away; closed surfaces become one BREP solid (STEP, Representation BREP), open ones (duct walls) stay surfaces |
| Stitched (one solid per component) | File > Export, split and stitch | each component's patches stitched into its own watertight BREP (an open shell for a duct); not intersected |
| Untrimmed Surfaces | File > Export, STEP / IGES | each component's surface patches, no topology |

**Units.** OpenVSP's CAD length unit does not scale the model, it declares the unit of its numbers (a pod 10 long
written as LEN_M reads as 10 m, with the default LEN_FT as 10 ft). OpenVSP Model builds in metres, so the files are
always declared in metres and the CAD program converts on import. The plain export reads the unit from OpenVSP's
STEP/IGES settings, the stitched and trimmed ones from its Surface Intersection settings; the driver sets both.

**Test** (`validation/cad_export/test_cad.py`, results in `validation/cad_export/results.txt`): the hobby rocket's body as a
Fuselage (24 sections, 6 mm slices), its +y fin as a Wing, and a waisted tube (inner radius 60 mm narrowing to 51 mm,
6 mm wall, 300 mm long) cut by an inlet plane turned 15 degrees about y and an outlet plane turned 25 degrees about
z, as two Ducts (outer and inner wall). Each CAD file is read back with OpenCASCADE (gmsh); CAD -> implicit is nTop's
own field at points on the CAD faces, implicit -> CAD the distance from nTop mesh vertices to the nearest CAD face.

| | Trimmed STEP | Stitched STEP | Untrimmed STEP | Trimmed IGES |
|---|---|---|---|---|
| solids / faces | 1 / 22 | 2 / 30 | 0 / 30 | 0 / 22 |
| body, CAD -> implicit p99 / max (mm) | 2.94 / 2.94 | 2.94 / 2.94 | 2.94 / 2.94 | 2.94 / 2.94 |
| fin, CAD -> implicit p99 / max (mm) | 1.33 / 2.94 | 1.33 / 2.94 | 1.31 / 2.93 | 1.32 / 2.94 |
| fin, implicit -> CAD p99 / max (mm) | 0.43 / 0.49 | 0.43 / 0.49 | 0.43 / 0.49 | 0.43 / 0.49 |
| duct, CAD -> implicit mean / max (mm) | 0.004 / 0.13 | 0.004 / 0.13 | 0.004 / 0.13 | 0.004 / 0.13 |
| duct, implicit -> CAD mean / max (mm) | inner wall dropped | 0.004 / 0.14 | 0.004 / 0.14 | inner wall dropped |

- The CAD is OpenVSP's surface: the stitched body solid is 2442.5 cm3 against OpenVSP's own 2438.1 cm3 (CompGeom).
  What separates it from the implicit is the OpenVSP fit: the body is 0.7% under the exact 2460.6 cm3 (the model
  check's mean of 0.13 mm over its surface) and its flat base sits 2.9 mm short of the body's (a blunt end is placed
  midway to the next, empty, slice: half the 6 mm Slice Spacing; finer slices bring it closer).
- Trimmed drops the duct's inner wall: Surface Intersection treats it as inside the outer wall's tube and trims it
  away. For both walls of a tube, use Stitched or Untrimmed, or export them separately.
- The four exports took 4 s each; the duct block's 24 cuts about 70 s in nTop.
- **Degenerate geometry** (nTop's import warning). Where OpenVSP closes a surface to a point, its patches have a
  control-net row collapsed to that point (the surface's derivative vanishes along it): the tip of a pointed nose, the
  centre of a flat end cap, and the leading and trailing edges of a wing's root and tip caps. Parasolid (nTop's CAD
  kernel) reports these as degenerate geometry: 2 on the rocket body (nose tip, base centre), 4 on the fin. The open
  duct walls have none. They are in every export method. nTop still reads the files: Import Part, then CAD Body to
  Implicit, gave the stitched body back within 0.05 mm (mean) of the original implicit and the fin within 0.5 mm
  (`validation/cad_export/ntop_round_trip.py`; singular points located with
  `validation/cad_export/singular_points.py`).
- **Planar End Caps** (CAD Export option, STEP; `vsp_ntop/step_caps.py`) edits the STEP file: each flat cap's patches
  become one ADVANCED_FACE on a PLANE bounded by the cap's rim edges (shared with the side surfaces, kept as they are),
  with holes where another part passes through the cap (trimmed export). OpenVSP's cap patches are not oriented
  consistently with each other (two of them run a shared edge the same way), so each rim loop takes its direction
  from a side face that shares one of its edges, and the plane's normal from that loop. Everything else in the file
  is kept (shells, names, units); OpenVSP's own file stays as `<Case Name>.openvsp.stp`. On the test, all three caps
  (body base, fin root and tip) were rebuilt in both the stitched and trimmed files; only the nose tip's collapsed
  patches are left, and the deviations above are unchanged (base face area 2290.0 mm2 as before, = pi 27^2).
  Rebuilding the file with gmsh/OpenCASCADE instead worked too, but split each open duct wall into one sheet body
  per face and dropped OpenVSP's surface names.
- **Trimmed export in nTop.** OpenCASCADE reads the trimmed STEP solid within 0.03 mm of the bodies, but nTop's CAD
  Body to Implicit of it was 2.8 mm mean and 20 mm at the 99th percentile off the original body (with planar caps
  nTop split it into more bodies). The stitched export round-trips cleanly (0.05 mm), so use Stitched for nTop.

## Lofted CAD: STEP straight from the slices

OpenVSP's fuselage skinning can't follow a flattened section. OpenVSP 3.56 sets a section's surface slope only at
four points (top, bottom and the two sides; the angles' "Set" flags can't be turned off) and interpolates it in
between. On a flattened nose the real slope drops from the side value to the top or bottom value almost at once. The
diagonals come out about 3 degrees too steep at every section, which shades as a ripple. More sections, more
contour points, exact or least-squares angles, lower tangent strength, and OpenVSP's own Fit Model (Levenberg-
Marquardt on the angles: 9 minutes for the nose alone, rms slope error 1.38 to 1.32 degrees) don't remove it.

**Lofted CAD Export (VSP)** skips OpenVSP for the CAD (`vsp_ntop/loft.py`, `vsp_ntop/step_write.py`, driver
`run_loft.py`; numpy only). It takes the same OpenVSP Fuselage blocks and fits each body's nTop slices directly:

- **Surface.** One bicubic B-spline surface through every slice (about 400 at Slice Spacing 0.0025), least squares
  with a light curvature penalty. Around the body it is periodic, Points Around (48) control points. Each contour is
  split where it crosses the horizontal and vertical lines through its centre (right, bottom, left, top), and each
  quarter is sampled evenly by arc length (splitting at the widest point made the seams wander on flat sides, visibly
  wavy in nTop). Along the body its parameter is the distance travelled by the sections' centres and equivalent
  radii, so a rounded nose that closes up at almost constant x is still followed. There is one control row per two
  slices (Slices per Control Point).
- **Ends.** A pointed end is a pole (the first row of control points is the tip). A blunt end carries the last slice
  out to the end plane (placed as OpenVSP Model places it) and is closed by a flat face.
- **Sharp edges.** Where the contours turn a corner at the sides, top or bottom (more than 30 degrees, median over
  the middle 80% of slices), that split point's knot is made triple, so the surface keeps the edge (a chine). A side
  with a chine is split at its widest point, which is the chine.
- **CAD.** The periodic surface is cut exactly (knot insertion) into clamped patches: four at the split points, or
  Faces Around 8, 12, 16 ... (each quarter cut evenly; a multiple of 4 dividing Points Around). Each is a B-spline
  face; neighbours share their edge curves exactly and join as smoothly as the uncut surface. One closed solid per
  body (the side faces, plus a flat cap per blunt end), all bodies in one STEP file, AP214, millimetres.
  `docs/STEP_FORMAT.md` walks through the format with records from these files.

Tests (`validation/loft_cad/test_loft.py`, results in `validation/loft_cad/results.txt`). Each body is one nTop
expression, also evaluated in numpy as the exact field. The lofted surface and OpenVSP's (24 sections, default
skinning) are scored on distance |f|/|grad f| and on the angle between the surface normal and grad f, over the body
between 2% and 98% of its length:

| Body | What it tests | Loft: mean / max distance, rms normal error | OpenVSP |
|---|---|---|---|
| sides | flattened elliptic nose, drooped, top-only boattail, flat base | 0.001 / 0.05 mm, 0.30 deg | 0.038 / 0.36 mm, 1.25 deg |
| canopy | round fuselage with a canopy smoothly blended on (concave sections) | 0.001 / 0.03 mm, 0.17 deg | 0.023 / 0.33 mm, 0.58 deg |
| blob | rounded-square sections on a curved centreline, poles at both ends | 0.000 / 0.003 mm, 0.20 deg | 0.030 / 0.31 mm, 0.71 deg |
| chine | lens sections, sharp chines the whole length | 0.000 / 0.001 mm, 0.01 deg | 0.016 / 0.24 mm, 3.8 deg (0.85 clear of the chines) |

- Without the sharp-edge knots, the chine body was rounded by up to 0.20 mm along its edges.
- On the sides body's nose, where OpenVSP had a +-0.2 mm wave between sections on the diagonals, the loft is a smooth
  offset of at most 0.013 mm.
- The worst remaining normal errors are on the sides body's own creases (x = 150 and 450 mm, where its exact surface
  has a slope break) and at its flattened rounded tip (4.7 degrees rms in the first 2%, against OpenVSP's 11).

Every file reads in OpenCASCADE as one closed solid, and the lofted surface's points lie on the CAD faces to 1e-13 mm.
nTop imports them with no warnings (pointed tips included). CAD Body to Implicit of each sits 0.013 to 0.020 mm mean
and at most 0.073 mm from the exact surface (OpenVSP's stitched file of the sides body: 0.041 mean, 0.354 max), also
with 16 faces round each body. End to end in nTop (two bodies in one Fuselage block into the Lofted CAD block): 27 s.

Limits: one closed contour per slice (the largest is used if a slice has several), bodies sliced along x (no wings or
ducts yet). Only the four split points can hold a sharp edge; an edge elsewhere on the section is rounded over about
one control-point spacing, as is a slope break along the body. Files are about 1.2 MB per body at the defaults.
OpenCASCADE's own volume of these solids reads 0.26% low on the sides body (1511.1 cm3), while integrating the
surface directly gives 1515.2 against the exact 1515.1.

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
  warning). `<case>.ntvsp.convergence.csv` lists every case, and the polar has a `Converged` column. On 53 earlier
  runs, every one whose numbers were sound passes (the X-59 baseline: 60 of 60 cases); it fails the NaN cases above
  and relaxed wakes still moving 15 to 290% per iteration (sound ones moved under 6%). VSPAERO's wake span efficiency
  `Ew` is infinite in every supersonic case (it computes no wake forces there) and is not checked.

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
the total angle of attack, so sideslip is included; no axial force. It is meant to go with Exposed Strip wings.

On the cruciform rocket, with Exposed Strip + Carry-over (TR 1307) fins and Body Lift (slender body and crossflow), the
results against inviscid SU2 at Mach 2 to 5 are:
- CN within −8% to +4%;
- the centre of pressure within 4 mm at three of five conditions, and 25 and 38 mm aft at Mach 4 and 5, alpha 4.

OpenVSP's tessellation-based Planar Slice was not usable for S(x): it was 80% off near a nose. OpenVSP's end caps
are several rings in one plane, sweeping across the cap. A cap at the front counts as a point; one further aft (a
flat base) is left out, so the body does not seem to shrink before its base.

**Potential Lift: Slender Body (Section Shape)** replaces the area with each section's crossflow apparent mass, which
is what slender-body theory actually uses:

$$\frac{d\mathbf{C}}{dx}=\frac{1}{S_{ref}}\,\frac{d(\mathbf{M}\,\hat{\mathbf e})}{dx}\,\sin 2\alpha_t\cos\frac{\alpha_t}{2}$$

Here $\hat{\mathbf e}$ is the crossflow direction in the y-z plane. $\mathbf M(x)$ is the section's 2×2 apparent-mass
tensor per unit fluid density, $M_{jk}=-\oint\phi_j\,n_k\,ds$, where $\phi_j$ is the exterior potential flow of the
section moving with unit velocity along $j$. The driver solves that flow on each OpenVSP ring with constant-strength
source panels, Richardson-extrapolated from the ring's edges and each edge split in two (the panels converge to first
order).

Checks:
- a circle gives its area to 0.02%, so round bodies are unchanged;
- an ellipse gives πa² and πb² across its axes to 0.2%;
- a flat plate gives πw²/4.

A section wide across the crossflow then lifts more. One not symmetric about the crossflow also gets a side component.
Several fuselages add their tensors, with no interference between them. The crossflow term still uses the radius of
the circle of equal area.

On the trapezoid-section rocket of `validation/noncircular`, the section's vertical-crossflow apparent mass is 1.37
times its area. That moves the forebody's lift and centre of pressure towards SU2's, but SU2's nose lifts 1.7 times
the round one's, so part of the gap is beyond slender-body theory (the report has the numbers).

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
method of the `ntop_newtonian` project, built here (`build/newtonian.py`, `build/viscous.py`) and
evaluated once per flight condition through nTop's list processing.

* **Modified Newtonian (VSP)**: Cp = Cp_max max(-n.V, 0)^2 with n the unit gradient of the body's field, so any
  implicit works; each force and moment component is one Surface Integral. Cp_max is the stagnation pressure behind
  a normal shock above Mach 1 (Rayleigh pitot) and the isentropic one below it (both 1.839 at Mach 1), 2
  (Classical), or **User Defined**: the **Cp Max Field** input, a scalar field, so Cp_max may differ over the body
  (one value on the fins, another on the nose); a constant field is one number everywhere. No self-shadowing: correct
  for convex bodies, a part in another's shadow is over-loaded.
* **Modified Newtonian on Mesh (VSP)**: the same pressure on the implicit **Body** (its unit gradient at each face
  centroid is the normal; the same Cp_max choices), integrated over a **Mesh** of it (Mesh from Implicit Body, or a
  remeshed or imported closed surface mesh) instead of Surface Integrals. The driver `vsp_ntop.run_newtonian_mesh`
  lists the faces' centroids, nTop evaluates the Body's gradient there (and the Cp Max Field, only when User
  Defined), and the driver sums -Cp S_f and -Cp (r_f - c) x S_f over the faces (S_f the outward vector area) at every
  condition. A mesh whose faces point inward overall is turned round. The pressure is the implicit's, the integral
  the mesh's. The log reports how well the Body's normals agree with the Mesh's faces; under half the area within 25
  degrees stops the driver, since a Mesh that is not where the Body is (other length units, a moved or scaled body)
  would otherwise give wrong forces silently. An OpenVSP Model Mesh passes (on the hobby rocket about 16% of its area
  is off: the fins' extensions, kept as separate shells inside the body) and gives CN within 1 to 4% and CA about 10%
  of nTop's own mesh. It writes the same coefficient table.
* **Tangent Wedge on Mesh (VSP)**: the same pipeline (the Body's normals, the Mesh's integrals) with a
  local-inclination pressure law (`vsp_ntop/local_methods.py`, after hypanel's `methods.py`). Each face's impact
  angle is $\delta$, with $\sin\delta=-\hat n\cdot\hat V$ (positive windward). Supersonic only (Mach Min defaults above 1); $\gamma=1.4$.
  - *Windward*, **Impact Method**: **Tangent Wedge**, the exact two-dimensional oblique shock that turns the flow
    through $\delta$ (the weak solution, from Thompson's closed form of the cubic in $\tan\beta$):
    $$C_p=\frac{2}{\gamma M^2}\left(\frac{p_2}{p_1}-1\right),\qquad
    \frac{p_2}{p_1}=1+\frac{2\gamma}{\gamma+1}\left(M^2\sin^2\beta-1\right)$$
    or **Tangent Cone**, the surface pressure of the exact cone of half angle $\delta$ (Taylor-Maccoll, weak shock;
    a port of hypanel's `conical`, tabulated once per Mach number and interpolated in $\sqrt{C_p}$), the usual choice
    for bodies of revolution and noses. Above the detachment angle $\delta_{max}(M)$ (the wedge's or the cone's),
    $C_p$ rises with $\sin^2\delta$ from its detachment value to $C_{p,max}$ at 90° (hypanel's bridge): the one use of
    **Stagnation Pressure** and the **Cp Max Field** here.
  - *Leeward*, **Expansion**: **Prandtl-Meyer**, the isentropic expansion from the freestream with
    $\nu(M_2)=\nu(M)+|\delta|$, limited at vacuum; or **Newtonian Shadow**, Cp 0.
  - *Separated faces:* with Prandtl-Meyer, faces turned more than 30° away blend (smoothstep, complete at 60°) to
    **Separated Faces**: Freestream (Cp 0, the default, so a flat base adds no base drag) or Gabeaud's base pressure
    $-1/M^2$.

  Checks:
  - against oblique-shock tables, $p_2/p_1$ to four digits at Mach 2 and 3, and the detachment angles to 0.01°;
  - $\nu(2)=26.380°$;
  - in the linear limit, windward plus leeward faces of a thin plate give its $4\alpha/\beta$ to 0.5%;
  - tangent cone: the cone table against an independent Taylor-Maccoll integration in $\theta$ (`tests/`) to 1e-6 in
    Cp, and a 15° cone's axial force at Mach 3 and 8 to 1%; detachment at Mach 2 is 40.7°.

  Faces are independent: no shock-expansion marching, no interference between parts.

  **Use it for thin surfaces, not bodies.** On a slender body each face of a constant section sees $\pm\alpha\cos\phi$,
  and a local two-dimensional law turns that into linear normal force all along the tube. Real crossflow round a
  constant section carries almost none. On the hobby rocket's body it gave CN 1.3 to 2.7 where SU2 has 0.09 to 1.0.

  Results on the fins of `validation/rasaero`, against SU2's fin loads:
  - within 3% at Mach 3 to 5 and alpha 4, 10% high at Mach 5 and alpha 8, 18% high at Mach 2 and alpha 2;
  - VSPAERO's exposed strips are 4 to 10% low at the same conditions.

  With Body Lift (slender + crossflow) and no carry-over (tangent-wedge fins already carry about SU2's fin load):
  - CN within −5% to +1% of SU2;
  - centre of pressure −2 to +44 mm.

  The fins' pressure drag is physical: about 0.22 in CA at Mach 2 from the 11° bevels, the fins' share of what Wave
  Drag gives the whole rocket. At Mach 1.5 the fins' leading edges are subsonic ($\beta m<1$), and the local law is
  poor (CN 50 to 60% off against RASAero II). 330 conditions on 540,000 faces take about 40 s of driver time, plus
  nTop's mesh.
* **Viscous Correction (VSP)**: skin friction at the inviscid edge state, hypanel's `viscous.correct`: edge pressure
  from the Modified Newtonian Cp of each condition (or the freestream, Cp = 0, for flat-plate friction), edge total
  pressure behind a normal shock above Mach 1 and isentropic below it, total temperature conserved; Blasius with
  Eckert's reference temperature (laminar) and White's van Driest II (turbulent), sudden transition at a running-length
  Reynolds number, adiabatic or fixed wall; shear along the surface flow, faded out between 30 and 60 degrees of
  leeward inclination. The freestream pressure and temperature of each condition come from the flight condition.
* **Modified Newtonian Pressure (VSP)** shows one flight condition case's static pressure on the body, p_inf + q Cp in
  pascals, with p_inf and q at the case's altitude (Case Number as in Surface Pressure (VSP)).

Axes: x along **Aft Axis** (the model's +x), z toward **Up Axis**, y = z x x: CA, CY, CN along them, moments about
them (the axes of VSPAERO too). Each writes a coefficient table, `<Flight Condition>_<Case Name>.ntvsp.coefficients.csv`,
with Export Text, for **Coefficients Contribution (VSP)**. Accuracy is nTop's Surface Integral (about 1e-4 on smooth
bodies, 0.1 to 1% across sharp edges, **Integration Tolerance**); validation against exact results and wind-tunnel
data, beside hypanel's implementation, is in `validation/newtonian/REPORT.md`: sphere CD within 8e-4, 15 degree cone
CA and CN within 0.002 and 0.010 from 0 to 60 degrees, the Apollo capsule's pressures at hypanel's RMS, skin friction
within 0.5% of hypanel's on a plate and a cone.

**Run time.** Each condition costs six Surface Integrals for each of Newtonian and friction, about 1 to 3 s each on
simple bodies, more on a full aircraft and much more on thin parts (fins); use Mach Min and Max, and few angles, where
Newtonian applies. Modified Newtonian on Mesh costs one mesh, one evaluation of the gradient at the face centroids
and a sum per condition. On the hobby rocket of `validation/rasaero` (1.2 m, four 4 mm fins), four Newtonian
runs (the whole rocket at 611,000 faces, the body, the fins, and the rocket again) over 330 conditions took 99 s in
all; Surface Integrals on the body alone took 750 s for 4 conditions. On the body at Mach 3 (alpha 0 and 8, beta 0
and 3 degrees) a 4 mm Mesh from Implicit Body (118,000 faces) agrees with the Surface Integral block within 0.5% in
CN, CY and the moments and to four digits in CA, and gives exact zeros where the Surface Integral's noise gives 7e-4.

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
20.94, 20.98, 21.58, 21.84 and 21.95 m. The starter's own AeroSandbox estimate is 20.90 m; with the CG at 18.27 m
both give a static margin of about 30% of the mean aerodynamic chord subsonic. The block also gives each case's
center of pressure, x_cp = x_ref - CMy cref / CFz, which runs away as the normal force goes to zero.

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
(`dev/inventory.txt`) and is a block input, out of scope, or not applicable (with the reason).

**VSPAERO**

| OpenVSP option | Here |
|---|---|
| Alpha, Beta, Mach; Re_cref, Vinf, Rho (from the standard atmosphere at each flight point); Sref, bref, cref; Xcg, Ycg, Zcg | Flight Condition |
| Vref, Machref (manual reference velocity) | not exposed: Vref = Vinf |
| NCPU, fixed wake, wake iterations, wake nodes, far field, ground effect, freezing, implicit wake, wake relaxation, convergence factors, core size, multipole, preconditioners, quad tree, `-save`/`-restart` | Solver Settings |
| Thick and thin geometry sets | each component's Set (Thick, Thin or None) |
| Find body wakes | Analysis: Body Wakes |
| Wing trailing-edge wakes | automatic: OpenVSP marks every wing trailing edge |
| X-Z symmetry; 2D FEM and Tecplot files; Sref/cref from the model; CG calculation; stall model; propellers, rotors, actuator disks, noise; stability types and unsteady time stepping; Cp slicer; control surfaces | out of scope |
| NRef (refinement levels), cull fraction, degenerate-region merging, intersection cleaning | OpenVSP's defaults |
| `-fem`, `-interrogate` with survey points | not implemented |

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

Inside every custom block, each top-level block is a variable: `ntopfile.as_variables` wraps it in a named
`core.var<type>` when the block is saved. A variable checks its block's output type, so the type must be nTop's own:
Export Text gives a `text_file`, and a block run once per list item a `list<...>` (`ntopfile.list_outputs`). A wrong
type shows up only at run time, as "Unknown block" or "requires X, but received X List".

## Surface pressures

Analyses that compute a surface pressure can write it for nTop, behind an **Export Pressures** option (Off by
default):

| Analysis | Points | Written |
| --- | --- | --- |
| VSPAERO Analysis (VSP) | panel centroids, nodes, or both (the option's choice) | Thick set: static pressure; Thin set: the jump across the sheet |
| Modified Newtonian on Mesh (VSP), Tangent Wedge on Mesh (VSP) | the Mesh's face centroids | static pressure |

Files are `<Case Name>.ntvsp.pressure.<Case>.csv` (static pressure p = p_inf + q Cp) and
`<Case Name>.ntvsp.pressure_jump.<Case>.csv` (Delta p = q Delta Cp), `x, y, z, value` in metres and pascals, one per
flight condition case (the Case column of `<Flight Condition>.ntvsp.flight.csv`: every flight point x sideslip x angle
of attack), with p_inf and q at that case's altitude. Each run first removes the pressure files of an earlier run of
the same Case Name. Read them with **Surface Pressure (VSP)** (Results, Case Number) or, for a jump, **Surface
Pressure Jump (VSP)**, which requires a **Reference Pressure**: a vortex-lattice sheet carries only the jump across
it, so the map shows Reference Pressure + Delta p. Both give a Scalar Point Map (Field from Point Map turns it into a
field). Modified Newtonian (VSP) has no points: its pressure is the field of **Modified Newtonian Pressure (VSP)**.
Body Lift computes sectional loads, not a surface pressure, so it exports none.

On a body-plus-fins VSPAERO run (13,138 thick and 2,584 thin panels) every exported value is within 2e-8 q of
p_inf + q Cp from VSPAERO's `.adb`.

## Tests

```bash
dev/venv/Scripts/python.exe -m pytest
```

`tests/` builds one small notebook (a sphere, a thin flat plate and a 15 degree cone; Mach 3 and 8, alpha 0, 10 and
30 degrees), runs it once with ntopcl and checks against closed-form answers:

- Modified Newtonian, Modified Newtonian on Mesh and Tangent Wedge on Mesh, and an Aero Database of each: the
  Newtonian sphere's drag (every Cp_max choice), the plate's Newtonian and tangent-wedge forces (exact oblique shock,
  Prandtl-Meyer), the cone's tangent-cone axial force (an independent Taylor-Maccoll integration);
- the pressure export at every face centroid, that only runs with Export Pressures on write pressure files, both
  pressure readers and Modified Newtonian Pressure.

About 80 s. Needs a signed-in nTop and a Python with this package, numpy, pytest and `implicit_io`; the blocks call
the same Python unless `NTOP_PYTHON_PATH` is set. ntopcl runs at below-normal priority on half the cores
(`NTOPCL_MAX_THREADS`) and is stopped above 8 GB (`NTOPCL_MEMORY_CAP_GB`).

## Known limits

* Runs with the same Working Directory and Case Name overwrite each other (the DOE runs one design at a time).
* A component block's body must have one closed contour per cut (the largest is used when there are several, with
  a note in the log); a fuselage is sliced along +x and a wing's span axis must be normal to x.
* The fuselage fit follows the body between sections with OpenVSP's own skinning, so abrupt steps (an inlet face)
  are smoothed: 0.43 m at the X-59's nacelle inlet.
* VSPAERO's exported pressures are read from its binary `.adb`, checked against its Tecplot output. Integrated over
  the surface they do not give VSPAERO's forces (it computes those from its vortex edges): use them for viewing.
* nTop pairs lists item by item: a Wing block's Span Axes must have one vector per body, and Aero Coefficients' lists
  one value per flight condition ("Input list length mismatch" otherwise).
* Modified Newtonian and the Viscous Correction cost Surface Integrals per condition (see above). The on-mesh friction
  sums of `ntop_newtonian` were not carried over: a list of fields evaluated at a list of points is a list of lists,
  which nTop does not allow, so they cannot run over a flight condition.
