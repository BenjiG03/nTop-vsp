# STEP files, as nTop_vsp writes them

A guide to the STEP format for reading the files from `vsp_ntop/step_write.py` (the lofted CAD) and OpenVSP's
exports. Examples are real records from `runs/loft_cad/sides_loft.stp`.

## What STEP is

STEP is ISO 10303, the standard for exchanging product data between CAD systems. A `.stp` / `.step` file is plain
text in the "Part 21" syntax (ISO 10303-21). What the records mean is set by an **application protocol** (AP):
AP203 (`CONFIG_CONTROL_DESIGN`, used by OpenVSP), AP214 (`AUTOMOTIVE_DESIGN`, used here) or AP242. For a solid made
of B-spline faces they hold the same geometry and topology; only the wrapper around it differs a little.

IGES is the older format (fixed 80-column records, mostly loose surfaces). STEP carries the **topology** too:
how faces join into a closed solid. That's why a STEP solid imports as one watertight body.

## Anatomy of a file

```
ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('nTop_vsp lofted body: B-spline faces fitted to nTop slices'),'2;1');
FILE_NAME('sides_loft.stp','2026-11-21T00:59:12',(''),(''),'','nTop_vsp step_write','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.04017800193690645,0.0,-11.824249569057999));
#2=VERTEX_POINT('',#1);
...
ENDSEC;
END-ISO-10303-21;
```

* **HEADER:** a description, the file name, a timestamp, the writing program, and the schema (which AP).
* **DATA:** one record per line here (a record may span lines; it ends with `;`). Each record is
  `#n=ENTITY_TYPE(arguments);`. The number `#n` is only an identifier: records point to each other with `#n`, and
  the order in the file doesn't matter.

Argument syntax:

| Written | Meaning |
|---|---|
| `'text'` | string (a quote inside is doubled: `'it''s'`) |
| `1.5`, `1.E-04`, `0.0` | real; it must contain a decimal point (`1.` not `1`) |
| `3` | integer |
| `.T.` `.F.` | true / false |
| `.MILLI.`, `.UNSPECIFIED.` | enumeration value |
| `#12` | reference to record 12 |
| `( ... )` | list; lists nest (a surface's control net is a list of rows) |
| `$` | not given (optional argument left out) |
| `*` | derived: the value comes from elsewhere (used in `ORIENTED_EDGE`) |
| `(A()B()C())` | a **complex instance**: one record that is several entity types at once, listed alphabetically. Units are written this way |

## The layers of a part

Reading from the top down, a part file is a chain of references:

```
SHAPE_DEFINITION_REPRESENTATION ── links the product to its shape
 ├─ PRODUCT_DEFINITION_SHAPE ── PRODUCT_DEFINITION ── PRODUCT_DEFINITION_FORMATION ── PRODUCT ('Sides')
 └─ ADVANCED_BREP_SHAPE_REPRESENTATION ('Sides')
     ├─ context: units (mm, radian, steradian) and the length tolerance
     ├─ AXIS2_PLACEMENT_3D: the part's origin
     └─ MANIFOLD_SOLID_BREP ('Sides')                    ┐
         └─ CLOSED_SHELL                                 │
             └─ ADVANCED_FACE ×5                         │ topology:
                 ├─ FACE_OUTER_BOUND ── EDGE_LOOP        │ what is connected
                 │    └─ ORIENTED_EDGE ── EDGE_CURVE     │ to what
                 │          ├─ VERTEX_POINT ×2           ┘
                 │          └─ B_SPLINE_CURVE_WITH_KNOTS ┐ geometry:
                 └─ B_SPLINE_SURFACE_WITH_KNOTS or PLANE │ where it is
                       └─ CARTESIAN_POINT ...            ┘
```

### 1. Product wrapper and units

From the end of `sides_loft.stp`:

```
#12936=PRODUCT('Sides','Sides','',(#12935));
#12940=PRODUCT_DEFINITION('design','',#12938,#12939);
#12941=PRODUCT_DEFINITION_SHAPE('','',#12940);
#12942=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#12945=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),#12942,'distance_accuracy_value','closure tolerance');
#12946=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#12945))
        GLOBAL_UNIT_ASSIGNED_CONTEXT((#12942,#12943,#12944))REPRESENTATION_CONTEXT('Context3D','3D'));
#12947=ADVANCED_BREP_SHAPE_REPRESENTATION('Sides',(#12928,#12932),#12946);
#12948=SHAPE_DEFINITION_REPRESENTATION(#12941,#12947);
```

* The **unit** belongs to the representation context, and every coordinate in the file is in it. Here it is
  `SI_UNIT(.MILLI.,.METRE.)`, millimetres. OpenVSP writes `SI_UNIT(.UNSET.,.METRE.)` (metres) and leaves its
  numbers in model units. That's why its CAD length unit only *labels* the numbers rather than scaling them.
* The **uncertainty** (1e-4 mm here) tells the reader how far apart two points can be and still count as the same.

### 2. Topology: what is connected to what

```
#3888=ORIENTED_EDGE('',*,*,#212,.T.);
#3889=ORIENTED_EDGE('',*,*,#835,.T.);
#3890=ORIENTED_EDGE('',*,*,#414,.F.);
#3891=EDGE_LOOP('',(#3888,#3889,#3890));
#3892=FACE_OUTER_BOUND('',#3891,.T.);
#3893=ADVANCED_FACE('Sides side 1',(#3892),#3887,.T.);
#12927=CLOSED_SHELL('',(#3893,#6900,#9907,#12914,#12926));
#12928=MANIFOLD_SOLID_BREP('Sides',#12927);
```

* `MANIFOLD_SOLID_BREP`: a solid, given by its boundary ("B-rep" = boundary representation).
* `CLOSED_SHELL`: the faces that enclose it. In a closed shell, every edge is used by exactly two faces, once in
  each direction.
* `ADVANCED_FACE(name, bounds, surface, same_sense)`: a piece of `surface`, cut out by its bounds. With
  `same_sense = .T.`, the face's outward normal is the surface's own normal.
* `FACE_OUTER_BOUND` → `EDGE_LOOP`: the boundary, as a list of `ORIENTED_EDGE`s in order. The loop runs
  **anticlockwise seen from outside** (from the tip of the outward normal). Holes would be `FACE_BOUND`s running the
  other way.
* `ORIENTED_EDGE('',*,*,edge,.T./.F.)`: an edge used forwards (`.T.`) or backwards (`.F.`). The two `*`s are its
  start and end, derived from the edge and the flag.
* `EDGE_CURVE(name, start vertex, end vertex, curve, same_sense)`: a stretch of a curve between two vertices.
  Faces that meet along an edge refer to the same `EDGE_CURVE`, which is what makes the solid watertight.
* `VERTEX_POINT(name, point)`: a corner.

The face above, side 1 (from the right-most line round to the bottom centreline), has three edges:
* `#212` forwards: the seam along the right side, from the nose tip to the base;
* `#835` forwards: its stretch of the base outline;
* `#414` backwards: the seam along the bottom, from the base back to the tip.

There is no fourth edge at the nose. The surface's whole first row of control points is the tip, so the face's
boundary there has shrunk to a point (a **pole**). All four side faces meet at that one tip vertex.

### 3. Geometry: where things are

```
#1=CARTESIAN_POINT('',(0.04017800193690645,0.0,-11.824249569057999));
#211=B_SPLINE_CURVE_WITH_KNOTS('',3,(#11,#12,...,#210),.UNSPECIFIED.,.F.,.F.,(4,1,1,...,1,4),
                               (0.,0.00502,...,1.),.UNSPECIFIED.);
#3887=B_SPLINE_SURFACE_WITH_KNOTS('',3,3,((#887,...,#901),(#902,...),...),.UNSPECIFIED.,.F.,.F.,.F.,
                                  (4,1,...,1,4),(4,1,1,1,1,1,1,1,1,1,1,1,4),(0.,...,1.),(0.,0.0208,...,0.25),
                                  .UNSPECIFIED.);
#12919=PLANE('',#12918);
```

`B_SPLINE_SURFACE_WITH_KNOTS` has 13 arguments:

1. name;
2–3. degree in u and in v: 3 and 3, bicubic;
4. control points: a list of rows. Each row is one value of u and runs along v. Here a row runs a quarter of the
   way round the body and the rows march from nose to tail;
5. surface form (`.UNSPECIFIED.` = general);
6–7. whether u and v are closed;
8. whether it self-intersects;
9–10. knot multiplicities in u and in v;
11–12. the distinct knot values in u and in v;
13. knot type.

**Knots, briefly.** A B-spline is a chain of polynomial pieces. The knots are the parameter values where the pieces
join; each piece is shaped by its `degree + 1` nearest control points. The multiplicities say how often each knot
repeats.
* A knot repeated `degree + 1` times (`4` for a cubic) at both ends **clamps** the spline: it starts exactly on its
  first control point and ends exactly on its last.
* Repeating an interior knot `degree` times makes the curve pass through a control point there, which is how the
  loft is cut into patches.
* Interior knots of multiplicity 1 keep a cubic smooth up to curvature (C2).
* The number of control points is always `sum(multiplicities) − degree − 1`.

`B_SPLINE_CURVE_WITH_KNOTS` is the one-direction version: the seams (curves along the body) and the base outline
pieces. Each edge's curve is exactly the boundary row or column of the faces either side. That's why the two faces
meeting there match exactly, not just within a tolerance.

`PLANE` is a plane defined by an `AXIS2_PLACEMENT_3D` (origin, normal, in-plane x direction). The flat base cap lies
on one.

## How a lofted body is laid out

```
              seam w=0 (right)           ┌────────── base cap (PLANE), at a blunt end
 tip ●═══════════════════════════════════╡
 (pole)       side 1   (right → bottom)  │   each side face: one bicubic B-spline patch,
       ═══════════════════════════════════   u nose → tail, w a quarter of the way round
              seam w=1/4 (bottom)        │
              side 2 … 4                 │   seams and base pieces: B-spline curve edges,
                                         ┘   each shared by the two faces it joins
```

`vsp_ntop/loft.py` fits one smooth surface all the way round the body (periodic in w). It then cuts that surface
exactly into the four patches by inserting knots, so the patches join with the same smoothness as the uncut
surface. With Faces Around 8, 12, 16 ... the layout is the same with more, narrower side faces and seams. Four faces rather than one also avoids a single face that wraps round and meets itself along a seam, which
some CAD readers handle badly.

At a pointed end, all four patches meet at the pole vertex. At a blunt end, each patch's last row of control points
is a boundary curve in the end plane. Those four curves form the loop of the flat cap. The cap's loop runs the other
way from the side faces' use of the same edges, so each edge is used once forwards and once backwards, as a closed
shell requires.

**Orientation rule used.** For the side faces, `dS/du × dS/dw` points out of the body; `write_loft` reverses w
if it doesn't. Each face loop then runs anticlockwise in the (u, w) parameter plane, which is anticlockwise about
the outward normal. So all faces have `same_sense .T.`.

## Reading a file yourself

* Search for `MANIFOLD_SOLID_BREP` (solids), `CLOSED_SHELL` / `OPEN_SHELL` (watertight or not) and `ADVANCED_FACE`
  (faces, with their names).
* Follow the `#n` references downward from a face to see its surface and edges.
* Units: search for `SI_UNIT`.
* `validation/cad_export/occ_check.py file.stp` reads a file with OpenCASCADE (through gmsh) and lists its solids,
  faces, edges, volume and area.

## What OpenVSP writes, for comparison

| | OpenVSP (stitched export) | nTop_vsp loft |
|---|---|---|
| Schema | AP203 `CONFIG_CONTROL_DESIGN` | AP214 `AUTOMOTIVE_DESIGN` |
| Unit | metres | millimetres |
| Faces | one bicubic patch per OpenVSP surface strip, each spanning the sections, plus cap patches | four B-spline faces round the body plus planar caps |
| Edges | `SURFACE_CURVE`s with `PCURVE`s (the curve also given in each face's (u, w) parameters) | plain 3D B-spline curves; readers compute the parameter-space curves themselves |
| Flat ends | B-spline patches collapsing to a point at the centre (nTop flags these as degenerate) | one `PLANE` face (OpenVSP's files get the same with the CAD Export block's Planar End Caps option) |
