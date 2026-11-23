"""Read a STEP file with OpenCASCADE (gmsh) and report what it finds: solids, shells, faces, volume, area, and
optionally the distance from given points to the CAD surface.

    dev/venv/Scripts/python.exe validation/cad_export/occ_check.py file.stp [points.npy]
"""
import sys
import numpy as np
import gmsh

path = sys.argv[1]
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("Geometry.OCCFixDegenerated", 0)
gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 0)
gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 0)
gmsh.option.setNumber("Geometry.OCCSewFaces", 0)
ents = gmsh.model.occ.importShapes(path, highestDimOnly=False)
gmsh.model.occ.synchronize()
vols = gmsh.model.getEntities(3)
faces = gmsh.model.getEntities(2)
print(f"{path.split('/')[-1]}: solids {len(vols)}, faces {len(faces)}, edges {len(gmsh.model.getEntities(1))}, "
      f"vertices {len(gmsh.model.getEntities(0))}")
for v in vols:
    print(f"   solid volume {gmsh.model.occ.getMass(3, v[1]) / 1e3:.3f} cm3")
area = sum(gmsh.model.occ.getMass(2, f[1]) for f in faces)
print(f"   surface area {area / 1e2:.3f} cm2")
for f in faces:
    print(f"   face {f[1]} {gmsh.model.getType(2, f[1])} area {gmsh.model.occ.getMass(2, f[1]) / 1e2:.3f} cm2")
if len(sys.argv) > 2:
    pts = np.load(sys.argv[2])
    best = np.full(len(pts), np.inf)
    for f in faces:
        c = np.array(gmsh.model.getClosestPoint(2, f[1], pts.ravel())[0]).reshape(-1, 3)
        best = np.minimum(best, np.linalg.norm(c - pts, axis=1))
    print(f"   points to CAD: mean {best.mean():.2e} max {best.max():.2e} mm (n={len(pts)})")
gmsh.finalize()
