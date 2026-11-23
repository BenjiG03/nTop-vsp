import sys, gmsh, numpy as np
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.occ.importShapes(sys.argv[1], highestDimOnly=False); gmsh.model.occ.synchronize()
N = 120
for _, s in gmsh.model.getEntities(2):
    lo, hi = gmsh.model.getParametrizationBounds(2, s)
    us = np.linspace(lo[0], hi[0], N); vs = np.linspace(lo[1], hi[1], N)
    uv = np.array([(u, v) for u in us for v in vs])
    d = np.array(gmsh.model.getDerivative(2, s, uv.reshape(-1).tolist())).reshape(-1, 6)
    du = np.linalg.norm(d[:, :3], axis=1); dv = np.linalg.norm(d[:, 3:], axis=1)
    cr = np.linalg.norm(np.cross(d[:, :3], d[:, 3:]), axis=1) / np.maximum(du * dv, 1e-300)
    ru, rv = du / np.median(du), dv / np.median(dv)
    bad = (ru < 1e-3) | (rv < 1e-3) | (cr < 1e-3)
    if bad.any():
        where = uv[bad]
        p = np.array(gmsh.model.getValue(2, s, where[:1].reshape(-1).tolist()))
        print(f"face {s}: {bad.sum()} singular grid points, min |du|/med {ru.min():.2e} min |dv|/med {rv.min():.2e} min sin(angle) {cr.min():.2e}; u range {where[:,0].min():.3f}-{where[:,0].max():.3f} (of {lo[0]:.3f}-{hi[0]:.3f}) v range {where[:,1].min():.3f}-{where[:,1].max():.3f} (of {lo[1]:.3f}-{hi[1]:.3f}) at {np.round(p, 2)}")
gmsh.finalize()
