"""nTop round trip: import a STEP (Import Part), CAD bodies -> implicit, evaluate at the original nTop mesh vertices."""
import os, sys
import numpy as np
ROOT = r"C:/Users/Benji/Documents/nTop_playground/nTop_vsp"
sys.path.insert(0, ROOT + "/build"); sys.path.insert(0, ROOT + "/validation/rasaero"); sys.path.insert(0, ROOT + "/validation/cad_export")
import rocket as rk
from ntopfile import Graph, Ref
W = os.environ.get("RT_DIR", ROOT + "/runs/cad_export/ntimport"); os.makedirs(W, exist_ok=True)
f, idx, pts_csv = sys.argv[1], int(sys.argv[2]), sys.argv[3]
g = Graph()
part = g.block("import_part<file_path,bool,cad_healing_enum,bool>[2.0.0]", "part",
               [g.file_path(f), g.boolean(False), g.enum("cad_healing_enum", int(os.environ.get("HEAL", "0"))), g.boolean(False)], "Part")
bodies = g.block("get_cad_bodies<part>", "list<brep>", [part], "Bodies")
body = g.block("core.list_element<list_interface,integer>", "brep", [bodies, g.integer(idx)], "Body") if os.environ.get("ELEM") else Ref(bodies, (idx,))
imp = g.block("brep_to_implicit<brep,real>[1.2.0]", "implicit", [body, g.real(float(os.environ.get("TOL", "0.0001")), {"length": 1})], "Implicit")
one = g.leaf("real_field", {"expression": "1", "units": {}})
pm = g.block("import_real_point_map<file_path,unit_length_enum,real_field>", "real_point_map", [g.file_path(pts_csv), g.unit_length("mm"), one], "Pts")
vals = g.block("evaluate_field<real_field,list<point>>", "list<real>", [imp, Ref(pm, ("points",))], "Vals")
m = g.block("point_map<list<point>,list<real>>", "real_point_map", [Ref(pm, ("points",)), vals], "M")
outp = f"{W}/{os.path.basename(f)}.{idx}.vals.csv"
out = g.block("export_point_map<file_path,real_point_map,unit_length_enum>", "point_map_file", [g.file_path(outp), m, g.unit_length("mm")], "Out")
res, log, dt = rk.run(g, [], out, "point_map_file", f"{W}/{os.path.basename(f)}.{idx}.rt.ntop", timeout=1200)
print("==", os.path.basename(f), idx, round(dt), "s", [l[:200] for l in log.splitlines() if "[E]" in l or "[W]" in l][:5])
if os.path.exists(outp):
    a = np.abs(np.loadtxt(outp, delimiter=",", ndmin=2)[:, 3]) * 1000
    print(f"   |field| at the original body's vertices: mean {a.mean():.3f} p99 {np.percentile(a, 99):.3f} max {a.max():.3f} mm (n={len(a)})")
