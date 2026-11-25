"""The E2E cases, shared by the tests and the local notebook generator."""
import numpy as np

MACHS = (3.0, 8.0)
ALPHAS = (0.0, 10.0, 30.0)
SREF = 0.01  # m^2, the plate's planform
CREF = 0.1  # m
SPHERE_R = 0.05  # m
PLATE = (0.1, 0.1, 0.002)  # m: chord (x), span (y), thickness (z)
CONE_HALF_ANGLE = 15.0  # deg
CONE_LENGTH = 0.1  # m
MESH_TOL = 0.5e-3  # m (finer takes many GB)
COARSE_TOL = 2e-3  # m, for the pressure export
USER_CP_MAX = 1.5
DEG = np.radians(1.0)

# (case name, block kind, body, Stagnation Pressure index, extra choice inputs: name -> index)
CASES = (
    ("si_sphere_classical", "si", "sphere", 1, {}),
    ("mesh_sphere_modified", "mesh", "sphere", 0, {}),
    ("mesh_sphere_export", "mesh", "sphere_coarse", 0, {"Export Pressures": 1}),
    ("mesh_sphere_user", "mesh", "sphere", 2, {}),
    ("mesh_plate_classical", "mesh", "plate", 1, {}),
    ("tw_plate", "tw", "plate", 0, {}),
    ("tw_cone", "tw", "cone", 0, {"Impact Method": 1}),
)
READ_CASE = ("mesh_sphere_export", 2)
READ_ROWS = (0, 1, 2)
FAKE_JUMP = ((0.0, 0.0, 0.0, 0.0), (0.01, 0.0, 0.0, 100.0), (0.0, 0.02, 0.0, -50.0))  # x, y, z (m), Delta p (Pa)
REFERENCE_PRESSURE = 1000.0  # Pa
FIELD_CASE = 2  # Mach 3, alpha 10
FIELD_POINTS_MM = ((-50.0, 0.0, 0.0), (0.0, 0.0, -50.0), (0.0, 0.0, 50.0), (0.0, -50.0, 0.0))
