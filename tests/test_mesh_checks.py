"""Checks of the mesh driver that need no nTop."""
from __future__ import annotations

import numpy as np
import pytest

from vsp_ntop.common import DriverError
from vsp_ntop.run_newtonian_mesh import unit_normals


def _cube_faces():
    """Vector areas of the six faces of a unit cube, centred at the origin (outward)."""
    return np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def test_matching_body_normals_pass():
    s = _cube_faces()
    n = unit_normals(s * 3.0, s, lambda msg: None)  # gradients of any length
    assert np.allclose(n, s)


def test_body_normals_not_matching_the_mesh_are_refused():
    """The Body evaluated somewhere else than the Mesh (wrong units, moved body): all normals point one way."""
    s = _cube_faces()
    gradients = np.tile([1.0, 0.2, 0.1], (6, 1))
    with pytest.raises(DriverError, match="do not match the Mesh"):
        unit_normals(gradients, s, lambda msg: None)
