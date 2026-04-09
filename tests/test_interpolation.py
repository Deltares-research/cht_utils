"""Tests for cht_utils.interpolation."""

import numpy as np
import pytest

from cht_utils.interpolation import interp2, interp2_bilinear, interp3


class TestInterp2:
    """Tests for regular grid interpolation."""

    def setup_method(self):
        self.x0 = np.array([0.0, 1.0, 2.0, 3.0])
        self.y0 = np.array([0.0, 1.0, 2.0])
        # z = x + y
        self.z0 = np.array([[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]], dtype=float)

    def test_exact_grid_points(self):
        x1 = np.array([0.0, 1.0, 2.0])
        y1 = np.array([0.0, 1.0, 2.0])
        z1 = interp2(self.x0, self.y0, self.z0, x1, y1)
        np.testing.assert_allclose(z1, [0.0, 2.0, 4.0])

    def test_interpolated_points(self):
        x1 = np.array([0.5])
        y1 = np.array([0.5])
        z1 = interp2(self.x0, self.y0, self.z0, x1, y1)
        np.testing.assert_allclose(z1, [1.0])

    def test_2d_target(self):
        x1 = np.array([[0.5, 1.5], [0.5, 1.5]])
        y1 = np.array([[0.5, 0.5], [1.5, 1.5]])
        z1 = interp2(self.x0, self.y0, self.z0, x1, y1)
        assert z1.shape == (2, 2)
        np.testing.assert_allclose(z1, [[1.0, 2.0], [2.0, 3.0]])

    def test_out_of_bounds_returns_nan(self):
        x1 = np.array([5.0])
        y1 = np.array([5.0])
        z1 = interp2(self.x0, self.y0, self.z0, x1, y1)
        assert np.isnan(z1[0])

    def test_nearest_method(self):
        x1 = np.array([0.3])
        y1 = np.array([0.3])
        z1 = interp2(self.x0, self.y0, self.z0, x1, y1, method="nearest")
        np.testing.assert_allclose(z1, [0.0])


class TestInterp2Bilinear:
    """Tests for bilinear interpolation."""

    def test_center_of_cell(self):
        xp = np.array([0.0, 1.0, 2.0])
        yp = np.array([0.0, 1.0, 2.0])
        zp = np.array([[0, 1, 2], [1, 2, 3], [2, 3, 4]], dtype=float)
        z = interp2_bilinear(xp, yp, zp, np.array([0.5]), np.array([0.5]))
        np.testing.assert_allclose(z, [1.0])

    def test_exact_grid_point(self):
        xp = np.array([0.0, 1.0])
        yp = np.array([0.0, 1.0])
        zp = np.array([[0, 1], [2, 3]], dtype=float)
        z = interp2_bilinear(xp, yp, zp, np.array([0.0]), np.array([0.0]))
        np.testing.assert_allclose(z, [0.0])


class TestInterp3:
    """Tests for unstructured interpolation."""

    def test_scattered_to_grid(self):
        # Source: scattered points on z = x + y
        x0 = np.array([0.0, 1.0, 0.0, 1.0])
        y0 = np.array([0.0, 0.0, 1.0, 1.0])
        z0 = x0 + y0

        x1 = np.array([0.5])
        y1 = np.array([0.5])
        z1 = interp3(x0, y0, z0, x1, y1)
        np.testing.assert_allclose(z1, [1.0], atol=0.01)

    def test_2d_target(self):
        x0 = np.array([0, 1, 0, 1, 0.5])
        y0 = np.array([0, 0, 1, 1, 0.5])
        z0 = x0 * 2 + y0 * 3

        x1 = np.array([[0.25, 0.75], [0.25, 0.75]])
        y1 = np.array([[0.25, 0.25], [0.75, 0.75]])
        z1 = interp3(x0, y0, z0, x1, y1)
        assert z1.shape == (2, 2)
