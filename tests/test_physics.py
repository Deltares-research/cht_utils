"""Tests for cht_utils.physics subpackage."""

import numpy as np
import pytest

from cht_utils.physics import (
    disper,
    disper_fentonmckee,
    deshoal,
    wavecelerity,
    split_waves_guza,
)
from cht_utils.physics.runup_vo21 import runup_vo21


class TestDisper:
    def test_known_values(self):
        """Deep water: k ~ omega^2 / g."""
        omega = 2 * np.pi / 10  # 10s wave
        h = 100.0  # deep water
        k = disper(omega, h)
        # Deep water approximation: k = omega^2 / g
        k_deep = omega**2 / 9.81
        np.testing.assert_allclose(k, k_deep, rtol=0.01)

    def test_shallow_water(self):
        """Shallow water: k ~ omega / sqrt(g*h)."""
        omega = 2 * np.pi / 100  # very long wave
        h = 1.0
        k = disper(omega, h)
        k_shallow = omega / np.sqrt(9.81 * h)
        np.testing.assert_allclose(k, k_shallow, rtol=0.05)

    def test_returns_positive_for_positive_omega(self):
        k = disper(1.0, 10.0)
        assert np.all(np.array(k) > 0)

    def test_list_inputs(self):
        k = disper([0.5, 1.0], [10.0, 10.0])
        assert len(k) == 2


class TestDisperFentonMcKee:
    def test_deep_water(self):
        sigma = 2 * np.pi / 10
        d = 100.0
        k = disper_fentonmckee(sigma, d)
        k_deep = sigma**2 / 9.81
        np.testing.assert_allclose(k, k_deep, rtol=0.05)

    def test_array_input(self):
        sigma = np.array([0.5, 1.0])
        d = np.array([10.0, 20.0])
        k = disper_fentonmckee(sigma, d)
        assert k.shape == (2,)


class TestDeshoal:
    def test_same_depth(self):
        """No shoaling if depths are equal."""
        hm0 = 1.0
        result = deshoal(hm0, 10.0, 20.0, 20.0)
        np.testing.assert_allclose(result, hm0, rtol=1e-10)

    def test_different_depths(self):
        """De-shoaling adjusts height based on group velocity ratio."""
        hm0 = 1.0
        result = deshoal(hm0, 10.0, 5.0, 20.0)
        # Result should differ from input when depths differ
        assert result != hm0
        assert result > 0


class TestWavecelerity:
    def test_positive_celerity(self):
        cg = wavecelerity(10.0, 20.0)
        assert cg > 0

    def test_deep_water_cg(self):
        """Deep water group velocity ~ g*T/(4*pi)."""
        Tp = 10.0
        d = 200.0
        cg = wavecelerity(Tp, d)
        cg_deep = 9.81 * Tp / (4 * np.pi)
        np.testing.assert_allclose(cg, cg_deep, rtol=0.01)


class TestSplitWavesGuza:
    def test_returns_expected_columns(self):
        import pandas as pd

        n = 100
        t = np.linspace(0, 10, n)
        df = pd.DataFrame({
            "z": 0.5 * np.sin(2 * np.pi * t / 5) + 1.0,
            "u": 0.1 * np.sin(2 * np.pi * t / 5),
        })
        zb = -1.0
        result = split_waves_guza(df, zb)
        assert "zin" in result.columns
        assert "zout" in result.columns
        assert "uin" in result.columns
        assert "uout" in result.columns
        assert len(result) == n


class TestRunupVo21:
    def test_basic_computation(self):
        """Smoke test: runup should be positive for reasonable inputs."""
        hm0 = np.array([1.0, 2.0])
        tp = 10.0
        ztide = 0.0
        sl1 = np.array([0.05, 0.05])
        sl2 = np.array([0.1, 0.1])

        r = runup_vo21(hm0, tp, ztide, sl1, sl2, "slope")
        r.compute_r2p(np.array([30.0, 30.0]), np.array([0.0, 0.0]))

        assert hasattr(r, "r2p")
        assert np.all(r.r2p > 0)
        # Larger waves should give larger runup
        assert r.r2p[1] > r.r2p[0]

    def test_setup_component(self):
        hm0 = np.array([1.5])
        tp = 8.0
        sl1 = np.array([0.04])
        sl2 = np.array([0.08])

        r = runup_vo21(hm0, tp, 0.0, sl1, sl2, "slope")
        r.compute_setup(np.array([25.0]), np.array([0.0]))
        assert r.setup > 0
