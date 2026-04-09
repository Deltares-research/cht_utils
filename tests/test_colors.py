"""Tests for cht_utils.colors."""

import matplotlib as mpl
import numpy as np

from cht_utils.colors import cm2png, read_color_maps, read_colormap, rgb2hex


class TestRgb2Hex:
    def test_black(self):
        assert rgb2hex((0, 0, 0)) == "000000"

    def test_white(self):
        assert rgb2hex((255, 255, 255)) == "ffffff"

    def test_red(self):
        assert rgb2hex((255, 0, 0)) == "ff0000"

    def test_mixed(self):
        assert rgb2hex((16, 32, 64)) == "102040"


class TestReadColormap:
    def test_reads_rgb_file(self, tmp_path):
        cmap_file = tmp_path / "test.txt"
        cmap_file.write_text("0.0 0.0 1.0\n0.5 0.5 0.5\n1.0 0.0 0.0\n")
        result = read_colormap(str(cmap_file))
        assert result.shape == (3, 3)
        np.testing.assert_allclose(result[0], [0.0, 0.0, 1.0])
        np.testing.assert_allclose(result[2], [1.0, 0.0, 0.0])


class TestReadColorMaps:
    def test_registers_colormaps(self, tmp_path):
        cmap_file = tmp_path / "mymap.txt"
        cmap_file.write_text("0.0 0.0 1.0\n0.5 0.5 0.5\n1.0 0.0 0.0\n")
        result = read_color_maps(str(tmp_path))
        assert "mymap" in result


class TestCm2Png:
    def test_creates_horizontal_png(self, tmp_path):
        out = tmp_path / "cbar.png"
        cmap = mpl.colormaps["viridis"]
        cm2png(cmap, str(out), orientation="horizontal", vmin=0, vmax=10)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_creates_vertical_png(self, tmp_path):
        out = tmp_path / "cbar_v.png"
        cmap = mpl.colormaps["viridis"]
        cm2png(cmap, str(out), orientation="vertical")
        assert out.exists()
