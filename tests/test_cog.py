"""Tests for cht_utils.cog subpackage.

These are import and basic logic tests. Full COG conversion tests
require rasterio/rioxarray with actual raster data.
"""


class TestImports:
    def test_import_geotiff_to_cog(self):
        from cht_utils.cog import geotiff_to_cog, is_cog

        assert callable(geotiff_to_cog)
        assert callable(is_cog)

    def test_import_netcdf_to_cog(self):
        from cht_utils.cog import netcdf_to_cog

        assert callable(netcdf_to_cog)

    def test_import_xyz_to_cog(self):
        from cht_utils.cog import xyz_to_cog

        assert callable(xyz_to_cog)


class TestIsCog:
    def test_nonexistent_file(self):
        from cht_utils.cog import is_cog

        assert is_cog("nonexistent_file.tif") is False
