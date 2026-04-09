"""Tests for cht_utils.probabilistic subpackage."""

import pytest


class TestImports:
    def test_import_prob_floodmaps(self):
        from cht_utils.probabilistic import prob_floodmaps
        assert callable(prob_floodmaps)

    def test_import_merge_nc_his(self):
        from cht_utils.probabilistic import merge_nc_his
        assert callable(merge_nc_his)

    def test_import_merge_nc_map(self):
        from cht_utils.probabilistic import merge_nc_map
        assert callable(merge_nc_map)
