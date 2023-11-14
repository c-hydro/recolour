"""
Library Features:

Name:          lib_interface_cci
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from pynetcf.time_series import GriddedNcOrthoMultiTs
from pygeogrids.netcdf import load_grid
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to interface cci time-series
class CCITs(GriddedNcOrthoMultiTs):

    def __init__(self, ts_path, grid_path=None):
        if grid_path is None:
            grid_path = os.path.join(ts_path, "grid.nc")

        grid = load_grid(grid_path)
        super(CCITs, self).__init__(ts_path, grid)
# -------------------------------------------------------------------------------------
