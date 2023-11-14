"""
Library Features:

Name:          lib_data_io_gpi
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
from pygeogrids.grids import BasicGrid
from netCDF4 import Dataset
# ----------------------------------------------------------------------------------------------------------------------


file_grid = Dataset(grid_path)

domain_lon_1d = file_grid['longitude'][:]
domain_lat_1d = file_grid['latitude'][:]
domain_lon_2d, domain_lat_2d = np.meshgrid(domain_lon_1d, domain_lat_1d)

# *************************************************************************
# The grid which is created is, as always, of the proprietary types
# defined by TUW as BasicGrid or CellGrid
# this is what it is returned
domain_grid = BasicGrid(domain_lon_2d.flatten(), domain_lat_2d.flatten()).to_cell_grid(cellsize=5.)
# *************************************************************************