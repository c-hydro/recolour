"""
Library Features:

Name:          lib_utils_debug
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231212'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
from repurpose.resample import resample_to_grid

from lib_info_args import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debugging
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert cell to grid
def convert_cell_to_grid(obj_cell, var_name_data,
                         var_name_geo_x='longitude', var_name_geo_y='latitude',
                         search_rad=50000, fill_values=np.nan, min_neighbours=1, neighbours=1):

    var_data_1d = obj_cell[var_name_data].values
    var_x_1d = obj_cell[var_name_geo_x].values
    var_y_1d = obj_cell[var_name_geo_y].values
    var_x_2d, var_y_2d = np.meshgrid(np.unique(var_x_1d), np.unique(var_y_1d))

    var_data_obj = resample_to_grid(
        {'data': var_data_1d},
        var_x_1d, var_y_1d, var_x_2d, var_y_2d,
        search_rad=search_rad, fill_values=fill_values,
        min_neighbours=min_neighbours, neighbours=neighbours)
    var_data_2d = np.flipud(var_data_obj['data'])

    return var_data_2d, var_x_2d, var_y_2d
# ----------------------------------------------------------------------------------------------------------------------
