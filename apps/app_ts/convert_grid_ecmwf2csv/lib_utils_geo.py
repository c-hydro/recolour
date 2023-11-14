"""
Library Features:

Name:          lib_utils_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.5.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import pyresample
import numpy as np

from lib_info_args import logger_name
from lib_info_args import (geo_coord_name_x, geo_coord_name_y,
                           geo_var_name_x, geo_var_name_y,
                           geo_dim_name_x, geo_dim_name_y)

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------------
# method to get value using coordinates x,y
def get_grid_value_from_xy(geo_grid_da, geo_point_x, geo_point_y, select_method='nearest'):

    geo_selector = {geo_var_name_x: geo_point_x, geo_var_name_y: geo_point_y, 'method': select_method}
    geo_grid_value = geo_grid_da.sel(**geo_selector)

    geo_grid_x = float(geo_grid_value[geo_var_name_x].values.ravel()[0])
    geo_grid_y = float(geo_grid_value[geo_var_name_y].values.ravel()[0])
    geo_grid_value = float(geo_grid_value.values.ravel()[0])

    return geo_grid_value, geo_grid_x, geo_grid_y
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get indexes using a spatial window
def get_idx_by_win(geo_idx_x, geo_index_y, geo_grid_x=None, geo_grid_y=None, geo_win_x=1, geo_win_y=1):

    top, bottom = geo_index_y + geo_win_y, geo_index_y - geo_win_y
    left, right = geo_idx_x - geo_win_x, geo_idx_x + geo_win_x

    geo_idx_left_right = np.arange(left, right + 1, 1)
    geo_idx_bottom_top = np.arange(bottom, top + 1, 1)

    geo_idx_list, geo_coord_list = [], []
    for ix in geo_idx_left_right:
        for iy in geo_idx_bottom_top:
            geo_idx_step = [ix, iy]
            if geo_idx_step not in geo_idx_list:
                if (geo_grid_x is not None) and (geo_grid_y is not None):
                    geo_coord_x, geo_coord_y = geo_grid_x[ix, iy], geo_grid_y[ix, iy]
                    geo_coord_step = [geo_coord_x, geo_coord_y]
                    geo_coord_list.append(geo_coord_step)

                geo_idx_list.append(geo_idx_step)

    return geo_idx_list, geo_coord_list
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get indexes using coordinates x,y
def get_grid_idx_from_xy(geo_point_x, geo_point_y, geo_grid_x, geo_grid_y,
                         geo_radius_influence=50000, geo_neighbours=1):

    geo_grid_obj = pyresample.geometry.GridDefinition(lats=geo_grid_y, lons=geo_grid_x)

    array_point_x, array_point_y = np.array([geo_point_x]), np.array([geo_point_y])

    geo_swath_obj = pyresample.geometry.SwathDefinition(lons=array_point_x, lats=array_point_y)

    # Determine nearest (w.r.t. great circle distance) neighbour in the grid.
    _, _, array_index_1d, array_distance = pyresample.kd_tree.get_neighbour_info(
        source_geo_def=geo_grid_obj, target_geo_def=geo_swath_obj, radius_of_influence=geo_radius_influence,
        neighbours=geo_neighbours)

    array_index_2d = np.unravel_index(array_index_1d, geo_grid_obj.shape)

    geo_index_continuous = int(array_index_1d[0])
    geo_distance = float(array_distance[0])
    geo_index_x = int(array_index_2d[0][0])
    geo_index_y = int(array_index_2d[1][0])

    geo_check_x = geo_grid_x[geo_index_x, geo_index_y]
    geo_check_y = geo_grid_y[geo_index_x, geo_index_y]

    return geo_index_continuous, geo_index_x, geo_index_y, geo_distance
# ----------------------------------------------------------------------------------------------------------------------
