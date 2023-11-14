"""
Library Features:

Name:          lib_utils_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220420'
Version:       '1.5.0'
"""
#######################################################################################
# Libraries
import logging
import pyresample
import numpy as np

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to get value using coordinates x,y
def get_grid_value_from_xy(geo_grid_da, geo_point_x, geo_point_y, select_method='nearest'):

    geo_grid_value = geo_grid_da.sel(Longitude=geo_point_x, Latitude=geo_point_y, method=select_method)

    geo_grid_x = float(geo_grid_value['Longitude'].values.ravel()[0])
    geo_grid_y = float(geo_grid_value['Latitude'].values.ravel()[0])
    geo_grid_value = float(geo_grid_value.values.ravel()[0])

    return geo_grid_value, geo_grid_x, geo_grid_y
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get indexes using a spatial window
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
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get indexes using coordinates x,y
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

    return geo_index_continuous, geo_index_x, geo_index_y, geo_distance
# -------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Method to convert curve number to s (vmax)
def convert_cn2s(data_cn, data_terrain):

    data_s = (1000.0 / data_cn - 10) * 25.4
    data_s[data_cn <= 0] = np.nan
    data_s[data_cn > 100] = np.nan

    data_s[(data_terrain >= 0) & (data_s < 1.0)] = 1.0

    data_s[data_s < 0] = 0.0

    data_s[data_terrain < 0] = np.nan

    data_s[0, :] = np.nan
    data_s[-1, :] = np.nan
    data_s[:, 0] = np.nan
    data_s[:, -1] = np.nan

    # Debug
    # plt.figure()
    # plt.imshow(data_s)
    # plt.colorbar()
    # plt.show()

    return data_s
# ------------------------------------------------------------------------------------
