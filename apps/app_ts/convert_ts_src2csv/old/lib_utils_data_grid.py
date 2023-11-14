"""
Library Features:

Name:          lib_utils_data_grid_sm
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import numpy as np
import pandas as pd

from lib_data_io_nc import read_file_nc
from lib_data_io_binary import read_file_binary

from lib_utils_obj import create_darray_2d
from lib_utils_geo import convert_cn2s
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to join data point
def join_data_point(point_time, point_data, point_collection):

    if point_collection is None:
        point_collection = pd.DataFrame(data=point_data, index=[point_time])
    else:
        # Create tmp collection
        tmp_collection = pd.DataFrame(data=point_data, index=[point_time])
        # append new line to dataframe
        point_collection = pd.concat([point_collection, tmp_collection])

    return point_collection
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to extract grid values to point values
def extract_data_grid2point(
        da_sm, dframe_points, method_spatial_operation='average',
        method_spatial_mask=True, value_spatial_mask=1):

    values_sm = da_sm.values.ravel()

    point_data = {}
    for point_id, point_row in dframe_points.iterrows():
        point_name = point_row['point_name']
        point_idx_1d_list = point_row['point_idx_1d']
        point_cnet_list = point_row['point_cnet']

        if not isinstance(point_idx_1d_list, list):
            point_idx_1d_list = [point_idx_1d_list]
        if not isinstance(point_cnet_list, list):
            point_cnet_list = [point_cnet_list]

        point_sm_list_unmasked, point_sm_list_masked = [], []
        for point_idx_1d_step, point_cnet_step in zip(point_idx_1d_list, point_cnet_list):
            if method_spatial_mask:
                if point_cnet_step != value_spatial_mask:
                    point_sm_step = values_sm[point_idx_1d_step]
                    point_sm_list_unmasked.append(point_sm_step)
                else:
                    point_sm_step = values_sm[point_idx_1d_step]
                    point_sm_list_masked.append(point_sm_step)

        if method_spatial_operation == 'average':
            point_sm_array = np.array(point_sm_list_unmasked, dtype=float)
            point_sm_cmp = float(np.mean(point_sm_array))
        else:
            log_stream.error(' ===> Spatial data operation "' + method_spatial_operation + '" is not supported')
            raise NotImplementedError('Case not implemented yet. Only "average" method is available.')

        point_data[point_name] = point_sm_cmp

    return point_data
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get data binary
def get_data_binary(file_name, da_geo, da_cn, da_cnet=None, mask_cnet=True, mask_limits=True,
                    geo_x_tag='Longitude', geo_y_tag='Latitude',
                    value_cnet_mask=1, value_sm_mask=-1,
                    value_sm_min=0.0, value_sm_max=1.0, var_sm_scale_factor=100.0):

    geo_x_1d = da_geo[geo_x_tag]
    geo_y_1d = da_geo[geo_y_tag]
    geo_values = da_geo.values
    cn_values = da_cn.values
    cnet_values = da_cnet.values

    vtot_values = read_file_binary(file_name, data_geo=geo_values)
    vmax_values = convert_cn2s(cn_values, geo_values)
    sm_values = vtot_values / vmax_values

    if mask_cnet:
        sm_values[cnet_values == value_cnet_mask] = value_sm_mask

    if mask_limits:
        sm_values[sm_values < value_sm_min] = np.nan
        sm_values[sm_values > value_sm_max] = np.nan

    sm_values = sm_values * var_sm_scale_factor

    da_sm = create_darray_2d(sm_values, geo_x_1d, geo_y_1d,
                             coord_name_x=geo_x_tag, coord_name_y=geo_y_tag,
                             dim_name_x=geo_x_tag, dim_name_y=geo_y_tag)

    return da_sm
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get data netcdf
def get_data_nc(file_name, da_geo, da_cn, da_cnet=None, mask_cnet=True, mask_limits=True,
                geo_x_tag='west_east', geo_y_tag='south_north',
                value_cnet_mask=1, value_sm_mask=-1,
                value_sm_min=0, value_sm_max=1, var_sm_scale_factor=100.0):

    geo_x_1d = da_geo[geo_x_tag].values
    geo_y_1d = da_geo[geo_y_tag].values
    # geo_values = da_geo.values
    # cn_values = da_cn.values
    cnet_values = da_cnet.values

    sm_values, geo_x_values, geo_y_values = read_file_nc(file_name)

    if mask_cnet:
        sm_values[cnet_values == value_cnet_mask] = value_sm_mask

    if mask_limits:
        sm_values[sm_values < value_sm_min] = np.nan
        sm_values[sm_values > value_sm_max] = np.nan

    sm_values = sm_values * var_sm_scale_factor

    da_sm = create_darray_2d(sm_values, geo_x_1d, geo_y_1d,
                             coord_name_x='west_east', coord_name_y='south_north',
                             dim_name_x='west_east', dim_name_y='south_north')

    return da_sm

# -------------------------------------------------------------------------------------
