"""
Library Features:

Name:          lib_data_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230522'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd
import warnings

from pytesmo.grid.resample import resample_to_grid

from lib_utils_io import create_darray_2d

# Logging
warnings.simplefilter("ignore", UserWarning)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get data
def get_data(dframe_obj_in, variable_name, variable_geo_x='lon', variable_geo_y='lat'):

    data_values = dframe_obj_in[variable_name].values
    data_geo_x = dframe_obj_in[variable_geo_x].values
    data_geo_y = dframe_obj_in[variable_geo_y].values
    data_indexes = dframe_obj_in.index.values

    data_obj = {variable_name: data_values, variable_geo_x: data_geo_x, variable_geo_y: data_geo_y}

    dframe_obj_out = pd.DataFrame(index=data_indexes, data=data_obj)

    return dframe_obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter data
def filter_data(dframe_obj_in,
                variable_name, variable_geo_x='lon', variable_geo_y='lat',
                min_value=0.0, max_value=None,
                no_data_value=np.nan, no_data_skip_rows=True):

    dframe_obj_tmp = dframe_obj_in.copy()
    if min_value is not None:
        dframe_obj_tmp[dframe_obj_tmp[variable_name] < min_value] = no_data_value
    if max_value is not None:
        dframe_obj_tmp[dframe_obj_tmp > max_value] = no_data_value

    if no_data_skip_rows:
        if np.isnan(no_data_value):
            dframe_obj_out = dframe_obj_tmp.dropna()
        else:
            dframe_obj_out = dframe_obj_tmp.drop(dframe_obj_tmp[dframe_obj_tmp[variable_name] == no_data_value].index)
    else:
        dframe_obj_out = dframe_obj_tmp.copy()

    return dframe_obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to mask data
def mask_data(da_obj_in, da_reference, mask_value_min=0, mask_value_max=None, mask_no_data=np.nan,
              var_name_data='variable', var_name_geo_x='longitude', var_name_geo_y='latitude',
              coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
              ):

    data_values = da_obj_in.values
    geo_x_values = da_obj_in[var_name_geo_x].values
    geo_y_values = da_obj_in[var_name_geo_y].values
    mask_values = da_reference.values

    if mask_value_min is not None:
        data_values[mask_values < mask_value_min] = mask_no_data
    if mask_value_max is not None:
        data_values[mask_values > mask_value_max] = mask_no_data

    # method to create data array
    da_obj_out = create_darray_2d(
        data_values, geo_x_values, geo_y_values, name=var_name_data,
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    return da_obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample args
def resample_data_args(
        grid_dframe, resampling_max_distance=18000,
        resampling_grid_resolution=None, resampling_min_neighbours=1, resampling_neighbours=8):

    grid_attrs = grid_dframe.attrs
    grid_bbox = grid_attrs['bbox']
    grid_bb_right, grid_bb_left = grid_attrs['bb_right'], grid_attrs['bb_left']
    grid_bb_bottom, grid_bb_top = grid_attrs['bb_bottom'], grid_attrs['bb_top']
    if resampling_grid_resolution is None:
        grid_res_geo_x, grid_res_geo_y = grid_attrs['res_lon'], grid_attrs['res_lat']
    else:
        grid_res_geo_x = grid_res_geo_y = resampling_grid_resolution

    # organize resample argument(s)
    resample_kwargs = {
        'geo_bbox': grid_bbox,
        'geo_x_right': grid_bb_right, 'geo_x_left': grid_bb_left,
        'geo_y_lower': grid_bb_bottom, 'geo_y_upper': grid_bb_top,
        'geo_x_res': grid_res_geo_x, 'geo_y_res': grid_res_geo_y,
        'resampling_max_distance': resampling_max_distance,
        'resampling_min_neighbours': resampling_min_neighbours,
        'resampling_neighbours': resampling_neighbours
    }

    return resample_kwargs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data
def resample_data_fx(file_dframe,
                     geo_x_left, geo_x_right, geo_y_lower, geo_y_upper, geo_x_res=0.01, geo_y_res=0.01,
                     var_name_data='surface_soil_moisture', var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
                     resampling_max_distance=18000, resampling_min_neighbours=1, resampling_neighbours=8,
                     **kwargs):

    geo_x_out = np.arange(geo_x_left, geo_x_right, geo_x_res)
    geo_y_out = np.arange(geo_y_lower, geo_y_upper, geo_y_res)

    data_in = file_dframe[var_name_data].values
    geo_x_in = file_dframe[var_name_geo_x].values
    geo_y_in = file_dframe[var_name_geo_y].values

    # make 2d grid out the 1D grid spacings
    geo_x_grid, geo_y_grid = np.meshgrid(geo_x_out, geo_y_out)

    # method to resample data to grid
    data_masked = resample_to_grid(
        {var_name_data: data_in},
        geo_x_in, geo_y_in, geo_x_grid, geo_y_grid, search_rad=resampling_max_distance,
        min_neighbours=resampling_min_neighbours, neighbours=resampling_neighbours)
    data_grid = data_masked[var_name_data].data

    # check south-north array orientation
    geo_y_upper, geo_y_lower = geo_y_grid[0, 0], geo_y_grid[-1, 0]
    if geo_y_upper < geo_y_lower:
        geo_y_grid = np.flipud(geo_y_grid)
        data_grid = np.flipud(data_grid)

    # method to create data array
    data_out = create_darray_2d(
        data_grid, geo_x_grid[0, :], geo_y_grid[:, 0], name=var_name_data,
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    return data_out
# ----------------------------------------------------------------------------------------------------------------------
