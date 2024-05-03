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

from copy import deepcopy

from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel
from pytesmo.grid.resample import resample_to_grid
from pyresample.geometry import GridDefinition
from pyresample.kd_tree import resample_nearest, resample_gauss, resample_custom

from pytesmo.time_series.filters import exp_filter

from lib_utils_io import create_darray_2d

# debugging
import matplotlib.pylab as plt

# logging
warnings.simplefilter("ignore", UserWarning)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute data swi
def compute_data_swi(ts_data, ts_var="swi", ts_method="exp_filter", ts_ctime=6):

    # get julian dates
    ts_jd = ts_data.index.to_julian_date().values
    # drop nan values
    ts_data.dropna(inplace=True)
    # get values
    ts_values_raw = ts_data.values.astype(np.float64)

    # iterate over characteristic times
    if ts_method == "exp_filter":
        ts_values_filtered = exp_filter(ts_values_raw, ts_jd, ctime=ts_ctime)
    else:
        logging.error(" ===> Method '" + ts_method + "' is not supported yet")
        raise NotImplementedError("Case not implemented yet")
    # define data output
    ts_data_out = pd.Series(ts_values_filtered, index=ts_data.index, name=ts_var)

    return ts_data_out

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to add data to the dataframe
def add_data(dframe_obj_in, time_reference='time_reference', time_data='time', var_time_delta='time_delta'):

    # add time delta to the dataframe
    time_reference = dframe_obj_in[time_reference].values
    time_data = dframe_obj_in[time_data].values

    time_delta = pd.DatetimeIndex(time_data) - pd.DatetimeIndex(time_reference)
    hours_delta = np.round(time_delta.total_seconds().values / 3600, 2)

    dframe_obj_in[var_time_delta] = hours_delta

    return dframe_obj_in
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get data
def get_data(dframe_obj_in, variable_name, variable_geo_x='longitude', variable_geo_y='latitude'):

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
# method to interpolate data
def interpolate_data(obj_da_in, obj_da_ref,
                     var_name_data='surface_soil_moisture', var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
                     interpolating_active=True,
                     interpolating_method='nn', interpolating_max_distance=18000,
                     interpolating_neighbours=8, interpolating_fill_value=np.nan,
                     ):

    # get data in
    values_in = obj_da_in.values
    geo_x_arr_in = obj_da_in[var_name_geo_x].values
    geo_y_arr_in = obj_da_in[var_name_geo_y].values
    geo_x_grid_in, geo_y_grid_in = np.meshgrid(geo_x_arr_in, geo_y_arr_in)

    # get data reference
    geo_x_arr_ref = obj_da_ref[var_name_geo_x].values
    geo_y_arr_ref = obj_da_ref[var_name_geo_y].values
    geo_x_grid_ref, geo_y_grid_ref = np.meshgrid(geo_x_arr_ref, geo_y_arr_ref)

    if (((geo_x_grid_in.shape[0] != geo_x_grid_ref.shape[0]) or (geo_y_grid_in.shape[1] != geo_y_grid_ref.shape[1]))
            or interpolating_active):

        geo_grid_in = GridDefinition(lons=geo_x_grid_in, lats=geo_y_grid_in)
        geo_grid_ref = GridDefinition(lons=geo_x_grid_ref, lats=geo_y_grid_ref)

        if interpolating_method == 'nn':
            values_out_tmp = resample_nearest(
                geo_grid_in, values_in, geo_grid_ref,
                radius_of_influence=interpolating_max_distance,
                fill_value=interpolating_fill_value)

        elif interpolating_method == 'gauss':
            values_out_tmp = resample_gauss(
                geo_grid_in, values_in, geo_grid_ref,
                radius_of_influence=interpolating_max_distance,
                neighbours=resampling_neighbours, sigmas=250000,
                fill_value=interpolating_fill_value)

        elif interpolating_method == 'idw':
            weight_fx = lambda r: 1 / r ** 2
            values_out_tmp = resample_custom(
                geo_grid_in, values_in, geo_grid_ref,
                radius_of_influence=interpolating_max_distance, neighbours=interpolating_neighbours,
                weight_funcs=weight_fx,
                fill_value=interpolating_fill_value)
        else:
            logging.error(' ===> Interpolating method "' + interpolating_method + '" is not available')
            raise NotImplemented('Case not implemented yet')

        if interpolating_fill_value is None:
            values_out_resampled = values_out_tmp.data
        else:
            values_out_resampled = deepcopy(values_out_tmp)

        obj_da_out = create_darray_2d(
            values_out_resampled, geo_x_grid_ref[0, :], geo_y_grid_ref[:, 0], name=var_name_data,
            coord_name_x=coord_name_x, coord_name_y=coord_name_y,
            dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    else:
        obj_da_out = deepcopy(obj_da_in)

    return obj_da_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample args
def resample_data_args(
        grid_dframe, resampling_max_distance=18000,
        resampling_grid_resolution=None, resampling_min_neighbours=1, resampling_neighbours=8,
        resampling_kernel_active=False,
        resampling_kernel_method='box', resampling_kernel_width=15,
        resampling_kernel_mode='center', resampling_kernel_stddev=4):

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
        'resampling_neighbours': resampling_neighbours,
        'filtering_active': resampling_kernel_active,
        'filtering_method': resampling_kernel_method,
        'filtering_width': resampling_kernel_width,
        'filtering_mode': resampling_kernel_mode,
        'filtering_stddev': resampling_kernel_stddev,
    }

    return resample_kwargs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data
def resample_data_fx(file_dframe,
                     geo_x_left, geo_x_right, geo_y_lower, geo_y_upper,
                     geo_x_res=0.01, geo_y_res=0.01, geo_no_data=np.nan,
                     var_name_data='surface_soil_moisture', var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
                     resampling_max_distance=18000, resampling_min_neighbours=1, resampling_neighbours=8,
                     filtering_active=True,
                     filtering_method='box', filtering_width=15, filtering_mode='center', filtering_stddev=4,
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

    filter_grid = None
    if filtering_active:
        filter_masked = resample_to_grid(
            {var_name_data: data_in},
            geo_x_in, geo_y_in, geo_x_grid, geo_y_grid, search_rad=resampling_max_distance / 5,
            min_neighbours=resampling_min_neighbours, neighbours=resampling_neighbours)
        filter_grid = filter_masked[var_name_data].data

    # check south-north array orientation
    geo_y_upper, geo_y_lower = geo_y_grid[0, 0], geo_y_grid[-1, 0]
    if geo_y_upper < geo_y_lower:
        geo_y_grid = np.flipud(geo_y_grid)
        data_grid = np.flipud(data_grid)
        if filtering_active:
            filter_grid = np.flipud(filter_grid)

    if filtering_active:
        if filtering_method == 'gaussian':
            # filter values and create a filtered grid
            obj_kernel = Gaussian2DKernel(
                x_stddev=filtering_stddev, mode=filtering_mode)
        elif filtering_method == 'box':
            obj_kernel = Box2DKernel(width=filtering_width, mode=Box2DKernel)
        else:
            raise NotImplementedError('Filtering method "' + filtering_method + '" not implemented yet')

        filter_kernel = convolve(filter_grid, obj_kernel)
        data_grid[filter_kernel == 0] = geo_no_data

    ''' debug
    plt.figure()
    plt.imshow(data_grid)
    plt.colorbar()
    plt.clim(0, 100)
    plt.show()
    '''

    # method to create data array
    data_out = create_darray_2d(
        data_grid, geo_x_grid[0, :], geo_y_grid[:, 0], name=var_name_data,
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    return data_out
# ----------------------------------------------------------------------------------------------------------------------
