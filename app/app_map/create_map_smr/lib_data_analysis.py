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

from lib_utils_geo import resample_points_to_grid

from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel
from pyresample.geometry import GridDefinition
from pyresample.kd_tree import resample_nearest, resample_gauss, resample_custom

from lib_utils_io import create_darray

from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)

# debugging
import matplotlib.pylab as plt

# logging
warnings.simplefilter("ignore", UserWarning)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select data
def select_data(dframe_obj_in,
                var_name='layer_1_scaled',
                min_value=0.0, max_value=100.0,
                no_data_value=np.nan, no_data_skip_rows=True):

    if no_data_value is None:
        no_data_value = np.nan

    dframe_obj_tmp = dframe_obj_in.copy()
    if min_value is not None:
        dframe_obj_tmp[dframe_obj_tmp[var_name] < min_value] = no_data_value
    if max_value is not None:
        dframe_obj_tmp[dframe_obj_tmp[var_name] > max_value] = no_data_value

    if no_data_skip_rows:
        if np.isnan(no_data_value):
            dframe_obj_out = dframe_obj_tmp.dropna()
        else:
            dframe_obj_out = dframe_obj_tmp.drop(dframe_obj_tmp[dframe_obj_tmp[var_name] == no_data_value].index)
    else:
        dframe_obj_out = dframe_obj_tmp.copy()

    return dframe_obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to adapt data
def adapt_data(obj_data_src, var_name_geo_x='longitude', var_name_geo_y='latitude', var_name_time='time'):

    # adapt geographical values x
    if var_name_geo_x in list(obj_data_src.keys()):
        geo_x_values_dst = obj_data_src[var_name_geo_x].values
        obj_data_src.pop(var_name_geo_x)
    else:
        alg_logger.error(' ===> Geographical variable "' + var_name_geo_x + '" is not available in the datasets')
        raise RuntimeError('Geographical variable is needed by the method')
    # adapt geographical values y
    if var_name_geo_y in list(obj_data_src.keys()):
        geo_y_values_dst = obj_data_src[var_name_geo_y].values
        obj_data_src.pop(var_name_geo_y)
    else:
        alg_logger.error(' ===> Geographical variable "' + var_name_geo_y + '" is not available in the datasets')
        raise RuntimeError('Geographical variable is needed by the method')
    # adapt values time
    if var_name_time in list(obj_data_src.keys()):
        time_values_dst = obj_data_src[var_name_time].values
        obj_data_src.pop(var_name_time)
    else:
        alg_logger.error(' ===> Time variable "' + var_name_time + '" is not available in the datasets')
        raise RuntimeError('Time variable is needed by the method')

    # adapt variable data
    obj_data_dst = {}
    for var_name, var_data_src in obj_data_src.items():
        var_values_src = var_data_src.values
        obj_data_dst[var_name] = var_values_src

    return obj_data_dst, geo_x_values_dst, geo_y_values_dst, time_values_dst
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data
def resample_data(obj_data_src, geo_x_values_src, geo_y_values_src, geo_x_values_dst, geo_y_values_dst, **kwargs):

    # iterate over variable(s)
    obj_data_dst = {}
    for var_name, var_values_src in obj_data_src.items():
        if var_values_src is not None:

            if var_name in kwargs:
                var_settings = kwargs.get(var_name)
            else:
                var_settings = {}

            var_values_dst, _, _ = resample_points_to_grid(
                var_values_src, geo_x_values_src, geo_y_values_src, geo_x_values_dst, geo_y_values_dst, **var_settings)

            var_values_dst = np.flipud(var_values_dst)

            ''' debug
            plt.figure()
            plt.imshow(var_values_dst)
            plt.colorbar()
            plt.show()
            plt.figure()
            plt.imshow(geo_y_values_dst)
            plt.colorbar()
            plt.show()
            '''

            obj_data_dst[var_name] = var_values_dst
        else:
            alg_logger.warning(' ===> Data "' + var_name + '" is not available in the datasets')
            obj_data_dst[var_name] = None

    return obj_data_dst
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter data
def filter_data(obj_data_src, filter_type='gauss', gauss_st_dev=3, gauss_mode='center', box_width=3, **kwargs):

    # iterate over variable(s)
    obj_data_dst = {}
    for var_name, var_values_src in obj_data_src.items():

        # box values and create a filtered grid
        var_values_tmp = deepcopy(var_values_src)

        # check filter is active
        if filter_type is not None:
            # check filter type
            if filter_type == 'gauss':

                obj_kernel_gauss = Gaussian2DKernel(x_stddev=gauss_st_dev, mode=gauss_mode)
                var_values_dst = convolve(var_values_tmp, obj_kernel_gauss)

            elif filter_type == 'box':

                obj_kernel_box = Box2DKernel(box_width)
                var_values_dst = convolve(var_values_tmp, obj_kernel_box)

            else:
                alg_logger.error(' ===> Filter type "' + filter_type + '" is not available')
                raise NotImplemented('Case not implemented yet')

        else:
            # skip filter
            var_values_dst = deepcopy(var_values_tmp)

        # store filtered data
        obj_data_dst[var_name] = var_values_dst

    return obj_data_dst

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to mask data
def mask_data(obj_data_src, var_values_ref, no_data_ref=-9999):

    # iterate over variable(s)
    obj_data_dst = {}
    for var_name, var_values_src in obj_data_src.items():
        if var_values_src is not None:
            var_values_src[var_values_ref == no_data_ref] = np.nan

            ''' debug
            plt.figure()
            plt.imshow(var_values_src)
            plt.colorbar()
            plt.show()
            '''
            obj_data_dst[var_name] = var_values_src
        else:
            alg_logger.warning(' ===> Data "' + var_name + '" is not available in the datasets')
            obj_data_dst[var_name] = None

    return obj_data_dst

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

        obj_da_out = create_darray(
            values_out_resampled, geo_x_grid_ref[0, :], geo_y_grid_ref[:, 0], name=var_name_data,
            coord_name_x=coord_name_x, coord_name_y=coord_name_y,
            dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    else:
        obj_da_out = deepcopy(obj_da_in)

    return obj_da_out
# ----------------------------------------------------------------------------------------------------------------------
