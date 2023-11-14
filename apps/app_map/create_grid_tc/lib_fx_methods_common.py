"""
Library Features:

Name:          lib_fx_methods_common
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import os
import numpy as np

from copy import deepcopy

from pyresample.geometry import GridDefinition
from pyresample.kd_tree import resample_nearest, resample_gauss, resample_custom
from repurpose.resample import resample_to_grid

from astropy.convolution import convolve, Gaussian2DKernel

from lib_info_args import logger_name
from lib_utils_io import create_darray_2d
from lib_utils_generic import pad_list

# set logger warnings
warnings.filterwarnings("ignore")
# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to crop data
def crop_data(obj_in,
              geo_x_left, geo_x_right, geo_y_lower, geo_y_upper,
              var_name_geo_x='longitude', var_name_geo_y='latitude'):

    obj_out, mask_lon, mask_lat = {}, None, None
    for var_key, var_da_in in obj_in.items():

        if (mask_lon is None) or (mask_lat is None):
            values_lon, values_lat = var_da_in[var_name_geo_x], var_da_in[var_name_geo_y]
            mask_lon = (values_lon >= geo_x_left) & (values_lon <= geo_x_right)
            mask_lat = (values_lat >= geo_y_lower) & (values_lat <= geo_y_upper)

        var_da_out = var_da_in.where(mask_lon & mask_lat, drop=True)

        ''' debug
        plt.figure(); plt.imshow(var_da_in.values); plt.colorbar()
        plt.figure(); plt.imshow(var_da_out.values); plt.colorbar()
        plt.show()
        '''

        obj_out[var_key] = var_da_out

    return obj_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data
def resample_data(
        obj_in,
        geo_x_out, geo_y_out,
        var_mapping_in=None, var_mapping_out=None,
        var_name_geo_x='longitude', var_name_geo_y='latitude',
        coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
        resampling_max_distance=18000, resampling_neighbours=8, resampling_min_neighbours=1,
        resampling_fill_value=np.nan,
        **kwargs):

    # define geographical information
    if (geo_x_out.ndim == 1) and (geo_y_out.ndim == 1):
        geo_x_grid_out, geo_y_grid_out = np.meshgrid(geo_x_out, geo_y_out)
    elif (geo_x_out.ndim == 2) and (geo_y_out.ndim == 2):
        geo_x_grid_out, geo_y_grid_out = deepcopy(geo_x_out), deepcopy(geo_y_out)
    else:
        alg_logger.error(' ===> Geographical information format is not supported')
        raise NotImplemented('Case not implemented yet')

    # check south-north array orientation
    geo_y_upper_out, geo_y_lower_out = geo_y_grid_out[0, 0], geo_y_grid_out[-1, 0]
    if geo_y_upper_out < geo_y_lower_out:
        geo_y_grid_out = np.flipud(geo_y_grid_out)

    geo_grid_out = GridDefinition(lons=geo_x_grid_out, lats=geo_y_grid_out)

    # iterate over varieble(s)
    obj_out, mask_lon, mask_lat = {}, None, None
    for var_key_in, var_key_out in zip(var_mapping_in.values(), var_mapping_out.values()):

        # check variable availability
        if var_key_in not in list(obj_in.columns):
            alg_logger.error(' ===> Variable "' + var_key_in + '" is not available in the source DataFrame')
            raise RuntimeError('Variable is needed by the algorithm. Check your settings and procedures')

        # get data
        values_arr_in = obj_in[var_key_in].values
        geo_x_arr_in, geo_y_arr_in = obj_in[var_name_geo_x].values, obj_in[var_name_geo_y].values
        # get nan(s)
        nan_idx_tmp = np.argwhere(np.isnan(values_arr_in))[:, 0]

        # extent to grid without nan(s)
        values_arr_tmp = np.delete(values_arr_in, nan_idx_tmp)
        geo_x_arr_tmp = np.delete(geo_x_arr_in, nan_idx_tmp)
        geo_y_arr_tmp = np.delete(geo_y_arr_in, nan_idx_tmp)

        # check finite values
        flag_any_finite = np.any(values_arr_tmp)
        if flag_any_finite:

            # resample point 2 map
            values_obj = resample_to_grid(
                {var_key_in: values_arr_tmp},
                geo_x_arr_tmp, geo_y_arr_tmp, geo_x_grid_out, geo_y_grid_out,
                search_rad=resampling_max_distance, fill_values=resampling_fill_value,
                min_neighbours=resampling_min_neighbours, neighbours=resampling_neighbours)

            if resampling_fill_value is None:
                values_grid_out = values_obj[var_key_in].data
            else:
                values_grid_out = values_obj[var_key_in]

            ''' # debug
            plt.figure()
            plt.imshow(values_grid_out); plt.colorbar(); plt.clim([0, 3]); plt.show()
            '''

        else:
            values_grid_out = None
            alg_logger.warning(' ===> Resample method was skipped; all values are defined by Nan(s)')

        # method to create data array
        if values_grid_out is not None:
            var_da_out = create_darray_2d(
                values_grid_out, geo_x_grid_out[0, :], geo_y_grid_out[:, 0], name=var_key_out,
                coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                dim_name_x=dim_name_x, dim_name_y=dim_name_y)
        else:
            var_da_out = None

        obj_out[var_key_out] = var_da_out

    return obj_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter data
def filter_data(obj_in,
                filter_type='gauss', filter_pixels_std=4, filter_mode='center',
                var_mapping_in=None, var_mapping_out=None,
                var_name_geo_x='longitude', var_name_geo_y='latitude',
                coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
                ):

    obj_out = deepcopy(obj_in)
    for var_id, (var_key_in, var_key_out) in enumerate(zip(var_mapping_in.values(), var_mapping_out.values())):

        # check variable availability
        if var_key_in in list(obj_in.keys()):

            # check variable dataset (not defined by NoneType)
            if obj_in[var_key_in] is not None:

                obj_da = obj_in[var_key_in]

                values_tmp = deepcopy(obj_da.values)
                geo_x_arr, geo_y_arr = obj_da[var_name_geo_x].values, obj_da[var_name_geo_y].values

                if filter_type == 'gauss':
                    obj_kernel = Gaussian2DKernel(x_stddev=filter_pixels_std, mode=filter_mode)
                    values_out = convolve(values_tmp, obj_kernel)
                else:
                    alg_logger.error(' ===> The filter mode "' + filter_type + '" is not supported')
                    raise NotImplemented('Case not implemented yet')

                # method to create data array
                var_da_out = create_darray_2d(
                    values_out, geo_x_arr, geo_y_arr, name=var_key_out,
                    coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                    dim_name_x=dim_name_x, dim_name_y=dim_name_y)

                ''' debug
                plt.figure(); plt.imshow(obj_in[var_key_in].values); plt.colorbar()
                plt.figure(); plt.imshow(var_da_out.values); plt.colorbar()
                plt.show()
                '''
            else:
                alg_logger.warning(' ===> Filtered datasets for variable "' +
                                   var_key_in + '" is not available. Datasets will be defined by NoneType')
                var_da_out = None

        else:
            alg_logger.warning(' ===> Variable "' + var_key_in + '" is not available in the source DataFrame')
            var_da_out = None

        obj_out[var_key_out] = var_da_out

    return obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to mask data
def mask_data(obj_in, reference_da,
              var_mapping_in=None, var_mapping_out=None,
              mask_value_min=0, mask_value_max=None, mask_no_data=np.nan,
              var_name_geo_x='longitude', var_name_geo_y='latitude',
              coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
              ):

    if not isinstance(mask_value_min, list):
        mask_value_min = [mask_value_min]
    if not isinstance(mask_value_max, list):
        mask_value_max = [mask_value_max]
    if not isinstance(mask_no_data, list):
        mask_no_data = [mask_no_data]

    var_keys_in, var_keys_out = list(var_mapping_in.keys()), list(var_mapping_out.keys())
    assert var_keys_in == var_keys_out, 'Input and output variable(s) must have the same length'

    mask_value_min = pad_list(var_keys_in, mask_value_min, mask_value_min[0])
    mask_value_max = pad_list(var_keys_in, mask_value_max, mask_value_max[0])
    mask_no_data = pad_list(var_keys_in, mask_no_data, mask_no_data[0])

    obj_out = deepcopy(obj_in)
    for var_id, (var_key_in, var_key_out) in enumerate(zip(var_mapping_in.values(), var_mapping_out.values())):

        # check variable availability
        if var_key_in in list(obj_in.keys()):

            # check variable dataset (not defined by NoneType)
            if obj_in[var_key_in] is not None:

                obj_da = obj_in[var_key_in]

                values_out = deepcopy(obj_da.values)
                geo_x_arr, geo_y_arr = obj_da[var_name_geo_x].values, obj_da[var_name_geo_y].values
                mask_values = reference_da.values

                value_min, value_max, value_no_data = mask_value_min[var_id], mask_value_max[var_id], mask_no_data[var_id]

                if value_min is not None:
                    values_out[mask_values < value_min] = value_no_data
                if value_max is not None:
                    values_out[mask_values > value_max] = value_no_data

                # method to create data array
                var_da_out = create_darray_2d(
                    values_out, geo_x_arr, geo_y_arr, name=var_key_out,
                    coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                    dim_name_x=dim_name_x, dim_name_y=dim_name_y)

                ''' debug
                plt.figure(); plt.imshow(obj_in[var_key_in].values); plt.colorbar()
                plt.figure(); plt.imshow(var_da_out.values); plt.colorbar()
                plt.show()
                '''
            else:
                alg_logger.warning(' ===> Masked datasets for variable "' +
                                   var_key_in + '" is not available. Datasets will be defined by NoneType')
                var_da_out = None

        else:
            alg_logger.warning(' ===> Variable "' + var_key_in + '" is not available in the source DataFrame')
            var_da_out = None

        obj_out[var_key_out] = var_da_out

    return obj_out

# ----------------------------------------------------------------------------------------------------------------------

