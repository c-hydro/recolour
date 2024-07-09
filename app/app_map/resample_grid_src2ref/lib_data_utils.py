"""
Library Features:

Name:          lib_data_utils
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240708'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np

from copy import deepcopy

from scipy.interpolate import griddata
from pyresample.geometry import GridDefinition
from pyresample.kd_tree import resample_nearest
from repurpose.resample import resample_to_grid

from lib_data_tiff import organize_file_tiff, write_file_tiff
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)

# debugging
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to apply regrid method
def apply_method_resample(file_lons_src, file_lats_src, file_data_src, file_lons_ref, file_lats_ref,
                          resample_mode='interpolate',
                          resample_method='nearest', resample_max_distance=25000, resample_fill_value=0,
                          resample_min_neighbours=1, resample_neighbours=4, **kwargs):

    # iterate over the source data
    file_data_dst = []
    for i_data, data_src in enumerate(file_data_src):

        data_tmp = deepcopy(data_src)

        if resample_mode == 'resample_grid_to_grid':

            # resample grid source to grid reference
            data_dst = resample_grid_src2dst(
                file_lons_src, file_lats_src, data_tmp, file_lons_ref, file_lats_ref,
                interp_method=resample_method, interpolating_max_distance=resample_max_distance,
                interpolating_fill_value=resample_fill_value)

        elif resample_mode == 'resample_points_to_grid':

            # resample points source (starting from grid source) to grid reference
            file_lons_1d, file_lats_1d = file_lons_src.ravel(), file_lats_src.ravel()

            data_dst, _, _ = resample_points_src2dst(
                data_tmp, file_lons_1d, file_lats_1d, file_lons_ref, file_lats_ref,
                search_rad=resample_max_distance, fill_values=resample_fill_value,
                min_neighbours=resample_min_neighbours, neighbours=resample_neighbours)

        elif resample_mode == 'interpolate':
            data_dst = interpolate_src2dst(file_lons_src, file_lats_src, data_tmp, file_lons_ref, file_lats_ref,
                                           interp_method=resample_method, no_data=resample_fill_value)
        else:
            log_stream.error(' ===> Regrid method "' + resample_mode + '" is not expected by the algorithm')
            raise NotImplementedError('Method not implemented yet')

        ''' debug
        plt.figure()
        plt.imshow(data_src)
        plt.colorbar()
        plt.clim(0, 100)
        plt.figure()
        plt.imshow(data_dst)
        plt.colorbar()
        plt.clim(0, 100)
        plt.show()
        '''

        # store resampled data
        file_data_dst.append(data_dst)

    return file_data_dst
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to interpolate grid data to a reference dataset
def interpolate_src2dst(lons_in, lats_in, values_in, lons_out, lats_out, interp_method='nearest', no_data=np.nan):

    values_out = griddata((lons_in.ravel(), lats_in.ravel()), values_in.ravel(),
                         (lons_out, lats_out), method=interp_method, fill_value=no_data)

    return values_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample points to grid
def resample_points_src2dst(var_data_in_1d, var_geox_1d_in, var_geoy_1d_in,
                            var_geox_2d_out, var_geoy_2d_out,
                            var_mask_2d_out=None,
                            search_rad=50000, fill_values=np.nan,
                            min_neighbours=1, neighbours=4, **kwargs):

    values_obj = resample_to_grid(
        {'data': var_data_in_1d},
        var_geox_1d_in, var_geoy_1d_in, var_geox_2d_out, var_geoy_2d_out,
        search_rad=search_rad, fill_values=fill_values,
        min_neighbours=min_neighbours, neighbours=neighbours)

    var_data_out = values_obj['data']
    if var_mask_2d_out is not None:
        var_data_out[var_mask_2d_out == 0] = fill_values

    return var_data_out, var_geox_2d_out, var_geoy_2d_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample grid data to a reference dataset (version with pyresample)
def resample_grid_src2dst(lons_in, lats_in, values_in, lons_out, lats_out, interp_method='nearest',
                          interpolating_max_distance=10000, interpolating_fill_value=0):

    geo_grid_in = GridDefinition(lons=lons_in, lats=lats_in)
    geo_grid_ref = GridDefinition(lons=lons_out, lats=lats_out)

    if interp_method == 'nearest':
        values_out = resample_nearest(
            geo_grid_in, values_in, geo_grid_ref,
            radius_of_influence=interpolating_max_distance,
            fill_value=interpolating_fill_value)
    else:
        log_stream.error(' ===> Interpolation method is not allowed')
        raise NotImplementedError('Method not implemented yet')

    ''' start debug
    file_height_in, file_width_in, file_geo_transform_in, file_geo_epsg_in = organize_file_tiff(values_in, lons_in, lats_in)
    write_file_tiff('file_data_in.tif', values_in, file_width_in, file_height_in, file_geo_transform_in, file_geo_epsg_in)

    file_height_out, file_width_out, file_geo_transform_out, file_geo_epsg_out = organize_file_tiff(values_out, lons_out, lats_out)
    write_file_tiff('file_data_out.tif', values_out, file_width_out, file_height_out, file_geo_transform_out, file_geo_epsg_out)

    plt.figure()
    plt.imshow(values_in)
    plt.colorbar()
    plt.clim(0, 1)
    plt.figure()
    plt.imshow(values_out)
    plt.colorbar()
    plt.clim(0, 1)
    plt.show()
    '''

    return values_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def apply_method_mask(file_data_src, file_data_ref, mask_undefined_value=0, mask_fill_value=np.nan, **kwargs):

    # iterate over the source data
    file_data_dst = []
    for i_data, data_src in enumerate(file_data_src):
        data_tmp = deepcopy(data_src)
        data_dst = mask_src2ref(data_tmp, file_data_ref,
                                mask_value=mask_undefined_value, fill_value=mask_fill_value)
        # store masked data
        file_data_dst.append(data_dst)

        ''' debug
        plt.figure()
        plt.imshow(data_src)
        plt.colorbar()
        plt.clim(0, 100)
        plt.figure()
        plt.imshow(data_dst)
        plt.colorbar()
        plt.clim(0, 100)
        plt.show()
        '''

    return file_data_dst
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to mask source data over reference dataset
def mask_src2ref(values_data, values_ref, mask_value=0, fill_value=np.nan):
    values_data[values_ref == mask_value] = fill_value
    return values_data
# ----------------------------------------------------------------------------------------------------------------------
