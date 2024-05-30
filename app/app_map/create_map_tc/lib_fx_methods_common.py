"""
Library Features:

Name:          lib_fx_methods_common
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240507'
Version:       '1.6.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import numpy as np
from scipy.interpolate import griddata

from copy import deepcopy

from pyresample.geometry import GridDefinition
from repurpose.resample import resample_to_grid
from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel

from lib_info_args import logger_name
from lib_utils_io import create_darray_2d
from lib_utils_generic import pad_list

# set logger warnings
warnings.filterwarnings("ignore")
# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
from lib_utils_plot import plot_data_2d
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to update data using masked area
def update_data_by_masked_area(
        obj_in, geo_obj, mask_obj, removed_obj, cell_obj,
        var_mapping_in=None, var_mapping_out=None,
        var_name_ref='soil_moisture_ref', var_name_k1='soil_moisture_k1', var_name_k2='soil_moisture_k2',
        var_name_geo_x='longitude', var_name_geo_y='latitude',
        coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
        max_distance=25000, min_neighbours=1, neighbours=8,
        gauss_stddev_ref=3, gauss_mode_ref='center',
        gauss_stddev_k1=4, gauss_mode_k1='center', box_width=2, gauss_stddev_k2=4, gauss_mode_k2='center',
        value_min=0.1, value_max=0.95,
        debug=False, **kwargs):

    # check obj_in
    for var_key, var_name in var_mapping_in.items():
        if var_name not in list(obj_in.keys()):
            alg_logger.error(' ===> Variable "' + var_name + '" is not available in the source DataFrame')
            raise RuntimeError('Variable is needed by the algorithm. Check your settings and procedures')

    # create obj_out
    obj_out = deepcopy(obj_in)
    
    # check obj_out dataset (and exit if is NoneType)
    if cell_obj is None or (cell_obj[var_name_geo_x] is None or cell_obj[var_name_geo_y] is None):
        alg_logger.warning(' ===> Resampled datasets for variable '
                           '"soil_moisture_resampled" is not available. Datasets will be defined by NoneType ')
        return obj_out

    obj_out.pop('soil_moisture_resampled', None)
    obj_out.pop('soil_moisture_filtered', None)

    # geo data reference
    geo_x_1d = geo_obj[var_name_geo_x].values
    geo_y_1d = geo_obj[var_name_geo_y].values
    geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
    geo_values = geo_obj.values
    # mask data reference
    mask_values = mask_obj.values
    removed_values = removed_obj.values
    # cell data reference
    var_geo_x = cell_obj[var_name_geo_x].values
    var_geo_y = cell_obj[var_name_geo_y].values
    var_cell_reference = cell_obj[var_name_ref].values
    var_cell_k1 = cell_obj[var_name_k1].values
    var_cell_k2 = cell_obj[var_name_k2].values

    # resample data reference
    var_obj_reference = resample_to_grid(
        {'data': var_cell_reference},
        var_geo_x, var_geo_y, geo_x_2d, geo_y_2d,
        search_rad=max_distance, fill_values=np.nan,
        min_neighbours=min_neighbours, neighbours=neighbours)
    var_data_reference = var_obj_reference['data']
    var_data_reference[geo_values == 0] = np.nan
    # resample data k1
    var_obj_k1 = resample_to_grid(
        {'data': var_cell_k1},
        var_geo_x, var_geo_y, geo_x_2d, geo_y_2d,
        search_rad=max_distance, fill_values=np.nan,
        min_neighbours=min_neighbours, neighbours=neighbours)
    var_data_k1 = var_obj_k1['data']
    var_data_k1[geo_values == 0] = np.nan

    # filter values and create a filtered grid
    obj_kernel = Gaussian2DKernel(
        x_stddev=gauss_stddev_ref, mode=gauss_mode_ref)
    var_data_gauss_ref = convolve(var_data_reference, obj_kernel)
    # filter values and create a filtered grid
    obj_kernel_k1 = Gaussian2DKernel(
        x_stddev=gauss_stddev_k1, mode=gauss_mode_k1)
    var_data_gauss_k1 = convolve(var_data_k1, obj_kernel)

    # get data resampled and filtered
    var_data_resampled = obj_in['soil_moisture_resampled'].values
    var_data_filtered = obj_in['soil_moisture_filtered'].values

    # apply reference data to resampled data
    alg_logger.info(' ------> Mask resampled data over masked areas ... ')
    tmp_data_step_1 = deepcopy(var_data_resampled)
    tmp_data_step_1[mask_values > 0] = -1
    idx_nd_1 = np.argwhere(tmp_data_step_1 == -1)
    tmp_data_step_1[idx_nd_1[:, 0], idx_nd_1[:, 1]] = var_data_reference[idx_nd_1[:, 0], idx_nd_1[:, 1]]
    alg_logger.info(' ------> Mask resampled data over masked areas ... DONE')

    # apply filtered data to resampled data
    alg_logger.info(' ------> Apply filtered data to resampled data over masked areas ... ')
    tmp_data_step_2 = deepcopy(tmp_data_step_1)
    tmp_data_step_2[np.isnan(var_data_reference)] = -1
    idx_nd_2 = np.argwhere(tmp_data_step_2 == -1)
    tmp_data_step_2[idx_nd_2[:, 0], idx_nd_2[:, 1]] = var_data_gauss_ref[idx_nd_2[:, 0], idx_nd_2[:, 1]]
    alg_logger.info(' ------> Apply filtered data to resampled data over masked areas ... DONE')

    # apply domain mask data to resampled data
    alg_logger.info(' ------> Mask resampled data over terrain areas ... ')
    tmp_data_step_3 = deepcopy(tmp_data_step_2)
    tmp_data_step_3[geo_values < 1] = np.nan
    tmp_data_step_3[removed_values > 0] = np.nan
    alg_logger.info(' ------> Mask resampled data over terrain areas ... DONE')

    # filtered data
    tmp_data_filtered = deepcopy(var_data_filtered)
    tmp_data_filtered[idx_nd_1[:, 0], idx_nd_1[:, 1]] = var_data_gauss_ref[idx_nd_1[:, 0], idx_nd_1[:, 1]]
    # filter values and create a filtered grid
    obj_kernel_filtered = Gaussian2DKernel(
        x_stddev=3, mode='center')
    tmp_data_filtered = convolve(tmp_data_filtered, obj_kernel_filtered)
    tmp_data_filtered[geo_values < 1] = np.nan
    tmp_data_filtered[removed_values > 0] = np.nan

    # apply value min and max to resampled data
    alg_logger.info(' ------> Fill resampled data out of the variable limits ... ')
    tmp_data_step_4 = deepcopy(tmp_data_step_3)
    tmp_data_step_4[tmp_data_step_4 < value_min] = np.nan
    tmp_data_step_4[tmp_data_step_4 > value_max] = np.nan
    idx_nd_4_nan = np.argwhere(np.isnan(tmp_data_step_4))
    tmp_data_step_4[idx_nd_4_nan[:, 0], idx_nd_4_nan[:, 1]] = tmp_data_filtered[idx_nd_4_nan[:, 0], idx_nd_4_nan[:, 1]]
    alg_logger.info(' ------> Fill resampled data out of the variable limits ... DONE')

    # box values and create a filtered grid
    tmp_data_step_5 = deepcopy(tmp_data_step_4)
    obj_kernel_box = Box2DKernel(box_width)
    tmp_data_box = convolve(tmp_data_step_5, obj_kernel_box)
    tmp_data_box[geo_values < 1] = np.nan
    tmp_data_box[removed_values > 0] = np.nan

    # active debug
    if debug:
        plot_data_2d(removed_values)
        plot_data_2d(tmp_data_step_3)
        plot_data_2d(tmp_data_step_4)
        plot_data_2d(tmp_data_filtered)
        plot_data_2d(geo_values)
        plot_data_2d(mask_values)
        plot_data_2d(var_data_reference)
        plot_data_2d(var_data_gauss_ref)
        plot_data_2d(var_data_gauss_k1)
        plot_data_2d(tmp_data_step_1)
        plot_data_2d(tmp_data_step_2)
        plot_data_2d(tmp_data_step_3)
        plot_data_2d(tmp_data_step_4)
        plot_data_2d(tmp_data_box)
        plot_data_2d(tmp_data_filtered)

    # method to create data array resampled
    var_da_out_resampled = create_darray_2d(
        tmp_data_box, geo_x_1d, geo_y_1d, name='soil_moisture_resampled',
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    # method to create data array filtered
    var_da_out_filtered = create_darray_2d(
        tmp_data_filtered, geo_x_1d, geo_y_1d, name='soil_moisture_filtered',
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    # update obj_out
    obj_out['soil_moisture_resampled'] = var_da_out_resampled
    obj_out['soil_moisture_filtered'] = var_da_out_filtered

    # check obj_out
    for var_key, var_name in var_mapping_out.items():
        if var_name not in list(obj_out.keys()):
            alg_logger.error(' ===> Variable "' + var_name + '" is not available in the destination DataFrame')
            raise RuntimeError('Variable is needed by the algorithm. Check your settings and procedures')

    return obj_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to interpolate points to grid using gridded data
def interp_data(lons_in, lats_in, values_in, lons_out, lats_out, nodata=-9999, interp_method='nearest'):

    values_out = griddata((lons_in.ravel(), lats_in.ravel()), values_in.ravel(),
                          (lons_out, lats_out), method=interp_method,
                          fill_value=nodata)
    return values_out
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
# method to get parameter value
def resample_parameter(var_name, var_parameters=None, var_value_default=None):

    var_value_param = deepcopy(var_value_default)
    if var_parameters is None:
        pass
    else:
        if isinstance(var_parameters, dict):
            if var_name in list(var_parameters.keys()):
                var_value_param = var_parameters[var_name]
            else:
                logger_name.warning(' ===> Variable "' + var_name +
                                    '" is not available in the removing artifacts dictionary')
        else:
            logger_name.error(' ===> Removing artifacts must be defined by dictionary')
            raise NotImplemented('Case not implemented yet')

    return var_value_param
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data (proxy
def resample_data(var_mapping_in, var_mapping_out,
                  resampling_max_distance=None, resampling_min_neighbours=None, resampling_neighbours=None,
                  resampling_remove_artifacts=None, resampling_datasets_artifacts=None,  **kwargs):

    # iterate over variable(s)
    var_obj = {}
    for (var_key_in, var_name_in), (var_key_out, var_name_out) in zip(var_mapping_in.items(), var_mapping_out.items()):

        # info variable resampling start
        alg_logger.info(' ------> Variable "' + var_name_in + '" ... ')

        # get resampling max distance
        value_max_distance = resample_parameter(var_name=var_key_in, var_parameters=resampling_max_distance,
                                                var_value_default=1000)
        # get resampling min neighbours
        value_min_neighbours = resample_parameter(var_name=var_key_in, var_parameters=resampling_min_neighbours,
                                                  var_value_default=1)
        # get resampling neighbours
        value_neighbours = resample_parameter(var_name=var_key_in, var_parameters=resampling_neighbours,
                                              var_value_default=8)
        # get resampling remove artifacts
        value_remove_artifacts = resample_parameter(var_name=var_key_in, var_parameters=resampling_remove_artifacts,
                                                    var_value_default=False)

        value_datasets_artifacts = resample_parameter(var_name=var_key_in, var_parameters=resampling_datasets_artifacts,
                                                      var_value_default=None)

        # resample data
        if value_remove_artifacts:
            var_da = resample_data_remove_artefacts(
                resampling_max_distance=value_max_distance,
                resampling_neighbours=value_neighbours, resampling_min_neighbours=value_min_neighbours,
                var_name_data_in=var_name_in, var_name_data_out=var_name_out,
                resampling_datasets_artifacts=value_datasets_artifacts,
                **kwargs)
        else:
            var_da = resample_data_generic(
                resampling_max_distance=value_max_distance,
                resampling_neighbours=value_neighbours, resampling_min_neighbours=value_min_neighbours,
                var_name_data_in=var_name_in, var_name_data_out=var_name_out,
                **kwargs)

        # store data array
        var_obj[var_name_out] = var_da

        alg_logger.info(' ------> Variable "' + var_name_in + '" ... DONE')

    return var_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data (remove artefacts)
def resample_data_remove_artefacts(
        var_data_obj, var_name_data_in, var_name_data_out,
        geo_mask_out, geo_x_out, geo_y_out,
        var_name_geo_x='longitude', var_name_geo_y='latitude',
        coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
        resampling_max_distance=1000, resampling_neighbours=4, resampling_min_neighbours=1,
        resampling_datasets_artifacts=None, resampling_distance_artifacts=25000,
        resampling_fill_value=np.nan,
        resampling_gauss_stddev_partial=1, resampling_gauss_stddev_global=8,
        resampling_gauss_mode_partial='center', resampling_gauss_mode_global='center',
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

    # iterate over variable(s)
    var_cell_list = np.unique(var_data_obj['cell'].values).tolist()

    # iterate over cell(s)
    alg_logger.info(' -------> Resample cell(s) ... ')
    values_grid_merged, values_grid_raw, idx_finite_merged = None, None, None
    for var_cell_value in var_cell_list:

        # info cell resampling start
        alg_logger.info(' --------> Cell Group "' + str(var_cell_value) + '" ... ')

        # select cell object
        obj_cell = var_data_obj.loc[var_data_obj['cell'] == var_cell_value]

        # check variable availability
        if var_name_data_in not in list(obj_cell.columns):
            alg_logger.error(' ===> Variable "' + var_name_data_in + '" is not available in the source DataFrame')
            raise RuntimeError('Variable is needed by the algorithm. Check your settings and procedures')

        # get data
        values_arr_in = obj_cell[var_name_data_in].values
        geo_x_arr_in, geo_y_arr_in = obj_cell[var_name_geo_x].values, obj_cell[var_name_geo_y].values
        # get nan(s)
        nan_idx_tmp = np.argwhere(np.isnan(values_arr_in))[:, 0]

        # extent to grid without nan(s)
        values_arr_tmp = np.delete(values_arr_in, nan_idx_tmp)
        geo_x_arr_tmp = np.delete(geo_x_arr_in, nan_idx_tmp)
        geo_y_arr_tmp = np.delete(geo_y_arr_in, nan_idx_tmp)

        # check finite values
        flag_any_finite = np.any(values_arr_tmp)
        if flag_any_finite:

            # initialize values grid merged
            if values_grid_merged is None:
                values_grid_merged = np.zeros(shape=(geo_y_grid_out.shape[0], geo_x_grid_out.shape[1]))
                values_grid_merged[:] = np.nan

            # initialize values grid raw
            if values_grid_raw is None:
                values_grid_raw = np.zeros(shape=(geo_y_grid_out.shape[0], geo_x_grid_out.shape[1]))
                values_grid_raw[:] = np.nan

            # resample point 2 map
            values_obj = resample_to_grid(
                {var_name_data_in: values_arr_tmp},
                geo_x_arr_tmp, geo_y_arr_tmp, geo_x_grid_out, geo_y_grid_out,
                search_rad=resampling_max_distance, fill_values=resampling_fill_value,
                min_neighbours=resampling_min_neighbours, neighbours=resampling_neighbours)

            # get values data
            if resampling_fill_value is None:
                values_grid_out = values_obj[var_name_data_in].data
            else:
                values_grid_out = values_obj[var_name_data_in]

            # raw values
            values_grid_raw = np.where(
                np.isnan(values_grid_raw), deepcopy(values_grid_out), deepcopy(values_grid_raw))

            # mask out values outside the mask
            values_grid_out[geo_mask_out == 0] = np.nan

            # filter values
            obj_kernel = Gaussian2DKernel(
                x_stddev=resampling_gauss_stddev_partial, mode=resampling_gauss_mode_partial)
            values_grid_tmp = convolve(values_grid_out, obj_kernel)

            # mask filtered values outside the mask
            values_grid_tmp[geo_mask_out == 0] = np.nan

            # compute finite indexes
            idx_finite_combined = np.argwhere((np.isfinite(values_grid_merged) & np.isfinite(values_grid_tmp)))
            if idx_finite_merged is None:
                idx_finite_merged = deepcopy(idx_finite_combined)
            else:
                idx_finite_merged = np.concatenate((idx_finite_merged, idx_finite_combined), axis=0)

            # merge values
            values_grid_merged = np.where(
                np.isnan(values_grid_merged), values_grid_out, deepcopy(values_grid_merged))

            # mask values with combined finite indexes
            values_grid_merged[idx_finite_combined[:, 0], idx_finite_combined[:, 1]] = np.nan

            ''' debug
            plt.figure(); plt.imshow(values_grid_out); plt.colorbar(); plt.clim([0, 1])
            plt.figure(); plt.imshow(values_grid_tmp); plt.colorbar(); plt.clim([0, 1])
            plt.figure(); plt.imshow(values_grid_merged); plt.colorbar(); plt.clim([0, 1])
            plt.show()
            '''

            # info cell resampling start
            alg_logger.info(' --------> Cell Group "' + str(var_cell_value) + '" ... DONE')

        else:
            alg_logger.warning(' ===> Resample method was skipped; all values are defined by Nan(s)')
            alg_logger.info(' --------> Cell Group "' + str(var_cell_value) + '" ... SKIPPED')

    # info cell resampling end
    alg_logger.info(' -------> Resample cell(s) ... DONE')

    values_grid_merged[geo_mask_out == 0] = np.nan
    values_grid_raw[geo_mask_out == 0] = np.nan

    ''' debug
    plt.figure(); plt.imshow(values_grid_merged); plt.colorbar(); plt.clim([0, 1])
    plt.figure(); plt.imshow(values_grid_raw); plt.colorbar(); plt.clim([0, 1])
    plt.show()
    '''

    # check grid merged availability
    alg_logger.info(' -------> Organize cell(s) ... ')
    if values_grid_merged is not None:

        # info remove artefacts start
        alg_logger.info(' --------> Remove artefacts ... ')

        # mask values with merged indexes (to remove spots)
        values_grid_merged[idx_finite_merged[:, 0], idx_finite_merged[:, 1]] = np.nan
        # values_grid_raw[idx_finite_merged[:, 0], idx_finite_merged[:, 1]] = np.nan

        values_grid_composite = np.zeros(shape=(geo_y_grid_out.shape[0], geo_x_grid_out.shape[1]))
        values_grid_composite[:] = np.nan

        # filter values and create a filtered grid
        obj_kernel = Gaussian2DKernel(
            x_stddev=resampling_gauss_stddev_global, mode=resampling_gauss_mode_global)
        # replace nan(s) with filtered values
        # values_grid_filtered = interpolate_replace_nans(values_grid_merged, obj_kernel)
        # values_grid_filtered = convolve(values_grid_merged, obj_kernel)
        values_grid_filtered = convolve(values_grid_raw, obj_kernel)

        # filled empty domain cell with filtered values
        values_grid_filled = np.where(np.isnan(values_grid_merged),
                                      values_grid_filtered, deepcopy(values_grid_merged))

        values_grid_composite[np.isfinite(values_grid_merged)] = values_grid_merged[np.isfinite(values_grid_merged)]
        values_grid_composite[np.isnan(values_grid_merged)] = values_grid_filled[np.isnan(values_grid_merged)]

        # mask filled values outside the mask
        values_grid_filled[geo_mask_out == 0] = np.nan
        values_grid_composite[geo_mask_out == 0] = np.nan

        # info remove artefacts end
        alg_logger.info(' --------> Remove artifacts ... DONE')

        # method to create data array
        var_da_out = create_darray_2d(
            values_grid_merged, geo_x_grid_out[0, :], geo_y_grid_out[:, 0], name=var_name_data_out,
            coord_name_x=coord_name_x, coord_name_y=coord_name_y,
            dim_name_x=dim_name_x, dim_name_y=dim_name_y)

        var_da_out.attrs = {'index': idx_finite_merged}

        # info cell merging end
        alg_logger.info(' -------> Organize cell(s) ... DONE')
    else:
        # info cell merging end
        alg_logger.info(' -------> Organize cell(s) ... SKIPPED. All values are defined by Nan(s)')
        var_da_out = None

    return var_da_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data (generic)
def resample_data_generic(
        var_data_obj, var_name_data_in, var_name_data_out,
        geo_mask_out, geo_x_out, geo_y_out,
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

    # geo_grid_out = GridDefinition(lons=geo_x_grid_out, lats=geo_y_grid_out)

    # check object format
    if var_data_obj is not None:

        # check variable availability
        if var_name_data_in not in list(var_data_obj.columns):
            alg_logger.error(' ===> Variable "' + var_name_data_in + '" is not available in the source DataFrame')
            raise RuntimeError('Variable is needed by the algorithm. Check your settings and procedures')

        # get data
        values_arr_in = var_data_obj[var_name_data_in].values
        geo_x_arr_in, geo_y_arr_in = var_data_obj[var_name_geo_x].values, var_data_obj[var_name_geo_y].values
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
                {var_name_data_in: values_arr_tmp},
                geo_x_arr_tmp, geo_y_arr_tmp, geo_x_grid_out, geo_y_grid_out,
                search_rad=resampling_max_distance, fill_values=resampling_fill_value,
                min_neighbours=resampling_min_neighbours, neighbours=resampling_neighbours)

            if resampling_fill_value is None:
                values_grid_out = values_obj[var_name_data_in].data
            else:
                values_grid_out = values_obj[var_name_data_in]

            # mask values outside the mask
            values_grid_out[geo_mask_out == 0] = np.nan

            ''' debug
            plt.figure(); plt.imshow(values_grid_out); plt.colorbar(); plt.clim([0, 1])
            plt.show()
            '''

        else:
            values_grid_out = None
            alg_logger.warning(' ===> Resample method was skipped; all values are defined by Nan(s)')

        # method to create data array
        if values_grid_out is not None:
            var_da_out = create_darray_2d(
                values_grid_out, geo_x_grid_out[0, :], geo_y_grid_out[:, 0], name=var_name_data_out,
                coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                dim_name_x=dim_name_x, dim_name_y=dim_name_y)
        else:
            var_da_out = None
    else:
        alg_logger.warning(' ===> Resampled datasets for variable "' + var_name_data_in +
                           '" is not available. Datasets will be defined by NoneType')
        var_da_out = None

    return var_da_out

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

