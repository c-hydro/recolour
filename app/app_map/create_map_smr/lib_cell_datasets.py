"""
Library Features:

Name:          lib_cell_datasets
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy

import pygeogrids.netcdf as grid2nc

from lib_cell_tools import extract_cell_data, extract_cell_attrs, extract_registry_info

from lib_utils_generic import remove_key_with_null_value, check_key_between_dicts
from lib_data_io_nc import read_file_cell, write_file_cell
from lib_fx_scale import compute_data_scaling
from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# message(s) suppressed in the console
pd.options.mode.chained_assignment = None
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get datasets cells
def get_datasets_cell(
        file_name, cell_name,
        list_variable_data_in=None, list_variable_data_out=None,
        list_variable_registry_in=None, list_variable_registry_out=None,
        index_name_data_in='time', index_name_registry_in='location_id', index_name_cell_in='cell',
        index_name_row_size_in='row_size',
        index_name_data_out='time', index_name_registry_out='location_id', index_name_cell_out='cell'):

    # if row_index is the first variable in the list, move it to the end
    if index_name_row_size_in == list_variable_registry_in[0]:
        row_size_tmp_in, row_size_tmp_out = list_variable_registry_in[0], list_variable_registry_out[0]
        list_variable_registry_in = list_variable_registry_in[1:] + [row_size_tmp_in]
        list_variable_registry_out = list_variable_registry_out[1:] + [row_size_tmp_out]

    # common list variable in/out
    list_variable_in = list_variable_data_in + list_variable_registry_in
    list_variable_out = list_variable_data_out + list_variable_registry_out

    # check file name availability
    collections_obj, registry_obj = None, None
    if os.path.exists(file_name):

        # get datasets in
        cell_datasets_in, attrs_datasets_in = read_file_cell(file_name, selected_variables=list_variable_in)

        # get datasets variable dimensions
        cell_data_dims = None
        for var_name in list_variable_data_in:
            if var_name in list(cell_datasets_in.keys()):
                if var_name != index_name_data_in:
                    cell_data_dims = cell_datasets_in[var_name].ndim
                    break
        if cell_data_dims is None:
            alg_logger.error(' ===> Variable data dimensions is not defined. Check the algortihm to solve the issue')
            raise RuntimeError('Variable data dimensions must be defined.')

        # iterate over variables
        cell_datasets_out, cell_attrs_out = {}, {}
        for var_name_in, var_name_out in zip(list_variable_in, list_variable_out):

            if var_name_in in list(cell_datasets_in.keys()):
                var_data_in = cell_datasets_in[var_name_in]
                attrs_data_in = attrs_datasets_in[var_name_in]

                if isinstance(var_data_in, np.ma.MaskedArray):
                    var_data_in = var_data_in.data

                if var_name_in == index_name_data_in:
                    var_data_in = pd.to_datetime(var_data_in).floor('S')

                if var_name_in == index_name_data_in:
                    var_data_in = pd.to_datetime(var_data_in).floor('S')

                # save variable(s) defined in list (or save all variable(s)
                if var_name_out not in cell_datasets_out:
                    cell_datasets_out[var_name_out] = {}
                    cell_datasets_out[var_name_out] = var_data_in
                elif var_name_out in cell_datasets_out:
                    var_data_tmp = cell_datasets_out[var_name_out]
                    var_data_tmp = np.concatenate([var_data_tmp, var_data_in])
                    cell_datasets_out[var_name_out] = var_data_tmp

                # save attributes defined in list (or save all attributes)
                if var_name_out not in cell_attrs_out:
                    cell_attrs_out[var_name_out] = {}
                    cell_attrs_out[var_name_out] = attrs_data_in

            else:
                alg_logger.warning(' ===> Variable "' + var_name_in + '" not available ih the cell datasets')

        collections_data, collections_index, collections_attrs = None, None, None
        for name_variable_data_out in list_variable_data_out:
            if name_variable_data_out == index_name_data_out:
                collections_index = cell_datasets_out[index_name_data_out]

            if collections_data is None:
                collections_data = {}
            collections_data[name_variable_data_out] = cell_datasets_out[name_variable_data_out]

            if collections_attrs is None:
                collections_attrs = {}
            collections_attrs[name_variable_data_out] = cell_attrs_out[name_variable_data_out]

        if collections_index is not None:
            collections_data[index_name_cell_out] = [int(cell_name)] * len(collections_index)
        else:
            alg_logger.warning(' ===> Variable "time" not available in the collections data')

        registry_data, registry_index, registry_attrs = None, None, None
        for name_variable_registry_out in list_variable_registry_out:
            if name_variable_registry_out == index_name_registry_out:
                registry_index = cell_datasets_out[index_name_registry_out]

            if registry_data is None:
                registry_data = {}
            registry_data[name_variable_registry_out] = cell_datasets_out[name_variable_registry_out]

            if registry_attrs is None:
                registry_attrs = {}
            registry_attrs[name_variable_registry_out] = cell_attrs_out[name_variable_registry_out]

        if registry_index is None:
            alg_logger.error(' ===> Variable "location_id" not available in the registry data')
            raise RuntimeError('Check your registry data object')

        if index_name_cell_out not in list(registry_data.keys()):
            registry_data[index_name_cell_out] = [int(cell_name)] * registry_index.shape[0]

        # organize collections
        if cell_data_dims == 1:
            # create object dataframe
            collections_obj = pd.DataFrame(data=collections_data, index=collections_index)
            for var_name in list(collections_obj.columns):
                if var_name in list(collections_attrs.keys()):
                    collections_obj[var_name].attrs = collections_attrs[var_name]

        elif cell_data_dims == 2:

            # get geo_x and geo_y variable(s)
            if index_name_registry_out in list(registry_data.keys()):
                var_location_id = registry_data[index_name_registry_out]
            else:
                alg_logger.error(' ===> Variable "location_id" not available in the registry data')
                raise RuntimeError('Check your registry data object')

            # get time variable(s)
            var_time = deepcopy(collections_index)

            # create object dataset
            collections_obj = xr.Dataset()
            for var_name in list_variable_data_out:

                # check variable name
                if var_name != index_name_data_out:

                    # get data
                    var_data = collections_data[var_name]
                    # create data array
                    var_da = xr.DataArray(
                        var_data,
                        dims=[index_name_registry_out, index_name_data_out],
                        coords={index_name_registry_out: (index_name_registry_out, var_location_id),
                                index_name_data_out: (index_name_data_out, var_time)})

                    # save data array
                    collections_obj[var_name] = var_da
            # save attributes
            collections_obj.attrs = cell_attrs_out

        else:
            # message if variable data dimensions is not supported
            alg_logger.error(' ===> Variable data dimensions is not supported.')
            raise RuntimeError('Check the algorithm to solve the issue.')

        # organize registry
        registry_obj = pd.DataFrame(data=registry_data, index=registry_index)
        for var_name in list(registry_obj.columns):
            if var_name in list(registry_attrs.keys()):
                registry_obj[var_name].attrs = registry_attrs[var_name]

    else:

        # info warning file not found
        alg_logger.warning(' ===> File cell "' + file_name + '" not found')

    return collections_obj, registry_obj

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extend gpi search
def extend_idx_search(gpi_in, lon_in, lat_in, gpi_list_in, grid_in, grid_ref, max_dist=np.Inf, max_gpi=5):

    gpi_ex_arr, distance_ex_arr = grid_ref.find_k_nearest_gpi(lon=lon_in, lat=lat_in, k=max_gpi, max_dist=max_dist)
    gpi_ex_list, distance_ex_list = gpi_ex_arr.tolist()[0], distance_ex_arr.tolist()[0]

    idx_out, gpi_out = None, None
    for gpi_ex, distance_ex in zip(gpi_ex_list, distance_ex_list):
        if gpi_ex != gpi_in:

            max_gpi = max_gpi + 10

            lon_ex, lat_ex = grid_ref.gpi2lonlat(gpi_ex)
            gpi_out_arr, distance_out_arr = grid_in.find_k_nearest_gpi(
                lon=lon_ex, lat=lat_ex, k=max_gpi, max_dist=max_dist)
            gpi_out_list, distance_out_list = gpi_out_arr.tolist()[0], distance_out_arr.tolist()[0]

            for gpi_out, distance_out in zip(gpi_out_list, distance_out_list):
                idx_out = np.argwhere(gpi_list_in == gpi_out)

                if idx_out.shape[0] != 0:
                    return idx_out.flatten()[0], gpi_out

    if idx_out.shape[0] == 0:
        idx_out, gpi_out = None, None

    if idx_out is None:
        alg_logger.warning(' ===> No gpi(s) found in extended distance methods')

    return idx_out, gpi_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to search idx lut between grids
def search_idx_between_grids(
        gpi_in, lon_in, lat_in, gpi_list_out, grid_obj_out, grid_obj_in=None, max_dist=25000,
        extend_dist=True, extend_gpi=10):

    # search idx lut between grids
    idx_out = np.argwhere(gpi_list_out == gpi_in)

    # check idx lut between grids
    if idx_out.shape[0] == 0:
        gpi_found = grid_obj_out.find_nearest_gpi(lon=lon_in, lat=lat_in, max_dist=max_dist)[0]
        lon_found, lat_found = grid_obj_out.gpi2lonlat(gpi_found)
        idx_out = np.argwhere(gpi_list_out == gpi_found)

        # check idx lut between grids
        if idx_out.shape[0] == 0:

            if extend_dist:
                alg_logger.info(' --------> Search gpi in extended distance ... ')
                idx_ex, gpi_ex = extend_idx_search(
                    gpi_found, lon_found, lat_found, gpi_list_out,
                    grid_obj_out, grid_obj_in,
                    max_dist=np.Inf, max_gpi=extend_gpi)

                if idx_ex is not None:
                    alg_logger.info(' ---------> Switch from "' + str(gpi_found) + '" to "' + str(gpi_ex) + '" gpi')
                    alg_logger.info(' --------> Search gpi in extended distance ... DONE')
                    return idx_ex
                else:
                    alg_logger.info(' --------> Search gpi in extended distance ... FAILED')
                    return None

            return None
        else:
            idx_out = idx_out.flatten()[0]

    elif idx_out.shape[0] == 1:
        idx_out = idx_out.flatten()[0]
    else:
        alg_logger.error(' ===> Multiple gpi(s) found for the same location')
        raise RuntimeError('Check your grid object')

    return idx_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get object 2d
def get_obj_2d(cell_data, cell_coord, cell_attrs, index_data, index_name_time='time'):

    select_data, select_coord, select_attrs = {}, {}, {}
    for cell_var, cell_values in cell_data.items():
        select_data[cell_var] = cell_values[index_data, :]

        if cell_var in list(cell_attrs.keys()):
            select_attrs[cell_var] = cell_attrs[cell_var]
        else:
            select_attrs[cell_var] = {}

    for cell_var, cell_values in cell_coord.items():
        select_coord[cell_var] = deepcopy(cell_values)

        if cell_var in list(cell_attrs.keys()):
            select_attrs[cell_var] = cell_attrs[cell_var]
        else:
            select_attrs[cell_var] = {}

    obj_dframe = pd.DataFrame(data=select_data, index=select_coord[index_name_time])
    obj_dframe = obj_dframe.sort_index()

    obj_attrs = deepcopy(select_attrs)

    return obj_dframe, obj_attrs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get object 1d
def get_obj_1d(cell_data, cell_coord, cell_attrs, row_start, row_end, index_name_time='time'):

    select_data, select_coord, select_attrs = {}, {}, {}
    for cell_var, cell_values in cell_data.items():
        tmp_values = cell_values[row_start:row_end]
        select_data[cell_var] = tmp_values

        if cell_var in list(cell_attrs.keys()):
            select_attrs[cell_var] = cell_attrs[cell_var]
        else:
            select_attrs[cell_var] = {}

    for cell_var, cell_values in cell_coord.items():
        select_coord[cell_var] = cell_values[row_start:row_end]

        if cell_var in list(cell_attrs.keys()):
            select_attrs[cell_var] = cell_attrs[cell_var]
        else:
            select_attrs[cell_var] = {}

    obj_dframe = pd.DataFrame(data=select_data, index=select_coord[index_name_time])
    obj_dframe = obj_dframe.sort_index()

    obj_attrs = deepcopy(select_attrs)

    return obj_dframe, obj_attrs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get object time-series
def get_obj_ts(var_name_list, obj_data, var_params=None):

    # iterate over variable(s) in the list
    var_values, var_time = None, None
    for var_name_step in var_name_list:

        # check object source and reference type
        if isinstance(obj_data, pd.DataFrame):

            if var_name_step in list(obj_data.columns):
                var_values = obj_data[var_name_step].values
                var_time = obj_data.index.values

                if var_params is not None:
                    if var_name_step in list(var_params.keys()):
                        var_min = var_params[var_name_step]['min_value']
                        var_max = var_params[var_name_step]['max_value']
                        var_scale_factor = var_params[var_name_step]['scale_factor']
                        var_no_data = var_params[var_name_step]['no_data']

                        if var_min is not None:
                            var_values[var_values < var_min] = np.nan
                        if var_max is not None:
                            var_values[var_values > var_max] = np.nan
                        if var_scale_factor is not None:
                            var_values = var_values * var_scale_factor
                        if var_no_data is not None:
                            var_values[var_values == var_no_data] = np.nan

                continue

        else:
            alg_logger.error(' ===> Data object must be a pandas DataFrame')
            raise RuntimeError('Check your data object')

    if var_values is None or var_time is None:
        alg_logger.error(' ===> Variable source data not available in the collections data')
        raise RuntimeError('Check your collections data object')

    return var_values, var_time
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select object values and time in a time-series
def select_obj_ts(time_obj_ref, data_obj_ref, method='nearest'):

    # apply method to select data and index
    if method == 'nearest':
        # get index for the nearest time
        var_idx_select = data_obj_ref.index.get_loc(time_obj_ref, method=method)
        # get data and time
        data_obj_select = data_obj_ref.iloc[var_idx_select]
        time_obj_select = data_obj_ref.index[var_idx_select]
    elif method == 'mean':
        # get data and time
        data_obj_select = data_obj_ref.mean()
        time_obj_select = data_obj_ref.index.mean()
    elif method == 'first':
        # get data and time
        data_obj_select = data_obj_ref.iloc[0]
        time_obj_select = data_obj_ref.index[0]
    elif method == 'last':
        # get data and time
        data_obj_select = data_obj_ref.iloc[-1]
        time_obj_select = data_obj_ref.index[-1]
    elif method == 'max':
        # get data and time
        data_obj_select = data_obj_ref.max()
        time_obj_select = data_obj_ref.index.max()
    elif method == 'min':
        # get data and time
        data_obj_select = data_obj_ref.min()
        time_obj_select = data_obj_ref.index.min()
    else:
        alg_logger.error(' ===> Select time-series method "' + method + '" not allowed')
        raise NotImplementedError('Case not implemented yet')

    data_obj_select = data_obj_select.values

    time_obj_delta = time_obj_select - time_obj_ref
    time_obj_diff = np.round(time_obj_delta.total_seconds() / 3600, 2)

    if data_obj_select.shape[0] == 1:
        data_obj_select = data_obj_select[0]
    else:
        alg_logger.warning(' ===> Select time-series returned with multiple values. Try to use the first value')
        data_obj_select = data_obj_select[0]

    return data_obj_select, time_obj_select, time_obj_diff
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to apply scaling datasets cell
def apply_scaling_datasets_cell(
        time_obj,
        collections_obj_ref_nrt, registry_obj_ref_nrt,
        collections_obj_ref_dr, registry_obj_ref_dr,
        collections_obj_other_dr, registry_obj_other_dr,
        grid_obj_ref_nrt, grid_obj_ref_dr, grid_obj_other_dr,
        var_params=None,
        max_dist=25000,
        fx_var_name_in=None, fx_var_name_out=None,
        fx_var_methods=None, fx_var_suffix_out="_scaled", fx_var_prefix_out=None,
        index_name_geo_x_ref_nrt='lon_ref_nrt', index_name_geo_y_ref_nrt='lat_ref_nrt',
        index_name_cell_ref_nrt='cell_ref_nrt', index_name_row_size_ref_nrt='row_size_ref_nrt',
        index_name_location_ref_nrt='location_id_ref_nrt', index_name_time_ref_nrt='time_ref_nrt',
        index_name_geo_x_ref_dr='lon_ref_dr', index_name_geo_y_ref_dr='lat_ref_dr',
        index_name_cell_ref_dr='cell_ref_dr', index_name_row_size_ref_dr='row_size_ref_dr',
        index_name_location_ref_dr='location_id_ref_dr', index_name_time_ref_dr='time_ref_dr',
        index_name_geo_x_other_dr='lon_other_dr', index_name_geo_y_other_dr='lat_other_dr',
        index_name_cell_other_dr='cell_other_dr', index_name_row_size_other_dr='row_size_other_dr',
        index_name_location_other_dr='location_id_other_dr', index_name_time_other_dr='time_other_dr',
        index_name_geo_x_ref_dst='lon_ref_dst', index_name_geo_y_ref_dst='lat_ref_dst',
        index_name_cell_ref_dst='cell_ref_dst', index_name_row_size_ref_dst='row_size_ref_dst',
        index_name_location_ref_dst='location_id_other_dr'):

    if fx_var_name_in is None:
        fx_var_name_in = {'layer_1': ["layer_1_ref_nrt", "layer_1_ref_dr", "layer_1_other_dr"]}
    if fx_var_name_out is None:
        fx_var_name_out = {'layer_1': ["layer_1_scaled"]}

    if not isinstance(fx_var_name_in, dict):
        alg_logger.error(' ===> Variable name for source data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    fx_var_name_in.pop('__comment__', None)
    if not isinstance(fx_var_name_out, dict):
        alg_logger.error(' ===> Variable name for destination data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    fx_var_name_out.pop('__comment__', None)
    if not isinstance(fx_var_methods, dict):
        alg_logger.error(' ===> Variable name for metrics data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    fx_var_methods.pop('__comment__', None)

    # check if lists variables and methods are equal
    if list(fx_var_name_in.keys()) != list(fx_var_name_out.keys()):
        alg_logger.error(' ===> Variable name for source and destination data must have the same keys')
        raise RuntimeError('Check your filter variable object')
    if list(fx_var_name_in.keys()) != list(fx_var_methods.keys()):
        alg_logger.error(' ===> Variable name for source data and methods data must have the same keys')
        raise RuntimeError('Check your filter variable object')
    if list(fx_var_name_out.keys()) != list(fx_var_methods.keys()):
        alg_logger.error(' ===> Variable name for destination data and methods data must have the same keys')
        raise RuntimeError('Check your filter variable object')

    # get registry reference nrt information
    alg_logger.info(' -------> Get source reference nrt registry ... ')

    # get reference nrt datasets
    (gpi_list_ref_nrt,
     lat_list_ref_nrt, lon_list_ref_nrt, cell_list_ref_nrt,
     row_size_n_list_ref_nrt, row_size_cumulative_list_ref_nrt) = extract_registry_info(
        registry_obj_ref_nrt,
        flag_type='registry_data',
        index_name_geo_x=index_name_geo_x_ref_nrt, index_name_geo_y=index_name_geo_y_ref_nrt,
        index_name_cell=index_name_cell_ref_nrt, index_name_row_size=index_name_row_size_ref_nrt,
        n_start=None, n_end=None)
    registry_obj_dst = deepcopy(registry_obj_ref_nrt)
    alg_logger.info(' -------> Get source reference nrt registry ... DONE')

    # get registry reference dr information
    alg_logger.info(' -------> Get source reference dr registry ... ')

    # get reference dr datasets
    (gpi_list_ref_dr,
     lat_list_ref_dr, lon_list_ref_dr, cell_list_ref_dr,
     row_size_n_list_ref_dr, row_size_cumulative_list_ref_dr) = extract_registry_info(
        registry_obj_ref_dr,
        flag_type='registry_data',
        index_name_geo_x=index_name_geo_x_ref_dr, index_name_geo_y=index_name_geo_y_ref_dr,
        index_name_cell=index_name_cell_ref_dr, index_name_row_size=index_name_row_size_ref_dr,
        n_start=None, n_end=None)

    alg_logger.info(' -------> Get source reference dr registry ... DONE')

    # get registry other information
    alg_logger.info(' -------> Get source other dr registry ... ')

    # get reference dr datasets
    (gpi_list_other_dr,
     lat_list_other_dr, lon_list_other_dr, cell_list_other_dr,
     row_size_n_list_other_dr, row_size_cumulative_list_other_dr) = extract_registry_info(
        registry_obj_other_dr,
        flag_type='registry_data',
        index_name_geo_x=index_name_geo_x_other_dr, index_name_geo_y=index_name_geo_y_other_dr,
        index_name_cell=index_name_cell_other_dr, index_name_row_size=index_name_row_size_other_dr,
        n_start=None, n_end=None)

    alg_logger.info(' -------> Get source other dr registry ... DONE')

    # define variable list expected
    var_list_expected = []
    for var_tmp in fx_var_name_in.values():
        if var_tmp is not None:
            var_list_expected.extend(var_tmp)

    # extract data cell and metrics datasets and attributes
    alg_logger.info(' -------> Get source reference nrt datasets, coord and attributes ... ')
    cell_data_ref_nrt, cell_coord_ref_nrt = extract_cell_data(collections_obj_ref_nrt, var_list_expected)
    cell_attrs_ref_nrt = extract_cell_attrs(collections_obj_ref_nrt)
    alg_logger.info(' -------> Get source reference nrt datasets, coord and attributes ... DONE')

    # extract data cell and metrics datasets and attributes
    alg_logger.info(' -------> Get source reference dr datasets, coord and attributes ... ')
    cell_data_ref_dr, cell_coord_ref_dr = extract_cell_data(collections_obj_ref_dr, var_list_expected)
    cell_attrs_ref_dr = extract_cell_attrs(collections_obj_ref_dr)
    alg_logger.info(' -------> Get source reference nrt datasets, coord and attributes ... DONE')

    # extract data cell and metrics datasets and attributes
    alg_logger.info(' -------> Get source other dr datasets, coord and attributes ... ')
    cell_data_other_dr, cell_coord_other_dr = extract_cell_data(collections_obj_other_dr, var_list_expected)
    cell_attrs_other_dr = extract_cell_attrs(collections_obj_other_dr)
    alg_logger.info(' -------> Get source other nrt datasets, coord and attributes ... DONE')

    # info start gpi(s) loop
    alg_logger.info(' -------> Compute variable(s) ... ')
    data_dst, time_step_dst, time_diff_dst = {}, {}, {}
    gpi_list_ref_dst, lon_list_ref_dst, lat_list_ref_dst, cell_list_ref_dst = [], [], [], []
    for idx_step_ref_nrt, (gpi_step_ref_nrt, lon_step_ref_nrt, lat_step_ref_nrt, cell_step_ref_nrt) in enumerate(
            zip(gpi_list_ref_nrt, lon_list_ref_nrt, lat_list_ref_nrt, cell_list_ref_nrt)):

        # method to search gpi between grids (reference nrt and reference dr)
        idx_step_ref_dr = search_idx_between_grids(
            gpi_step_ref_nrt, lon_step_ref_nrt, lat_step_ref_nrt, gpi_list_ref_dr,
            grid_obj_out=grid_obj_ref_dr, grid_obj_in=grid_obj_ref_nrt,
            max_dist=max_dist, extend_dist=True, extend_gpi=10)

        # method to search gpi between grids (reference nrt and other dr)
        idx_step_other_dr = search_idx_between_grids(
            gpi_step_ref_nrt, lon_step_ref_nrt, lat_step_ref_nrt, gpi_list_other_dr,
            grid_obj_out=grid_obj_other_dr, grid_obj_in=grid_obj_ref_nrt,
            max_dist=max_dist, extend_dist=True, extend_gpi=10)

        # check idx lut between grids
        if idx_step_ref_dr is not None and idx_step_other_dr is not None:

            # check data type for reference nrt (1d or 2d)
            if row_size_cumulative_list_ref_nrt is not None:

                # set indexes start and end
                row_size_start_ref_nrt = row_size_cumulative_list_ref_nrt[idx_step_ref_nrt]
                row_size_end_ref_nrt = row_size_cumulative_list_ref_nrt[idx_step_ref_nrt + 1]
                # get object 1d
                obj_ref_nrt, attrs_ref_nrt = get_obj_1d(
                    cell_data_ref_nrt, cell_coord_ref_nrt, cell_attrs_ref_nrt,
                    row_size_start_ref_nrt, row_size_end_ref_nrt,
                    index_name_time=index_name_time_ref_nrt)

            else:
                # message error for different format
                alg_logger.error(' ===> Data type for reference nrt must be 1d format')
                raise NotImplemented('Case not implemented yet')

            # check data type for reference dr (1d or 2d)
            if row_size_cumulative_list_ref_dr is not None:

                # set indexes start and end
                row_size_start_ref_dr = row_size_cumulative_list_ref_dr[idx_step_ref_dr]
                row_size_end_ref_dr = row_size_cumulative_list_ref_dr[idx_step_ref_dr + 1]
                # get object 1d
                obj_ref_dr, attrs_ref_dr = get_obj_1d(
                    cell_data_ref_dr, cell_coord_ref_dr, cell_attrs_ref_dr,
                    row_size_start_ref_dr, row_size_end_ref_dr,
                    index_name_time=index_name_time_ref_dr)

            else:
                # message error for different format
                alg_logger.error(' ===> Data type for reference nrt must be 1d format')
                raise NotImplemented('Case not implemented yet')

            # check data type for reference other (1d or 2d)
            if row_size_cumulative_list_other_dr is not None:
                # message error for different format
                alg_logger.error(' ===> Data type for other dr must be 2d format')
                raise NotImplemented('Case not implemented yet')
            else:
                # get object 2d
                obj_other_dr, attrs_other_dr = get_obj_2d(
                    cell_data_other_dr, cell_coord_other_dr, cell_attrs_other_dr,
                    idx_step_other_dr,
                    index_name_time=index_name_time_other_dr)

            # update registry information
            gpi_list_ref_dst.append(gpi_step_ref_nrt)
            lon_list_ref_dst.append(lon_step_ref_nrt)
            lat_list_ref_dst.append(lat_step_ref_nrt)
            cell_list_ref_dst.append(cell_step_ref_nrt)

        else:
            obj_ref_nrt, obj_ref_dr, obj_other_dr = None, None, None
            attrs_ref_nrt, attrs_ref_dr, attrs_other_dr = None, None, None
            alg_logger.warning(' ===> No gpi(s) found for the reference nrt location "' + str(gpi_step_ref_nrt) + '"')
            alg_logger.warning(' ===> Reference dr idx: "' + str(idx_step_ref_dr) + '" -- Other dr idx: "' +
                               str(idx_step_other_dr) + '"')

        # check datasets objects
        if obj_ref_nrt is not None and obj_ref_dr is not None and obj_other_dr is not None:

            # clean dictionary from NoneType
            fx_var_name_in = remove_key_with_null_value(fx_var_name_in)
            fx_var_name_out = remove_key_with_null_value(fx_var_name_out)
            fx_var_methods = remove_key_with_null_value(fx_var_methods)
            # check dictionary compatibility
            dict_other = [fx_var_name_out, fx_var_methods]
            check_key_between_dicts(dict_ref=fx_var_name_in, dict_other_list=dict_other)

            # iterate over variable(s) and method(s)
            for ((data_key_in, var_name_in), (data_key_out, var_name_out),
                 (fx_key, fx_signature)) in zip(fx_var_name_in.items(), fx_var_name_out.items(), fx_var_methods.items()):

                # check if variable key(s) are equal
                assert data_key_in == data_key_out == fx_key, ' ===> Variable key(s) must the same'

                # get source and reference variable(s)
                var_data_src, var_attrs_src, var_data_ref, var_attrs_ref = None, None, None, None
                var_time_src, var_time_ref = None, None
                collections_data = {'ref_nrt': None, 'ref_dr': None, 'other_dr': None}
                collections_time = {'ref_nrt': None, 'ref_dr': None, 'other_dr': None}

                # get values and time for reference nrt
                var_values_ref_nrt, var_time_ref_nrt = get_obj_ts(var_name_in, obj_ref_nrt, var_params)
                collections_data['ref_nrt'], collections_time['ref_nrt'] = var_values_ref_nrt, var_time_ref_nrt
                # get values and time for reference dr
                var_values_ref_dr, var_time_ref_dr = get_obj_ts(var_name_in, obj_ref_dr, var_params)
                collections_data['ref_dr'], collections_time['ref_dr'] = var_values_ref_dr, var_time_ref_dr
                # get values and time for other dr
                var_values_other_dr, var_time_other_dr = get_obj_ts(var_name_in, obj_other_dr, var_params)
                collections_data['other_dr'], collections_time['other_dr'] = var_values_other_dr, var_time_other_dr

                # check source and reference data
                if collections_data['ref_nrt'] is None:
                    alg_logger.error(' ===> Variable data reference nrt not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                if collections_data['ref_dr'] is None:
                    alg_logger.error(' ===> Variable data reference dr not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                if collections_data['other_dr'] is None:
                    alg_logger.error(' ===> Variable data other dr not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                # check source and reference time
                if collections_time['ref_nrt'] is None:
                    alg_logger.error(' ===> Variable time reference nrt not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                if collections_time['ref_dr'] is None:
                    alg_logger.error(' ===> Variable time reference dr not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                if collections_time['other_dr'] is None:
                    alg_logger.error(' ===> Variable time other dr not available in the collections data')
                    raise RuntimeError('Check your collections data object')

                # get fx signature
                if fx_signature is not None:
                    # get fx settings
                    fx_type, fx_args = fx_signature['type'], fx_signature['args']
                    # apply scaling method
                    var_data_dst, var_time_dst = compute_data_scaling(
                        collections_data, collections_time, fx_method=fx_type, fx_args=fx_args)
                    # check scaled datasets values
                    if (var_data_dst is not None) and (var_time_dst is not None):
                        # define variable dframe
                        var_dframe_dst = pd.DataFrame(data={var_name_out: var_data_dst}, index=var_time_dst)
                    else:
                        var_data_dst, var_time_dst, var_dframe_dst = None, None, None
                else:
                    var_data_dst, var_time_dst, var_dframe_dst = None, None, None

                # check scaled datasets values
                if var_dframe_dst is not None:

                    # method to select value and time in a time-series according to the method
                    var_value, var_time_step, var_time_diff = select_obj_ts(
                        time_obj, var_dframe_dst, method='nearest')

                    # initialize datasets objects
                    if var_name_out not in data_dst:
                        data_dst[var_name_out] = []
                        time_step_dst[var_name_out], time_diff_dst[var_name_out] = [], []
                    # update datasets objects
                    if var_value is not None:
                        # data
                        data_tmp = data_dst[var_name_out]
                        data_tmp.append(var_value)
                        data_dst[var_name_out] = data_tmp
                        # time step
                        data_tmp = time_step_dst[var_name_out]
                        data_tmp.append(var_time_step)
                        time_step_dst[var_name_out] = data_tmp
                        # time difference
                        data_tmp = time_diff_dst[var_name_out]
                        data_tmp.append(var_time_diff)
                        time_diff_dst[var_name_out] = data_tmp

                else:
                    # message warning if datasets variable is defined by NoneType
                    alg_logger.warning(' ===> Datasets variable "' + var_name_out + '" is defined by NoneType')

                    # initialize datasets objects
                    if var_name_out not in data_dst:
                        data_dst[var_name_out] = []
                        time_step_dst[var_name_out], time_diff_dst[var_name_out] = [], []

                    # data
                    data_tmp = data_dst[var_name_out]
                    data_tmp.append(np.nan)
                    data_dst[var_name_out] = data_tmp
                    # time step
                    data_tmp = time_step_dst[var_name_out]
                    data_tmp.append(time_obj)
                    time_step_dst[var_name_out] = data_tmp
                    # time difference
                    data_tmp = time_diff_dst[var_name_out]
                    data_tmp.append(np.nan)
                    time_diff_dst[var_name_out] = data_tmp

        else:
            # message warning if some or all datasets are defined by NoneType
            alg_logger.warning(' ===> Datasets object not available for the selected gpi')

    # info end gpi(s) loop
    alg_logger.info(' -------> Compute variable(s) ... DONE')

    # info start registry object
    alg_logger.info(' -------> Organize registry obj ... ')
    # define registry object
    data_registry_ref_dst = {
        index_name_geo_x_ref_dst: lon_list_ref_dst, index_name_geo_y_ref_dst: lat_list_ref_dst,
        index_name_location_ref_dst: gpi_list_ref_dst, index_name_cell_ref_dst: cell_list_ref_dst}
    # define registry dframe
    registry_obj_ref_dst = pd.DataFrame(data=data_registry_ref_dst, index=gpi_list_ref_dst)
    # info end registry object
    alg_logger.info(' -------> Organize registry obj ... DONE')

    # info start collections object
    alg_logger.info(' -------> Organize collections obj ... ')
    # define datasets object
    var_data_workspace, var_time_workspace, var_delta_workspace = {}, None, None
    for var_name_step in list(data_dst.keys()):

        var_data_arr = data_dst[var_name_step]

        if var_time_workspace is None:
            var_time_workspace = time_step_dst[var_name_step]
            var_data_workspace['time_scaled'] = var_time_workspace
        if var_delta_workspace is None:
            var_delta_workspace = time_diff_dst[var_name_step]
            var_data_workspace['delta_scaled'] = var_delta_workspace

        var_data_workspace[var_name_step] = var_data_arr
    # add other information
    var_data_workspace['cell_scaled'] = cell_list_ref_dst
    var_data_workspace['gpi_scaled'] = gpi_list_ref_dst
    var_data_workspace['lon_scaled'] = lon_list_ref_dst
    var_data_workspace['lat_scaled'] = lat_list_ref_dst

    # define collections dframe
    collections_obj_ref_dst = pd.DataFrame(data=var_data_workspace, index=gpi_list_ref_dst)

    # info end collections object
    alg_logger.info(' -------> Organize collections obj ... DONE')

    return collections_obj_ref_dst, registry_obj_ref_dst
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to dump datasets cell
def dump_datasets_cell(file_path,
                       data_obj, registry_obj,
                       data_attrs, registry_attrs,
                       list_variable_data_in, list_variable_data_out,
                       list_variable_registry_in, list_variable_registry_out):

    # change variable names data
    for var_name_in, var_name_out in zip(list_variable_data_in, list_variable_data_out):
        if isinstance(data_obj, pd.DataFrame):
            if var_name_in in list(data_obj.columns):
                data_obj.rename(columns={var_name_in: var_name_out}, inplace=True)
            else:
                alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the data object')
            if var_name_in in list(data_attrs.keys()):
                if var_name_out not in list(data_attrs.keys()):
                    data_attrs[var_name_out] = data_attrs[var_name_in]
                    data_attrs.pop(var_name_in)
        elif isinstance(data_obj, xr.Dataset):
            if var_name_in in list(data_obj.data_vars):
                data_obj = data_obj.rename({var_name_in: var_name_out})
            elif var_name_in in list(data_obj.coords):
                data_obj = data_obj.rename({var_name_in: var_name_out})
            else:
                alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the data object')
            if var_name_in in list(data_attrs.keys()):
                if var_name_out not in list(data_attrs.keys()):
                    data_attrs[var_name_out] = data_attrs[var_name_in]
                    data_attrs.pop(var_name_in)
        else:
            alg_logger.error(' ===> Data object must be a pandas DataFrame or a xarray Dataset')
            raise RuntimeError('Check your data object')

    # change variable names (mixed data and registry)
    for var_name_in, var_name_out in zip(list_variable_registry_in, list_variable_registry_out):
        if isinstance(data_obj, xr.Dataset):
            if var_name_in in list(data_obj.data_vars):
                data_obj = data_obj.rename({var_name_in: var_name_out})
            elif var_name_in in list(data_obj.coords):
                data_obj = data_obj.rename({var_name_in: var_name_out})

            if var_name_in in list(data_attrs.keys()):
                if var_name_out not in list(data_attrs.keys()):
                    data_attrs[var_name_out] = data_attrs[var_name_in]
                    data_attrs.pop(var_name_in)

    # change registry names (registry)
    for var_name_in, var_name_out in zip(list_variable_registry_in, list_variable_registry_out):
        if var_name_in in list(registry_obj.columns):
            registry_obj.rename(columns={var_name_in: var_name_out}, inplace=True)
        else:
            alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the registry object')

        if var_name_in in list(registry_attrs.keys()):
            if var_name_out not in list(registry_attrs.keys()):
                registry_attrs[var_name_out] = registry_attrs[var_name_in]
                registry_attrs.pop(var_name_in)

    # change registry names (mixed data and registry)
    for var_name_in, var_name_out in zip(list_variable_data_in, list_variable_data_out):
        if var_name_in in list(registry_obj.columns):
            registry_obj.rename(columns={var_name_in: var_name_out}, inplace=True)

        if var_name_in in list(registry_attrs.keys()):
            if var_name_out not in list(registry_attrs.keys()):
                registry_attrs[var_name_out] = registry_attrs[var_name_in]
                registry_attrs.pop(var_name_in)

    # drop variable not in the list (data and registry)
    if isinstance(data_obj, pd.DataFrame):
        for var_name in list(data_obj.columns):
            if var_name not in list_variable_data_out:
                data_obj.drop(var_name, axis=1, inplace=True)
                if var_name in list(data_attrs.keys()):
                    data_attrs.pop(var_name)
    elif isinstance(data_obj, xr.Dataset):
        for var_name in list(data_obj.data_vars):
            if var_name not in list_variable_data_out:
                data_obj = data_obj.drop_vars(var_name)
                if var_name in list(data_attrs.keys()):
                    data_attrs.pop(var_name)
    else:
        alg_logger.error(' ===> Data object must be a pandas DataFrame or a xarray Dataset')
        raise RuntimeError('Check your data object')

    for reg_name in list(registry_obj.columns):
        if reg_name not in list_variable_registry_out:
            registry_obj.drop(reg_name, axis=1, inplace=True)
            if reg_name in list(registry_attrs.keys()):
                registry_attrs.pop(reg_name)

    # dump destination datasets in netcdf format
    folder_name, file_name = os.path.split(file_path)
    os.makedirs(folder_name, exist_ok=True)

    # write datasets
    write_file_cell(file_path,
                    data_obj, registry_obj, data_attrs, registry_attrs,
                    file_tag_location='location_id')

# ----------------------------------------------------------------------------------------------------------------------
