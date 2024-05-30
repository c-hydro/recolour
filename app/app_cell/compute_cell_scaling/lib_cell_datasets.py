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
# method to apply scaling datasets cell
def apply_scaling_datasets_cell(
        time_obj,
        collections_obj_data_src, registry_obj_data_src,
        collections_obj_ref, registry_obj_ref,
        grid_obj_data, grid_obj_ref,
        max_dist=25000, var_min=0, var_max=100, var_no_data=np.nan,
        fx_var_name_in=None, fx_var_name_out=None,
        fx_var_methods=None, fx_var_suffix_out="_scaled", fx_var_prefix_out=None,
        index_name_geo_x_data='lon', index_name_geo_y_data='lat',
        index_name_cell_data='cell', index_name_row_size_data='row_size',
        index_name_location_data='location_id', index_name_time_data='time',
        index_name_geo_x_ref='lon', index_name_geo_y_ref='lat',
        index_name_cell_ref='cell', index_name_row_size_ref='row_size',
        index_name_location_ref='location_id', index_name_time_ref='time'):

    if fx_var_name_in is None:
        fx_var_name_in = {'layer_1': "layer_1"}
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

    # get registry source information
    alg_logger.info(' -------> Get registry source and destination datasets ... ')
    (gpi_list_data,
     lat_list_data, lon_list_data, cell_list_data,
     row_size_n_list_data, row_size_cumulative_list_data) = extract_registry_info(
        registry_obj_data_src,
        flag_type='registry_data',
        index_name_geo_x=index_name_geo_x_data, index_name_geo_y=index_name_geo_y_data,
        index_name_cell=index_name_cell_data, index_name_row_size=index_name_row_size_data,
        n_start=None, n_end=None)
    registry_obj_data_dst = deepcopy(registry_obj_data_src)
    alg_logger.info(' -------> Get registry source and destination datasets ... DONE')

    # get registry reference information
    alg_logger.info(' -------> Get registry reference datasets ... ')
    (gpi_list_ref,
     lat_list_ref, lon_list_ref, cell_list_ref,
     row_size_n_list_ref, row_size_cumulative_list_ref) = extract_registry_info(
        registry_obj_ref,
        flag_type='registry_reference',
        index_name_geo_x=index_name_geo_x_ref, index_name_geo_y=index_name_geo_y_ref,
        index_name_cell=index_name_cell_ref, index_name_row_size=index_name_row_size_ref,
        n_start=None, n_end=None)
    alg_logger.info(' -------> Get registry reference datasets ... DONE')

    # define variable list expected
    var_list_expected = []
    for var_tmp in fx_var_name_in.values():
        if var_tmp is not None:
            var_list_expected.extend(var_tmp)

    # extract data cell and metrics datasets and attributes
    alg_logger.info(' -------> Get source datasets, coord and attributes ... ')
    cell_data_src, cell_coord_src = extract_cell_data(collections_obj_data_src, var_list_expected)
    cell_attrs_src = extract_cell_attrs(collections_obj_data_src)
    alg_logger.info(' -------> Get source datasets, coord and attributes ... DONE')

    # extract reference cell and metrics datasets and attributes
    alg_logger.info(' -------> Get reference datasets, coord and attributes ... ')
    cell_data_ref, cell_coord_ref = extract_cell_data(collections_obj_ref, var_list_expected)
    cell_attrs_ref = extract_cell_attrs(collections_obj_ref)
    alg_logger.info(' -------> Get reference datasets, coord and attributes ... DONE')

    # iterate over gpi(s), lon(s) and lat(s)
    cell_data_dst, cell_coord_dst, cell_attrs_dst = {}, {}, {}
    dims_dst, row_size_dst = None, []
    attrs_dst = {}
    flag_not_set, flag_none = {}, {}
    idx_not_found = []
    # info start gpi(s) loop
    alg_logger.info(' -------> Compute datasets scaling ... ')
    for idx_step_data, (gpi_step_data, lon_step_data, lat_step_data) in enumerate(
            zip(gpi_list_data, lon_list_data, lat_list_data)):

        # search idx lut between grids
        idx_step_ref = np.argwhere(gpi_list_ref == gpi_step_data)

        # check idx lut between grids
        if idx_step_ref.shape[0] == 0:
            gpi_step_ref = grid_obj_ref.find_nearest_gpi(lon=lon_step_data, lat=lat_step_data, max_dist=max_dist)[0]
            idx_step_ref = np.argwhere(gpi_list_ref == gpi_step_ref)

            # check idx lut between grids
            if idx_step_ref.shape[0] == 0:
                alg_logger.warning(' ===> No gpi(s) found for the location "' + str(gpi_step_data) + '"')
                idx_not_found.append(gpi_step_data)
                continue

        elif idx_step_ref.shape[0] == 1:
            idx_step_ref = idx_step_ref.flatten()[0]
        else:
            alg_logger.error(' ===> Multiple gpi(s) found for the same location')
            raise RuntimeError('Check your grid object')

        # check data type (1d or 2d)
        if row_size_cumulative_list_ref is not None:

            # set indexes start and end
            row_size_start_ref = row_size_cumulative_list_ref[idx_step_ref]
            row_size_end_ref = row_size_cumulative_list_ref[idx_step_ref + 1]

            select_data_ref, select_coord_ref = {}, {}
            for cell_var, cell_values in cell_data_ref.items():
                tmp_values = cell_values[row_size_start_ref:row_size_end_ref]
                tmp_values[tmp_values < var_min] = np.nan
                tmp_values[tmp_values > var_max] = np.nan
                select_data_ref[cell_var] = tmp_values

            for cell_var, cell_values in cell_coord_ref.items():
                select_coord_ref[cell_var] = cell_values[row_size_start_ref:row_size_end_ref]

            obj_ref = pd.DataFrame(data=select_data_ref, index=select_coord_ref[index_name_time_ref])
            obj_ref = obj_ref.sort_index()

        else:

            select_data_ref, select_coord_ref = {}, {}
            for cell_var, cell_values in cell_data_ref.items():
                select_data_ref[cell_var] = cell_values[idx_step_data, :]
            for cell_var, cell_values in cell_coord_ref.items():
                select_coord_ref[cell_var] = deepcopy(cell_values)

            obj_ref = pd.DataFrame(data=select_data_ref, index=select_coord_ref[index_name_time_ref])
            obj_ref = obj_ref.sort_index()

        row_size_start_src, row_size_end_src = None, None
        if row_size_cumulative_list_data is not None:
            # set indexes start and end
            row_size_start_src = row_size_cumulative_list_data[idx_step_data]
            row_size_end_src = row_size_cumulative_list_data[idx_step_data + 1]

            select_data_value_src, select_coord_value_src = {}, {}
            select_data_attrs_src, select_coord_attrs_src = {}, {}
            for cell_var, cell_values in cell_data_src.items():
                tmp_values = cell_values[row_size_start_src:row_size_end_src]
                tmp_values[tmp_values < var_min] = np.nan
                tmp_values[tmp_values > var_max] = np.nan
                select_data_value_src[cell_var] = tmp_values

                if cell_var in list(cell_attrs_src.keys()):
                    select_data_attrs_src[cell_var] = cell_attrs_src[cell_var]
                else:
                    select_data_attrs_src[cell_var] = {}

            for cell_var, cell_values in cell_coord_src.items():
                select_coord_value_src[cell_var] = cell_values[row_size_start_src:row_size_end_src]

                if cell_var in list(cell_attrs_src.keys()):
                    select_coord_attrs_src[cell_var] = cell_attrs_src[cell_var]
                else:
                    select_coord_attrs_src[cell_var] = {}

            obj_src = pd.DataFrame(data=select_data_value_src, index=select_coord_value_src[index_name_time_data])
            obj_src = obj_src.sort_index()

            attrs_src = {**select_data_attrs_src, **select_coord_attrs_src}
            attrs_dst = deepcopy(attrs_src)

        else:

            select_data_value_src, select_coord_value_src = {}, {}
            select_data_attrs_src, select_coord_attrs_src = {}, {}
            for cell_var, cell_values in cell_data_src.items():
                select_data_value_src[cell_var] = cell_values[idx_step_data, :]

                if cell_var in list(cell_attrs_src.keys()):
                    select_data_attrs_src[cell_var] = cell_attrs_src[cell_var]
                else:
                    select_data_attrs_src[cell_var] = {}

            for cell_var, cell_values in cell_coord_src.items():
                select_coord_value_src[cell_var] = deepcopy(cell_values)

                if cell_var in list(cell_attrs_src.keys()):
                    select_coord_attrs_src[cell_var] = cell_attrs_src[cell_var]
                else:
                    select_coord_attrs_src[cell_var] = {}

            obj_src = pd.DataFrame(data=select_data_value_src, index=select_coord_ref[index_name_time_ref])
            obj_src = obj_src.sort_index()

            attrs_src = {**select_data_attrs_src, **select_coord_attrs_src}
            attrs_dst = deepcopy(attrs_src)

        # iterate over variable(s) and method(s)
        row_size_tmp = None
        for ((data_key_in, var_name_in), (data_key_out, var_name_out),
             (fx_key, fx_signature)) in zip(fx_var_name_in.items(), fx_var_name_out.items(), fx_var_methods.items()):

            # check if variable key(s) are equal
            assert data_key_in == data_key_out == fx_key, ' ===> Variable key(s) must the same'

            # get source and reference variable(s)
            var_data_src, var_attrs_src, var_data_ref, var_attrs_ref = None, None, None, None
            var_time_src, var_time_ref = None, None
            collections_data, collections_time = {'src': None, 'ref': None}, {'src': None, 'ref': None}

            # check variables defined in the list
            var_link, var_n_data = {}, None
            if var_name_in is not None:

                # iterate over variable(s) in the list
                for var_name_step in var_name_in:

                    # check object source and reference type
                    if isinstance(obj_src, pd.DataFrame) and isinstance(obj_ref, pd.DataFrame):

                        # get variable source data
                        if var_data_src is None:
                            if var_name_step in list(obj_src.columns):
                                var_data_src = obj_src[var_name_step].values
                                collections_data['src'] = var_data_src

                                var_link[var_name_out] = var_name_step
                                if var_name_step not in list(cell_data_dst.keys()):
                                    if var_name_step in list(cell_data_src.keys()):
                                        cell_data_dst[var_name_step] = cell_data_src[var_name_step]
                                        var_n_data = cell_data_src[var_name_step].shape[0]
                                    else:
                                        cell_data_dst[var_name_step] = None

                        # get variable reference time
                        if var_time_src is None:
                            var_time_src = obj_src.index.values
                            collections_time['src'] = var_time_src

                        # get variable reference data
                        if var_data_ref is None:
                            if var_name_step in list(obj_ref.columns):
                                var_data_ref = obj_ref[var_name_step].values
                                collections_data['ref'] = var_data_ref
                        # get variable reference time
                        if var_time_ref is None:
                            var_time_ref = obj_ref.index.values
                            collections_time['ref'] = var_time_ref
                    else:
                        alg_logger.error(' ===> Data object must be a pandas DataFrame')
                        raise RuntimeError('Check your data object')

                # check source and reference data
                if collections_data['src'] is None:
                    alg_logger.error(' ===> Variable source data not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                if collections_data['ref'] is None:
                    alg_logger.error(' ===> Variable reference data not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                # check source and reference time
                if collections_time['src'] is None:
                    alg_logger.error(' ===> Variable source time not available in the collections data')
                    raise RuntimeError('Check your collections data object')
                if collections_time['ref'] is None:
                    alg_logger.error(' ===> Variable reference time not available in the collections data')
                    raise RuntimeError('Check your collections data object')

                # get fx signature
                if fx_signature is not None:
                    # get fx settings
                    fx_type, fx_args = fx_signature['type'], fx_signature['args']

                    # apply scaling method
                    var_data_dst, var_time_dst = compute_data_scaling(
                        collections_data, collections_time, fx_method=fx_type, fx_args=fx_args)
                else:
                    var_data_dst, var_time_dst = None, None

                # check scaled datasets values
                if var_data_dst is not None:

                    # check scaled datasets format
                    if row_size_cumulative_list_data is not None:

                        time_n = var_time_dst.shape[0]

                        if var_name_out not in cell_data_dst:
                            var_data_tmp = np.zeros(shape=(var_n_data,))
                            var_data_tmp[:] = np.nan
                            cell_data_dst[var_name_out] = var_data_tmp
                            cell_data_dst[var_name_out][row_size_start_src:row_size_end_src] = var_data_dst
                        else:
                            cell_data_dst[var_name_out][row_size_start_src:row_size_end_src] = var_data_dst

                        if index_name_time_data not in cell_coord_dst:

                            var_data_type = np.dtype('datetime64[s]')
                            var_data_tmp = np.empty((var_n_data,), dtype=var_data_type)
                            var_data_tmp[:] = np.datetime64('NaT')
                            cell_coord_dst[index_name_time_data] = var_data_tmp
                            cell_coord_dst[index_name_time_data][row_size_start_src:row_size_end_src] = var_time_dst
                        else:
                            cell_coord_dst[index_name_time_data][row_size_start_src:row_size_end_src] = var_time_dst

                        # check datasets
                        check_data = pd.DataFrame(
                            index=cell_coord_dst[index_name_time_data],
                            data={'data': cell_data_dst[var_name_out], 'time': cell_coord_dst[index_name_time_data]})

                        check_data = check_data.dropna(axis=0)
                        check_n_data = check_data['data'].shape[0]
                        check_n_time = check_data['time'].shape[0]

                        if check_n_data != check_n_time:
                            alg_logger.warning(' ===> Data and time size not matching')

                        if row_size_tmp is None:
                            row_size_dst.append(time_n)
                            row_size_tmp = True

                        if dims_dst is None:
                            dims_dst = cell_data_dst[var_name_out].ndim

                        if var_name_out in list(var_link.keys()):
                            var_name_in = var_link[var_name_out]
                            if var_name_in in list(attrs_dst.keys()):
                                attrs_var = attrs_dst[var_name_in]
                                if 'name' in attrs_var:
                                    var_tag_in = attrs_var['name']
                                    var_tag_out = 'variable_scaled'
                                    if fx_var_prefix_out is not None:
                                        var_tag_out = fx_var_prefix_out + var_tag_in
                                    if fx_var_suffix_out is not None:
                                        var_tag_out = var_tag_in + fx_var_suffix_out

                                    attrs_dst[var_name_out] = {'name': var_tag_out}

                    else:

                        time_n = cell_coord_src[index_name_time_data].shape[0]
                        location_n = cell_coord_src[index_name_location_data].shape[0]

                        if var_name_out not in cell_data_dst:
                            var_data_tmp = np.zeros(shape=(location_n, time_n))
                            var_data_tmp[:, :] = np.nan
                            cell_data_dst[var_name_out] = var_data_tmp

                        if index_name_time_data not in cell_coord_dst:
                            cell_coord_dst[index_name_time_data] = var_time_dst
                        if index_name_location_data not in cell_coord_dst:
                            cell_coord_dst[index_name_location_data] = cell_coord_src[index_name_location_data]

                        if row_size_tmp is None:
                            row_size_dst.append(time_n)
                            row_size_tmp = True

                        if var_data_dst is not None:
                            cell_data_dst[var_name_out][idx_step_data, :] = var_data_dst

                        if dims_dst is None:
                            dims_dst = cell_data_dst[var_name_out].ndim

                        if var_name_out in list(var_link.keys()):
                            var_name_in = var_link[var_name_out]
                            if var_name_in in list(attrs_dst.keys()):
                                attrs_var = attrs_dst[var_name_in]
                                if 'name' in attrs_var:
                                    var_tag_in = attrs_var['name']
                                    var_tag_out = 'variable_scaled'
                                    if fx_var_prefix_out is not None:
                                        var_tag_out = fx_var_prefix_out + var_tag_in
                                    if fx_var_suffix_out is not None:
                                        var_tag_out = var_tag_in + fx_var_suffix_out

                                    attrs_dst[var_name_out] = {'name': var_tag_out}

                else:
                    if data_key_in not in list(flag_none.keys()):
                        flag_none[data_key_in] = True
                        alg_logger.warning(' ===> Variable "' + data_key_in + '" defined by NoneType object')

            else:
                if data_key_in not in list(flag_not_set.keys()):
                    flag_not_set[data_key_in] = True
                    alg_logger.warning(' ===> Variable "' + data_key_in + '" not set for scaling method')

    # info end gpi(s) loop
    alg_logger.info(' -------> Compute datasets scaling ... DONE')

    # remove gpi(s) not found from the datasets and registry
    alg_logger.info(' -------> Remove dataset empty gpi(s) ... ')
    if idx_not_found:

        # define registry object
        data_registry_tmp = {
            index_name_geo_x_data: lon_list_data, index_name_geo_y_data: lat_list_data,
            index_name_location_data: gpi_list_data,
            index_name_cell_data: cell_list_data, index_name_row_size_data: row_size_n_list_data}
        obj_registry_tmp = pd.DataFrame(data=data_registry_tmp, index=gpi_list_data)
        obj_registry_tmp = obj_registry_tmp.drop(index=idx_not_found)

        # get registry reference information
        (gpi_list_data,
         lat_list_data, lon_list_data, cell_list_data,
         row_size_n_list_data, row_size_cumulative_list_data) = extract_registry_info(
            obj_registry_tmp,
            flag_type='registry_reference',
            index_name_geo_x=index_name_geo_x_data, index_name_geo_y=index_name_geo_y_data,
            index_name_cell=index_name_cell_data, index_name_row_size=index_name_row_size_data,
            n_start=None, n_end=None)
        # info end remove gpi(s) not found
        alg_logger.info(' -------> Remove dataset empty gpi(s) ... DONE')

    else:
        # info end remove gpi(s) not found
        alg_logger.info(' -------> Remove dataset empty gpi(s) ... SKIPPED. All gpi(s) found')

    # organize collections
    alg_logger.info(' -------> Organize dataset collections ... ')
    if dims_dst == 1:

        # get object data
        obj_data = deepcopy(cell_data_dst)
        obj_attrs = deepcopy(attrs_dst)
        obj_time = deepcopy(cell_coord_dst[index_name_time_data])

        # check time not available
        idx_no_time = np.argwhere(np.isnat(obj_time))[:, 0]
        tmp_time = deepcopy(cell_coord_src[index_name_time_data])
        obj_time[idx_no_time] = deepcopy(tmp_time[idx_no_time])
        idx_test = np.argwhere(np.isnat(obj_time))[:, 0]

        # check data size
        n_data = obj_time.shape[0]
        for data_key, data_values in obj_data.items():
            if data_values.shape[0] > n_data:
                obj_data[data_key] = data_values[:n_data]
            elif data_values.shape[0] == n_data:
                pass
            else:
                alg_logger.error(' ===> Data size not matching')
                raise RuntimeError('Check your data size')
        obj_data[index_name_time_data] = obj_time

        # define registry object
        obj_registry = {index_name_geo_x_data: lon_list_data, index_name_geo_y_data: lat_list_data,
                        index_name_location_data: gpi_list_data,
                        index_name_cell_data: cell_list_data, index_name_row_size_data: row_size_n_list_data}

        if row_size_n_list_data.tolist() != row_size_dst:
            alg_logger.error(' ===> Row size list not matching')
            raise RuntimeError('Check your row size list')

        registry_obj_data_dst = pd.DataFrame(data=obj_registry, index=gpi_list_data)

        # define collections data
        collections_obj_data_dst = pd.DataFrame(data=obj_data, index=obj_time)
        collections_obj_data_dst = collections_obj_data_dst.fillna(var_no_data)

        # define collections attributes
        collections_attrs_dst, registry_attrs_dst = {}, {}
        for var_name in list(obj_data.keys()):
            if var_name in list(obj_attrs.keys()):
                collections_attrs_dst[var_name] = obj_attrs[var_name]

        for var_name in [index_name_time_data, index_name_location_data]:
            if var_name == index_name_time_data:
                if var_name in list(obj_attrs.keys()):
                    registry_attrs_dst[var_name] = obj_attrs[index_name_time_data]
                else:
                    registry_attrs_dst[var_name]= {}
            elif var_name == index_name_location_data:
                if var_name in list(obj_attrs.keys()):
                    registry_attrs_dst[var_name] = obj_attrs[index_name_location_data]
                else:
                    registry_attrs_dst[var_name] = {}

    elif dims_dst == 2:

        obj_data = deepcopy(cell_data_dst)
        obj_attrs = deepcopy(attrs_dst)
        coord_time = deepcopy(cell_coord_dst[index_name_time_data])
        coord_location = deepcopy(cell_coord_dst[index_name_location_data])

        collections_obj_data_dst = xr.Dataset()
        collections_attrs_dst, registry_attrs_dst = {}, {}
        for var_name in list(obj_data.keys()):
            if var_name not in [index_name_time_data, index_name_location_data]:
                collections_obj_data_dst[var_name] = xr.DataArray(
                    obj_data[var_name],
                    dims=[index_name_location_data, index_name_time_data],
                    coords={
                        index_name_location_data: (index_name_location_data, coord_location),
                        index_name_time_data: (index_name_time_data, coord_time)
                    }
                )
                collections_attrs_dst[var_name] = obj_attrs[var_name]

        for var_name in [index_name_time_data, index_name_location_data]:
            if var_name == index_name_time_data:
                if var_name in list(obj_attrs.keys()):
                    registry_attrs_dst[var_name] = obj_attrs[index_name_time_data]
                else:
                    registry_attrs_dst[var_name] = {}
            elif var_name == index_name_location_data:
                if var_name in list(obj_attrs.keys()):
                    registry_attrs_dst[var_name] = obj_attrs[index_name_location_data]
                else:
                    registry_attrs_dst[var_name] = {}

    else:
        alg_logger.error(' ===> Dimension not allowed for collections object')
        raise RuntimeError('Check dimensions of your collections object')

    # info end organize collections
    alg_logger.info(' -------> Organize dataset collections ... DONE')

    return collections_obj_data_dst, registry_obj_data_dst, collections_attrs_dst, registry_attrs_dst
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


# ----------------------------------------------------------------------------------------------------------------------
# method to dump datasets grid
def dump_datasets_grid(file_path, file_grid_obj, file_grid_reference=None):

    # check file availability
    if not os.path.exists(file_path):

        # update gpis and active gpis (if needed)
        # file_gpis = file_grid_reference['gpi'].values
        # file_grid_obj.gpis = file_gpis
        # file_grid_obj.activegpis = file_gpis

        # make grid folder
        folder_name, file_name = os.path.split(file_path)
        os.makedirs(folder_name, exist_ok=True)
        # save grid information in file
        grid2nc.save_grid(file_path, file_grid_obj)
# ----------------------------------------------------------------------------------------------------------------------
