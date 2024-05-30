"""
Library Features:

Name:          lib_utils_datasets
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

from lib_data_io_nc import read_file_cell, write_file_cell
from lib_data_analysis import compute_data_metrics
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
        index_name_data='time', index_name_registry='location_id', index_name_cell='cell',
        index_name_row_size='row_size'):

    # if row_index is the first variable in the list, move it to the end
    if index_name_row_size == list_variable_registry_in[0]:
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
                if var_name != index_name_data:
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

                if var_name_in == index_name_data:
                    var_data_in = pd.to_datetime(var_data_in).floor('S')

                if var_name_in == index_name_data:
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
            if name_variable_data_out == index_name_data:
                collections_index = cell_datasets_out[index_name_data]

            if collections_data is None:
                collections_data = {}
            collections_data[name_variable_data_out] = cell_datasets_out[name_variable_data_out]

            if collections_attrs is None:
                collections_attrs = {}
            collections_attrs[name_variable_data_out] = cell_attrs_out[name_variable_data_out]

        collections_data[index_name_cell] = [int(cell_name)] * len(collections_index)

        registry_data, registry_index, registry_attrs = None, None, None
        for name_variable_registry_out in list_variable_registry_out:
            if name_variable_registry_out == index_name_registry:
                registry_index = cell_datasets_out[index_name_registry]

            if registry_data is None:
                registry_data = {}
            registry_data[name_variable_registry_out] = cell_datasets_out[name_variable_registry_out]

            if registry_attrs is None:
                registry_attrs = {}
            registry_attrs[name_variable_registry_out] = cell_attrs_out[name_variable_registry_out]

        if index_name_cell not in list(registry_data.keys()):
            registry_data[index_name_cell] = [int(cell_name)] * registry_index.shape[0]

        # organize collections
        if cell_data_dims == 1:
            # create object dataframe
            collections_obj = pd.DataFrame(data=collections_data, index=collections_index)
            for var_name in list(collections_obj.columns):
                if var_name in list(collections_attrs.keys()):
                    collections_obj[var_name].attrs = collections_attrs[var_name]

        elif cell_data_dims == 2:

            # get geo_x and geo_y variable(s)
            if index_name_registry in list(registry_data.keys()):
                var_location_id = registry_data[index_name_registry]
            else:
                alg_logger.error(' ===> Variable "location_id" not available in the registry data')
                raise RuntimeError('Check your registry data object')

            # get time variable(s)
            var_time = deepcopy(collections_index)

            # create object dataset
            collections_obj = xr.Dataset()
            for var_name in list_variable_data_out:

                # check variable name
                if var_name != index_name_data:

                    # get data
                    var_data = collections_data[var_name]
                    # create data array
                    var_da = xr.DataArray(
                        var_data,
                        dims=[index_name_registry, index_name_data],
                        coords={index_name_registry: (index_name_registry, var_location_id),
                                index_name_data: (index_name_data, var_time)})

                    # save data array
                    collections_obj[var_name] = var_da

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
# method to apply metrics datasets cell
def apply_metrics_datasets_cell(
        time_obj, collections_obj_in, registry_obj_in,
        fx_var_name_in=None, fx_var_name_out=None,
        fx_var_methods=None):

    if fx_var_name_in is None:
        fx_var_name_in = {'sm_data': "sm"}
    if fx_var_name_out is None:
        fx_var_name_out = {'sm_data': ["sm_mean", "sm_min", "sm_max", "sm_std"]}

    if not isinstance(fx_var_name_in, dict):
        alg_logger.error(' ===> Variable name for source data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    if not isinstance(fx_var_name_out, dict):
        alg_logger.error(' ===> Variable name for destination data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    if not isinstance(fx_var_methods, dict):
        alg_logger.error(' ===> Variable name for metrics data must be a dictionary')
        raise RuntimeError('Check your filter variable object')

    # get registry info
    reg_index_list = registry_obj_in.index.values
    reg_latitude_list, reg_latitude_attrs = registry_obj_in['lat'].values, registry_obj_in['lat'].attrs
    reg_longitude_list, reg_longitude_attrs = registry_obj_in['lon'].values, registry_obj_in['lon'].attrs
    reg_cell_list, reg_cell_attrs = registry_obj_in['cell'].values, registry_obj_in['cell'].attrs
    if 'row_size' in list(registry_obj_in.columns):
        reg_row_size_list, reg_row_size_attrs = registry_obj_in['row_size'].values, registry_obj_in['row_size'].attrs
        # organize row size
        reg_row_size_cumulative = np.append(0, np.cumsum(reg_row_size_list))
    else:
        reg_row_size_list, reg_row_size_attrs, reg_row_size_cumulative = None, None, None
        alg_logger.warning(' ===> Variable "row_size" not available in the registry data. Set as NoneType object')

    # get collections attributes
    attrs_obj = collections_obj_in.attrs

    # iterate over registry
    collections_obj_stack, collections_time_stack = None, None
    collections_obj_out, registry_obj_out, collections_obj_attrs, registry_obj_attrs = None, None, None, None
    reg_id_out, reg_loc_out, reg_lon_out, reg_lat_out, reg_cell_out, reg_row_size_out = [], [], [], [], [], []
    reg_time_start_out, reg_time_end_out, reg_time_reference_out = [], [], []
    message_lock_variable_out, message_lock_method_out, message_lock_method_null = False, False, False
    for reg_id, (reg_loc, reg_lon, reg_lat, reg_cell) in enumerate(
            zip(reg_index_list, reg_longitude_list, reg_latitude_list, reg_cell_list)):

        # check collections object type
        if isinstance(collections_obj_in, pd.DataFrame):

            # set indexes start and end
            reg_row_start, reg_row_end = reg_row_size_cumulative[reg_id], reg_row_size_cumulative[reg_id + 1]
            # select collections by rows
            collections_selected = collections_obj_in.iloc[reg_row_start:reg_row_end]
            # sort collections by time index
            collections_sorted = collections_selected.sort_index()
            # get row size
            reg_row_size = collections_sorted.shape[0]
            # get time start, end and reference
            reg_time_start, reg_time_end = collections_sorted.index[0], collections_sorted.index[-1]
            reg_time_reference = time_obj

        elif isinstance(collections_obj_in, xr.Dataset):

            # organize collections by stack
            if collections_obj_stack is None:
                collections_obj_stack, collections_time_stack = {}, None
                for var_name in list(collections_obj_in.variables):
                    if var_name != 'time':
                        var_obj_stack = collections_obj_in[var_name].values
                        collections_obj_stack[var_name] = var_obj_stack
                    else:
                        collections_time_stack = collections_obj_in[var_name].values

            if collections_time_stack is None:
                alg_logger.error(' ===> Time variable not available in the collections data')
                raise RuntimeError('Check your collections data object')

            # select collections by rows
            collections_data = {}
            for var_name in list(collections_obj_stack.keys()):
                var_data = collections_obj_stack[var_name]

                if var_data.ndim == 2:
                    tmp_data = var_data[reg_id, :]
                    collections_data[var_name] = tmp_data

            collections_data['time'] = collections_time_stack
            collections_selected = pd.DataFrame(data=collections_data, index=collections_time_stack)

            # sort collections by time index
            collections_sorted = collections_selected.sort_index()
            # get row size
            reg_row_size = collections_sorted.shape[0]
            # get time start, end and reference
            reg_time_start, reg_time_end = collections_sorted.index[0], collections_sorted.index[-1]
            reg_time_reference = time_obj

        else:
            alg_logger.error(' ===> Collections object must be a pandas DataFrame or a xarray Dataset')
            raise RuntimeError('Check your collections object')

        # compute data
        collections_data_out, collections_registry_out = None, None
        for key_var, var_name_in in fx_var_name_in.items():

            # get data in
            ts_data = collections_sorted[var_name_in]
            attrs_data = collections_obj_in[var_name_in].attrs

            # organize attrs
            if collections_obj_attrs is None:
                collections_obj_attrs = {}
            if var_name_in not in list(collections_obj_attrs.keys()):
                collections_obj_attrs[var_name_in] = attrs_data

            # check key in destination dictionary
            if key_var in list(fx_var_name_out.keys()):

                # check key in method dictionary
                if key_var in list(fx_var_methods.keys()):

                    # get method list
                    methods_obj = fx_var_methods[key_var]
                    # get variable out list
                    list_name_out = fx_var_name_out[key_var]

                    # check variable out list
                    if (list_name_out is not None) and (methods_obj is not None):

                        # set variable out as list
                        if not isinstance(list_name_out, list):
                            list_name_out = [list_name_out]

                        # workspace collections - iterate over method(s)
                        for var_name_out, (method_name, method_args) in zip(list_name_out, methods_obj.items()):
                            # copy data source
                            ts_tmp = deepcopy(ts_data)
                            # apply methods
                            values_out = compute_data_metrics(
                                ts_tmp, fx_method=method_name, fx_args=method_args)
                            # collect data destination
                            if collections_data_out is None:
                                collections_data_out = pd.DataFrame(data={var_name_out: values_out}, index=[reg_loc])
                            else:
                                collections_data_out[var_name_out] = values_out
                    else:
                        # message if variable(s) or/and method(s) are defined by NoneType object
                        if not message_lock_method_null:
                            alg_logger.warning(' ===> Variable(s) and/or method(s) for variable "' + key_var +
                                               '" defined by NoneType object')
                            message_lock_method_null = True
                else:
                    # message if methods are not defined
                    if not message_lock_method_out:
                        alg_logger.warning(' ===> Method(s) for variable "' + key_var + '" not available')
                        message_lock_method_out = True
            else:
                # message if variable out are not defined
                if not message_lock_variable_out:
                    alg_logger.warning(' ===> Variable(s) for key "' + key_var + '" not available')
                    message_lock_variable_out = True

        # workspace registry
        reg_row_size_out.append(int(reg_row_size))
        reg_id_out.append(int(reg_id))
        reg_loc_out.append(int(reg_loc))
        reg_lon_out.append(reg_lon)
        reg_lat_out.append(reg_lat)
        reg_cell_out.append(int(reg_cell))
        reg_time_start_out.append(reg_time_start)
        reg_time_end_out.append(reg_time_end)
        reg_time_reference_out.append(reg_time_reference)

        # collections workspace out
        if collections_obj_out is None:
            collections_interface = deepcopy(collections_data_out)
            collections_interface.reset_index()
            collections_obj_out = collections_interface
        else:
            collections_interface = deepcopy(collections_data_out)
            collections_interface.reset_index()
            collections_obj_out = pd.concat([collections_obj_out, collections_interface], axis=0)

    # organize registry dataframe
    registry_obj_out = {'gpi': reg_loc_out, 'location_id': reg_loc_out, 'lon': reg_lon_out, 'lat': reg_lat_out,
                        'cell': reg_cell_out, 'name': reg_id_out, 'row_size': reg_row_size_out}
    registry_obj_out = pd.DataFrame(registry_obj_out, index=reg_loc_out)

    # organize registry attributes
    registry_obj_attrs = {
        'lat': reg_latitude_attrs, 'lon': reg_longitude_attrs, 'cell': reg_cell_attrs, 'row_size': reg_row_size_attrs}
    # organize collections attributes
    collections_obj_out.attrs = attrs_obj

    return collections_obj_out, registry_obj_out, collections_obj_attrs, registry_obj_attrs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to dump datasets cell
def dump_datasets_cell(file_path,
                       data_obj, registry_obj,
                       data_attrs, registry_attrs,
                       list_variable_data_in, list_variable_data_out,
                       list_variable_registry_in, list_variable_registry_out):

    # change variable names (data and registry)
    for var_name_in, var_name_out in zip(list_variable_data_in, list_variable_data_out):
        if var_name_in in list(data_obj.columns):
            data_obj.rename(columns={var_name_in: var_name_out}, inplace=True)
        else:
            alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the data object')
        if var_name_in in list(data_attrs.keys()):
            if var_name_out not in list(data_attrs.keys()):
                data_attrs[var_name_out] = data_attrs[var_name_in]
                data_attrs.pop(var_name_in)

    for var_name_in, var_name_out in zip(list_variable_registry_in, list_variable_registry_out):
        if var_name_in in list(registry_obj.columns):
            registry_obj.rename(columns={var_name_in: var_name_out}, inplace=True)
        else:
            alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the registry object')

        if var_name_in in list(registry_attrs.keys()):
            if var_name_out not in list(registry_attrs.keys()):
                registry_attrs[var_name_out] = registry_attrs[var_name_in]
                registry_attrs.pop(var_name_in)

    # drop variable not in the list (data and registry)
    for var_name in list(data_obj.columns):
        if var_name not in list_variable_data_out:
            data_obj.drop(var_name, axis=1, inplace=True)
            if var_name in list(data_attrs.keys()):
                data_attrs.pop(var_name)

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
                    data_obj, registry_obj, data_attrs, registry_attrs, file_tag_location='location_id')

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to dump datasets grid
def dump_datasets_grid(file_path, file_grid_obj, file_grid_reference):

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
