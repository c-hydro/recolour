"""
Library Features:

Name:          lib_utils_datasets
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240507'
Version:       '1.5.0'
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
import lib_data_fx

from lib_data_io_nc import read_file_cell, write_file_cell
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

        collections_data, collections_index, collections_attrs, collections_dims = None, None, None, None
        for name_variable_data_out in list_variable_data_out:
            if name_variable_data_out == index_name_data:
                collections_index = cell_datasets_out[index_name_data]

            if collections_data is None:
                collections_data = {}
            collections_data[name_variable_data_out] = cell_datasets_out[name_variable_data_out]
            if collections_dims is None:
                collections_dims = collections_data[name_variable_data_out].ndim
            if collections_attrs is None:
                collections_attrs = {}
            collections_attrs[name_variable_data_out] = cell_attrs_out[name_variable_data_out]

        if collections_dims == 1:
            collections_data[index_name_cell] = [int(cell_name)] * len(collections_index)
        elif collections_dims == 2:
            collections_data[index_name_registry] = cell_datasets_out[index_name_registry]
        else:
            alg_logger.error(' ===> Dimension not allowed for collections object')
            raise RuntimeError('Check dimensions of your collections object')

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
        if collections_dims == 1:
            collections_obj = pd.DataFrame(data=collections_data, index=collections_index)
            for var_name in list(collections_obj.columns):
                if var_name in list(collections_attrs.keys()):
                    collections_obj[var_name].attrs = collections_attrs[var_name]
        elif collections_dims == 2:

            if index_name_data in list(collections_data.keys()):
                values_data = collections_data[index_name_data]
            else:
                alg_logger.error(' ===> Index name "' + index_name_data + '" not available in collections data')
                raise RuntimeError('Check your collections data object')

            if index_name_registry in list(collections_data.keys()):
                values_registry = collections_data[index_name_registry]
            else:
                alg_logger.error(' ===> Index name "' + index_name_registry + '" not available in collections data')
                raise RuntimeError('Check your collections data object')

            collections_obj = xr.Dataset()
            for var_name in list(collections_data.keys()):
                if var_name not in [index_name_data, index_name_registry]:
                    collections_obj[var_name] = xr.DataArray(
                        collections_data[var_name],
                        dims=[index_name_registry, index_name_data],
                        coords={
                            index_name_registry: (index_name_registry, values_registry),
                            index_name_data: (index_name_data, values_data)
                        }
                    )
                    collections_obj[var_name].attrs = collections_attrs[var_name]
                elif var_name == index_name_data:
                    collections_obj[var_name].attrs = collections_attrs[var_name]

        else:
            alg_logger.error(' ===> Dimension not allowed for collections object')
            raise RuntimeError('Check dimensions of your collections object')

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
# method to apply root zone soil moisture datasets cell
def apply_rzsm_datasets_cell(
        time_obj, collections_obj_in, registry_obj_in,
        fx_var_name_in=None, fx_var_name_out=None, fx_var_methods=None,
        index_name_time='time', index_name_locations='location_id',
        data_var_name=None, data_var_mode=None, data_var_min_value=None, data_var_max_value=None,
        data_var_scale_factor=None, data_var_undef=None
):

    if fx_var_name_in is None:
        fx_var_name_in = {'layer_1': ['var40']}
    if fx_var_name_out is None:
        fx_var_name_out = {'layer_1': ["var_0_7"]}

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

    # copy data and registry source
    collections_obj_out = deepcopy(collections_obj_in)
    registry_obj_out = deepcopy(registry_obj_in)

    if index_name_time in list(collections_obj_in.coords):
        values_time = collections_obj_in[index_name_time].values
        # attrs_time = collections_obj_in[index_name_time].attrs
    else:
        alg_logger.error(' ===> Index name "' + index_name_time + '" not available in collections data')
        raise RuntimeError('Check your collections data object')

    if index_name_locations in list(collections_obj_in.coords):
        values_locations = collections_obj_in[index_name_locations].values
        # attrs_locations = collections_obj_in[index_name_locations].attrs
    else:
        alg_logger.error(' ===> Index name "' + index_name_locations + '" not available in collections data')
        raise RuntimeError('Check your collections data object')

    # compute data
    collections_obj_attrs, registry_obj_attrs = None, None
    message_lock_variable_out, message_lock_method_out, message_lock_method_null = False, False, False
    for key_var, list_name_in in fx_var_name_in.items():

        # set variable out as list
        if not isinstance(list_name_in, list):
            list_name_in = [list_name_in]

        # filter data using min, max, scale factor and undef
        for var_name_in in list_name_in:
            var_mode, var_min_value, var_max_value, var_scale_factor, var_undef = None, None, None, None, None
            if var_name_in in data_var_name:
                var_idx = data_var_name.index(var_name_in)
                var_mode = data_var_mode[var_idx]
                var_min_value, var_max_value = data_var_min_value[var_idx], data_var_max_value[var_idx]
                var_scale_factor, var_undef = data_var_scale_factor[var_idx], data_var_undef[var_idx]

            if (var_mode is not None) and (var_mode == 'reference'):
                if var_undef is not None:
                    collections_obj_in = collections_obj_in.where(
                        collections_obj_in[var_name_in] != var_undef, np.nan)
                if var_min_value is not None:
                    collections_obj_in = collections_obj_in.where(
                        collections_obj_in[var_name_in] > var_min_value, np.nan)
                if var_max_value is not None:
                    collections_obj_in = collections_obj_in.where(
                        collections_obj_in[var_name_in] < var_max_value, np.nan)

            if var_scale_factor is not None:
                collections_obj_in[var_name_in] = collections_obj_in[var_name_in] * var_scale_factor


        # get data in
        obj_data = collections_obj_in[list_name_in]

        # organize attrs
        if collections_obj_attrs is None:
            collections_obj_attrs = {}
        for var_name_in in list_name_in:
            if var_name_in in list(collections_obj_in.data_vars) + list(collections_obj_in.dims):
                if var_name_in not in list(collections_obj_attrs.keys()):
                    collections_obj_attrs[var_name_in] = collections_obj_in[var_name_in].attrs
            else:
                collections_obj_attrs[var_name_in] = {}

        # check key in destination dictionary
        if key_var in list(fx_var_name_out.keys()):

            # check key in method dictionary
            if key_var in list(fx_var_methods.keys()):

                # get method list
                methods_obj = fx_var_methods[key_var]
                # get variable out list
                list_name_out = fx_var_name_out[key_var]

                # check variable out and methods list
                if (list_name_out is not None) and (methods_obj is not None):

                    # set variable out as list
                    if not isinstance(list_name_out, list):
                        list_name_out = [list_name_out]

                    # workspace collections - iterate over method(s)
                    for var_name_out, (method_name, method_args) in zip(list_name_out, methods_obj.items()):

                        # copy data source
                        tmp_data = deepcopy(obj_data)

                        # set method obj for preparing variable(s)
                        if hasattr(lib_data_fx, 'prepare_data'):
                            var_layers_data = lib_data_fx.prepare_data(tmp_data, list_name_in)
                        else:
                            alg_logger.error(' ===> Method "prepare_data" not available in the library')
                            raise RuntimeError('Check your method object')

                        # set method obj for computing variable(s)
                        if hasattr(lib_data_fx, method_name):
                            compute_obj_fx = getattr(lib_data_fx, method_name)
                            values_out = compute_obj_fx(**var_layers_data, **method_args)
                        else:
                            alg_logger.error(' ===> Method "' + method_name + '" not available in the library')
                            raise RuntimeError('Check your method object')

                        # store values in datasets object
                        collections_obj_out[var_name_out] = xr.DataArray(
                            values_out,
                            dims=[index_name_locations, index_name_time],
                            coords={
                                index_name_locations: (index_name_locations, values_locations),
                                index_name_time: (index_name_time, values_time)
                            }
                        )

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

    # organize registry dataframe
    registry_obj_attrs = {
        'lat': reg_latitude_attrs, 'lon': reg_longitude_attrs, 'cell': reg_cell_attrs}

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

        if isinstance(data_obj, pd.DataFrame):
            if var_name_in in list(data_obj.columns):
                data_obj.rename(columns={var_name_in: var_name_out}, inplace=True)
            else:
                alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the data object')
        elif isinstance(data_obj, xr.Dataset):
            if var_name_in in list(data_obj.data_vars):
                data_obj = data_obj.rename({var_name_in: var_name_out})
            else:
                alg_logger.warning(' ===> Variable "' + var_name_in + '" not available in the data object')
        else:
            alg_logger.error(' ===> Data object type not allowed')
            raise RuntimeError('Check your data object type')
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
        alg_logger.error(' ===> Data object type not allowed')
        raise RuntimeError('Check your data object type')

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
