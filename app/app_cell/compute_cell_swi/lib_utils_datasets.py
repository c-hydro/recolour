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

from copy import deepcopy

from pytesmo.validation_framework.data_manager import get_result_names
from lib_data_io_nc import read_file_cell, write_file_cell
from lib_data_analysis import compute_data_swi
from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# message(s) suppressed in the console
pd.options.mode.chained_assignment = None
# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# method to get committed datasets cells
def get_datasets_committed_OLD(df_globals, field='committed_area'):

    comm = df_globals[field] == 1
    df_committed = df_globals[comm]
    df_committed['type'] = [field] * len(df_committed)

    return df_committed
# -----------------------------------------------------------------------------


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
    collections_dframe, registry_dframe = None, None
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
        collections_dframe = pd.DataFrame(data=collections_data, index=collections_index)
        for var_name in list(collections_dframe.columns):
            if var_name in list(collections_attrs.keys()):
                collections_dframe[var_name].attrs = collections_attrs[var_name]
        # organize registry
        registry_dframe = pd.DataFrame(data=registry_data, index=registry_index)
        for var_name in list(registry_dframe.columns):
            if var_name in list(registry_attrs.keys()):
                registry_dframe[var_name].attrs = registry_attrs[var_name]

    else:

        # info warning file not found
        alg_logger.warning(' ===> File cell "' + file_name + '" not found')

    return collections_dframe, registry_dframe

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter datasets cell
def filter_datasets_cell(
        time_obj, collections_obj_in, registry_obj_in,
        var_name_data_in=None, var_name_data_out=None,
        ctime_data=6):

    if var_name_data_out is None:
        var_name_data_out = {'swi_data': ['swi_values_t06']}
    if var_name_data_in is None:
        var_name_data_in = {'swi_data': "sm"}

    if not isinstance(var_name_data_in, dict):
        alg_logger.error(' ===> Variable name for source data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    if not isinstance(var_name_data_out, dict):
        alg_logger.error(' ===> Variable name for destination data must be a dictionary')
        raise RuntimeError('Check your filter variable object')
    if not isinstance(ctime_data, list):
        ctime_data = [ctime_data]

    # get registry info
    reg_index_in = registry_obj_in.index.values
    reg_latitude_in, reg_latitude_attrs = registry_obj_in['lat'].values, registry_obj_in['lat'].attrs
    reg_longitude_in, reg_longitude_attrs = registry_obj_in['lon'].values, registry_obj_in['lon'].attrs
    reg_cell_in, reg_cell_attrs = registry_obj_in['cell'].values, registry_obj_in['cell'].attrs
    reg_row_size_in, reg_row_size_attrs = registry_obj_in['row_size'].values, registry_obj_in['row_size'].attrs

    # organize row size
    row_size_cumulative_in = np.append(0, np.cumsum(reg_row_size_in))

    # get collections attributes
    attrs_obj_in = collections_obj_in.attrs

    # iterate over registry
    collections_obj_out, registry_obj_out, collections_obj_attrs, registry_obj_attrs = None, None, None, None
    reg_id_out, reg_loc_out, reg_lon_out, reg_lat_out, reg_cell_out, reg_row_size_out = [], [], [], [], [], []
    reg_time_start_out, reg_time_end_out, reg_time_reference_out = [], [], []
    for reg_id, (reg_loc, reg_lon, reg_lat, reg_cell) in enumerate(
            zip(reg_index_in, reg_longitude_in, reg_latitude_in, reg_cell_in)):

        # set indexes start and end
        row_size_start_in, row_size_end_in = row_size_cumulative_in[reg_id], row_size_cumulative_in[reg_id + 1]
        # select collections by rows
        collections_selected_in = collections_obj_in.iloc[row_size_start_in:row_size_end_in]
        # sort collections by time index
        collections_sorted_in = collections_selected_in.sort_index()

        # compute data
        df_data_out, attrs_data_out = pd.DataFrame(), {}
        for key_var_in, var_name_in in var_name_data_in.items():

            # get data in
            ts_data_in = collections_sorted_in[var_name_in]
            attrs_data_in = collections_obj_in[var_name_in].attrs

            # organize attrs
            if collections_obj_attrs is None:
                collections_obj_attrs = {}
            if var_name_in not in list(collections_obj_attrs.keys()):
                collections_obj_attrs[var_name_in] = attrs_data_in

            # compute data out
            if key_var_in in list(var_name_data_out.keys()):

                # get variable out list
                list_name_out = var_name_data_out[key_var_in]
                # check variable out list
                if list_name_out is not None:

                    if not isinstance(list_name_out, list):
                        list_name_out = [list_name_out]

                    # iterate to compute swi output
                    for var_name_out, ctime_value in zip(list_name_out, ctime_data):

                        ts_data_tmp = deepcopy(ts_data_in)
                        ts_data_out = compute_data_swi(ts_data_tmp, ts_var=var_name_out, ts_ctime=ctime_value)

                        df_data_out[var_name_out] = ts_data_out

                else:
                    # save data in the common dataframe
                    df_data_out[var_name_in] = ts_data_in
            else:
                # save data in the common dataframe
                df_data_out[var_name_in] = ts_data_in

        # get time start, end and reference
        time_start_out, time_end_out, time_reference_out = df_data_out.index[0], df_data_out.index[-1], time_obj

        # collections workspace out
        if collections_obj_out is None:
            collections_interface = df_data_out
            collections_interface.reset_index()
            collections_obj_out = collections_interface
        else:
            collections_interface = df_data_out
            collections_interface.reset_index()
            collections_obj_out = pd.concat([collections_obj_out, collections_interface], axis=0)

        # workspace registry info out
        reg_row_size_out.append(int(df_data_out.shape[0]))
        reg_id_out.append(int(reg_id))
        reg_loc_out.append(int(reg_loc))
        reg_lon_out.append(reg_lon)
        reg_lat_out.append(reg_lat)
        reg_cell_out.append(int(reg_cell))
        # workspace time info out
        reg_time_start_out.append(time_start_out)
        reg_time_end_out.append(time_end_out)
        reg_time_reference_out.append(time_reference_out)

    # organize registry dataframe
    registry_obj_out = {'gpi': reg_loc_out, 'location_id': reg_loc_out, 'lon': reg_lon_out, 'lat': reg_lat_out,
                        'cell': reg_cell_out, 'time_start': reg_time_start_out, 'time_end': reg_time_end_out,
                        'time_reference': reg_time_reference_out, 'name': reg_id_out, 'row_size': reg_row_size_out}
    registry_obj_out = pd.DataFrame(registry_obj_out, index=reg_loc_out)

    registry_obj_attrs = {
        'lat': reg_latitude_attrs, 'lon': reg_longitude_attrs, 'cell': reg_cell_attrs, 'row_size': reg_row_size_attrs}

    # organize collections dataframe
    collections_obj_out.reset_index(inplace=True)
    collections_obj_out.drop('index', axis=1, inplace=True)
    collections_obj_out.attrs = attrs_obj_in

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


# -----------------------------------------------------------------------------
# method to convert datasets obj
def convert_datasets_obj_OLD(data_dict):
    data_dframe= pd.DataFrame.from_dict(data_dict)
    data_dframe['type'] = ['global'] * len(data_dframe)
    return data_dframe
# -----------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get dataset mode(s)
def get_dataset_modes_OLD(dset_obj, dset_mode='reference'):
    dset_name = []
    for dset_key, dset_fields in dset_obj.items():
        if 'type' in list(dset_fields.keys()):
            tmp_mode = dset_fields['type']
            if tmp_mode == dset_mode:
                dset_name.append(dset_key)
    if dset_name.__len__() == 1:
        dset_name = dset_name[0]
    return dset_name
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get dataset name(s)
def get_dataset_names_OLD(ref_key, datasets):
    """
    Get dataset names in correct order as used in the validation framework
        -) reference dataset = ref
        -) first other dataset = k1
        -) second other dataset = k2
    This is important to correctly iterate through the H-SAF metrics and to
    save each metric with the name of the used datasets

    Parameters
    ----------
    ref_key: basestring
        Name of the reference dataset
    datasets: dict
        Dictionary of dictionaries as provided to the validation framework
        in order to perform the validation process.

    Returns
    -------
    dataset_names: list
        List of the dataset names in correct order

    """
    ds_dict = {}
    for ds in datasets.keys():
        ds_dict[ds] = datasets[ds]['columns']
    ds_names = get_result_names(ds_dict, ref_key, n=3)
    dataset_names = []
    for name in ds_names[0]:
        dataset_names.append(name[0])

    return dataset_names

# ----------------------------------------------------------------------------------------------------------------------
