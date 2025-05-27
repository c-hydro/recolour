"""
Library Features:

Name:          lib_utils_datasets
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""


# -----------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy

from pytesmo.validation_framework.data_manager import get_result_names

from lib_info_args import logger_name
from lib_utils_io import create_darray_2d
from lib_utils_grid import read_grid_file
from lib_data_io_nc import read_file_cell, read_file_collection, write_file_collection
from lib_data_io_pickle import write_file_obj, read_file_obj

# set logger
alg_logger = logging.getLogger(logger_name)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# default data type(s)
data_type_expected = ['ASCAT', 'ECMWF', 'HMC']
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to get committed datasets cells
def get_datasets_committed(df_globals, field='committed_area'):

    comm = df_globals[field] == 1
    df_committed = df_globals[comm]
    df_committed['type'] = [field] * len(df_committed)

    return df_committed
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to organize datasets cells
def organize_datasets_cell(
        file_name, file_time, cell_name,
        time_frequency_left='D', time_frequency_right='D', time_window_left=5, time_window_right=5,
        list_variable_data_in=None, list_variable_data_out=None,
        list_variable_registry_in=None, list_variable_registry_out=None,
        index_name_data='time', index_name_registry='location', index_name_cell='cell',
        index_name_row_size='row_size'):

    # info datasets start
    alg_logger.info(' ------> Get datasets ... ')

    # if row_index is the first variable in the list, move it to the end
    if index_name_row_size == list_variable_registry_in[0]:
        row_size_tmp_in, row_size_tmp_out = list_variable_registry_in[0], list_variable_registry_out[0]
        list_variable_registry_in = list_variable_registry_in[1:] + [row_size_tmp_in]
        list_variable_registry_out = list_variable_registry_out[1:] + [row_size_tmp_out]

    # common list variable in/out
    list_variable_in = list_variable_data_in + list_variable_registry_in
    list_variable_out = list_variable_data_out + list_variable_registry_out

    # check file name availability
    file_time_mask = None, None
    collections_dset, registry_dframe = None, None
    if os.path.exists(file_name):

        # get datasets in
        cell_datasets_in, attrs_datasets_in = read_file_cell(file_name, selected_variables=list_variable_in)

        # iterate over variables
        cell_datasets_out = {}
        for var_name_in, var_name_out in zip(list_variable_in, list_variable_out):

            if var_name_in in list(cell_datasets_in.keys()):
                var_data_in = cell_datasets_in[var_name_in]

                if isinstance(var_data_in, np.ma.MaskedArray):
                    var_data_in = var_data_in.data

                if var_name_in == index_name_data:
                    var_data_in = pd.to_datetime(var_data_in).floor('S')
                    file_time_from = pd.date_range(end=file_time, periods=time_window_left, freq=time_frequency_left)
                    file_time_to = pd.date_range(start=file_time, periods=time_window_right, freq=time_frequency_right)
                    file_time_mask = (var_data_in >= file_time_from[0]) & (var_data_in <= file_time_to[-1])

                # Save variable(s) defined in list (or save all variable(s)
                if var_name_out not in cell_datasets_out:
                    cell_datasets_out[var_name_out] = {}
                    cell_datasets_out[var_name_out] = var_data_in
                elif var_name_out in cell_datasets_out:
                    var_data_tmp = cell_datasets_out[var_name_out]
                    var_data_tmp = np.concatenate([var_data_tmp, var_data_in])
                    cell_datasets_out[var_name_out] = var_data_tmp

            else:

                alg_logger.warning(' ===> Variable "' + var_name_in + '" not available ih the cell datasets')

        # check time mask availability
        if file_time_mask is None:
            alg_logger.error(' ===> Time mask is not defined')
            raise NotImplementedError('Case not implemented yet')

        # data collections
        select_index_out = None
        collections_data, collections_index = None, None
        for name_variable_data_out in list_variable_data_out:
            if name_variable_data_out == index_name_data:
                tmp_index_out = cell_datasets_out[index_name_data]
                select_index_out = tmp_index_out[file_time_mask]
                collections_index = select_index_out
            else:
                if collections_data is None:
                    collections_data = {}

                tmp_datasets_out = cell_datasets_out[name_variable_data_out]
                if tmp_datasets_out.ndim == 1:
                    select_datasets_out = tmp_datasets_out[file_time_mask]
                elif tmp_datasets_out.ndim == 2:
                    select_datasets_out = tmp_datasets_out[:, file_time_mask]
                else:
                    alg_logger.error(' ===> Dimension "' + str(tmp_datasets_out.ndim) + '" not allowed')
                    raise NotImplementedError('Case not implemented yet')
                collections_data[name_variable_data_out] = select_datasets_out

        # registry collections
        registry_data, registry_index = None, None
        for name_variable_registry_out in list_variable_registry_out:
            if name_variable_registry_out == index_name_registry:
                registry_index = cell_datasets_out[index_name_registry]
            else:
                if registry_data is None:
                    registry_data = {}
                if name_variable_registry_out in cell_datasets_out:
                    registry_data[name_variable_registry_out] = cell_datasets_out[name_variable_registry_out]
                else:
                    if name_variable_registry_out == index_name_row_size:
                        pass
                    else:
                        alg_logger.error(' ===> Variable "' + name_variable_registry_out +
                                         '" not available in the cell datasets')
                        raise RuntimeError('Variable must be included in the cell datasets')

        if index_name_cell not in list(registry_data.keys()):
            registry_data[index_name_cell] = [int(cell_name)] * registry_index.shape[0]
        if index_name_row_size not in list(registry_data.keys()):
            registry_data[index_name_row_size] = [int(select_index_out.shape[0])] * registry_index.shape[0]

        # organize in data object
        collections_dset = xr.Dataset()
        for name_variable_data_out in list_variable_data_out:
            if name_variable_data_out != index_name_data:
                variable_values = collections_data[name_variable_data_out]
                variable_da = create_darray_2d(
                    variable_values,  collections_index, registry_index,
                    coord_name_x='time', coord_name_y='gpi',
                    dim_name_x='time', dim_name_y='gpi')
                collections_dset[name_variable_data_out] = variable_da

        collections_dset.attrs = attrs_datasets_in

        # drop duplicates
        _, index = np.unique(collections_dset['time'], return_index=True)
        collections_dset = collections_dset.isel(time=index)

        # organize registry object
        registry_dframe = pd.DataFrame(data=registry_data, index=registry_index)

        # info datasets start
        alg_logger.info(' ------> Get datasets ... DONE')

    else:

        # info datasets end
        alg_logger.warning(' ===> File cell "' + file_name + '" not found')
        alg_logger.info(' ------> Get datasets ... SKIPPED')

    return collections_dset, registry_dframe

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to add datasets grid
def organize_datasets_points(
        time_obj, collections_obj, registry_obj,
        collections_pivot='sm_values', collections_fx='nearest'):

    # info datasets start
    alg_logger.info(' ------> Convert datasets  ... ')

    # get registry info
    registry_index = registry_obj.index
    registry_latitude = registry_obj['latitude']
    registry_longitude = registry_obj['longitude']
    registry_cell = registry_obj['cell']
    registry_row_size = registry_obj['row_size']

    # organize time start and end
    time_start = time_obj.replace(hour=0, minute=0, second=0)
    time_end = time_obj.replace(hour=23, minute=59, second=59)

    # get collections attributes
    attrs_obj = collections_obj.attrs

    # iterate over registry
    alg_logger.info(' -------> Get points ... ')
    collections_workspace = None
    for reg_id, (reg_loc, reg_lon, reg_lat, reg_cell) in enumerate(zip(
            registry_index, registry_longitude, registry_latitude, registry_cell)):

        # select collections by gpis
        collections_masked = collections_obj.loc[dict(gpi=reg_loc, time=slice(time_start, time_end))]

        # check collections availability
        if collections_masked.dims['time'] > 0:

            collections_tmp = deepcopy(collections_masked)

            if collections_fx == 'mean':
                collections_dset = collections_tmp.mean()
                collections_time = collections_tmp['time'].mean().values
            elif collections_fx == 'max':
                collections_dset = collections_tmp.max()
                collections_time = collections_tmp[collections_pivot].idxmax(dim='time').values
            elif collections_fx == 'nearest':
                collections_dset = collections_tmp.sel(time=time_obj, method='nearest')
                collections_time = collections_dset['time'].values
            elif collections_fx == 'min':
                collections_dset = collections_tmp.min()
                collections_time = collections_tmp[collections_pivot].idxmin(dim='time').values
            elif collections_fx == 'first':
                collections_dset = collections_tmp.sel(time=collections_tmp['time'][0])
                collections_time = collections_dset['time'].values
            elif collections_fx == 'last':
                collections_dset = collections_tmp.sel(time=collections_tmp['time'][-1])
                collections_time = collections_dset['time'].values
            else:
                alg_logger.error(' ===> Collections function "' + collections_fx + '" not allowed')
                raise NotImplementedError('Case not implemented yet')

            # organize collections data
            collections_vars = list(collections_dset.data_vars)

            collections_data = {}
            for var in collections_vars:
                values = collections_dset[var].values
                collections_data[var] = values

            collections_index = pd.DatetimeIndex([collections_time])

            collections_df = pd.DataFrame(data=collections_data, index=collections_index)
            collections_df['time'] = collections_time
            collections_df['fx'] = collections_fx

            # organize collections registry
            collections_df['gpi'] = np.array(reg_loc, dtype=int)
            collections_df['location'] = np.array(reg_loc, dtype=int)
            collections_df['longitude'] = reg_lon
            collections_df['latitude'] = reg_lat
            collections_df['cell'] = np.array(reg_cell, dtype=int)
            collections_df['time_start'] = time_start
            collections_df['time_end'] = time_end
            collections_df['time_reference'] = time_obj
            collections_df['pivot'] = collections_pivot
            collections_df.name = reg_id

            # merge collections
            if collections_workspace is None:
                collections_interface = deepcopy(collections_df)
                collections_interface.reset_index()
                collections_workspace = collections_interface
            else:
                collections_interface = deepcopy(collections_df)
                collections_interface.reset_index()
                collections_workspace = pd.concat([collections_workspace, collections_interface], axis=0)

        else:
            alg_logger.warning(' ===> Collections for location "' + str(reg_loc) + '" are not available')

    alg_logger.info(' -------> Get points ... DONE')

    # organize workspace (transpose the dataframe)
    alg_logger.info(' -------> Organize points... ')
    if collections_workspace is not None:
        collections_workspace.reset_index(inplace=True)
        collections_workspace.drop('index', inplace=True, axis=1)
        collections_workspace.set_index('gpi', inplace=True)
        collections_workspace.index.name = 'gpi'
        collections_workspace.attrs = attrs_obj
        alg_logger.info(' -------> Organize points ... DONE')
    else:
        alg_logger.warning(' ===> Points datasets are not available.')
        alg_logger.info(' -------> Organize points ... SKIPPED. Collections is defined by NoneType')

    # info datasets end
    alg_logger.info(' ------> Convert datasets ... DONE')

    return collections_workspace
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to convert datasets obj
def convert_datasets_obj(data_dict):
    data_dframe= pd.DataFrame.from_dict(data_dict)
    data_dframe['type'] = ['global'] * len(data_dframe)
    return data_dframe
# -----------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get dataset mode(s)
def get_dataset_modes(dset_obj, dset_mode='reference'):
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
def get_dataset_names(ref_key, datasets):
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
