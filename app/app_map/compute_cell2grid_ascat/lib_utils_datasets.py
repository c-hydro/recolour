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

from copy import deepcopy

from pytesmo.validation_framework.data_manager import get_result_names

from lib_info_args import logger_name
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
        file_name, cell_name,
        list_variable_data_in=None, list_variable_data_out=None,
        list_variable_registry_in=None, list_variable_registry_out=None,
        index_name_data='time', index_name_registry='location', index_name_cell='cell'):

    # info datasets start
    alg_logger.info(' ------> Get datasets cell ... ')

    # common list variable in/out
    list_variable_in = list_variable_data_in + list_variable_registry_in
    list_variable_out = list_variable_data_out + list_variable_registry_out

    # check file name availability
    collections_dframe, registry_dframe = None, None
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

        collections_data, collections_index = None, None
        for name_variable_data_out in list_variable_data_out:
            if name_variable_data_out == index_name_data:
                collections_index = cell_datasets_out[index_name_data]
            else:
                if collections_data is None:
                    collections_data = {}
                collections_data[name_variable_data_out] = cell_datasets_out[name_variable_data_out]
        collections_data[index_name_cell] = [int(cell_name)] * len(collections_index)

        registry_data, registry_index = None, None
        for name_variable_registry_out in list_variable_registry_out:
            if name_variable_registry_out == index_name_registry:
                registry_index = cell_datasets_out[index_name_registry]
            else:
                if registry_data is None:
                    registry_data = {}
                registry_data[name_variable_registry_out] = cell_datasets_out[name_variable_registry_out]

        # organize in dataframes
        collections_dframe = pd.DataFrame(data=collections_data, index=collections_index)
        collections_dframe.attrs = attrs_datasets_in

        registry_dframe = pd.DataFrame(data=registry_data, index=registry_index)

        # info datasets start
        alg_logger.info(' ------> Get datasets cell ... DONE')

    else:

        # info datasets end
        alg_logger.warning(' ===> File cell "' + file_name + '" not found')
        alg_logger.info(' ------> Get datasets cell ... SKIPPED')

    return collections_dframe, registry_dframe

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to add datasets grid
def organize_datasets_points(
        time_obj, collections_obj, registry_obj, data_type='ASCAT', grid_dframe=None,
        collections_fx='mean'):

    # info datasets start
    alg_logger.info(' ------> Convert datasets cell to points ... ')

    # get registry info
    registry_index = registry_obj.index
    registry_latitude = registry_obj['latitude']
    registry_longitude = registry_obj['longitude']
    registry_row_size = registry_obj['row_size']

    # organize time start and end
    time_start = time_obj.replace(hour=0, minute=0, second=0)
    time_end = time_obj.replace(hour=23, minute=59, second=59)
    # organize row size
    row_size_cumulative = np.append(0, np.cumsum(registry_row_size))

    # get collections attributes
    attrs_obj = collections_obj.attrs

    # iterate over registry
    alg_logger.info(' -------> Get collections ... ')
    collections_workspace = None
    for reg_id, (reg_loc, reg_lon, reg_lat) in enumerate(zip(registry_index, registry_longitude, registry_latitude)):

        # set indexes start and end
        row_size_start, row_size_end = row_size_cumulative[reg_id], row_size_cumulative[reg_id + 1]

        # select collections by rows
        collections_selected = collections_obj.iloc[row_size_start:row_size_end]
        # sort collections by time index
        collections_sorted = collections_selected.sort_index()

        # mask collections by time period
        time_mask = (collections_sorted.index > time_start) & (collections_sorted.index <= time_end)
        collections_masked = collections_sorted.loc[time_mask]

        # check collections availability
        if collections_masked.shape[0] > 0:

            collections_tmp = deepcopy(collections_masked)
            if collections_masked.shape[0] > 1:

                if collections_fx == 'mean':

                    collections_defined = collections_tmp.mean()
                    collections_defined.columns = list(collections_tmp.columns)
                    collections_index = collections_tmp.index.mean()
                    collections_defined['time'] = collections_index
                    collections_defined['fx'] = collections_fx

                elif collections_fx == 'max':

                    collections_defined = collections_tmp.max()
                    collections_defined.columns = list(collections_tmp.columns)
                    collections_index = collections_tmp.index.max()
                    collections_defined['time'] = collections_index
                    collections_defined['fx'] = collections_fx

                elif collections_fx == 'min':

                    collections_defined = collections_tmp.min()
                    collections_defined.columns = list(collections_tmp.columns)
                    collections_index = collections_tmp.index.min()
                    collections_defined['time'] = collections_index
                    collections_defined['fx'] = collections_fx

                elif collections_fx == 'first':

                    collections_defined = collections_tmp.iloc[0]
                    collections_defined.columns = list(collections_tmp.columns)
                    collections_index = collections_tmp.index[0]
                    collections_defined['time'] = collections_index
                    collections_defined['fx'] = collections_fx

                elif collections_fx == 'last':

                    collections_defined = collections_tmp.iloc[-1]
                    collections_defined.columns = list(collections_tmp.columns)
                    collections_index = collections_tmp.index[-1]
                    collections_defined['time'] = collections_index
                    collections_defined['fx'] = collections_fx

                else:
                    alg_logger.error(' ===> Collections function "' + collections_fx + '" not allowed')
                    raise NotImplementedError('Case not implemented yet')
            else:
                # add time info
                collections_defined = collections_tmp.max()
                collections_defined.columns = list(collections_tmp.columns)
                collections_index = collections_tmp.index[0]
                collections_defined['time'] = collections_index
                collections_defined['fx'] = 'single'

            # add registry info
            collections_defined['location'] = np.array(reg_loc, dtype=int)
            collections_defined['longitude'] = reg_lon
            collections_defined['latitude'] = reg_lat
            collections_defined['time_start'] = time_start
            collections_defined['time_end'] = time_end
            collections_defined['time_reference'] = time_obj
            collections_defined.name = reg_id

            if collections_workspace is None:
                collections_interface = collections_defined.to_frame()
                collections_interface.reset_index()
                collections_workspace = collections_interface
            else:
                collections_interface = collections_defined.to_frame()
                collections_interface.reset_index()
                collections_workspace = pd.concat([collections_workspace, collections_interface], axis=1)

        else:
            alg_logger.warning(' ===> Collections for location "' + str(reg_loc) + '" are not available')

    alg_logger.info(' -------> Get collections ... DONE')

    # organize workspace (transpose the dataframe)
    alg_logger.info(' -------> Order collections ... ')
    collections_workspace = collections_workspace.T
    collections_workspace.attrs = attrs_obj
    alg_logger.info(' -------> Order collections ... DONE')

    # info datasets end
    alg_logger.info(' ------> Convert datasets cell to points ... DONE')

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
