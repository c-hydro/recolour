"""
Library Features:

Name:          lib_cell_tools
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy

from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# message(s) suppressed in the console
pd.options.mode.chained_assignment = None
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extract cell data
def extract_cell_data(collections_cell, var_list_expected=None):

    if isinstance(collections_cell, pd.DataFrame):
        var_list = list(collections_cell.columns)
        coord_list = collections_cell.index.name

        tmp_list = deepcopy(var_list)
        if coord_list is None:
            coord_list = []
            for var_name in var_list:
                if 'time' in var_name:
                    coord_list.append(var_name)
                    tmp_list.remove(var_name)
            var_list = deepcopy(tmp_list)

        if var_list_expected is None:
            var_list_expected = deepcopy(var_list)

    elif isinstance(collections_cell, xr.Dataset):
        var_list = list(collections_cell.data_vars)
        coord_list = list(collections_cell.coords)

        if var_list_expected is None:
            var_list_expected = deepcopy(var_list)
    else:
        alg_logger.error(' ===> Collections object must be a pandas DataFrame or a xarray Dataset')
        raise RuntimeError('Check your collections object')

    var_data_dict, coord_data_dict = None, None
    for var_name in var_list:

        if var_name in var_list_expected:
            var_data = collections_cell[var_name].values

            if var_data_dict is None:
                var_data_dict = {}
            var_data_dict[var_name] = var_data

    if coord_list is not None:
        for coord_name in coord_list:
            coord_data = collections_cell[coord_name].values

            if coord_data_dict is None:
                coord_data_dict = {}
            coord_data_dict[coord_name] = coord_data

    return var_data_dict, coord_data_dict
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extract cell attributes
def extract_cell_attrs(collections_cell):

    if isinstance(collections_cell, pd.DataFrame):
        var_list = list(collections_cell.columns)
    elif isinstance(collections_cell, xr.Dataset):
        var_list = list(collections_cell.data_vars)
    else:
        alg_logger.error(' ===> Collections object must be a pandas DataFrame or a xarray Dataset')
        raise RuntimeError('Check your collections object')

    global_attrs_dict = collections_cell.attrs

    var_attrs_dict = None
    for var_name in var_list:
        var_attrs = collections_cell[var_name].attrs

        if (not var_attrs) and (var_name in list(global_attrs_dict.keys())):
            pass
        else:
            if var_attrs_dict is None:
                var_attrs_dict = {}
            var_attrs_dict[var_name] = var_attrs

    if global_attrs_dict is None:
        global_attrs_dict = {}
    if var_attrs_dict is None:
        var_attrs_dict = {}

    collections_attrs_dict = {**global_attrs_dict, **var_attrs_dict}

    return collections_attrs_dict
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extract registry information
def extract_registry_info(
        registry_obj, index_name_geo_x='lon', index_name_geo_y='lat',
        index_name_cell='cell', index_name_row_size='row_size', n_start=None, n_end=None,
        flag_type='registry'):

    # get registry reference information
    index_list = registry_obj.index.values
    lat_list, lat_attrs = registry_obj[index_name_geo_y].values, registry_obj[index_name_geo_y].attrs
    lon_list, lon_attrs = registry_obj[index_name_geo_x].values, registry_obj[index_name_geo_x].attrs
    cell_list, cell_attrs = registry_obj[index_name_cell].values, registry_obj[index_name_cell].attrs
    if index_name_row_size in list(registry_obj.columns):
        row_size_n_list = registry_obj[index_name_row_size].values
        row_size_n_attrs = registry_obj[index_name_row_size].attrs
        # organize row size
        row_size_cumulative_list = np.append(0, np.cumsum(row_size_n_list))
    else:
        row_size_n_list, row_size_n_attrs, row_size_cumulative_list = None, None, None
        alg_logger.warning(' ===> Variable "row_size" not available in the "' +
                           flag_type + '" data. Set as NoneType object')

    # select data for debugging mode (subset of data)
    if (n_start is not None) and (n_end is not None):
        alg_logger.warning(' ===> Selecting data from index ' + str(n_start) + ' to index ' + str(n_end) + ' ...')
        index_list = index_list[n_start:n_end]
        lat_list = lat_list[n_start:n_end]
        lon_list = lon_list[n_start:n_end]
        cell_list = cell_list[n_start:n_end]
        if row_size_n_list is not None:
            row_size_n_list = row_size_n_list[n_start:n_end]
            row_size_cumulative_list = row_size_cumulative_list[n_start:n_end + 1]
        alg_logger.warning(' ===> Selecting data from index ' +
                           str(n_start) + ' to index ' + str(n_end) + ' ... DONE. Usage of data for debugging mode')

    return index_list, lat_list, lon_list, cell_list, row_size_n_list, row_size_cumulative_list
# ----------------------------------------------------------------------------------------------------------------------
