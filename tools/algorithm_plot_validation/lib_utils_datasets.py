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
import netCDF4

from copy import deepcopy

from pytesmo.validation_framework.data_manager import get_result_names

from lib_utils_grid import read_grid_file
from lib_data_io_nc import read_file_cell, read_file_collection, write_file_collection
from lib_data_io_pickle import write_file_obj, read_file_obj
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# default data type(s)
data_type_expected = ['ASCAT', 'RZSM', 'ECMWF', 'HMC']
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
        cell_list=None, cell_digits=4, cell_grid=None,
        list_variable_in=None, list_variable_out=None,
        folder_name_datasets='', file_name_datasets='{cell}.nc',
        variable_not_found='warning', active_empty_cell=False, active_fill_cell=True):

    # info script start
    logging.info(' ----> Get cell datasets ... ')

    # define variables (for undefined case)
    if list_variable_in is None:
        list_variable_in = ['gpi', 'lon', 'lat', '{:}_R']
    if list_variable_out is None:
        list_variable_out = ['gpi', 'lon', 'lat', '{:}_R']

    # Iterate over cell(s)
    var_list = None
    cell_datasets_out = None
    for cell_n in cell_list:

        # format cell
        cell_n = int(cell_n)
        cell_string = str(cell_n).zfill(cell_digits)

        # info start cell
        logging.info(' ----> Cell "' + cell_string + '" ... ')

        file_path_datasets = os.path.join(folder_name_datasets, file_name_datasets)
        file_path_datasets = file_path_datasets.format(cell=cell_string)

        if os.path.exists(file_path_datasets) or active_fill_cell:
            # init datasets out
            if cell_datasets_out is None:
                cell_datasets_out = {}

            # get datasets in
            if os.path.exists(file_path_datasets):
                cell_datasets_in = read_file_cell(file_path_datasets, file_variables=list_variable_in)
                if var_list is None:
                    var_list = list(cell_datasets_in.keys())
            else:
                if var_list is not None:
                    gpis, lons, lats = cell_grid.grid_points_for_cell(cell_n)
                    cell_n = gpis.__len__()
                    cell_datasets_tmp = {}
                    for var_name in var_list:
                        if var_name == 'lon':
                            cell_datasets_tmp[var_name] = lons.data
                        elif var_name == 'lat':
                            cell_datasets_tmp[var_name] = lats.data
                        elif var_name == 'gpi':
                            cell_datasets_tmp[var_name] = gpis.data
                        elif var_name == 'ALL_n_obs':
                            empty_arr = np.zeros(cell_n)
                            cell_datasets_tmp[var_name] = empty_arr
                        else:
                            empty_arr = np.zeros(cell_n)
                            empty_arr[:] = np.nan
                            cell_datasets_tmp[var_name] = empty_arr
                    cell_datasets_in = deepcopy(cell_datasets_tmp)
                else:
                    logging.warning(' ===> Variable list is not defined. Cell datasets wont be filled by zeros')
                    logging.info(' ----> Cell "' + cell_string + '" ... SKIPPED')
                    continue

            flag_n_obs = True
            if 'ALL_n_obs' in list(cell_datasets_in.keys()):
                n_obs = cell_datasets_in['ALL_n_obs']
                if np.all(n_obs == 0):
                    logging.warning(' ===> Cell "' + cell_string +
                                    '" is empty. All gpi(s) have observations equal to zero')
                    if active_empty_cell:
                        flag_n_obs = False

                    if active_fill_cell:

                        if var_list is not None:
                            gpis, lons, lats = cell_grid.grid_points_for_cell(cell_n)
                            cell_n = gpis.__len__()

                            cell_datasets_tmp = {}
                            for var_name in var_list:
                                if var_name == 'lon':
                                    cell_datasets_tmp[var_name] = lons.data
                                elif var_name == 'lat':
                                    cell_datasets_tmp[var_name] = lats.data
                                elif var_name == 'gpi':
                                    cell_datasets_tmp[var_name] = gpis.data
                                elif var_name == 'ALL_n_obs':
                                    empty_arr = np.zeros(cell_n)
                                    cell_datasets_tmp[var_name] = empty_arr
                                else:
                                    empty_arr = np.zeros(cell_n)
                                    empty_arr[:] = np.nan
                                    cell_datasets_tmp[var_name] = empty_arr
                            cell_datasets_in = deepcopy(cell_datasets_tmp)
                        else:
                            logging.warning(' ===> Variable list is not defined. Cell datasets wont be filled by zeros')
                            logging.info(' ----> Cell "' + cell_string + '" ... SKIPPED')
                            continue

            # check if samples are empty or not
            if flag_n_obs:
                # iterate over variables
                var_n = None
                for var_name_in, var_name_out in zip(list_variable_in, list_variable_out):

                    if var_name_in in list(cell_datasets_in.keys()):
                        var_data_in = cell_datasets_in[var_name_in]

                        if isinstance(var_data_in, np.ma.MaskedArray):
                            var_data_in = var_data_in.data

                        if var_n is None:
                            var_n = var_data_in.shape[0]

                        # Save variable(s) defined in list (or save all variable(s)
                        if var_name_out not in cell_datasets_out:
                            cell_datasets_out[var_name_out] = {}
                            cell_datasets_out[var_name_out] = var_data_in
                        elif var_name_out in cell_datasets_out:
                            var_data_tmp = cell_datasets_out[var_name_out]
                            var_data_tmp = np.concatenate([var_data_tmp, var_data_in])
                            cell_datasets_out[var_name_out] = var_data_tmp

                    else:
                        if variable_not_found == 'warning':
                            logging.warning(' ===> Variable "' + var_name_in + '" not found in file "' +
                                            file_path_datasets + '"')
                            logging.warning(' ===> Variable "' + var_name_out + '" not saved in datasets')
                        elif variable_not_found == 'error':
                            logging.error(' ===> Variable "' + var_name_in + '" not found in file "' +
                                          file_path_datasets + '"')
                            raise RuntimeError('Variable "' + var_name_out + '" not saved in datasets')

                # Save extra variable(s)
                var_cell_out = [int(cell_string)] * var_n
                if 'cell' not in cell_datasets_out:
                    cell_datasets_out['cell'] = {}
                    cell_datasets_out['cell'] = var_cell_out
                elif 'cell' in cell_datasets_out:
                    var_cell_tmp = cell_datasets_out['cell']
                    var_cell_tmp = np.concatenate([var_cell_tmp, var_cell_out])
                    cell_datasets_out['cell'] = var_cell_tmp

                # info end cell
                logging.info(' ----> Cell "' + cell_string + '" ... DONE')

            else:
                # info end cell
                logging.warning(' ===> File cell "' + file_path_datasets + '" samples are always zeros')
                logging.info(' ----> Cell "' + cell_string + '" ... SKIPPED')

        else:
            # info end cell
            logging.warning(' ===> File cell "' + file_path_datasets + '" not found')
            logging.info(' ----> Cell "' + cell_string + '" ... SKIPPED')

    if cell_datasets_out is None:
        logging.error(' ===> Cell datasets is defined by NoneType. All cells are not found or empty')
        raise RuntimeError('Cell datasets is empty. Check your datasets path(s)')

    return cell_datasets_out

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to add datasets grid
def organize_datasets_grid(data_obj, data_type='ASCAT',
                           file_name_grid='TUW_WARP5_grid_info_2_3.nc'):

    # get reference grid
    if os.path.exists(file_name_grid):
        grid_dframe = read_grid_file(file_name_grid)
    else:
        logging.error(' ===> File grid "' + file_name_grid + '" is not found!')
        raise FileNotFoundError('File must be available to correctly run the algorithm')

    # convert in upper case (just to avoid errors in upper/lower case of digits)
    data_type = data_type.upper()
    if data_type not in data_type_expected:
        logging.error(' ===> Datasets type is not included in expected type(s)')
        raise NotImplemented('Case not implemented yet')

    if ('gpi' in data_obj) and ('gpi' in grid_dframe) and ('committed_area' in grid_dframe):

        gpi_grid = grid_dframe['gpi'].values
        car_grid = grid_dframe['committed_area'].values
        land_grid = grid_dframe['land_flag'].values
        gpi_data = data_obj['gpi']

        if data_type == 'ASCAT':
            idx_data = np.apply_along_axis(lambda f: gpi_grid.searchsorted(f), 0, gpi_data)
            car_data = car_grid[idx_data]
            land_data = land_grid[idx_data]
            data_obj['committed_area'] = car_data
            data_obj['land_flag'] = land_data
        elif data_type == 'RZSM' or data_type == 'ECMWF':
            car_data = np.ones(shape=[gpi_data.shape[0]]) # if using np.zeros(), also use car_data[:] = 1
            data_obj['committed_area'] = car_data
            land_data = np.ones(shape=[gpi_data.shape[0]]) # if using np.zeros(), also use land_data[:] = 1
            data_obj['land_flag'] = land_data
        elif data_type == 'HMC':
            car_data = np.ones(shape=[gpi_data.shape[0]]) # if using np.zeros(), also use car_data[:] = 1
            data_obj['committed_area'] = car_data
            land_data = np.ones(shape=[gpi_data.shape[0]]) # if using np.zeros(), also use land_data[:] = 1
            data_obj['land_flag'] = land_data
        else:
            logging.error(' ===> Datasets type is not expected ("ASCAT" or "ECMWF/RZSM" or "HMC") ')
            raise NotImplemented('Case not implemented yet')

    else:
        logging.error(' ===> Grid fields are not available in the grid file "' + file_name_grid + '" ')
        raise RuntimeError('Fields must be defined in the grid obj')

    return data_obj
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
