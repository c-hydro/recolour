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
# method to fill a given list of missing cells using default values
# only cells with at least min_gpi_land land points are filled
def fill_datasets_cell(
        cell_datasets_obj,
        cell_missing_list,
        cell_grid=None,
        info_grid=None,
        variable_cell='cell',
        variable_gpi='gpi',
        variable_lon='lon',
        variable_lat='lat',
        variable_carea='ca',
        variable_label_carea='label_ca',
        variable_label_type='label_type',
        min_gpi_land=1000,
        default_values=None):

    logging.info(' ----> Fill missing cells by list ... ')

    if cell_datasets_obj is None:
        cell_datasets_obj = pd.DataFrame()

    info_grid_idx = info_grid.set_index(variable_gpi)

    # -----------------------------------------------------------------------------
    # create default values automatically if not defined
    if default_values is None:

        default_values = {}

        for var_name in cell_datasets_obj.columns:

            if var_name in [
                variable_gpi,
                variable_lon,
                variable_lat,
                variable_carea,
                variable_label_carea,
                variable_label_type,
                variable_cell
            ]:
                continue

            var_dtype = cell_datasets_obj[var_name].dtype

            if np.issubdtype(var_dtype, np.number):
                default_values[var_name] = np.nan
            else:
                default_values[var_name] = None

        logging.info(' -----> Default values created automatically')
    # -----------------------------------------------------------------------------

    datasets_missing = []

    for cell_n in cell_missing_list:

        cell_n = int(cell_n)

        logging.info(' -----> Cell "{:}" ... '.format(cell_n))

        # select grid points for this cell
        mask = cell_grid.arrcell == cell_n

        gpis_n = cell_grid.activegpis[mask]
        lons_n = cell_grid.activearrlon[mask]
        lats_n = cell_grid.activearrlat[mask]

        if gpis_n.size == 0:
            logging.warning(
                ' ===> Cell "{:}" has no grid points. SKIPPED'.format(cell_n)
            )
            continue

        # select info only for available gpis
        info_n = info_grid_idx.loc[gpis_n]

        # land / committed-area flag
        ca_n = info_n[variable_carea].values

        # check minimum land gpis
        idx_land = ~np.isnan(ca_n)
        n_gpi_land = np.sum(idx_land)

        if n_gpi_land < min_gpi_land:

            logging.warning(
                ' ===> Cell "{:}" has only {:} land gpis. '
                'Minimum is {:}. SKIPPED'.format(
                    cell_n, n_gpi_land, min_gpi_land
                )
            )

            continue

        # keep only land gpis
        gpis_n = gpis_n[idx_land]
        lons_n = lons_n[idx_land]
        lats_n = lats_n[idx_land]
        ca_n = ca_n[idx_land]

        data_tmp = {
            variable_gpi: gpis_n,
            variable_lon: lons_n,
            variable_lat: lats_n,
            variable_carea: ca_n,
            variable_label_carea: np.where(
                ca_n == 1,
                'committed_area',
                'no_committed_area'
            ),
            variable_label_type: np.full(gpis_n.shape, 'global'),
            variable_cell: np.full(gpis_n.shape, cell_n)
        }

        # fill remaining variables using defaults
        for var_name, var_value in default_values.items():

            if var_name in data_tmp:
                continue

            if np.isscalar(var_value) or var_value is None:
                data_tmp[var_name] = np.full(gpis_n.shape, var_value)
            else:
                data_tmp[var_name] = np.asarray(var_value)

        datasets_missing.append(pd.DataFrame(data_tmp))

        logging.info(
            ' -----> Cell "{:}" ... DONE. Filled {:} land gpis'.format(
                cell_n, gpis_n.size
            )
        )

    if datasets_missing:

        datasets_missing = pd.concat(
            datasets_missing,
            ignore_index=True
        )

        cell_datasets_obj = pd.concat(
            [cell_datasets_obj, datasets_missing],
            ignore_index=True
        )

        cell_datasets_obj = cell_datasets_obj.sort_values(
            by=[variable_cell, variable_gpi]
        ).reset_index(drop=True)

    logging.info(' ----> Fill missing cells by list ... DONE')

    return cell_datasets_obj
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to filter datasets cell
def filter_datasets_cell(cell_list_in=None, cell_digits=4,
                         folder_name_datasets='', file_name_datasets='{cell}.nc'):

        # info script start
        logging.info(' ----> Filter cell datasets (select/missing)... ')

        # iterate over cell(s)
        cell_list_select, cell_list_missing = [], []
        for cell_n in cell_list_in:

            # format cell
            cell_n = int(cell_n)
            cell_string = str(cell_n).zfill(cell_digits)

            # info start cell
            logging.info(' -----> Cell "' + cell_string + '" ... ')

            # compose file name
            file_path_datasets = os.path.join(folder_name_datasets, file_name_datasets)
            file_path_datasets = file_path_datasets.format(cell=cell_string)

            # search file name (if exists or not)
            logging.info(' ------> File "' + file_path_datasets + '" ...')
            if os.path.exists(file_path_datasets):
                logging.info(' ------> File "' + file_path_datasets + '" ... EXISTS. Cell selected.')
                cell_list_select.append(cell_n)
            else:
                logging.info(' ------> File "' + file_path_datasets + '" ... DOES NOT EXISTS. Cell removed')
                cell_list_missing.append(cell_n)

            logging.info(' -----> Cell "' + cell_string + '" ... DONE')

        # info script end
        logging.info(' ----> Filter cell datasets (select/missing) ... DONE')

        return cell_list_select, cell_list_missing
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to organize datasets cells
def organize_datasets_cell(
        cell_list=None, cell_digits=4, cell_grid=None, info_grid=None,
        list_variable_in=None, list_variable_out=None,
        folder_name_datasets='', file_name_datasets='{cell}.nc',
        fill_empty_data=True):

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
        logging.info(' -----> Cell "' + cell_string + '" ... ')

        # define folder and file name
        file_path_datasets = os.path.join(folder_name_datasets, file_name_datasets)
        file_path_datasets = file_path_datasets.format(cell=cell_string)

        # get gpis, lons and lats
        mask = cell_grid.arrcell == cell_n
        cell_n = cell_grid.activearrcell[mask]
        gpis_n = cell_grid.activegpis[mask]
        lons_n = cell_grid.activearrlon[mask]
        lats_n = cell_grid.activearrlat[mask]

        cell_test = cell_grid.gpi2cell(gpis_n)

        # select rows using gpis_n = info_grid["gpi"], info_grid["gpi"]
        info_n = info_grid[info_grid["gpi"].isin(gpis_n)]

        # info start read
        logging.info(' ------> Read datasets ... ')
        # check if file exists or not
        if os.path.exists(file_path_datasets):

            # get datasets in
            cell_datasets_in = read_file_cell(file_path_datasets, expected_variables=list_variable_in)
            if var_list is None:
                var_list = list(cell_datasets_in.keys())

            # check datasets in
            if cell_datasets_in is None or not cell_datasets_in:

                # info end read - failed
                logging.warning(' ===> Cell datasets are defined by NoneType')
                logging.info(' ------> Read datasets ... FAILED')
                logging.info(' -----> Cell "' + cell_string + '" ... SKIPPED')

            else:

                # info end read - done
                logging.info(' ------> Read datasets ... DONE')

                # info start organize
                logging.info(' ------> Organize datasets ... ')
                # iterate over variables
                logging.info(' -------> Get data ... ')
                cell_datasets_tmp = {}
                for var_name_in, var_name_out in zip(list_variable_in, list_variable_out):

                    # check variable in datasets in
                    if var_name_in in list(cell_datasets_in.keys()):

                        # get data
                        var_data_in = cell_datasets_in[var_name_in]

                        # Save variable(s) defined in list (or save all variable(s)
                        if var_name_out not in cell_datasets_tmp:
                            cell_datasets_tmp[var_name_out] = {}
                            cell_datasets_tmp[var_name_out] = var_data_in

                # get gpi values
                gpis_tmp = cell_datasets_tmp['gpi']
                # select committed area using gpi
                ca_tmp = info_grid.set_index("gpi").loc[gpis_tmp, "ca"].values
                # save committed area in temporary datasets
                cell_datasets_tmp['ca'] = ca_tmp

                # create committed area labels
                cell_datasets_tmp['label_ca'] = np.where(ca_tmp == 1,'committed_area','no_committed_area')
                # create type label
                cell_datasets_tmp['label_type'] = np.full(len(ca_tmp),'global')

                logging.info(' -------> Get data ... DONE')

                # fill data
                logging.info(' -------> Fill data ... ')
                if fill_empty_data:

                    # existing gpis in dictionary
                    gpis_dict = cell_datasets_tmp['gpi']

                    # find missing gpis
                    mask_missing = ~np.isin(gpis_n, gpis_dict)

                    gpis_missing = gpis_n[mask_missing]
                    lons_missing = lons_n[mask_missing]
                    lats_missing = lats_n[mask_missing]

                    # get ca for missing gpis
                    ca_missing = (
                        info_grid
                        .set_index("gpi")
                        .loc[gpis_missing, "ca"]
                        .values
                    )

                    # append mandatory variables
                    cell_datasets_tmp['gpi'] = np.concatenate([cell_datasets_tmp['gpi'], gpis_missing])
                    cell_datasets_tmp['lon'] = np.concatenate([cell_datasets_tmp['lon'], lons_missing])
                    cell_datasets_tmp['lat'] = np.concatenate([cell_datasets_tmp['lat'], lats_missing])

                    # append ca and labels
                    cell_datasets_tmp['ca'] = np.concatenate([cell_datasets_tmp['ca'], ca_missing])

                    label_ca_missing = np.where(
                        ca_missing == 1,
                        'committed_area',
                        'no_committed_area'
                    )

                    cell_datasets_tmp['label_ca'] = np.concatenate([
                        cell_datasets_tmp['label_ca'],
                        label_ca_missing
                    ])

                    cell_datasets_tmp['label_type'] = np.concatenate([
                        cell_datasets_tmp['label_type'],
                        np.full(gpis_missing.shape, 'global')
                    ])

                    # fill all other variables with nan
                    for key, values in list(cell_datasets_tmp.items()):
                        if key in ['gpi', 'lon', 'lat', 'ca', 'label_ca', 'label_type']:
                            continue

                        nan_values = np.full(gpis_missing.shape, np.nan)
                        cell_datasets_tmp[key] = np.concatenate([values, nan_values])

                    # sorting indices using gpi
                    idx_sort = np.argsort(cell_datasets_tmp['gpi'])

                    # reorder all dictionary arrays
                    for key in cell_datasets_tmp.keys():
                        cell_datasets_tmp[key] = np.asarray(cell_datasets_tmp[key])[idx_sort]

                    var_n = len(cell_datasets_tmp['gpi'])
                    cell_datasets_tmp['cell'] = np.full(var_n, int(cell_string))

                    logging.info(' -------> Fill data ... DONE')

                else:
                    var_n = len(cell_datasets_tmp['gpi'])
                    cell_datasets_tmp['cell'] = np.full(var_n, int(cell_string))
                    logging.info(' -------> Fill data ... SKIPPED')

                # merge data
                logging.info(' -------> Merge data ... ')
                if cell_datasets_out is None:
                    cell_datasets_out = pd.DataFrame(cell_datasets_tmp)
                else:
                    cell_datasets_tmp = pd.DataFrame(cell_datasets_tmp)
                    cell_datasets_step = pd.DataFrame(cell_datasets_out)
                    cell_datasets_out = pd.concat([cell_datasets_step, cell_datasets_tmp], ignore_index=True)

                logging.info(' -------> Merge data ... DONE')

                # info end organize
                logging.info(' ------> Organize datasets ... DONE')

                # info end cell
                logging.info(' -----> Cell "' + cell_string + '" ... DONE')

        else:
            # info end cell
            logging.warning(' ===> File cell "' + file_path_datasets + '" not found')
            logging.info(' ------> Read datasets ... FAILED')
            logging.info(' -----> Cell "' + cell_string + '" ... SKIPPED')

    # check if datasets out is None or not defined
    if cell_datasets_out is None or cell_datasets_out.empty:
        logging.error(' ===> Cell datasets is defined by NoneType. All cells are not found or empty')
        raise RuntimeError('Cell datasets is empty. Check your datasets path(s)')

    # info script end
    logging.info(' ----> Get cell datasets ... DONE')

    return cell_datasets_out

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to add datasets grid
def organize_datasets_grid(
        data_obj,
        data_type='ASCAT',
        file_name_grid='TUW_WARP5_grid_info_2_3.nc'):

    # info script start
    logging.info(' ----> Organize cell datasets ... ')

    # get reference grid
    if os.path.exists(file_name_grid):
        grid_dframe = read_grid_file(file_name_grid)
    else:
        logging.error(f' ===> File grid "{file_name_grid}" is not found!')
        raise FileNotFoundError(
            'File must be available to correctly run the algorithm'
        )

    # upper case
    data_type = data_type.upper()

    if data_type not in data_type_expected:
        logging.error(' ===> Datasets type is not included in expected type(s)')
        raise NotImplementedError('Case not implemented yet')

    # check fields
    required_grid_fields = ['gpi', 'committed_area', 'land_flag']
    if not all(field in grid_dframe.columns for field in required_grid_fields):
        logging.error(f' ===> Grid fields are not available in "{file_name_grid}"')
        raise RuntimeError('Fields must be defined in the grid obj')
    if 'gpi' not in data_obj.columns:
        logging.error(' ===> "gpi" field not available in data object')
        raise RuntimeError('Field "gpi" is mandatory')

    # ASCAT
    if data_type == 'ASCAT':
        grid_info = grid_dframe[['gpi', 'committed_area', 'land_flag']]
        data_obj = data_obj.merge(grid_info,on='gpi', how='left')
    # ECMWF / RZSM / HMC
    elif data_type in ['RZSM', 'ECMWF', 'HMC']:
        data_obj['committed_area'] = 1
        data_obj['land_flag'] = 1
    else:
        logging.error(
            ' ===> Datasets type is not expected '
            '("ASCAT" or "ECMWF/RZSM" or "HMC")'
        )
        raise NotImplementedError('Case not implemented yet')

    # info script end
    logging.info(' ----> Organize cell datasets ... DONE')

    return data_obj
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to convert datasets obj
def convert_datasets_obj(data_obj):

    # check if dataframe
    if isinstance(data_obj, pd.DataFrame):
        data_dframe = data_obj.copy()
    # convert dictionary
    elif isinstance(data_obj, dict):
        data_dframe = pd.DataFrame.from_dict(data_obj)
    else:
        logging.error(' ===> Datasets object format is not supported')
        raise TypeError('Expected pandas.DataFrame or dict')

    # add type column if missing
    if 'type' not in data_dframe.columns:
        data_dframe['type'] = 'global'

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
