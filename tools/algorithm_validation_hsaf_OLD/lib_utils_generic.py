"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import pickle
import os
import numpy as np

from netCDF4 import Dataset
from pytesmo.validation_framework.data_manager import get_result_names
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to read data obj
def read_obj(filename):
    if os.path.exists(filename):
        data = pickle.load(open(filename, "rb"))
    else:
        data = None
    return data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to write data obj
def write_obj(filename, data):
    if os.path.exists(filename):
        os.remove(filename)
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get grid reference
def get_grid_reference(dset_obj, dset_key_root='ref',
                       dset_key_path='path_static', dset_key_file='file_grid'):

    if dset_key_root in list(dset_obj.keys()):
        dset_fields = dset_obj[dset_key_root]

        if dset_key_path in list(dset_fields.keys()):
            dset_path_grid = dset_fields[dset_key_path]
        else:
            logging.error(' ===> The field key "' + dset_key_path + ' is not available in the datasets object')
            raise RuntimeError('The field is needed by the algorithm and must be defined in the datasets obj')

        if dset_key_file in list(dset_fields.keys()):
            dset_file_grid = dset_fields[dset_key_file]
        else:
            logging.error(' ===> The field key "' + dset_key_path + ' is not available in the datasets object')
            raise RuntimeError('The field is needed by the algorithm and must be defined in the datasets obj')
    else:
        logging.error(' ===> The datasets key "' + dset_key_root + '" is not available in the datasets collection')
        raise RuntimeError('The field is needed by the algorithm and must be defined in the datasets collection')

    return dset_path_grid, dset_file_grid
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to define grid cells
def get_grid_cells(cell_start=0, cell_end=2566, cells_list=None,
                   path_grid='', file_grid='TUW_WARP5_grid_info_2_2.nc'):

    # info start
    logging.info(' ---> Compute cells list ... ')

    # check if cell list is defined or not
    if cells_list is None:

        # use grid file
        logging.info(' ----> Use information defined in the grid files ... ')
        file_path = os.path.join(path_grid, file_grid)
        if os.path.exists(file_path):

            file_handle = Dataset(os.path.join(path_grid, file_grid), mode='r')
            cell = file_handle['cell'][:]
            gpi = file_handle['gpi'][:]

            land = None
            if 'land_flag' in list(file_handle.variables):
                land = file_handle['land_flag'][:]
            file_handle.close()

            # get cells
            if land is not None:
                cells = np.unique(cell[land == 1])
            else:
                cells = np.unique(cell)

            # set idx start and end
            try:
                idx_start = np.where(cells == cell_start)[0][0]
            except BaseException as base_exp:
                logging.warning(' ===> Idx start is not available in the cells obj. Start is set to 0.')
                logging.warning(' Warning "' + str(base_exp) + '" found')
                idx_start = 0
            try:
                idx_end = np.where(cells == cell_end)[0][0] + 1
            except BaseException as base_exp:
                logging.warning(' ===> Idx end is not available in the cells obj. End is set to cell maximum length.')
                logging.warning(' Warning "' + str(base_exp) + '" found')
                idx_end = cells.shape[0]

            # select cells
            cells_obj = cells[idx_start:idx_end]
            cells_array = cells_obj.data
            # select gpis
            if land is not None:
                gpis = gpi[land == 1]
            else:
                gpis = gpi

            logging.info(' ----> Use information defined in the grid files ... DONE')

        else:

            logging.warning(' ===> Open grid file "' + file_path + '" ... FILE NOT FOUND')
            logging.info(' ----> Use information defined in the grid files ... FAILED')
            logging.info(' ----> Use information defined by "cell_start" and "cell_end"')

            # grid file is not defined
            cells_array = range(cell_start, cell_end)
            gpis = None

        # transform array to list
        cells_list = cells_array.tolist()

    else:
        # use cell list defined in the settings file
        logging.info(' ----> Use information defined in the settings file ... ')
        if not isinstance(cells_list, list):
            cells_list = [cells_list]
        gpis = None
        logging.info(' ----> Use information defined in the settings file ... DONE')

    # info end
    logging.info(' ---> Compute cells list ... DONE')

    return cells_list, gpis

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to slice list
def slice_list(elements, step):
    return [elements[x:x + step] for x in list(range(0, len(elements), step))]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to make folder
def make_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to check bit setting
def get_bit(a, bit_pos):
    """
    Returns 1 or 0 if bit is set or not.

    Parameters
    ----------
    a : int or numpy.ndarray
      Input array.
    bit_pos : int
      Bit position. First bit position is right.

    Returns
    -------
    b : numpy.ndarray
      1 if bit is set and 0 if not.
    """
    return np.clip(np.bitwise_and(a, 2 ** (bit_pos-1)), 0, 1)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get dataset mode(s)
def get_dataset_modes(dset_obj, dset_mode='ref'):
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

    # WARNING: get_result_names (pytesmo.validation_framework.data_manager)
    # DOES NOT provide the datasets' names in the correct order

    ds_names = get_result_names(ds_dict, ref_key, n=len(datasets))
    dataset_names = [ref_key]
    for name in ds_names[0]:
        if name[0] != ref_key: dataset_names.append(name[0])

    return dataset_names

# ----------------------------------------------------------------------------------------------------------------------
