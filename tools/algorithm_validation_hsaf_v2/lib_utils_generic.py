"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260506'
Version:       '1.1.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import pickle
import os
import errno
import json
import shutil
import numpy as np

from pathlib import Path
from netCDF4 import Dataset
from pytesmo.validation_framework.data_manager import get_result_names
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to dump file workspace
def dump_file_workspace(file_name, datasets, attributes, results):

    print('save file pickle')

    workspace = {
        'datasets': datasets,
        'attributes': attributes,
        'results': results
    }

    with open(file_name, 'wb') as file_handle:
        pickle.dump(workspace, file_handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f'pickle file saved: {file_name}')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to open log file
def manage_file_workspace(file_name_tmpl, file_gpi, file_cell, fill_cell=4, file_update=True):

    file_format_tmpl = {'cell': str(file_cell).zfill(fill_cell), 'gpi': str(file_gpi)}
    file_string = file_name_tmpl.format(**file_format_tmpl)
    file_obj = Path(file_string)

    if file_obj.parent == Path("."):
        file_name = file_obj.name
        folder_name = None
        path_name = file_name
    else:
        folder_name = file_obj.parent
        file_name = file_obj.name
        path_name = os.path.join(folder_name, file_name)

    # create folder (if needed)
    if folder_name is not None:
        os.makedirs(folder_name, exist_ok=True)

    return path_name

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
        file_path = os.path.join(path_grid, file_grid)

        logging.info(f' ----> Use information defined in the grid files {file_path} ... ')
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
            if cell_start is not None:
                try:
                    idx_start = np.where(cells == cell_start)[0][0]
                except BaseException as base_exp:
                    logging.warning(' ===> Idx start is not available in the cells obj. Start is set to 0.')
                    logging.warning(' Warning "' + str(base_exp) + '" found')
                    idx_start = 0
            else:
                logging.warning(
                    ' ===> Idx start equal to cell start selected from grid file. Cell start is null in settings file')
                idx_start = 0

            if cell_end is not None:
                try:
                    idx_end = np.where(cells == cell_end)[0][0] + 1
                except BaseException as base_exp:
                    logging.warning(' ===> Idx end is not available in the cells obj. End is set to cell maximum length.')
                    logging.warning(' Warning "' + str(base_exp) + '" found')
                    idx_end = cells.shape[0]

            else:
                logging.warning(
                    ' ===> Idx end equal to cell start selected from grid file. Cell end is null in settings file')
                idx_end = cell_end

            # select cells
            cells_obj = cells[idx_start:idx_end]
            cells_array = cells_obj.data
            # select gpis
            if land is not None:
                gpis = gpi[land == 1]
            else:
                gpis = gpi

            logging.info(f' ----> Use information defined in the grid files {file_path} ... DONE')

        else:

            logging.warning(' ===> Open grid file "' + file_path + '" ... FILE NOT FOUND')
            logging.info(' ----> Use information defined in the grid files ... FAILED')
            logging.info(' ----> Use information defined by "cell_start" and "cell_end"')

            if (cell_start is None) or (cell_end is None):
                logging.error(' ===> Variable(s) "cell_start" and "cell_end" must be defined by integer')
                raise RuntimeError('Cells are defined by NoneType')

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
# method to clean folder
def clean_folder(path):

    for name in os.listdir(path):
        file_path = os.path.join(path, name)

        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except FileNotFoundError:
            pass

        except OSError as e:
            if e.errno == errno.EBUSY and name.startswith(".nfs"):
                logging.warning(f" ===> Skipping busy NFS file: {file_path}")
            else:
                logging.error(f" ===> Failed removing {file_path}: {e}")
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
