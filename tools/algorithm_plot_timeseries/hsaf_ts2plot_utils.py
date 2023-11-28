# ----------------------------------------------------------------------------------------------------------------------
# Libraries
import numpy as np
import json
import os

from datetime import datetime
from dateutil.relativedelta import relativedelta

from netCDF4 import Dataset

from pygeogrids.grids import BasicGrid
from pytesmo.validation_framework.data_manager import get_result_names
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Method to compute dates group by month and year
def compute_times(startDate, endDate, split='month'):

    cur_date = start = datetime.strptime(startDate, '%Y-%m-%d').date()
    end = datetime.strptime(endDate, '%Y-%m-%d').date()

    dates_start = []
    dates_end = []
    if split:
        dates_start.append(start)
        while cur_date < end:
            #print(cur_date)
            if split == 'month':
                cur_date += relativedelta(months=1)
            elif split == 'year':
                cur_date += relativedelta(years=1)
            dates_end.append(cur_date)
            dates_start.append(cur_date)

        del dates_start[-1]
    else:
        dates_start.append(start)
        dates_end.append(end)

    return dates_start, dates_end
# ----------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Method to open netcdf file
def open_file(filename, filevars=None):

    fileobj = Dataset(filename, 'r')

    fileresult = {}
    for filevar in fileobj.variables:

        filevar = str(filevar)

        if filevars is not None:
            if filevar in filevars:
                filedata = fileobj[filevar][:]
            else:
                filedata = None
        else:
            filedata = fileobj[filevar][:]

        if filedata is not None:
            fileresult[filevar] = filedata

    fileobj.close()

    return fileresult
# -----------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to create folder (and check if folder exists)
def create_folder(sPathName=None, sPathDelimiter=None):

    from os import makedirs
    from os.path import exists

    if sPathName:
        if sPathDelimiter:
            sPathNameSel = sPathName.split(sPathDelimiter)[0]
        else:
            sPathNameSel = sPathName

        if not exists(sPathNameSel):
            makedirs(sPathNameSel)
        else:
            pass
    else:
        pass
# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to load grid risico datasets
def load_grid_risico(filename):

    data_grid = Dataset(filename)

    lon_1d = data_grid['longitude'][:]
    lat_1d = data_grid['latitude'][:]
    lon, lat = np.meshgrid(lon_1d, lat_1d)

    return BasicGrid(lon.flatten(), lat.flatten()).to_cell_grid(cellsize=5.)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to open log file
def open_log_file(idx, cell, path):

    filename = 'log_group_' + str(idx) + '_cell_' + str(cell) + '.txt'

    file_path = os.path.join(path, filename)

    make_folder(path)
    delete_file(file_path)

    file_handle = open(file_path, 'w')

    return file_handle

    #logging.basicConfig(filename=file_path, level=logging.INFO)

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to close log file
def close_log_file(file_handle):
    file_handle.close()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to slice list
def list_slice(S, step):
    return [S[x:x + step] for x in xrange(0, len(S), step)]
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to delete old file
def delete_file(file):

    if os.path.isfile(file):
        os.remove(file)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# function to make folder
def make_folder(path):

    if not os.path.exists(path):
        os.makedirs(path)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to save group order
def write_group_cells(groups, path):

    filename = 'log_list_group.txt'
    file = os.path.join(path, filename)

    make_folder(path)
    delete_file(file)

    file_obj = open(file, 'w')

    for idx, group in enumerate(groups):
        file_obj.write('group: ' + str(idx) + ' cell: ' + str(group) + ' \n')

    file_obj.close()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to divide cells in subgroups
def get_group_cells(cells, mp_cpu=1):

    n_cells = len(cells)
    n_groups = n_cells/mp_cpu

    if n_groups < 1:
        n_groups = 1

    n_elem = int(np.floor(n_cells/n_groups))

    groups = list_slice(cells, n_elem)

    return groups

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to get grid cells
def get_grid_cells(cell_start=0, cell_end=2566, path_grid='', file_grid='TUW_WARP5_grid_info_2_2.nc', cell_list=None):

    if os.path.exists(os.path.join(path_grid, file_grid)):

        ncFile = Dataset(os.path.join(path_grid, file_grid), mode='r')
        cell = ncFile['cell'][:]
        gpi = ncFile['gpi'][:]
        land = ncFile['land_flag'][:]
        ncFile.close()

        cells = np.unique(cell[land == 1])

        try:
            idx_start = np.where(cells == cell_start)[0][0]
        except BaseException:
            idx_start = 0

        try:
            idx_end = np.where(cells == cell_end)[0][0] + 1
        except BaseException:
            idx_end = cells.shape[0]

        cells = cells[idx_start:idx_end]

        gpis = gpi[land == 1]

    else:

        cells = range(cell_start, cell_end)
        gpis = None

    if cell_list:

        cell_num = map(int, cell_list)
        cells = np.asarray(cell_num)

    return cells, gpis

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Method to parser configuration file
def parser_config_file(file):

    with open(file) as json_config_file:
        config = json.load(json_config_file)

    return config

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to check bit setting
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
# Method to get dataset name(s)
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
