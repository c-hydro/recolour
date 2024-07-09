"""
Library Features:

Name:          lib_utils_grid
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import netCDF4
import pandas as pd
import numpy as np

from copy import deepcopy
from pygeogrids.grids import BasicGrid

from lib_data_io_geo import read_file_raster
from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to join gpi between grid and data
def find_grid_idx_to_data(gpi_grid, gpi_data):
    idx_data = np.apply_along_axis(lambda f: gpi_grid.searchsorted(f), 0, gpi_data)
    return idx_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read domain file
def read_domain_file(file_name):
    grid_da = read_file_raster(file_name, output_format='data_array')
    return grid_da
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to read ascat committed area
def read_grid_file(file_name):

    grid_obj = {}
    with netCDF4.Dataset(file_name) as file_handle:
        grid_obj['lon'] = file_handle.variables['lon'][:]
        grid_obj['lat'] = file_handle.variables['lat'][:]
        grid_obj['gpi'] = file_handle.variables['gpi'][:]

        grid_dim = file_handle.variables['gpi'][:].__len__()
        # get committed area data
        if 'committed_area' in list(file_handle.variables):
            grid_obj['committed_area'] = file_handle.variables['committed_area'][:]
        else:
            alg_logger.warning(' ===> Variable "committed_area" not available. All values will be equal to 1')
            grid_committed_area_default = np.zeros(shape=grid_dim)
            grid_committed_area_default[:] = 1
            grid_obj['committed_area'] = grid_committed_area_default
        # get land flag data
        if 'land_flag' in list(file_handle.variables):
            grid_obj['land_flag'] = file_handle.variables['land_flag'][:]
        else:
            alg_logger.warning(' ===> Variable "land_flag" not available. All values will be equal to 1')
            grid_land_flag_default = np.zeros(shape=grid_dim)
            grid_land_flag_default[:] = 1
            grid_obj['land_flag'] = grid_land_flag_default
        # create mask data
        grid_obj['mask'] = (grid_obj['committed_area'] == 0) & (grid_obj['land_flag'] == 1)

    # check masked array in the grid obj
    grid_tmp = {}
    for grid_key, grid_values in grid_obj.items():
        if isinstance(grid_values, np.ma.MaskedArray):
            grid_values = grid_values.data
            grid_tmp[grid_key] = grid_values
        else:
            grid_tmp[grid_key] = grid_values
    grid_obj = grid_tmp.copy()

    grid_dframe = pd.DataFrame(grid_obj)

    return grid_dframe
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get grid reference
def get_grid_reference(dset_obj, dset_key_root='reference',
                       dset_key_path='path_static', dset_key_file='file_grid'):

    if dset_key_root in list(dset_obj.keys()):
        dset_fields = dset_obj[dset_key_root]

        if dset_key_path in list(dset_fields.keys()):
            dset_path_grid = dset_fields[dset_key_path]
        else:
            alg_logger.error(' ===> The field key "' + dset_key_path + ' is not available in the datasets object')
            raise RuntimeError('The field is needed by the algorithm and must be defined in the datasets obj')

        if dset_key_file in list(dset_fields.keys()):
            dset_file_grid = dset_fields[dset_key_file]
        else:
            alg_logger.error(' ===> The field key "' + dset_key_path + ' is not available in the datasets object')
            raise RuntimeError('The field is needed by the algorithm and must be defined in the datasets obj')
    else:
        alg_logger.error(' ===> The datasets key "' + dset_key_root + '" is not available in the datasets collection')
        raise RuntimeError('The field is needed by the algorithm and must be defined in the datasets collection')

    return dset_path_grid, dset_file_grid
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get grid obj
def get_grid_obj(grid_reference):

    # get lons, lats and gpis
    lat_1d, lon_1d, gpi_1d = grid_reference['lat'], grid_reference['lon'], grid_reference['gpi']
    
    land_flag, committed_area, mask = None, None, None
    if 'land_flag' in list(grid_reference.keys()):
        land_flag = grid_reference['land_flag'].values
    if 'committed_area' in list(grid_reference.keys()):
        committed_area = grid_reference['committed_area'].values
    if 'mask' in list(grid_reference.keys()):
        mask = grid_reference['mask'].values
    
    # create grid object
    grid_obj = BasicGrid(
        lon=lon_1d, lat=lat_1d, gpis=gpi_1d).to_cell_grid(cellsize=5.)

    if land_flag is not None:
        grid_obj.land_flag = land_flag
    if committed_area is not None:
        grid_obj.committed_area = committed_area
    if mask is not None:
        grid_obj.mask = mask

    return grid_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to define grid cells
def get_grid_cells(cell_start=0, cell_end=2566, cells_list=None, cell_grid_reference=None):

    # info start
    alg_logger.info(' -----> Compute cells list ... ')

    # check if cell list is defined or not
    if cells_list is None:

        # use grid file
        alg_logger.info(' ------> Use information defined in the grid files ... ')
        if cell_grid_reference is not None:

            cells = cell_grid_reference.activearrcell[:]
            gpis = cell_grid_reference.activegpis[:]

            land_flag, committed_area, mask = None, None, None
            if hasattr(cell_grid_reference, 'land_flag'):
                land_flag = cell_grid_reference.land_flag[:]
            if hasattr(cell_grid_reference, 'committed_area'):
                committed_area = cell_grid_reference.committed_area[:]
            if hasattr(cell_grid_reference, 'mask'):
                mask = cell_grid_reference.mask[:]

            # filter land flag
            if land_flag is not None:
                cells = np.unique(cells[land_flag == 1])

            # set idx start and end
            try:
                idx_start = np.where(cells == cell_start)[0][0]
            except BaseException as base_exp:
                alg_logger.warning(' ===> Idx start is not available in the cells obj. Start is set to 0.')
                alg_logger.warning(' Warning "' + str(base_exp) + '" found')
                idx_start = 0
            try:
                idx_end = np.where(cells == cell_end)[0][0] + 1
            except BaseException as base_exp:
                alg_logger.warning(' ===> Idx end is not available in the cells obj. End is set to cell maximum length.')
                alg_logger.warning(' Warning "' + str(base_exp) + '" found')
                idx_end = cells.shape[0]

            # select cells
            cells_obj = cells[idx_start:idx_end]
            cells_array = cells_obj.data

            # select gpis
            if land_flag is not None:
                gpis_list = gpis[land_flag == 1]
            else:
                gpis_list = deepcopy(gpis)

            alg_logger.info(' ------> Use information defined in the grid files ... DONE')

        else:

            alg_logger.info(' ------> Use information defined in the grid files ... FAILED')
            alg_logger.info(' ------> Use information defined by "cell_start" and "cell_end"')

            # grid file is not defined
            cells_array = range(cell_start, cell_end)
            gpis_list = None

        # transform array to list
        cells_list = cells_array.tolist()

    else:
        # use cell list defined in the settings file
        alg_logger.info(' ------> Use information defined in the settings file ... ')
        if not isinstance(cells_list, list):
            cells_list = [cells_list]
        gpis_list = None
        alg_logger.info(' ------> Use information defined in the settings file ... DONE')

    # info end
    alg_logger.info(' -----> Compute cells list ... DONE')

    return cells_list, gpis_list

# ----------------------------------------------------------------------------------------------------------------------
