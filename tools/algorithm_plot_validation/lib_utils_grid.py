"""
Library Features:

Name:          lib_utils_grid
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# -----------------------------------------------------------------------------
# libraries
import logging
import os
import netCDF4
import time
import pandas as pd
import numpy as np
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to join gpi between grid and data
def find_grid_idx_to_data(gpi_grid, gpi_data):
    idx_data = np.apply_along_axis(lambda f: gpi_grid.searchsorted(f), 0, gpi_data)
    return idx_data
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Method to read ascat committed area
def read_grid_file(file_name):

    grid_obj = {}
    with netCDF4.Dataset(file_name) as file_handle:
        grid_obj['lon'] = file_handle.variables['lon'][:]
        grid_obj['lat'] = file_handle.variables['lat'][:]
        grid_obj['gpi'] = file_handle.variables['gpi'][:]
        gpi_data = grid_obj['gpi']

        # check if committed_area and land_flag are present in the grid
        # otherwise mask completely/put equal to global data
        # committed_area: 0 = not_committed, 1 = committed
        # land_flag: 0 = no_land, 1 = land
        for f in ['committed_area', 'land_flag']:
            try:
                grid_obj[f] = file_handle.variables[f][:]
            except KeyError as e:
                print(f'KeyError: {e}',
                      f'\nGrid does not contain {e} information: {e} is taken equal to global area')
                grid_obj[f] = np.zeros(shape=[gpi_data.shape[0]])

        grid_obj['land_flag'][:] = 1
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
# -----------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get grid reference
def get_grid_reference(dset_obj, dset_key_root='reference',
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

            file_handle = netCDF4.Dataset(os.path.join(path_grid, file_grid), mode='r')
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

