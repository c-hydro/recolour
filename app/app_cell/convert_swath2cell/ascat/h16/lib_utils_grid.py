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
        grid_obj['committed_area'] = file_handle.variables['committed_area'][:]
        grid_obj['land_flag'] = file_handle.variables['land_flag'][:]
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


# -------------------------------------------------------------------------------------
# method to create subgrid using a bbox
def subgrid4bbox(grid, min_lon, min_lat, max_lon, max_lat):
    """
    Select a spatial subset for the grid by bound box corner points

    Parameters
    ----------
    grid: BasicGrid or CellGrid
        Grid object to trim.
    min_lon: float
        Lower left corner longitude
    min_lat: float
        Lower left corner latitude
    max_lon: float
        Upper right corner longitude
    max_lat: float
        Upper right corner latitude

    Returns
    -------
    subgrid: BasicGrid or CellGrid
        Subset of the input grid.
    """
    gpis, lons, lats, _ = grid.get_grid_points()
    assert len(gpis) == len(lats) == len(lons)
    bbox_gpis = gpis[
        np.where(
            (lons <= max_lon)
            & (lons >= min_lon)
            & (lats <= max_lat)
            & (lats >= min_lat)
        )
    ]

    return grid.subgrid_from_gpis(bbox_gpis)
# -------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Method to open netcdf grid file
def read_grid_file_OLD(filename):

    fileobj = netCDF4.Dataset(filename, 'r')

    fileresult = {}
    for filevar in fileobj.variables:

        filevar = str(filevar)
        filedata = fileobj[filevar][:]

        if isinstance(filedata, np.ma.MaskedArray):
            filedata = filedata.data

        fileresult[filevar] = filedata

    fileobj.close()

    if 'land_flag' in fileresult:
        land_id = np.where(fileresult['land_flag'][:] == 1)

        for var in fileresult:
            data = fileresult[var][land_id]
            fileresult[var] = data

    return fileresult
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

