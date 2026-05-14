"""
Library Features:

Name:          lib_utils_grid
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260512'
Version:       '1.1.0'
"""

# -----------------------------------------------------------------------------
# libraries
import logging
import os
import netCDF4
import time
import pandas as pd
import numpy as np

from copy import deepcopy

import pygeogrids.grids as grids

from lib_data_io_pickle import read_file_obj, write_file_obj

# set default values
committed_area_default, land_data_default = 1, 1
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
        for f in ['committed_area', 'land_flag']:
            if f in file_handle.variables:
                grid_obj[f] = file_handle.variables[f][:]
            else:
                logging.warning(' ===> Variable "' + f + '" is not available in the grid file, Use default value')
                grid_obj[f] = np.zeros(shape=[gpi_data.shape[0]])
                if f == 'committed_area':
                    grid_obj[f][:] = committed_area_default
                if f == 'land_flag':
                    grid_obj[f][:] = land_data_default

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
def get_grid_objects(cell_start=0, cell_end=2566, cells_list=None, land_filter=True,
                     path_grid='', file_grid='TUW_WARP5_grid_info_2_2.nc',
                     reset_anc=True,
                     path_anc=None, file_anc='TUW_WARP5_grid_info_2_2.workspace'):

    # get grid objects - start
    logging.info(' ---> Get grid objects ... ')

    # manage paths source and ancillary
    file_path = (
        os.path.join(path_grid, file_grid)
        if path_grid is not None else file_grid
    )

    path_anc = (
        os.path.join(path_anc, file_anc)
        if path_anc is not None else file_anc
    )
    if path_anc is not None:
        os.makedirs(os.path.dirname(path_anc), exist_ok=True)

    # remove ancillary file (if needed)
    if reset_anc:
        if os.path.exists(path_anc):
            os.remove(path_anc)

    # get grid files - start
    logging.info(f' ----> Get information from grid file {file_path} ... ')
    if not os.path.exists(path_anc):
        if os.path.exists(file_path):

            file_handle = netCDF4.Dataset(os.path.join(path_grid, file_grid), mode='r')
            cell = file_handle['cell'][:]
            gpi = file_handle['gpi'][:]
            lons = file_handle['lon'][:]
            lats = file_handle['lat'][:]
            ca = file_handle['committed_area'][:]

            # apply filtering or not
            if land_filter:
                idx = np.where(file_handle.variables['land_flag'][:] == 1)[0]
            else:
                idx = slice(None)

            # create grid
            grid = grids.CellGrid(lons[idx], lats[idx], cell[idx], gpi[idx])

            # selected gpis and committed are
            gpi = grid.activegpis.data
            ca = ca[idx]

            # store data in workspace files
            collections_anc = {'grid': grid, 'gpi': gpi, 'committed_area': ca}
            write_file_obj(path_anc, collections_anc)

            # get grid files - end (done)
            logging.info(f' ----> Get information from grid file {file_path} ... DONE')
        else:

            # get grid files - end (failed)
            logging.error(' ===> Open grid file "' + file_path + '" ... FILE NOT FOUND')
            logging.info(f' ----> Get information from grid file {file_path} ... FAILED')
            raise RuntimeError('File grid is mandatory. Check your settings.')

    else:

        # get grid files - end (previously saved)
        logging.info(f' ----> Get information from grid file {file_path} ... PREVIOUSLY SAVED')
        collections_anc = read_file_obj(path_anc)
        grid, gpi, ca = collections_anc['grid'], collections_anc['gpi'], collections_anc['committed_area']

    # get grid objects - end
    logging.info(' ---> Get grid objects ... DONE')

    # compute grid objects - start
    logging.info(' ---> Compute grid objects ... ')
    if cells_list is None:

        # info
        logging.info(' ::: List of cell: NOT DEFINED')
        logging.info(' ::: Domain: GLOBAL')

        cell = grid.arrcell
        active_cells = np.unique(cell)

        # set idx start and end
        if cell_start is None:
            idx_start = 0
        else:
            idx_start = np.where(active_cells == cell_start)[0][0]

        if cell_end is None:
            idx_end = active_cells.shape[0]
        else:
            idx_end = np.where(active_cells == cell_end)[0][0] + 1

        # select cells
        cells_select = active_cells[idx_start:idx_end]
        cells_array = cells_select.data

        # transform array to list
        cells_list = cells_array.tolist()

        # grid
        grid_selected = deepcopy(grid)
        # info
        info_selected = pd.DataFrame({"gpi": gpi, "ca": ca, 'cell': cell})

    else:

        # info
        logging.info(' ::: List of cell: DEFINED')

        # check if cell_list or cell_start/cell_end are defined
        if cells_list is None:
            if cell_start is None or cell_end is None:
                logging.error(' ===> The cell_start and cell_end must be defined in the datasets collection')
                raise RuntimeError('The cell_start and cell_end must be defined in the datasets collection')
            else:
                # grid file is not defined
                cells_array = range(cell_start, cell_end)
                # transform array to list
                cells_list = cells_array.tolist()

            logging.info(f' ::: Domain: CELLS RANGE :: {cells_list}')

        else:
            logging.info(f' ::: Domain: CELLS LIST :: {cells_list}')

        # mask over the grid cells
        mask = np.isin(grid.arrcell, cells_list)

        # create filtered grid
        grid_selected = grids.CellGrid(
            grid.activearrlon[mask],
            grid.activearrlat[mask],
            grid.arrcell[mask],
            grid.activegpis[mask]
        )

        # create filtered gpis
        gpi = grid.activegpis[mask]
        cell = grid.arrcell[mask]
        ca = ca[mask]

        # info
        info_selected = pd.DataFrame({"gpi": gpi, "ca": ca, 'cell': cell})

    # compute grid objects - end
    logging.info(' ---> Compute grid objects ... DONE')

    return cells_list, info_selected, grid_selected

# ----------------------------------------------------------------------------------------------------------------------

