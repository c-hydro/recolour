"""
Library Features:

Name:           lib_data_io_cell
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org),
                Martina Natali (martina01.natali@edu.unife.it)
Date:           '20230711'
Version:        '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import xarray as xr

from repurpose.resample import resample_to_grid
from pygeogrids.grids import BasicGrid

from lib_info_args import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get grid data
def get_grid_data(dframe_cell, mask_grid, geox_grid, geoy_grid, no_data=-9999.0,
                  var_name='soil_moisture_ref', geox_name='longitude', geoy_name='latitude',
                  resampling_radius=10000, resampling_min_neighbours=1, resampling_neighbours=8):

    var_data_grid, idx_data_grid = None, None
    if dframe_cell is not None:
        var_data_arr = dframe_cell[var_name].values[:, 0]
        var_geox_arr = dframe_cell[geox_name].values
        var_geoy_arr = dframe_cell[geoy_name].values

        var_obj = resample_to_grid(
            {'data': var_data_arr},
            var_geox_arr, var_geoy_arr, geox_grid, geoy_grid,
            search_rad=resampling_radius, fill_values=np.nan,
            min_neighbours=resampling_min_neighbours, neighbours=resampling_neighbours)
        var_data_grid = var_obj['data']
        var_data_grid[mask_grid == 0] = np.nan

        idx_data_grid = np.argwhere(np.isfinite(var_data_grid))

    return var_data_grid, idx_data_grid
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to find grid cells
def find_grid_cells(cells_list=None, cell_start=None, cell_end=None):

    if cells_list is None and ((cell_start is not None) and (cell_end is not None)):
        cells_array = range(cell_start, cell_end)
        cells_list = cells_array.tolist()
    elif cells_list is not None:
        pass
    else:
        alg_logger.error(' ===> Cells data format is not supported')
        raise NotImplemented('Case not implemented yet')
    return cells_list
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create grid obj
def create_grid_obj(cell, longitude, latitude, gpis, luts):
    grid_obj = {'longitude': longitude, 'latitude': latitude,
                'cell': cell, 'gpis': gpis, 'luts': luts}
    return grid_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to find grid fields
def find_grid_fields(grid_ref, grid_k1, grid_k2, max_distance_k1=25000, max_distance_k2=25000):

    lut_idxs_cmp_k1 = grid_ref.calc_lut(other=grid_k1, max_dist=max_distance_k1)
    lut_idxs_cmp_k2 = grid_ref.calc_lut(other=grid_k2, max_dist=max_distance_k2)

    lut_idxs_ref, lut_idxs_k1, lut_idxs_k2 = [], [], []
    for el_idx, (el_value_k1, el_value_k2) in enumerate(zip(lut_idxs_cmp_k1, lut_idxs_cmp_k2)):
        if (el_value_k1 > 0) and (el_value_k2 > 0):
            lut_idxs_ref.append(el_idx)
            lut_idxs_k1.append(el_value_k1)
            lut_idxs_k2.append(el_value_k2)

    if lut_idxs_cmp_k1.__len__() > lut_idxs_ref.__len__():
        alg_logger.warning(' ===> Grid indexes computed and filtered of no data have different lengths. Try to set'
                           'a different maximum distance for k1 to avoid the undefined point between the grids')
    if lut_idxs_cmp_k2.__len__() > lut_idxs_ref.__len__():
        alg_logger.warning(' ===> Grid indexes computed and filtered of no data have different lengths. Try to set'
                           'a different maximum distance for k2 to avoid the undefined point between the grids')

    # dataset k1
    cell_k1 = grid_k1.activearrcell[lut_idxs_k1]
    lons_k1 = grid_k1.activearrlon[lut_idxs_k1]
    lats_k1 = grid_k1.activearrlat[lut_idxs_k1]
    gpis_k1 = grid_k1.activegpis[lut_idxs_k1]
    grid_k1 = create_grid_obj(cell_k1, lons_k1, lats_k1, gpis_k1, lut_idxs_k1)

    # dataset k2
    cell_k2 = grid_k2.activearrcell[lut_idxs_k2]
    lons_k2 = grid_k2.activearrlon[lut_idxs_k2]
    lats_k2 = grid_k2.activearrlat[lut_idxs_k2]
    gpis_k2 = grid_k2.activegpis[lut_idxs_k2]
    grid_k2 = create_grid_obj(cell_k2, lons_k2, lats_k2, gpis_k2, lut_idxs_k2)

    # dataset ref
    cell_ref = grid_ref.activearrcell[lut_idxs_ref]
    lons_ref = grid_ref.activearrlon[lut_idxs_ref]
    lats_ref = grid_ref.activearrlat[lut_idxs_ref]
    gpis_ref = grid_ref.activegpis[lut_idxs_ref]
    grid_ref = create_grid_obj(cell_ref, lons_ref, lats_ref, gpis_ref, lut_idxs_ref)

    return grid_ref, grid_k1, grid_k2
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
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
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create grid cell
def create_grid_cell(cell_lons, cell_lats):
    grid_obj = BasicGrid(cell_lons, cell_lats).to_cell_grid(cellsize=5.)
    return grid_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read grid cell (based on BasicGrid)
def read_grid_cell(grid_path, file_grid='grid.nc'):

    if grid_path is None:
        grid_folder, grid_file = os.path.dirname(__file__), file_grid
        grid_path = os.path.join(grid_folder, grid_file)

    if not os.path.exists(grid_path):
        alg_logger.error(' ===> Grid file "' + grid_path + '" is not available')
        raise RuntimeError('Grid file is needed by the procedure')

    grid_handle = xr.open_dataset(grid_path)

    lon_1d = grid_handle['lon'].values
    lat_1d = grid_handle['lat'].values

    lon_max = np.nanmax(lon_1d)
    if lon_max > 180:
        lon_1d = (lon_1d - 180)

    grid_obj = BasicGrid(lon_1d, lat_1d).to_cell_grid(cellsize=5.)

    return grid_obj
# ----------------------------------------------------------------------------------------------------------------------
