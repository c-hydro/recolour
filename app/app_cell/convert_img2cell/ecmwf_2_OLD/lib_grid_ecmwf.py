"""
Library Features:

Name:           lib_grid_ecmwf
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

from pygeogrids.grids import BasicGrid
from netCDF4 import Dataset

# debugging
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


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


# ----------------------------------------------------------------------------------------------------------------------
# method to create grid
def cell_grid(grid_path, file_grid='grid.nc'):

    if grid_path is None:
        grid_folder, grid_file = os.path.dirname(__file__), file_grid
        grid_path = os.path.join(grid_folder, grid_file)

    if not os.path.exists(grid_path):
        logging.error(' ===> Grid file "' + grid_path + '" is not available')
        raise RuntimeError('Grid file is needed by the procedure')

    data_grid = Dataset(grid_path)
    if 'lon' in list(data_grid.variables):
        domain_lon_1d = data_grid['lon'][:]
    elif 'longitude' in list(data_grid.variables):
        domain_lon_1d = data_grid['longitude'][:]
    else:
        logging.error(' ===> Longitude variable in grid file is not defined by expected tags (lon, longitude)')
        raise RuntimeError('Check the grid file to define the tag of longitude variable')

    if 'lat' in list(data_grid.variables):
        domain_lat_1d = data_grid['lat'][:]
    elif 'latitude' in list(data_grid.variables):
        domain_lat_1d = data_grid['latitude'][:]
    else:
        logging.error(' ===> Latitude variable in grid file is not defined by expected tags (lat, latitude)')
        raise RuntimeError('Check the grid file to define the tag of latitude variable')

    domain_lon_max = np.nanmax(domain_lon_1d)
    if domain_lon_max > 180:
        domain_lon_1d = (domain_lon_1d - 180) # % 360 - 180

    domain_lon_2d, domain_lat_2d = np.meshgrid(domain_lon_1d, domain_lat_1d)

    ''' debug
    plt.figure(); plt.imshow(domain_lon_2d); plt.colorbar()
    plt.figure(); plt.imshow(domain_lat_2d); plt.colorbar()
    plt.show()
    '''

    domain_grid = BasicGrid(domain_lon_2d.flatten(), domain_lat_2d.flatten()).to_cell_grid(cellsize=5.)

    return domain_grid
# ----------------------------------------------------------------------------------------------------------------------
