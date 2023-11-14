"""
Library Features:

Name:           lib_grid_hmc
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org),
                Martina Natali (martina01.natali@edu.unife.it)
Date:           '20230711'
Version:        '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
from pygeogrids.grids import BasicGrid
from netCDF4 import Dataset
import os
# -------------------------------------------------------------------------------------


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


# -------------------------------------------------------------------------------------
# method to define grids
def grids(only_land=False, grid_path=None):
    """
    Create global 0.25 DEG grids (origin in bottom left)

    Parameters
    ---------
    only_land : bool, optional (default: False)
        Uses the land mask to reduce the GLDAS 0.25DEG land grid to land points
        only.

    Returns
    --------
    grid : pygeogrids.CellGrid
        Either a land grid or a global grid
    """

    if grid_path is None:
        grid_folder, grid_file = os.path.dirname(__file__), 'hmc_grid.nc'
        grid_path = os.path.join(grid_folder, grid_file)

    if not os.path.exists(grid_path):
        logging.error(' ===> Grid file "' + grid_path + '" is not available')
        raise RuntimeError('Grid file is needed by the procedure')

    file_grid = Dataset(grid_path)

    domain_lon_1d = file_grid['longitude'][:]
    domain_lat_1d = file_grid['latitude'][:]
    domain_lon_2d, domain_lat_2d = np.meshgrid(domain_lon_1d, domain_lat_1d)

    # *************************************************************************
    # The grid which is created is, as always, of the proprietary types
    # defined by TUW as BasicGrid or CellGrid
    # this is what it is returned
    domain_grid = BasicGrid(domain_lon_2d.flatten(), domain_lat_2d.flatten()).to_cell_grid(cellsize=5.)

    if only_land:
        ds = Dataset(
            os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                "hmc_grid.nc",
            )
        )
        land_lats_1d = ds.variables["latitude_mask"][:]
        land_mask_1d = ds.variables["mask"][:].flatten().filled() == 0.0
        dlat = domain_lat_1d.size - land_lats_1d.size

        land_mask = np.concatenate((np.ones(dlat * domain_lon_1d.size), land_mask_1d))
        land_points = np.ma.masked_array(
            domain_grid.get_grid_points()[0], land_mask
        )

        land_grid = domain_grid.subgrid_from_gpis(
            land_points[~land_points.mask].filled()
        )
        return land_grid
    else:
        return domain_grid
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to configure cell grid
def cell_grid(grid_path=None):
    """Alias to create a global 0.25 DEG grid without gaps w. 5 DEG cells """
    return grids(only_land=False, grid_path=grid_path)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to configure land grid
def land_grid(grid_path=None):
    """Alias to create a global 0.25 DEG grid over land only w. 5 DEG cells """
    return grids(only_land=False, grid_path=grid_path)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to load grid
def load_grid(land_points=True, bbox=None, grid_path=None):
    """
    Load hmc grid.

    Parameters
    ----------
    land_points : bool, optional (default: True)
        Reshuffle only land points
    bbox : tuple, optional (default: True)
        (min_lat, min_lon, max_lat, max_lon)
        Bounding box to limit reshuffling to.
    """
    if land_points:
        subgrid = land_grid(grid_path)
        if bbox is not None:
            subgrid = subgrid4bbox(subgrid, *bbox)
    else:
        if bbox is not None:
            subgrid = subgrid4bbox(cell_grid(grid_path), *bbox)
        else:
            subgrid = None

    return subgrid
# -------------------------------------------------------------------------------------
