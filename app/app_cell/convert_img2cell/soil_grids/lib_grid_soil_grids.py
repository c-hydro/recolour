"""
Library Features:

Name:           lib_grid_soilgrids
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260521'
Version:        '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import os
import logging
import numpy as np

from netCDF4 import Dataset
from pygeogrids.grids import BasicGrid
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
def grids(
        only_land=False,
        grid_path=None,
        cellsize=5.0):
    """
    Create SoilGrids reference grid.

    Parameters
    ----------
    only_land : bool, optional (default: False)
        Uses the land mask to reduce the grid to land points only.

    Returns
    -------
    grid : pygeogrids.CellGrid
        Either a land grid or a global grid.
    """

    if grid_path is None:

        grid_folder = os.path.dirname(__file__)
        grid_file = 'soilgrids_grid.nc'

        grid_path = os.path.join(
            grid_folder,
            grid_file
        )

    if not os.path.exists(grid_path):
        logging.error(' ===> Grid file "' + grid_path + '" is not available')
        raise RuntimeError('Grid file is needed by the procedure')

    file_grid = Dataset(grid_path)

    domain_lon_1d = file_grid['longitude'][:]
    domain_lat_1d = file_grid['latitude'][:]

    domain_lon_2d, domain_lat_2d = np.meshgrid(
        domain_lon_1d,
        domain_lat_1d
    )

    domain_grid = BasicGrid(
        domain_lon_2d.flatten(),
        domain_lat_2d.flatten()
    ).to_cell_grid(
        cellsize=cellsize
    )

    if only_land:

        if 'mask' in file_grid.variables:

            land_mask_2d = file_grid.variables['mask'][:]
            land_mask_1d = land_mask_2d.flatten() == 1

            land_points = np.ma.masked_array(
                domain_grid.get_grid_points()[0],
                np.logical_not(land_mask_1d)
            )

            land_grid = domain_grid.subgrid_from_gpis(
                land_points[~land_points.mask].filled()
            )

            file_grid.close()

            return land_grid

        else:

            logging.warning(
                ' ===> Land mask not available in grid file. '
                'Return complete grid.'
            )

            file_grid.close()

            return domain_grid

    else:

        file_grid.close()

        return domain_grid

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to configure cell grid
def cell_grid(
        grid_path=None,
        cellsize=5.0):
    """
    Alias to create SoilGrids cell grid.
    """

    return grids(
        only_land=False,
        grid_path=grid_path,
        cellsize=cellsize
    )

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to configure land grid
def land_grid(
        grid_path=None,
        cellsize=5.0):
    """
    Alias to create SoilGrids land grid.
    """

    return grids(
        only_land=True,
        grid_path=grid_path,
        cellsize=cellsize
    )

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to load grid
def load_grid(
        land_points=True,
        bbox=None,
        grid_path=None,
        cellsize=5.0):
    """
    Load SoilGrids grid.

    Parameters
    ----------
    land_points : bool, optional (default: True)
        Use only land points.

    bbox : tuple, optional (default: None)
        (min_lon, min_lat, max_lon, max_lat)

    Returns
    -------
    subgrid : CellGrid
    """

    if land_points:

        subgrid = land_grid(
            grid_path=grid_path,
            cellsize=cellsize
        )

        if bbox is not None:

            subgrid = subgrid4bbox(
                subgrid,
                *bbox
            )

    else:

        if bbox is not None:

            subgrid = subgrid4bbox(
                cell_grid(
                    grid_path=grid_path,
                    cellsize=cellsize
                ),
                *bbox
            )

        else:

            subgrid = cell_grid(
                grid_path=grid_path,
                cellsize=cellsize
            )

    return subgrid

# -------------------------------------------------------------------------------------