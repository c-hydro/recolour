"""
Library Features:

Name:           lib_utils_geo
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
from pygeogrids.grids import BasicGrid

from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to regrid domain
def regrid_domain(lons, lats, mask, resolution_km=9.0):

    valid_lons = lons[mask].ravel().astype(np.float64)
    valid_lats = lats[mask].ravel().astype(np.float64)

    resolution_deg = resolution_km / 111.0

    lon_min = np.nanmin(valid_lons)
    lon_max = np.nanmax(valid_lons)
    lat_min = np.nanmin(valid_lats)
    lat_max = np.nanmax(valid_lats)

    reg_lons_1d = np.arange(
        lon_min,
        lon_max + resolution_deg,
        resolution_deg,
        dtype=np.float64
    )

    reg_lats_1d = np.arange(
        lat_min,
        lat_max + resolution_deg,
        resolution_deg,
        dtype=np.float64
    )

    reg_lons_2d, reg_lats_2d = np.meshgrid(reg_lons_1d, reg_lats_1d)

    source_grid = BasicGrid(
        lon=valid_lons,
        lat=valid_lats
    )

    nearest_gpi, distance = source_grid.find_nearest_gpi(
        reg_lons_2d.ravel(),
        reg_lats_2d.ravel(),
        max_dist=resolution_km * 1000.0
    )

    nearest_gpi = np.asarray(nearest_gpi)
    distance = np.asarray(distance)

    reg_mask = (
        (nearest_gpi >= 0) &
        np.isfinite(distance)
    ).reshape(reg_lons_2d.shape)

    return reg_lons_2d, reg_lats_2d, reg_mask
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to map points to grid indices
def map_points_to_grid_indices(src_lons, src_lats, grid_lons, grid_lats, domain_mask):
    nrows, ncols = domain_mask.shape

    if ncols < 2 or nrows < 2:
        raise RuntimeError("Grid must have at least 2 rows and 2 columns")

    x0 = float(grid_lons[0, 0])
    y0 = float(grid_lats[0, 0])

    dx = float(grid_lons[0, 1] - grid_lons[0, 0])
    dy = float(grid_lats[1, 0] - grid_lats[0, 0])

    cols = np.rint((src_lons - x0) / dx).astype(np.int64)
    rows = np.rint((src_lats - y0) / dy).astype(np.int64)

    inside = (
        (rows >= 0) & (rows < nrows) &
        (cols >= 0) & (cols < ncols)
    )

    rows = rows[inside]
    cols = cols[inside]

    inside_domain = domain_mask[rows, cols]
    rows = rows[inside_domain]
    cols = cols[inside_domain]

    return rows, cols
# ----------------------------------------------------------------------------------------------------------------------

