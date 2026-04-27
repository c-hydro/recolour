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

from config_utils import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
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

