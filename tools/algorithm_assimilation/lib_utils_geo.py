# -------------------------------------------------------------------------------------
# Libraries
import logging
import numpy as np
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create index map
def create_map_idx(grid_values):

    grid_rows = grid_values.shape[0]
    grid_cols = grid_values.shape[1]
    grid_n = grid_rows * grid_cols
    grid_idxs = np.arange(grid_n).reshape(grid_rows, grid_cols)

    return grid_idxs
# -------------------------------------------------------------------------------------
