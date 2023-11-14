"""
Library Features:

Name:          lib_data_io_binary
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import struct
import numpy as np

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read file binary
def read_file_binary(file_name, data_geo, scale_factor=10000):

    file_handle = open(file_name, 'rb')

    data_n = data_geo.shape[0] * data_geo.shape[1]
    data_format = 'i' * data_n

    file_obj = file_handle.read(-1)
    data_array = struct.unpack(data_format, file_obj)

    data_grid = np.reshape(data_array, (data_geo.shape[0], data_geo.shape[1]), order='F')

    data_grid = np.float32(data_grid / scale_factor)

    data_grid[data_geo < 0] = np.nan
    data_grid[0, :] = np.nan
    data_grid[-1, :] = np.nan
    data_grid[:, 0] = np.nan
    data_grid[:, -1] = np.nan

    # Debug
    # plt.figure()
    # plt.imshow(data_grid)
    # plt.colorbar()
    # plt.show()

    return data_grid

# -------------------------------------------------------------------------------------
