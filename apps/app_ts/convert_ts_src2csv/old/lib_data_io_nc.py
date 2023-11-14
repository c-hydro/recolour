"""
Library Features:

Name:          lib_data_io_nc
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import numpy as np
import xarray as xr

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read mat obj
def read_file_nc(file_name, var_name='SM', geo_x_name='Longitude', geo_y_name='Latitude'):

    dset_file = xr.open_dataset(file_name)

    var_data = dset_file[var_name].values
    geo_x = dset_file[geo_x_name].values
    geo_y = dset_file[geo_y_name].values

    geo_y_upper, geo_y_lower = geo_y[0, 0], geo_y[-1, 0]
    if geo_y_lower > geo_y_upper:
        geo_y = np.flipud(geo_y)
        var_data = np.flipud(var_data)

    return var_data, geo_x, geo_y
# -------------------------------------------------------------------------------------

