"""
Library Features:

Name:          lib_utils_data_point
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import numpy as np
import pandas as pd

from copy import deepcopy

from lib_data_io_mat import read_file_mat
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to join series
def join_data_series(series_sm, dframe_sm_in=None):

    if dframe_sm_in is None:
        dframe_sm_out = series_sm.to_frame()
    else:
        dframe_sm_tmp = series_sm.to_frame()
        dframe_sm_out = dframe_sm_in.merge(dframe_sm_tmp, left_index=True, right_index=True)

    return dframe_sm_out
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to concatenate series
def concatenate_data_series(series_sm_step, series_sm_concat=None):
    if series_sm_concat is None:
        series_sm_concat = deepcopy(series_sm_step)
    else:
        series_sm_concat = pd.concat([series_sm_concat, series_sm_step])
    return series_sm_concat
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get data matlab
def get_data_mat(file_name, var_name='sm', var_value_min=0.0, var_value_max=100.0, var_scale_factor=1.0):

    file_obj = read_file_mat(file_name)

    file_datetime_idx = pd.DatetimeIndex(file_obj['a1sDateVet'])
    file_values = file_obj['a1dVWC'].flatten()

    file_values[file_values < float(var_value_min)] = np.nan
    file_values[file_values > float(var_value_max)] = np.nan
    file_values = file_values * float(var_scale_factor)

    series_sm = pd.Series(data=file_values, index=file_datetime_idx, name=var_name)

    return series_sm

# -------------------------------------------------------------------------------------
