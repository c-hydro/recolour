"""
Library Features:

Name:          lib_fx_swi
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240523'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging

import numpy as np
import pandas as pd

from pytesmo.time_series.filters import exp_filter

from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute data swi
def compute_data_swi(ts_data, ts_var="swi", ts_method="exp_filter", ts_ctime=6):

    # get julian dates
    ts_jd = ts_data.index.to_julian_date().values
    # drop nan values
    ts_data.dropna(inplace=True)
    # get values
    ts_values_raw = ts_data.values.astype(np.float64)

    # iterate over characteristic times
    if ts_method == "exp_filter":
        ts_values_filtered = exp_filter(ts_values_raw, ts_jd, ctime=ts_ctime)
    else:
        alg_logger.error(" ===> Method '" + ts_method + "' is not supported yet")
        raise NotImplementedError("Case not implemented yet")
    # define data output
    ts_data_out = pd.Series(ts_values_filtered, index=ts_data.index, name=ts_var)

    return ts_data_out

# ----------------------------------------------------------------------------------------------------------------------
