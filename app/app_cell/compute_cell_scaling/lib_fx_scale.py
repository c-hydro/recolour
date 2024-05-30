"""
Library Features:

Name:          lib_fx_scale
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings

import numpy as np
import pandas as pd

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pytesmo.scaling import get_scaling_function, get_scaling_method_lut

from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute data scaling
def compute_data_scaling(fx_data, fx_time, fx_method="norm_min_max", fx_args=None):

    # get fx handle
    fx_handle = get_scaling_function(fx_method)
    if fx_handle is None:
        alg_logger.error(" ===> Method '" + fx_method + "' is not supported yet")
        raise NotImplementedError("Case not implemented yet")

    # get ts n
    ts_n_src, ts_n_ref = fx_data['src'].shape[0], fx_data['ref'].shape[0]

    # get fx time
    time_src, time_ref = fx_time['src'], fx_time['ref']
    # get fx data
    data_src, data_ref = fx_data['src'], fx_data['ref']

    # organize dataframes
    ts_src = pd.DataFrame(data={'data': data_src}, index=time_src)
    ts_ref = pd.DataFrame(data={'data': data_ref}, index=time_ref)
    # remove nan(s)
    ts_src.dropna(inplace=True)
    ts_ref.dropna(inplace=True)

    # check if time-series are not empty
    if (not ts_src.empty) and (not ts_ref.empty):

        # cast data to float64
        values_src = ts_src.values.astype(np.float64).flatten()
        values_ref = ts_ref.values.astype(np.float64).flatten()
        values_args = {'ref': values_ref, 'src': values_src}

        # destination data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data_dst = fx_handle(**values_args)
        # destination time
        time_dst = ts_src.index.values

    else:
        data_dst, time_dst = None, None

    return data_dst, time_dst
# ----------------------------------------------------------------------------------------------------------------------
