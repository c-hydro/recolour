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

import lib_fx_method as library_method
from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute data scaling
def compute_data_scaling(fx_data, fx_time, fx_method="mean_std", fx_args=None):

    # get fx handle
    if hasattr(library_method, fx_method):
        fx_handle = getattr(library_method, fx_method)
    else:
        alg_logger.error(" ===> Method '" + fx_method + "' is not supported yet")
        raise NotImplementedError("Case not implemented yet")

    # get ts n
    ts_n_ref_nrt, ts_n_ref_dr = fx_data['ref_nrt'].shape[0], fx_data['ref_dr'].shape[0]
    ts_n_other_dr = fx_data['other_dr'].shape[0]

    # get fx time
    time_ref_nrt, time_ref_dr, time_other_dr = fx_time['ref_nrt'], fx_time['ref_dr'], fx_time['other_dr']
    # get fx data
    data_ref_nrt, data_ref_dr, data_other_dr = fx_data['ref_nrt'], fx_data['ref_dr'], fx_data['other_dr']

    # organize dataframes
    ts_ref_nrt = pd.DataFrame(data={'data': data_ref_nrt}, index=time_ref_nrt)
    ts_ref_dr = pd.DataFrame(data={'data': data_ref_dr}, index=time_ref_dr)
    ts_other_dr = pd.DataFrame(data={'data': data_other_dr}, index=time_other_dr)

    # remove nan(s)
    ts_ref_nrt.dropna(inplace=True)
    ts_ref_dr.dropna(inplace=True)
    ts_other_dr.dropna(inplace=True)

    # check if time-series are not empty
    if ts_ref_nrt.empty:
        alg_logger.warning(' ===> Time-series "ref_nrt" is empty. Scaling is not performed')
    if ts_ref_dr.empty:
        alg_logger.warning(' ===> Time-series "ref_dr" is empty. Scaling is not performed')
    if ts_other_dr.empty:
        alg_logger.warning(' ===> Time-series "other_dr" is empty. Scaling is not performed')

    # check if time-series are not empty
    if (not ts_ref_nrt.empty) and (not ts_ref_dr.empty) and (not ts_other_dr.empty):

        # cast data to float64
        values_ref_nrt = ts_ref_nrt.values.astype(np.float64).flatten()
        values_ref_dr = ts_ref_dr.values.astype(np.float64).flatten()
        values_other_dr = ts_other_dr.values.astype(np.float64).flatten()
        values_args = {'ref_nrt': values_ref_nrt, 'ref_dr': values_ref_dr, 'other_dr': values_other_dr}

        # destination data
        data_dst = fx_handle(**values_args)
        # destination time
        time_dst = ts_ref_nrt.index.values

    else:
        data_dst, time_dst = None, None

    return data_dst, time_dst
# ----------------------------------------------------------------------------------------------------------------------
