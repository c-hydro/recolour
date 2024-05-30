"""
Library Features:

Name:          lib_fx_metrics
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240523'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
f
import numpy as np

from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute data metrics
def compute_data_metrics(ts_data, fx_method="mean", fx_args=None, fx_decimal_round=2):

    if fx_args is None:
        fx_args = {}

    if fx_method == "mean":
        value = ts_data.mean()
    elif fx_method == "median":
        value = ts_data.median()
    elif fx_method == "std":
        value = ts_data.std()
    elif fx_method == "min":
        value = ts_data.min()
    elif fx_method == "max":
        value = ts_data.max()
    elif "perc" in fx_method:
        value = ts_data.quantile(fx_args['percentile'])
    else:
        alg_logger.error(" ===> Method '" + fx_method + "' is not supported yet")
        raise NotImplementedError("Case not implemented yet")

    value = np.round(value, decimals=fx_decimal_round)

    return value

# ----------------------------------------------------------------------------------------------------------------------
