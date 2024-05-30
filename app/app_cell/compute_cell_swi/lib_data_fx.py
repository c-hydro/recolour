"""
Library Features:

Name:          lib_data_fx
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240508'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd
from copy import deepcopy

from pytesmo.time_series.filters import exp_filter

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to prepare data
def prepare_data(var_layer_obj, var_name_list, var_layer_template='var_data_in'):
    var_layer_data = {}
    for i, var_name_i in enumerate(var_name_list):
        var_layer_i = var_layer_template.format(i + 1)
        var_layer_data[var_layer_i] = var_layer_obj[var_name_i].values
    return var_layer_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute soil-water-index
def compute_swi(var_data_in, var_name='swi', var_method='exp_filter', var_ctime=6, var_no_data=-999999, **kwargs):

    # get julian dates
    var_jd = var_data_in.index.to_julian_date().values
    # drop nan values
    var_data_in[var_data_in == var_no_data] = np.nan
    var_data_in.dropna(inplace=True)
    # get values
    var_values_raw = var_data_in.values.astype(np.float64)
    var_values_raw = np.squeeze(var_values_raw)

    # iterate over characteristic times
    if var_method == "exp_filter":
        var_values_swi = exp_filter(var_values_raw, var_jd, ctime=var_ctime)
    else:
        logging.error(" ===> Method '" + var_method + "' is not supported yet")
        raise NotImplementedError("Case not implemented yet")

    # define data output
    var_data_out = pd.Series(var_values_swi, index=var_data_in.index, name=var_name)

    return var_data_out

# ----------------------------------------------------------------------------------------------------------------------
