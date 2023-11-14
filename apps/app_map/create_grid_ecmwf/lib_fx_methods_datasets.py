"""
Library Features:

Name:          lib_fx_methods_datasets
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from copy import deepcopy

import numpy as np

# debug
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute rzsm_OLD
def compute_rzsm(var_data_k1, var_data_k2, var_data_k3, **kwargs):

    # check variable(s)
    if (var_data_k1 is None) or (var_data_k2 is None) or (var_data_k3 is None):
        logging.error(' ===> All/some variable(s) are defined by NoneType')
        raise RuntimeError('All variable(s) must be defined to compute weight(s)')

    # var no_data
    var_no_data = np.isnan(var_data_k1)

    # var_0_7 = var40
    var_rzsm_k1 = deepcopy(var_data_k1)
    # var_0_28=0.25*var40 + 0.75*var41
    var_rzsm_k2 = 0.25 * var_data_k1 + 0.75 * var_data_k2
    # var_0_100=0.07*var40 + 0.21*var41 + 0.72*var42
    var_rzsm_k3 = 0.07 * var_data_k1 + 0.21 * var_data_k2 + 0.72 * var_data_k3

    # mask no data
    var_rzsm_k1[var_no_data] = np.nan
    var_rzsm_k2[var_no_data] = np.nan
    var_rzsm_k3[var_no_data] = np.nan

    return var_rzsm_k1, var_rzsm_k2, var_rzsm_k3

# ----------------------------------------------------------------------------------------------------------------------
