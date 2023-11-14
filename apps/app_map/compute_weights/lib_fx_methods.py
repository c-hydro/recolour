"""
Library Features:

Name:          lib_fx_methods
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import numpy as np

# debug
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute weights using error variance (3 vars)
def compute_weights_k3(var_data_k1, var_data_k2, var_data_k3, **kwargs):

    if (var_data_k1 is None) or (var_data_k2 is None) or (var_data_k3 is None):
        logging.error(' ===> All/some variable(s) are defined by NoneType')
        raise RuntimeError('All variable(s) must be defined to compute weight(s)')

    # lower part of the weight(s) computation
    var_weight_lower = var_data_k1 * var_data_k2 + var_data_k1 * var_data_k3 + var_data_k2 * var_data_k3
    # upper part of the weight(s) computation
    var_weight_k1_upper = var_data_k2 * var_data_k3
    var_weight_k2_upper = var_data_k1 * var_data_k3
    var_weight_k3_upper = var_data_k1 * var_data_k2

    # compute weight(s)
    var_weights_k1 = var_weight_k1_upper / var_weight_lower
    var_weights_k2 = var_weight_k2_upper / var_weight_lower
    var_weights_k3 = var_weight_k3_upper / var_weight_lower

    # check weigh(s) minimum and maximum value(s)
    var_weight_min_k1 = np.nanmin(var_weights_k1)
    var_weight_max_k1 = np.nanmax(var_weights_k1)
    var_weight_min_k2 = np.nanmin(var_weights_k2)
    var_weight_max_k2 = np.nanmax(var_weights_k2)
    var_weight_min_k3 = np.nanmin(var_weights_k3)
    var_weight_max_k3 = np.nanmax(var_weights_k3)

    return var_weights_k1, var_weights_k2, var_weights_k3

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute weights using error variance (2 vars)
def compute_weights_k2(var_data_k1, var_data_k2, **kwargs):

    if (var_data_k1 is None) or (var_data_k2 is None):
        logging.error(' ===> All/some variable(s) are defined by NoneType')
        raise RuntimeError('All variable(s) must be defined to compute weight(s)')

    # lower part of the weight(s) computation
    var_weight_lower = var_data_k1 + var_data_k2

    # compute weight(s)
    var_weights_k1 = var_data_k2 / var_weight_lower
    var_weights_k2 = var_data_k1 / var_weight_lower

    return var_weights_k1, var_weights_k2

# ----------------------------------------------------------------------------------------------------------------------
