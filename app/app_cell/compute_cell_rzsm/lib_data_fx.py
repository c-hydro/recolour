"""
Library Features:

Name:          lib_data_fx
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240505'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
from copy import deepcopy

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to prepare data
def prepare_data(var_layer_obj, var_name_list, var_layer_template='var_data_k{:}'):
    var_layer_data = {}
    for i, var_name_i in enumerate(var_name_list):
        var_layer_i = var_layer_template.format(i + 1)
        var_layer_data[var_layer_i] = var_layer_obj[var_name_i].values
    return var_layer_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute rzsm layer_1
def compute_rzsm_layer_1(var_data_k1, **kwargs):

    # check variable(s)
    if var_data_k1 is None:
        logging.error(' ===> All/some variable(s) are defined by NoneType')
        raise RuntimeError('All variable(s) must be defined to compute weight(s)')

    # get no data value
    if 'no_data' in kwargs:
        no_data = kwargs['no_data']
    else:
        no_data = np.nan

    # var undefined
    var_nans = np.isnan(var_data_k1)

    # var_0_7 = var40
    var_rzsm = deepcopy(var_data_k1)

    # mask no data
    var_rzsm[var_nans] = no_data
    var_rzsm[var_data_k1 == no_data] = no_data

    return var_rzsm

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute rzsm layer_2
def compute_rzsm_layer_2(var_data_k1, var_data_k2, **kwargs):

    # check variable(s)
    if (var_data_k1 is None) or (var_data_k2 is None):
        logging.error(' ===> All/some variable(s) are defined by NoneType')
        raise RuntimeError('All variable(s) must be defined to compute weight(s)')

    # get no data value
    if 'no_data' in kwargs:
        no_data = kwargs['no_data']
    else:
        no_data = np.nan

    # var undefined
    var_nans = np.isnan(var_data_k1)

    # var_0_28 = 0.25*var40 + 0.75*var41
    var_rzsm = 0.25 * var_data_k1 + 0.75 * var_data_k2

    # mask no data
    var_rzsm[var_nans] = no_data
    var_rzsm[var_data_k1 == no_data] = no_data
    var_rzsm[var_data_k2 == no_data] = no_data

    return var_rzsm

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute rzsm layer_3
def compute_rzsm_layer_3(var_data_k1, var_data_k2, var_data_k3, **kwargs):

    # check variable(s)
    if (var_data_k1 is None) or (var_data_k2 is None) or (var_data_k3 is None):
        logging.error(' ===> All/some variable(s) are defined by NoneType')
        raise RuntimeError('All variable(s) must be defined to compute weight(s)')

    # get no data value
    if 'no_data' in kwargs:
        no_data = kwargs['no_data']
    else:
        no_data = np.nan

    # var undefined
    var_nans = np.isnan(var_data_k1)

    # var_0_100 = 0.07*var40 + 0.21*var41 + 0.72*var42
    var_rzsm = 0.07 * var_data_k1 + 0.21 * var_data_k2 + 0.72 * var_data_k3

    # mask no data
    var_rzsm[var_nans] = no_data
    var_rzsm[var_data_k1 == no_data] = no_data
    var_rzsm[var_data_k2 == no_data] = no_data
    var_rzsm[var_data_k3 == no_data] = no_data

    return var_rzsm

# ----------------------------------------------------------------------------------------------------------------------
