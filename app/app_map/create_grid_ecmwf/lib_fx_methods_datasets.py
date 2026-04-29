"""
Library Features:

Name:          lib_fx_methods_datasets
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260428'
Version:       '1.5.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings

from copy import deepcopy

import numpy as np

# debug
try:
    import matplotlib.pylab as plt
except ImportError:
    pass
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute root-zoom-soil-moisture profile
def compute_rzsm_profile(native_layers=None, weighted_layers=None, **kwargs):
    """
    Compute weighted SWI/RZSM layers from native soil layers.

    Example:
        native_layers = {
            "var40": (0, 7),
            "var41": (7, 28),
            "var42": (28, 100),
            "var43": (100, 289),
        }

        weighted_layers = {
            "var_0_10cm": (0, 10),
            "var_0_28cm": (0, 28),
            "var_0_100cm": (0, 100),
        }
    """

    if native_layers is None:
        logging.error(' ===> Native layers are not defined')
        raise RuntimeError('Native layers are needed to compute weighted layers')

    if weighted_layers is None:
        logging.error(' ===> Weighted layers are not defined')
        raise RuntimeError('Weighted layers are needed to compute weighted layers')

    data_out = {}
    for var_out, (target_top, target_bottom) in weighted_layers.items():

        weighted_sum = None
        total_weight = 0
        mask_no_data = None

        for var_name, (native_top, native_bottom) in native_layers.items():

            if var_name not in kwargs:
                logging.error(' ===> Variable "' + var_name + '" is not available in input datasets')
                raise RuntimeError('Variable is needed to compute weighted layers')

            var_data = kwargs[var_name]

            if var_data is None:
                logging.error(' ===> Variable "' + var_name + '" is defined by NoneType')
                raise RuntimeError('Variable is needed to compute weighted layers')

            overlap_top = max(target_top, native_top)
            overlap_bottom = min(target_bottom, native_bottom)
            overlap = overlap_bottom - overlap_top

            if overlap > 0:

                if weighted_sum is None:
                    weighted_sum = var_data * overlap
                    mask_no_data = np.isnan(var_data)
                else:
                    weighted_sum = weighted_sum + var_data * overlap
                    mask_no_data = mask_no_data | np.isnan(var_data)

                total_weight += overlap

        if total_weight == 0:
            logging.error(
                ' ===> No overlap found for weighted layer "' +
                var_out + '" from ' + str(target_top) + ' to ' + str(target_bottom) + ' cm'
            )
            raise RuntimeError('Layer overlap is not valid')

        var_weighted = weighted_sum / total_weight
        var_weighted[mask_no_data] = np.nan

        data_out[var_out] = var_weighted

    return data_out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute rzsm_OLD
def compute_rzsm_profile_OLD(var_data_k1, var_data_k2, var_data_k3,
                             var_name_k1='var_0_7', var_name_k2='var_0_28', var_name_k3='var_0_100',
                             **kwargs):

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

    var_out = {var_name_k1: var_data_k1, var_name_k2: var_data_k2, var_name_k3: var_data_k3}

    return var_out

# ----------------------------------------------------------------------------------------------------------------------
