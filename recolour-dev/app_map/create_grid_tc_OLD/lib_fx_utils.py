"""
Library Features:

Name:          lib_fx_utils
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
from copy import deepcopy

from lib_info_args import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get fx settings
def get_fx_settings(fx_settings,
                    tag_fx_active='active', tag_fx_name='name',
                    tag_fx_variables='variables', tag_fx_parameters='parameters'):
    if tag_fx_active in list(fx_settings.keys()):
        fx_active = fx_settings[tag_fx_active]
    else:
        alg_logger.error(' ===> Fx settings "' + tag_fx_active + '" is not available in the fx obj')
        raise RuntimeError('Fx settings is needed by the fx class to correctly run')
    if tag_fx_name in list(fx_settings.keys()):
        fx_name = fx_settings[tag_fx_name]
    else:
        alg_logger.error(' ===> Fx settings "' + tag_fx_name + '" is not available in the fx obj')
        raise RuntimeError('Fx settings is needed by the fx class to correctly run')
    if tag_fx_variables in list(fx_settings.keys()):
        fx_vars = fx_settings[tag_fx_variables]
    else:
        alg_logger.error(' ===> Fx settings "' + tag_fx_variables + '" is not available in the fx obj')
        raise RuntimeError('Fx name is needed by the fx class to correctly run')
    if tag_fx_parameters in list(fx_settings.keys()):
        fx_params = fx_settings[tag_fx_parameters]
    else:
        alg_logger.error(' ===> Fx settings "' + tag_fx_parameters + '" is not available in the fx obj')
        raise RuntimeError('Fx name is needed by the fx class to correctly run')
    return fx_active, fx_name, fx_vars, fx_params
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get fx method
def get_fx_method(fx_name, fx_methods):
    if hasattr(fx_methods, fx_name):
        fx_handle = getattr(fx_methods, fx_name)
    else:
        alg_logger.error(' ===> Function "' + fx_name + '" is not available')
        raise RuntimeError('Function must be included in the methods repository')
    return fx_handle
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remove_nans (using a reference column)
def remove_nans(obj_data, reference_column='soil_moisture_weighted', keep_finite=True):
    # info start
    alg_logger.info(' -------> Remove nans data by "' + reference_column + '" .... ')

    # activate filter
    if keep_finite:

        if reference_column in list(obj_data.columns):
            n_nans = np.argwhere(np.isnan(obj_data[reference_column].values))[:, 0].size
            if n_nans == 0:
                obj_data_filtered = deepcopy(obj_data)
                alg_logger.info(' -------> Remove nans data by "' + reference_column +
                                '" ... SKIPPED. Columns nans are equal to "' + str(n_nans) + '"')
            else:
                obj_data_filtered = obj_data[np.isfinite(obj_data[reference_column])]
                alg_logger.info(' -------> Remove nans data by "' + reference_column +
                                '" ... DONE. Columns nans are equal to "' + str(n_nans) + '"')
        else:
            obj_data_filtered = deepcopy(obj_data)
            alg_logger.warning(' ===> DataFrame column "' + reference_column + '" is not available')
            alg_logger.info(' -------> Remove nans data by "' + reference_column + '" ... FAILED')

    else:
        obj_data_filtered = deepcopy(obj_data)
        alg_logger.info(' -------> Remove nans data by "' + reference_column + '" ... SKIPPED. Filter is not activated')

    return obj_data_filtered
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to check data
def check_data(obj_data, thr_data='all'):
    var_check_list = []
    for var_key, var_da in obj_data.items():
        if var_da is None:
            var_check_step = False
        else:
            var_check_step = True
        var_check_list.append(var_check_step)

    if thr_data == 'all':
        flag_data = all(var_check_list)
    elif thr_data == 'any':
        flag_data = any(var_check_list)
    else:
        alg_logger.error(' ===> Check flag "' + thr_data + '" is not supported')
        raise NotImplemented('Case not implemented yet')
    return flag_data
# ----------------------------------------------------------------------------------------------------------------------
