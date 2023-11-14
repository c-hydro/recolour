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
import numpy as np

from copy import deepcopy

from lib_info_args import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to scale data
def scale_data(obj_data, variables_data=None, parameters_data=None,
               var_name_k1='soil_moisture_k1', var_name_k2='soil_moisture_k2',
               var_beta_k1='beta_k1', var_beta_k2='beta_k2',
               var_mean_ref='mean_ref', var_mean_k1='mean_k1', var_mean_k2='mean_k2', **kwargs):

    vars_name_expected = [var_name_k1, var_name_k2]
    pars_name_expected = [var_beta_k1, var_beta_k2, var_mean_ref, var_mean_k1, var_mean_k2]

    obj_name = list(obj_data.columns)
    for vars_name_step in vars_name_expected:
        if vars_name_step not in obj_name:
            alg_logger.error(' ===> Variable "' + vars_name_step + '" must be included in the source DataFrame')
            raise RuntimeError('Check the source DataFrame to correctly run the method')
    for pars_name_step in pars_name_expected:
        if pars_name_step not in obj_name:
            alg_logger.error(' ===> Parameter "' + pars_name_step + '" must be included in the source DataFrame')
            raise RuntimeError('Check the source DataFrame to correctly run the method')

    pars_tmp = {}
    for par_name_fx, par_name_obj in parameters_data.items():
        if par_name_obj in list(obj_data.columns):
            pars_tmp[par_name_fx] = obj_data[par_name_obj].values
        else:
            alg_logger.error(' ===> Parameter "' + par_name_obj + '" is not available in the dataframe obj')
            raise RuntimeError('Parameter is needed by the algorithm to properly work')

    vars_tmp_in = {}
    for var_name_fx, var_name_obj in variables_data['in'].items():
        if var_name_obj in list(obj_data.columns):
            vars_tmp_in[var_name_fx] = obj_data[var_name_obj].values
        else:
            alg_logger.error(' ===> Variable "' + var_name_obj + '" is not available in the dataframe obj')
            raise RuntimeError('Variable is needed by the algorithm to properly work')

    var_scale_k1 = (pars_tmp['var_beta_k1'] * (vars_tmp_in['var_data_k1'] - pars_tmp['var_mean_k1']) +
                    pars_tmp['var_mean_ref'])
    var_scale_k2 = (pars_tmp['var_beta_k2'] * (vars_tmp_in['var_data_k2'] - pars_tmp['var_mean_k2']) +
                    pars_tmp['var_mean_ref'])

    var_scale_k1[var_scale_k1 > 1] = 1
    var_scale_k2[var_scale_k2 > 1] = 1

    vars_tmp_out = {'var_data_k1_scaled': var_scale_k1, 'var_data_k2_scaled': var_scale_k2}

    for var_name_fx, var_name_obj in variables_data['out'].items():
        if var_name_fx in list(vars_tmp_out.keys()):
            obj_data[var_name_obj] = vars_tmp_out[var_name_fx]
        else:
            alg_logger.error(' ===> Variable "' + var_name_fx + '" is not available in the temporary obj')
            raise RuntimeError('Variable is needed by the algorithm to properly work')

    return obj_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to weigh data
def weigh_data(obj_data, variables_data=None, parameters_data=None,
               var_name_ref='soil_moisture_ref', var_name_composite='soil_moisture_composite',
               var_name_k1='soil_moisture_k1_scaled', var_name_k2='soil_moisture_k2_scaled',
               var_flag_composite='flags_composite',
               var_weight_ref='weights_ref', var_weight_k1='weights_k1', var_weight_k2='weights_k2',
               var_weight_ref_k1='weights_ref_k1', var_weight_ref_k2='weights_ref_k2',
               var_weight_k1_ref='weights_k1_ref', var_weight_k2_ref='weights_k2_ref',
               flag_ref_k1_k2=0, flag_ref=1, flag_ref_k1=2, flag_ref_k2=3, flag_k1=4, flag_k2=5,
               active_ref_k1=True, active_ref_k2=False, active_ref=True, active_k1=False, active_k2=False,
               **kwargs):

    vars_name_expected = [var_name_ref, var_name_k1, var_name_k2, var_name_composite, var_flag_composite]
    pars_name_expected = [var_weight_ref, var_weight_k1, var_weight_k2,
                          var_weight_ref_k1, var_weight_ref_k2, var_weight_k1_ref, var_weight_k2_ref]

    obj_name = list(obj_data.columns)
    for vars_name_step in vars_name_expected:
        if vars_name_step not in obj_name:
            alg_logger.error(' ===> Variable "' + vars_name_step + '" must be included in the source DataFrame')
            raise RuntimeError('Check the source DataFrame to correctly run the method')
    for pars_name_step in pars_name_expected:
        if pars_name_step not in obj_name:
            alg_logger.error(' ===> Parameter "' + pars_name_step + '" must be included in the source DataFrame')
            raise RuntimeError('Check the source DataFrame to correctly run the method')

    pars_tmp = {}
    for par_name_fx, par_name_obj in parameters_data.items():
        if par_name_obj in list(obj_data.columns):
            pars_tmp[par_name_fx] = obj_data[par_name_obj].values
        else:
            alg_logger.error(' ===> Parameter "' + par_name_obj + '" is not available in the dataframe obj')
            raise RuntimeError('Parameter is needed by the algorithm to properly work')

    vars_tmp_in = {}
    for var_name_fx, var_name_obj in variables_data['in'].items():
        if var_name_obj in list(obj_data.columns):
            vars_tmp_in[var_name_fx] = obj_data[var_name_obj].values
        else:
            alg_logger.error(' ===> Variable "' + var_name_obj + '" is not available in the dataframe obj')
            raise RuntimeError('Variable is needed by the algorithm to properly work')

    var_tmp_ref = vars_tmp_in['var_data_ref']
    var_tmp_k1, var_tmp_k2 = vars_tmp_in['var_data_k1_scaled'], vars_tmp_in['var_data_k2_scaled']
    par_tmp_ref = pars_tmp['var_weights_ref']
    par_tmp_k1, par_tmp_k2 = pars_tmp['var_weights_k1'], pars_tmp['var_weights_k2']
    par_tmp_ref_k1, par_tmp_k1_ref = pars_tmp['var_weights_ref_k1'], pars_tmp['var_weights_k1_ref']
    par_tmp_ref_k2, par_tmp_k2_ref = pars_tmp['var_weights_ref_k2'], pars_tmp['var_weights_k2_ref']

    var_tmp_composite, flag_tmp_composite = None, None
    if ('var_data_composite' in list(vars_tmp_in.keys())) and ('var_flags_composite' in list(vars_tmp_in.keys())):
        var_tmp_composite = vars_tmp_in['var_data_composite']
        flag_tmp_composite = vars_tmp_in['var_flags_composite']

    idxs_tmp_fk1_fk2 = np.argwhere(np.isfinite(var_tmp_k1) & np.isfinite(var_tmp_k2))[:, 0]
    idxs_tmp_nk1_nk2 = np.argwhere(np.isnan(var_tmp_k1) & np.isnan(var_tmp_k2))[:, 0]
    idxs_tmp_fk1_nk2 = np.argwhere(np.isfinite(var_tmp_k1) & np.isnan(var_tmp_k2))[:, 0]
    idxs_tmp_nk1_fk2 = np.argwhere(np.isnan(var_tmp_k1) & np.isfinite(var_tmp_k2))[:, 0]

    # case 1
    var_tmp_ref_fk1_fk2, var_tmp_k1_fk1_fk2, var_tmp_k2_fk1_fk2, var_tmp_flag_fk1_fk2 = apply_idxs(
        var_tmp_ref, var_tmp_k1, var_tmp_k2, idxs_tmp_fk1_fk2, flag_select=flag_ref_k1_k2)

    par_tmp_ref_fk1_fk2, par_tmp_k1_fk1_fk2, par_tmp_k2_fk1_fk2, _ = apply_idxs(
        par_tmp_ref, par_tmp_k1, par_tmp_k2, idxs_tmp_fk1_fk2, flag_select=flag_ref_k1_k2)

    var_tmp_weigh_ref_fk1_fk2 = var_tmp_ref_fk1_fk2 * par_tmp_ref_fk1_fk2
    var_tmp_weigh_k1_fk1_fk2 = var_tmp_k1_fk1_fk2 * par_tmp_k1_fk1_fk2
    var_tmp_weigh_k2_fk1_fk2 = var_tmp_k2_fk1_fk2 * par_tmp_k2_fk1_fk2

    var_def_weight_fk1_fk2 = var_tmp_weigh_ref_fk1_fk2 + var_tmp_weigh_k1_fk1_fk2 + var_tmp_weigh_k2_fk1_fk2

    # case 2
    var_tmp_ref_nk1_nk2, var_tmp_k1_nk1_nk2, var_tmp_k2_nk1_nk2, var_tmp_flag_nk1_nk2 = apply_idxs(
        var_tmp_ref, var_tmp_k1, var_tmp_k2, idxs_tmp_nk1_nk2, flag_select=flag_ref)

    var_def_weight_nk1_nk2 = deepcopy(var_tmp_ref_nk1_nk2)

    # case 3
    var_tmp_ref_fk1_nk2, var_tmp_k1_fk1_nk2, var_tmp_k2_fk1_nk2, var_tmp_flag_fk1_nk2 = apply_idxs(
        var_tmp_ref, var_tmp_k1, var_tmp_k2, idxs_tmp_fk1_nk2, flag_select=flag_ref_k1)

    par_tmp_ref_k1_fk1_nk2, par_tmp_k1_ref_fk1_nk2, _, _ = apply_idxs(
        par_tmp_ref_k1, par_tmp_k1_ref, None, idxs_tmp_fk1_nk2, flag_select=flag_ref_k1)

    var_tmp_weigh_ref_fk1_nk2 = var_tmp_ref_fk1_nk2 * par_tmp_ref_k1_fk1_nk2
    var_tmp_weigh_k1_fk1_nk2 = var_tmp_k1_fk1_nk2 * par_tmp_k1_ref_fk1_nk2

    var_def_weight_fk1_nk2 = var_tmp_weigh_ref_fk1_nk2 + var_tmp_weigh_k1_fk1_nk2

    # case 4
    var_tmp_ref_nk1_fk2, var_tmp_k1_nk1_fk2, var_tmp_k2_nk1_fk2, var_tmp_flag_nk1_fk2 = apply_idxs(
        var_tmp_ref, var_tmp_k1, var_tmp_k2, idxs_tmp_nk1_fk2, flag_select=flag_ref_k2)

    par_tmp_ref_k2_fk1_nk2, _, par_tmp_ref_k2_fk1_nk2, _ = apply_idxs(
        par_tmp_ref_k2, None, par_tmp_k2_ref, idxs_tmp_nk1_fk2, flag_select=flag_ref_k2)

    var_tmp_weigh_ref_nk1_fk2 = var_tmp_ref_nk1_fk2 * par_tmp_ref_k2_fk1_nk2
    var_tmp_weigh_k2_nk1_fk2 = var_tmp_k2_nk1_fk2 * par_tmp_ref_k2_fk1_nk2

    var_def_weight_nk1_fk2 = var_tmp_weigh_ref_nk1_fk2 + var_tmp_weigh_k2_nk1_fk2

    # define out variable(s)
    var_weight_flag = np.zeros(shape=[var_tmp_ref.shape[0]], dtype=int)
    var_weight_flag[:] = -1
    var_weight_flag[idxs_tmp_fk1_fk2] = var_tmp_flag_fk1_fk2[[idxs_tmp_fk1_fk2]]
    if active_ref:
        var_weight_flag[idxs_tmp_nk1_nk2] = var_tmp_flag_nk1_nk2[[idxs_tmp_nk1_nk2]]
    if active_ref_k1:
        var_weight_flag[idxs_tmp_fk1_nk2] = var_tmp_flag_fk1_nk2[[idxs_tmp_fk1_nk2]]
    if active_ref_k2:
        var_weight_flag[idxs_tmp_nk1_fk2] = var_tmp_flag_nk1_fk2[[idxs_tmp_nk1_fk2]]

    var_weight_def = np.zeros(shape=[var_tmp_ref.shape[0]])
    var_weight_def[:] = np.nan
    var_weight_def[idxs_tmp_fk1_fk2] = var_def_weight_fk1_fk2[[idxs_tmp_fk1_fk2]]
    if active_ref:
        var_weight_def[idxs_tmp_nk1_nk2] = var_def_weight_nk1_nk2[[idxs_tmp_nk1_nk2]]
    if active_ref_k1:
        var_weight_def[idxs_tmp_fk1_nk2] = var_def_weight_fk1_nk2[[idxs_tmp_fk1_nk2]]
    if active_ref_k2:
        var_weight_def[idxs_tmp_nk1_fk2] = var_def_weight_nk1_fk2[[idxs_tmp_nk1_fk2]]

    var_weight_ref = np.zeros(shape=[var_tmp_ref.shape[0]])
    var_weight_ref[:] = np.nan
    var_weight_ref[idxs_tmp_fk1_fk2] = var_tmp_weigh_ref_fk1_fk2[[idxs_tmp_fk1_fk2]]
    if active_ref:
        var_weight_ref[idxs_tmp_nk1_nk2] = var_def_weight_nk1_nk2[[idxs_tmp_nk1_nk2]]
    if active_ref_k1:
        var_weight_ref[idxs_tmp_fk1_nk2] = var_tmp_weigh_ref_fk1_nk2[[idxs_tmp_fk1_nk2]]
    if active_ref_k2:
        var_weight_ref[idxs_tmp_nk1_fk2] = var_tmp_weigh_ref_nk1_fk2[[idxs_tmp_nk1_fk2]]

    var_weight_k1 = np.zeros(shape=[var_tmp_ref.shape[0]])
    var_weight_k1[:] = np.nan
    var_weight_k1[idxs_tmp_fk1_fk2] = var_tmp_weigh_k1_fk1_fk2[[idxs_tmp_fk1_fk2]]
    var_weight_k1[idxs_tmp_fk1_nk2] = var_tmp_weigh_k1_fk1_nk2[[idxs_tmp_fk1_nk2]]

    var_weight_k2 = np.zeros(shape=[var_tmp_ref.shape[0]])
    var_weight_k2[:] = np.nan
    var_weight_k2[idxs_tmp_fk1_fk2] = var_tmp_weigh_k2_fk1_fk2[[idxs_tmp_fk1_fk2]]
    var_weight_k2[idxs_tmp_nk1_fk2] = var_tmp_weigh_k2_nk1_fk2[[idxs_tmp_nk1_fk2]]

    # case to fill data using only k1 or k2 value(s)
    if (var_tmp_composite is not None) and (flag_tmp_composite is not None):
        idxs_nans_def = np.squeeze(np.argwhere(np.isnan(var_weight_def))).tolist()
        if not isinstance(idxs_nans_def, list): idxs_nans_def = [idxs_nans_def]

        var_def_composite = var_tmp_composite[idxs_nans_def]
        flag_def_composite = flag_tmp_composite[idxs_nans_def]

        # case 5 (in the composite k1 flag == 1)
        idxs_k1_local = np.squeeze(np.argwhere(flag_def_composite == 1))
        idxs_k1_global = np.array(idxs_nans_def)[idxs_k1_local].tolist()

        if active_k1:
            var_weight_def[idxs_k1_global] = var_def_composite[idxs_k1_local]
            var_weight_flag[idxs_k1_global] = flag_k1

        # case 6 (in the composite k2 flag == 2)
        idxs_k2_local = np.squeeze(np.argwhere(flag_def_composite == 2))
        idxs_k2_global = np.array(idxs_nans_def)[idxs_k2_local].tolist()

        if active_k2:
            var_weight_def[idxs_k2_global] = var_def_composite[idxs_k2_local]
            var_weight_flag[idxs_k2_global] = flag_k2

    # check limits (after weights)
    var_weight_def[var_weight_def > 1] = 1
    var_weight_def[var_weight_def < 0] = 0

    idxs_nans_def = np.squeeze(np.argwhere(np.isnan(var_weight_def))).tolist()
    if not isinstance(idxs_nans_def, list):
        idxs_nans_def = [idxs_nans_def]
    idxs_nans_n = idxs_nans_def.__len__()
    if idxs_nans_def.__len__() > 0:
        alg_logger.warning(' ===> Dataset weighted have "' + str(idxs_nans_n) + '" nan values. ')

    vars_tmp_out = {
        'var_part_ref_weighted': var_tmp_ref,
        'var_part_k1_weighted': var_weight_k1, 'var_part_k2_weighted': var_weight_k2,
        'var_data_weighted': var_weight_def, 'var_data_flag': var_weight_flag}

    for var_name_fx, var_name_obj in variables_data['out'].items():
        if var_name_fx in list(vars_tmp_out.keys()):
            obj_data[var_name_obj] = vars_tmp_out[var_name_fx]
        else:
            alg_logger.error(' ===> Variable "' + var_name_fx + '" is not available in the temporary obj')
            raise RuntimeError('Variable is needed by the algorithm to properly work')

    return obj_data

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to apply idxs to array(s)
def apply_idxs(var_raw_ref, var_raw_k1, var_raw_k2, idxs_select, flag_select=0):

    var_sel_ref = np.zeros(shape=[var_raw_ref.shape[0]])
    var_sel_ref[:] = np.nan
    var_sel_ref[idxs_select] = var_raw_ref[idxs_select]

    var_sel_k1 = None
    if var_raw_k1 is not None:
        var_sel_k1 = np.zeros(shape=[var_raw_ref.shape[0]])
        var_sel_k1[:] = np.nan
        var_sel_k1[idxs_select] = var_raw_k1[idxs_select]

    var_sel_k2 = None
    if var_raw_k2 is not None:
        var_sel_k2 = np.zeros(shape=[var_raw_ref.shape[0]])
        var_sel_k2[:] = np.nan
        var_sel_k2[idxs_select] = var_raw_k2[idxs_select]

    var_sel_flag = np.zeros(shape=[var_raw_ref.shape[0]])
    var_sel_flag[:] = np.nan
    var_sel_flag[idxs_select] = int(flag_select)

    return var_sel_ref, var_sel_k1, var_sel_k2, var_sel_flag
# ----------------------------------------------------------------------------------------------------------------------
