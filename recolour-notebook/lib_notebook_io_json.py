"""
Library Features:

Name:          lib_notebook_io_json
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20200401'
Version:       '3.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import json
import os
import pandas as pd
import numpy as np

from copy import deepcopy
from datetime import datetime
from scipy import stats

import pytesmo.scaling as scaling
import pytesmo.metrics as metrics

from pytesmo.time_series.filters import exp_filter, boxcar_filter
# ----------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to add seasonal description
def add_dframe_seasons(dframe_data, column_season='season', lut_season=None):

    dframe_time = dframe_data.index

    grp_season = [lut_season.get(pd.Timestamp(t_stamp).month) for t_stamp in dframe_time]
    dframe_data[column_season] = grp_season

    return dframe_data

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to print time series metrics
def print_time_series_metrics(metrics_obj):
    for metrics_datasets, metrics_fields in metrics_obj.items():
        print(' --- Datasets: "' + metrics_datasets + '"')
        for metrics_name, metrics_value in metrics_fields.items():
            print(' ---- ' + metrics_name + ' = ' + metrics_value)
        print(' ')
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to compute time series metrics
def compute_time_series_metrics(ts_dframe_generic, ts_reference='sm_obs'):

    x_df = ts_dframe_generic[ts_reference].to_frame()

    ts_list = list(ts_dframe_generic.columns)
    if ts_reference in ts_list:
        ts_list.remove(ts_reference)
    else:
        raise RuntimeError('Reference field "' + ts_reference + '" must be included in the source DataFrame')

    if 'control' in ts_list:
        ts_list.remove('control')

    metrics_obj = {}
    for ts_name in ts_list:

        y_df = ts_dframe_generic[ts_name].to_frame()

        xy_df = x_df.join(y_df)
        xy_df = xy_df.dropna()

        x_values = xy_df[ts_reference].values
        y_values = xy_df[ts_name].values

        pearson_r, person_p = stats.pearsonr(x_values, y_values)
        pearson_r, person_p = '{:.2f}'.format(pearson_r), '{:.2e}'.format(person_p)
        spearman_rho, spearman_p = stats.spearmanr(x_values, y_values)
        spearman_rho, spearman_p = '{:.2f}'.format(spearman_rho), '{:.2e}'.format(spearman_p)
        kendall_tau, kendall_p = stats.kendalltau(x_values, y_values)
        kendall_tau, kendall_p = '{:.2f}'.format(kendall_tau), '{:.2e}'.format(kendall_p)
        rmsd = '{:.2f}'.format(metrics.rmsd(x_values, y_values))
        bias = '{:.2f}'.format(metrics.bias(x_values, y_values))
        nash_sutcliffe = '{:.2f}'.format(metrics.nash_sutcliffe(x_values, y_values))

        # Calculate correlation coefficients, RMSD, bias, Nash Sutcliffe
        metrics_ts = {'Pearson R': pearson_r, 'Pearson p': person_p,
                      'Spearman rho': spearman_rho, 'Spearman p': spearman_p,
                      'Kendall tau': kendall_tau, 'Kendall p': kendall_p,
                      'RMSD': rmsd,
                      'Bias': bias,
                      'Nash Sutcliffe': nash_sutcliffe}

        metrics_obj[ts_name] = metrics_ts

    return metrics_obj
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to scale time series
def scale_time_series(ts_dframe_generic, ts_reference='sm_obs', ts_scale_method='cdf_beta_match'):

    if 'season' in list(ts_dframe_generic.columns):
        ts_dframe_generic = ts_dframe_generic.drop(columns=['season'])
    if 'control' in list(ts_dframe_generic.columns):
        ts_dframe_generic = ts_dframe_generic.drop(columns=['control'])

    ts_fields = list(ts_dframe_generic.columns)
    if ts_reference in ts_fields:
        ts_fields.remove(ts_reference)
    else:
        raise RuntimeError('Reference field "' + ts_reference + '" must be included in the source DataFrame')

    ts_index = ts_dframe_generic.index.values
    ts_n = ts_index.__len__()
    ts_data = {ts_reference: ts_dframe_generic[ts_reference].values}
    ts_dframe_reference = pd.DataFrame(data=ts_data, index=ts_index)

    for ts_other in ts_fields:

        ts_dframe_ref = ts_dframe_generic[ts_reference]
        ts_dframe_other = ts_dframe_generic[ts_other]

        ts_series_ref = ts_dframe_ref.dropna()
        ts_dframe_ref = ts_series_ref.to_frame()

        ts_series_other = ts_dframe_other.dropna()
        ts_dframe_other = ts_series_other.to_frame()

        if ts_series_ref.__len__() > ts_series_other.__len__():
            ts_dframe_tmp = ts_dframe_other.join(ts_dframe_ref)
        elif ts_series_ref.__len__() < ts_series_other.__len__():
            ts_dframe_tmp = ts_dframe_ref.join(ts_dframe_other)
        else:
            ts_dframe_tmp = ts_dframe_ref.join(ts_dframe_other)

        ts_dframe_tmp = ts_dframe_tmp.dropna()

        tmp_fields = list(ts_dframe_tmp.columns)
        tmp_idx = tmp_fields.index(ts_reference)

        # scale datasets using a method defined in pytestmo library
        ts_dframe_scaled = scaling.scale(ts_dframe_tmp, method=ts_scale_method, reference_index=tmp_idx)
        tmp_dframe_scaled = ts_dframe_scaled[ts_other].to_frame()

        ts_dframe_reference = ts_dframe_reference.join(tmp_dframe_scaled)

    # add control series
    ts_control = np.zeros(shape=[ts_n])
    ts_control[:] = 0
    ts_dframe_control = pd.DataFrame(data={'control':ts_control}, index=ts_index)
    ts_dframe_reference = ts_dframe_reference.join(ts_dframe_control)

    return ts_dframe_reference
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to filter time series by season
def filter_time_series_by_season(ts_dframe_generic, season_filter='DJF,MAM,JJA,SON'):

    # season filter
    season_filter = season_filter.replace(',', '_')
    season_filter = season_filter.split(',')

    if not isinstance(season_filter, list):
        season_filter = [season_filter]

    # Get the seasonal dframe
    if season_filter == ['ALL']:
        ts_dframe_filtered = deepcopy(ts_dframe_generic)
    else:
        ts_dframe_filtered = ts_dframe_generic[ts_dframe_generic['season'].isin(season_filter)]

    return ts_dframe_filtered
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to join time series
def join_time_series(ts_dframe_product, ts_dframe_datasets):

    if 'season' in list(ts_dframe_datasets.columns):
        ts_dframe_datasets = ts_dframe_datasets.drop(columns=['season'])
    ts_dframe_common = ts_dframe_product.join(ts_dframe_datasets)

    return ts_dframe_common

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read time series datasets
def read_time_series_datasets(
        file_name, file_sep=',', file_time_format='%Y-%m-%d %H:%M',
        file_min_value=0, file_max_value=100,
        flag_season_lut=True,
        flag_ts_ecmwf_layer_1='sm_ecmwf_layer_0_7',
        flag_ts_ecmwf_layer_2='sm_ecmwf_layer_0_28',
        flag_ts_ecmwf_layer_3='sm_ecmwf_layer_0_100',
        flag_ts_hmc='sm_hmc', flag_ts_smap='sm_smap',
        flag_ts_control='control',
        flag_time='time'):

    # season filter
    if flag_season_lut:
        file_season_lut = {
            1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA',
            7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON', 12: 'DJF'}
    else:
        file_season_lut = {
            1: 'ALL', 2: 'ALL', 3: 'ALL', 4: 'ALL', 5: 'ALL', 6: 'ALL',
            7: 'ALL', 8: 'ALL', 9: 'ALL', 10: 'ALL', 11: 'ALL', 12: 'ALL'}

    with open(file_name, 'r') as file_handle:
        file_data = json.load(file_handle)

    if flag_time in list(file_data.keys()):
        time_str = file_data[flag_time]
    else:
        print(' ===> Flag "' + flag_time + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')

    if flag_ts_ecmwf_layer_1 in list(file_data.keys()):
        ts_ecmwf_layer_1_str = file_data[flag_ts_ecmwf_layer_1]
    else:
        print(' ===> Flag "' + flag_ts_ecmwf_layer_1 + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')
    if flag_ts_ecmwf_layer_2 in list(file_data.keys()):
        ts_ecmwf_layer_2_str = file_data[flag_ts_ecmwf_layer_2]
    else:
        print(' ===> Flag "' + flag_ts_ecmwf_layer_2 +
              '" not available in the source file. Initialize with Nans time-series')
        ts_ecmwf_layer_2_str = None
    if flag_ts_ecmwf_layer_3 in list(file_data.keys()):
        ts_ecmwf_layer_3_str = file_data[flag_ts_ecmwf_layer_3]
    else:
        print(' ===> Flag "' + flag_ts_ecmwf_layer_3 +
              '" not available in the source file. Initialize with Nans time-series')
        ts_ecmwf_layer_3_str = None

    if flag_ts_hmc in list(file_data.keys()):
        ts_hmc_str = file_data[flag_ts_hmc]
    else:
        print(' ===> Flag "' + flag_ts_hmc + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')

    if flag_ts_smap in list(file_data.keys()):
        ts_smap_str = file_data[flag_ts_smap]
    else:
        print(' ===> Flag "' + flag_ts_smap + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')

    time_list = time_str.split(file_sep)
    time_n = time_list.__len__()

    ts_control = np.zeros(shape=[time_n])
    ts_control[:] = 0

    ts_ecmwf_layer_1_list = ts_ecmwf_layer_1_str.split(file_sep)
    if ts_ecmwf_layer_2_str is None:
        ts_ecmwf_layer_2_list = [-9999] * time_n
    else:
        ts_ecmwf_layer_2_list = ts_ecmwf_layer_2_str.split(file_sep)
    if ts_ecmwf_layer_3_str is None:
        ts_ecmwf_layer_3_list = [-9999] * time_n
    else:
        ts_ecmwf_layer_3_list = ts_ecmwf_layer_3_str.split(file_sep)

    ts_hmc_list = ts_hmc_str.split(file_sep)
    ts_smap_list = ts_smap_str.split(file_sep)

    ts_time, ts_ecmwf_l1, ts_ecmwf_l2, ts_ecmwf_l3, ts_hmc, ts_smap = [], [], [], [], [], []
    for time_step, e_l1_step, e_l2_step, e_l3_step, hmc_step, smap_step in zip(
            time_list,
            ts_ecmwf_layer_1_list, ts_ecmwf_layer_2_list, ts_ecmwf_layer_3_list,
            ts_hmc_list, ts_smap_list):

        ts_time.append(datetime.strptime(time_step, file_time_format).date())
        ts_ecmwf_l1.append(float(e_l1_step))
        ts_ecmwf_l2.append(float(e_l2_step))
        ts_ecmwf_l3.append(float(e_l3_step))
        ts_hmc.append(float(hmc_step))
        ts_smap.append(float(smap_step))

    ts_data_complete = {flag_ts_ecmwf_layer_1: ts_ecmwf_l1,
                        flag_ts_ecmwf_layer_2: ts_ecmwf_l2,
                        flag_ts_ecmwf_layer_3: ts_ecmwf_l3,
                        flag_ts_hmc: ts_hmc,
                        flag_ts_smap: ts_smap,
                        flag_ts_control: ts_control
                        }

    ts_data_ecmwf_l1 = {flag_ts_ecmwf_layer_1: ts_ecmwf_l1}
    ts_data_ecmwf_l2 = {flag_ts_ecmwf_layer_2: ts_ecmwf_l2}
    ts_data_ecmwf_l3 = {flag_ts_ecmwf_layer_3: ts_ecmwf_l3}
    ts_data_hmc = {flag_ts_hmc: ts_hmc}
    ts_data_smap = {flag_ts_smap: ts_smap}

    ts_index = pd.DatetimeIndex(ts_time)

    ts_dframe_complete = pd.DataFrame(data=ts_data_complete, index=ts_index)
    ts_dframe_ecmwf_l1 = pd.DataFrame(data=ts_data_ecmwf_l1, index=ts_index)
    ts_dframe_ecmwf_l2 = pd.DataFrame(data=ts_data_ecmwf_l2, index=ts_index)
    ts_dframe_ecmwf_l3 = pd.DataFrame(data=ts_data_ecmwf_l3, index=ts_index)
    ts_dframe_hmc = pd.DataFrame(data=ts_data_hmc, index=ts_index)
    ts_dframe_smap = pd.DataFrame(data=ts_data_smap, index=ts_index)

    if file_min_value is not None:
        ts_dframe_complete[ts_dframe_complete < file_min_value] = np.nan
        ts_dframe_ecmwf_l1[ts_dframe_ecmwf_l1 < file_min_value] = np.nan
        ts_dframe_ecmwf_l2[ts_dframe_ecmwf_l2 < file_min_value] = np.nan
        ts_dframe_ecmwf_l3[ts_dframe_ecmwf_l3 < file_min_value] = np.nan
        ts_dframe_hmc[ts_dframe_hmc < file_min_value] = np.nan
        ts_dframe_smap[ts_dframe_smap < file_min_value] = np.nan

    if file_max_value is not None:
        ts_dframe_complete[ts_dframe_complete > file_max_value] = np.nan
        ts_dframe_ecmwf_l1[ts_dframe_ecmwf_l1 > file_max_value] = np.nan
        ts_dframe_ecmwf_l2[ts_dframe_ecmwf_l2 > file_max_value] = np.nan
        ts_dframe_ecmwf_l3[ts_dframe_ecmwf_l3 > file_max_value] = np.nan
        ts_dframe_hmc[ts_dframe_hmc > file_max_value] = np.nan
        ts_dframe_smap[ts_dframe_smap > file_max_value] = np.nan

    ts_dframe_complete = add_dframe_seasons(ts_dframe_complete, column_season='season', lut_season=file_season_lut)
    ts_dframe_ecmwf_l1 = add_dframe_seasons(ts_dframe_ecmwf_l1, column_season='season', lut_season=file_season_lut)
    ts_dframe_ecmwf_l2 = add_dframe_seasons(ts_dframe_ecmwf_l2, column_season='season', lut_season=file_season_lut)
    ts_dframe_ecmwf_l3 = add_dframe_seasons(ts_dframe_ecmwf_l3, column_season='season', lut_season=file_season_lut)
    ts_dframe_hmc = add_dframe_seasons(ts_dframe_hmc, column_season='season', lut_season=file_season_lut)
    ts_dframe_smap = add_dframe_seasons(ts_dframe_smap, column_season='season', lut_season=file_season_lut)

    #ts_dframe_complete = ts_dframe_complete.dropna(how='all')
    #ts_dframe_ecmwf_l1 = ts_dframe_ecmwf_l1.dropna(how='any')
    #ts_dframe_ecmwf_l2 = ts_dframe_ecmwf_l2.dropna(how='any')
    #ts_dframe_ecmwf_l3 = ts_dframe_ecmwf_l3.dropna(how='any')
    #ts_dframe_hmc = ts_dframe_hmc.dropna(how='any')
    #ts_dframe_smap = ts_dframe_smap.dropna(how='any')

    return (ts_dframe_complete,
            ts_dframe_ecmwf_l1, ts_dframe_ecmwf_l2, ts_dframe_ecmwf_l3,
            ts_dframe_hmc, ts_dframe_smap)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read time series product
def read_time_series_product(
    file_name, file_sep=',', file_time_format='%Y-%m-%d %H:%M',
    file_min_value=0, file_max_value=100,
    flag_season_lut=True,
    flag_ts_obs='sm_obs', flag_ts_mod='sm_tc', flag_time='time'):
    
    # season filter
    if flag_season_lut:
        file_season_lut = {
            1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA',
            7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON', 12: 'DJF'}
    else:
        file_season_lut = {
            1: 'ALL', 2: 'ALL', 3: 'ALL', 4: 'ALL', 5: 'ALL', 6: 'ALL',
            7: 'ALL', 8: 'ALL', 9: 'ALL', 10: 'ALL', 11: 'ALL', 12: 'ALL'}

    with open(file_name, 'r') as file_handle:
        file_data = json.load(file_handle)

    if flag_time in list(file_data.keys()):
        time_str = file_data[flag_time]
    else:
        print(' ===> Flag "' + flag_time + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')
    if flag_ts_obs in list(file_data.keys()):
        ts_obs_str = file_data[flag_ts_obs]
    else:
        print(' ===> Flag "' + flag_ts_obs + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')
    if flag_ts_mod in list(file_data.keys()):
        ts_mod_str = file_data[flag_ts_mod]
    else:
        print(' ===> Flag "' + flag_ts_mod + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')

    time_list = time_str.split(file_sep)
    ts_obs_list = ts_obs_str.split(file_sep)
    ts_mod_list = ts_mod_str.split(file_sep)

    ts_time, ts_obs, ts_mod = [], [], []
    for time_step, obs_step, mod_step in zip(time_list, ts_obs_list, ts_mod_list):
        ts_time.append(datetime.strptime(time_step, file_time_format).date())
        ts_obs.append(float(obs_step))
        ts_mod.append(float(mod_step))

    ts_data_complete = {flag_ts_obs: ts_obs, flag_ts_mod: ts_mod}
    ts_data_obs = {flag_ts_obs: ts_obs}
    ts_data_mod = {flag_ts_mod: ts_mod}
    ts_index = pd.DatetimeIndex(ts_time)

    ts_dframe_complete = pd.DataFrame(data=ts_data_complete, index=ts_index)
    ts_dframe_obs = pd.DataFrame(data=ts_data_obs, index=ts_index)
    ts_dframe_mod = pd.DataFrame(data=ts_data_mod, index=ts_index)

    if file_min_value is not None:
        ts_dframe_complete[ts_dframe_complete < file_min_value] = np.nan
        ts_dframe_obs[ts_dframe_obs < file_min_value] = np.nan
        ts_dframe_mod[ts_dframe_mod < file_min_value] = np.nan

    if file_max_value is not None:
        ts_dframe_complete[ts_dframe_complete > file_max_value] = np.nan
        ts_dframe_obs[ts_dframe_obs > file_max_value] = np.nan
        ts_dframe_mod[ts_dframe_mod > file_max_value] = np.nan

    ts_dframe_complete = add_dframe_seasons(ts_dframe_complete, column_season='season', lut_season=file_season_lut)
    ts_dframe_obs = add_dframe_seasons(ts_dframe_obs, column_season='season', lut_season=file_season_lut)
    ts_dframe_mod = add_dframe_seasons(ts_dframe_mod, column_season='season', lut_season=file_season_lut)

    ts_dframe_complete = ts_dframe_complete.dropna(how='any')
    ts_dframe_obs = ts_dframe_obs.dropna(how='any')
    ts_dframe_mod = ts_dframe_mod.dropna(how='any')

    return ts_dframe_complete, ts_dframe_obs, ts_dframe_mod
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to write json time series
def write_time_series(file_name, file_dict_raw, file_indent=4, file_sep=','):

    file_dict_parser = {}
    for file_key, file_value in file_dict_raw.items():
        if isinstance(file_value, list):
            file_value = [str(i) for i in file_value]
            file_value = file_sep.join(file_value)
        elif isinstance(file_value, (int, float)):
            file_value = str(file_value)
        elif isinstance(file_value, str):
            pass
        else:
            print(' ===> Error in parsering json time series')
            raise RuntimeError('Parsering case not implemented yet')

        file_dict_parser[file_key] = file_value

    file_data = json.dumps(file_dict_parser, indent=file_indent, ensure_ascii=False, sort_keys=True)
    with open(file_name, "w", encoding='utf-8') as file_handle:
        file_handle.write(file_data)

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read file json
def read_file_json(file_name):

    env_ws = {}
    for env_item, env_value in os.environ.items():
        env_ws[env_item] = env_value

    with open(file_name, "r") as file_handle:
        json_block = []
        for file_row in file_handle:

            for env_key, env_value in env_ws.items():
                env_tag = '$' + env_key
                if env_tag in file_row:
                    env_value = env_value.strip("'\\'")
                    file_row = file_row.replace(env_tag, env_value)
                    file_row = file_row.replace('//', '/')

            # Add the line to our JSON block
            json_block.append(file_row)

            # Check whether we closed our JSON block
            if file_row.startswith('}'):
                # Do something with the JSON dictionary
                json_dict = json.loads(''.join(json_block))
                # Start a new block
                json_block = []

    return json_dict
# -------------------------------------------------------------------------------------
