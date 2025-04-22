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

from datetime import datetime

from lib_notebook_time_series import add_dframe_seasons, add_dframe_year
# ----------------------------------------------------------------------------------------------------------------------


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

    # message read time-series datasets start
    print(' ---> Read time-series datasets ... ')

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

    ts_dframe_complete = add_dframe_year(ts_dframe_complete, column_year='year')
    ts_dframe_ecmwf_l1 = add_dframe_year(ts_dframe_ecmwf_l1, column_year='year')
    ts_dframe_ecmwf_l2 = add_dframe_year(ts_dframe_ecmwf_l2, column_year='year')
    ts_dframe_ecmwf_l3 = add_dframe_year(ts_dframe_ecmwf_l3, column_year='year')
    ts_dframe_hmc = add_dframe_year(ts_dframe_hmc, column_year='year')
    ts_dframe_smap = add_dframe_year(ts_dframe_smap, column_year='year')

    # ts_dframe_complete = ts_dframe_complete.dropna(how='all')
    # ts_dframe_ecmwf_l1 = ts_dframe_ecmwf_l1.dropna(how='any')
    # ts_dframe_ecmwf_l2 = ts_dframe_ecmwf_l2.dropna(how='any')
    # ts_dframe_ecmwf_l3 = ts_dframe_ecmwf_l3.dropna(how='any')
    # ts_dframe_hmc = ts_dframe_hmc.dropna(how='any')
    # ts_dframe_smap = ts_dframe_smap.dropna(how='any')

    # message read time-series datasets end
    print(' ---> Read time-series datasets ... DONE')

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
        flag_ts_obs='sm_obs', flag_ts_mod='sm_tc', flag_time='time',
        mandatory_ts_obs=True, mandatory_ts_mod=False):

    # message read time-series product start
    print(' ---> Read time-series product ... ')

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
        time_list = time_str.split(file_sep)
    else:
        print(' ===> Flag "' + flag_time + '" not available in the source file')
        raise RuntimeError('Datasets is needed by the procedure')
    if flag_ts_obs in list(file_data.keys()):
        ts_obs_str = file_data[flag_ts_obs]
        ts_obs_list = ts_obs_str.split(file_sep)
    else:
        if mandatory_ts_obs:
            print(' ===> Flag "' + flag_ts_obs + '" not available in the source file')
            raise RuntimeError('Time-series is needed by the procedure')
        else:
            print(' ===> Flag "' + flag_ts_obs +
                  '" not available in the source file.'
                  ' Time-series is not mandatory and it will be initialized with "-9999.0" value')
            ts_obs_list = ['-9999.0'] * time_list.__len__()

    if flag_ts_mod in list(file_data.keys()):
        ts_mod_str = file_data[flag_ts_mod]
        ts_mod_list = ts_mod_str.split(file_sep)
    else:
        if mandatory_ts_mod:
            print(' ===> Flag "' + flag_ts_mod + '" not available in the source file')
            raise RuntimeError('Time-series is needed by the procedure')
        else:
            print(' ===> Flag "' + flag_ts_mod +
                  '" not available in the source file.'
                  ' Time-series is not mandatory and it will be initialized with "-9999.0" value')
            ts_mod_list = ['-9999.0'] * time_list.__len__()

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

    ts_dframe_complete = add_dframe_year(ts_dframe_complete, column_year='year')
    ts_dframe_obs = add_dframe_year(ts_dframe_obs, column_year='year')
    ts_dframe_mod = add_dframe_year(ts_dframe_mod, column_year='year')

    if mandatory_ts_obs and mandatory_ts_mod:
        ts_dframe_complete = ts_dframe_complete.dropna(how='any')

    ts_dframe_obs = ts_dframe_obs.dropna(how='any')
    if ts_dframe_obs.empty:
        print(' ===> Time series datasets "' + flag_ts_obs + '" is empty')
        if mandatory_ts_obs:
            raise RuntimeError('Error in reading time series datasets "' + flag_ts_obs + '"')

    ts_dframe_mod = ts_dframe_mod.dropna(how='any')
    if ts_dframe_mod.empty:
        print(' ===> Time series datasets "' + flag_ts_mod + '" is empty')
        if mandatory_ts_mod:
            raise RuntimeError('Error in reading time series datasets "' + flag_ts_mod + '"')

    # message read time-series product end
    print(' ---> Read time-series product ... DONE')

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

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file settings
def read_file_settings(file_name):
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

# ----------------------------------------------------------------------------------------------------------------------
