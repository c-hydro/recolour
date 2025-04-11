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
# ----------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to add seasonal description
def add_dframe_seasons(dframe_data, column_season='season', lut_season=None):

    dframe_time = dframe_data.index

    grp_season = [lut_season.get(pd.Timestamp(t_stamp).month) for t_stamp in dframe_time]
    dframe_data[column_season] = grp_season

    return dframe_data


# -------------------------------------------------------------------------------------






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
