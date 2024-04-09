"""
Library Features:

Name:          lib_data_statistics
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import os
import numpy as np
import pandas as pd

from copy import deepcopy

from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

logging.getLogger('pandas').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter dataframe nan(s) for stats
def filter_dataframe_nan(
        df_obj: pd.DataFrame,
        var_carea_name: str = 'committed_area', var_carea_value: int = 1,
        var_data_name: str = 'xy_pr'):

    sub_obj = ((df_obj[var_carea_name] == var_carea_value) & (np.isnan(df_obj[var_data_name])))
    value_filtered_nan = sub_obj.sum()

    return value_filtered_nan
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter dataframe data for stats
def filter_dataframe_data(
        df_obj: pd.DataFrame,
        var_carea_name: str = 'committed_area', var_carea_value: int = 1,
        var_data_name: str = 'xy_pr', var_stats_name: str = 'stats_pr',
        bnd_min_value: float or None = -1, bnd_max_value: float or None = 0.5, bnd_int_value: float = 3.5,
        bnd_min_extend: bool = True, bnd_max_extend: bool = False):

    if (bnd_min_value is not None) and (bnd_max_value is not None):
        if bnd_min_extend and (not bnd_max_extend):
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & \
                      (df_obj[var_data_name] >= bnd_min_value) & (df_obj[var_data_name] < bnd_max_value)
        elif (not bnd_min_extend) and (not bnd_max_extend):
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & \
                      (df_obj[var_data_name] > bnd_min_value) & (df_obj[var_data_name] < bnd_max_value)
        elif bnd_min_extend and bnd_max_extend:
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & \
                      (df_obj[var_data_name] >= bnd_min_value) & (df_obj[var_data_name] <= bnd_max_value)
        elif (not bnd_min_extend) and bnd_max_extend:
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & \
                      (df_obj[var_data_name] > bnd_min_value) & (df_obj[var_data_name] <= bnd_max_value)
        else:
            logging.error(' ===> DataFrame boundary extent (min and max case) are not expected')
            raise RuntimeError('Check the DataFrame extent to avoid errors like this')

    elif (bnd_min_value is None) and (bnd_max_value is not None):
        if bnd_max_extend:
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & (df_obj[var_data_name] <= bnd_max_value)
        elif not bnd_max_extend:
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & (df_obj[var_data_name] < bnd_max_value)
        else:
            logging.error(' ===> DataFrame boundary extent (only max case) are not expected')
            raise RuntimeError('Check the DataFrame extent to avoid errors like this')

    elif (bnd_min_value is not None) and (bnd_max_value is None):
        if bnd_min_extend:
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & (df_obj[var_data_name] >= bnd_max_value)
        elif not bnd_min_extend:
            sub_obj = (df_obj[var_carea_name] == var_carea_value) & (df_obj[var_data_name] > bnd_max_value)
        else:
            logging.error(' ===> DataFrame boundary conditions (only min case) are not expected')
            raise RuntimeError('Check the DataFrame conditions to avoid errors like this')

    elif (bnd_min_value is None) and (bnd_max_value is None):
        sub_obj = deepcopy(df_obj)
    else:
        logging.error(' ===> DataFrame boundary conditions are not expected. Both are ')
        raise RuntimeError('Check the DataFrame conditions to avoid errors like this')

    value_filtered_n = sub_obj.sum()
    df_obj.loc[sub_obj, var_stats_name] = bnd_int_value

    return df_obj, value_filtered_n
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute snr statistics
def compute_stats_snr(df_obj, variable_carea='committed_area',
                      variable_data='xyz_x_snr', variable_stats='stats_snr'):

    # initialized variable(s)
    num_comm, num_global = [], []
    df_obj[variable_stats] = np.nan

    # committed area case
    # apply limits
    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = None, 0, 3.5, False, False, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0, 3, 2.5, True, False, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 3, 6, 1.5, True, True, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 6, None, 0.5, False, False, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    carea_value = 1
    df_filtered_nan = filter_dataframe_nan(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data)
    num_comm.append(df_filtered_nan)

    # committed area case
    # apply limits
    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = None, 0, 7.5, False, False, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0, 3, 6.5, True, False, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 3, 6, 5.5, True, True, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 6, None, 4.5, False, False, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    carea_value = 0
    df_filtered_nan = filter_dataframe_nan(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data)
    num_global.append(df_filtered_nan)

    # compute percentage(s) of committed and global case(s)
    perc_comm = num_comm / (df_obj[variable_carea] == 1).sum() * 100.

    # ----------------------------------------------
    # debug
    # if no committed area mask is present, perc_global will be zero in all its entries
    # leading to errors
    # to account for that, for grids which have no committed area or land masks,
    # the global percentages are taken equal to committed percentages

    perc_global = num_global / (df_obj[variable_carea] == 0).sum() * 100.

    if np.any(np.isnan(perc_comm)):
        perc_comm = perc_global
        logging.warning(' ===> Committed case is taken equal to committed-area case.')
    elif np.any(np.isnan(perc_global)):
        perc_global = perc_comm
        logging.warning(' ===> Global case is taken equal to committed-area case.')
    # ----------------------------------------------

    '''
    if variable_datasets == 'ascat':
        num_non_comm_perc = num_non_comm / (df['committed_area'] == 0).sum() * 100.
    elif variable_datasets == 'ecmwf':
        num_non_comm_perc = num_non_comm
    else:
        raise IOError(' datasets name is not supported ("ascat" or "ecmwf")')
    '''
    sum_comm, sum_global = sum(perc_comm), sum(perc_global)
    # check percentage
    if sum_comm > 100:
        sum_str = "{:.1f}".format(sum_comm)
        logging.warning(' ===> Committed case ' + sum_str + ' [%] exceeded the 100 [%]. Check your stats.')
    if sum_global > 100:
        sum_str = "{:.1f}".format(sum_global)
        logging.warning(' ===> Global case ' + sum_str + ' [%] exceeded the 100 [%]. Check your stats.')

    return df_obj, perc_comm, perc_global

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute pearson statistics
def compute_stats_pearson(df_obj, variable_carea='committed_area',
                          variable_data='xy_pr', variable_stats='stats_pr'):

    # initialized variable(s)
    num_comm, num_global = [], []
    df_obj[variable_stats] = np.nan

    # committed area case
    # apply limits
    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = -1, 0.5, 3.5, True, False, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0.5, 0.65, 2.5, True, False, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0.65, 0.8, 1.5, True, False, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0.8, 1, 0.5, True, True, 1
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_comm.append(df_filtered_n)

    carea_value = 1
    df_filtered_nan = filter_dataframe_nan(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data)
    num_comm.append(df_filtered_nan)

    # global case
    # apply limits
    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = -1, 0.5, 7.5, True, False, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0.5, 0.65, 6.5, True, False, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0.65, 0.8, 5.5, True, False, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    bnd_min, bnd_max, bnd_value, extend_min, extend_max, carea_value = 0.8, 1, 4.5, True, True, 0
    df_obj, df_filtered_n = filter_dataframe_data(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data, var_stats_name=variable_stats,
        bnd_min_value=bnd_min, bnd_max_value=bnd_max, bnd_int_value=bnd_value,
        bnd_min_extend=extend_min, bnd_max_extend=extend_max)
    num_global.append(df_filtered_n)

    carea_value = 0
    df_filtered_nan = filter_dataframe_nan(
        df_obj,
        var_carea_name=variable_carea, var_carea_value=carea_value,
        var_data_name=variable_data)
    num_global.append(df_filtered_nan)

    # compute percentage(s) of committed and global case(s)
    perc_comm = num_comm / (df_obj[variable_carea] == 1).sum() * 100.

    # ----------------------------------------------
    # debug
    # if no committed area mask is present, perc_global will be zero in all its entries
    # leading to errors
    # to account for that, for grids which have no committed area or land masks,
    # the global percentages are taken equal to committed percentages

    perc_global = num_global / (df_obj[variable_carea] == 0).sum() * 100.

    if np.any(np.isnan(perc_comm)):
        perc_comm = perc_global
        logging.warning(' ===> Committed case is taken equal to committed-area case.')
    elif np.any(np.isnan(perc_global)):
        perc_global = perc_comm
        logging.warning(' ===> Global case is taken equal to committed-area case.')
    # ----------------------------------------------

    '''
    if variable_datasets == 'ascat':
        num_non_comm_perc = num_non_comm / (df['committed_area'] == 0).sum() * 100.
    elif variable_datasets == 'ecmwf':
        num_non_comm_perc = num_non_comm
    else:
        raise IOError(' datasets name is not supported ("ascat" or "ecmwf")')
    '''

    sum_comm, sum_global = sum(perc_comm), sum(perc_global)
    # check percentage
    if sum_comm > 100:
        sum_str = "{:.1f}".format(sum_comm)
        logging.warning(' ===> Committed case ' + sum_str + ' [%] exceeded the 100 [%]. Check your stats.')
    if sum_global > 100:
        sum_str = "{:.1f}".format(sum_global)
        logging.warning(' ===> Global case ' + sum_str + ' [%] exceeded the 100 [%]. Check your stats.')

    return df_obj, perc_comm, perc_global
# ----------------------------------------------------------------------------------------------------------------------
