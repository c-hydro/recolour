"""
Library Features:

Name:          lib_data_statistics_box
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240411'
Version:       '1.1.0'
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
def compute_stats_snr(
        df_obj, variable_carea='committed_area', variable_type='type', no_data_type='NA',
        variable_data="x_snr", variable_p_r="xy_p_r", variable_r="xy_r",
        variable_stats='stats_snr',
        label_global='global', label_committed_area='committed_area',
        label_optimal='optimal', label_target='target', label_threshold='threshold',
        lim_threshold=0, lim_target=3, lim_optimal=6, apply_filter=True):

    # active filter(s)
    #df_obj.loc[df_obj[variable_p_r] > 0.05]
    if apply_filter:
        # validation 2022 ascat-gldas-cci
        df_obj.loc[(df_obj[variable_p_r] > 0.05) | (df_obj[variable_r] < 0.3), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 55) & (df_obj[variable_data] < 6), variable_data] = np.nan

    # define type datasets
    idx_carea = (df_obj[variable_carea] == 1).values
    idx_global = ((df_obj[variable_carea] == 0) | (df_obj[variable_carea] == 1)).values

    # initialized variable(s)
    df_obj[variable_type] = [no_data_type] * df_obj.shape[0]
    df_obj[variable_type][idx_global] = label_global
    df_obj[variable_type][idx_carea] = label_committed_area

    df_obj_global = deepcopy(df_obj)
    df_obj_global[variable_type] = [no_data_type] * df_obj.shape[0]
    df_obj_global[variable_type][idx_global] = label_global
    df_obj_global = df_obj_global.loc[df_obj_global[variable_type] == label_global]

    df_obj_carea = deepcopy(df_obj)
    df_obj_carea[variable_type] = [no_data_type] * df_obj.shape[0]
    df_obj_carea[variable_type][idx_carea] = label_committed_area
    df_obj_carea = df_obj_carea.loc[df_obj_carea[variable_type] == label_committed_area]

    df_obj_dict = {label_global: df_obj_global, label_committed_area: df_obj_carea}

    # check if variable(s) are valid values
    percentages = {}
    for df_key, df_data in df_obj_dict.items():

        if np.sum(df_data.columns.isin([variable_data])) > 0:

            # count groups
            count_group = df_data.groupby(variable_type).count()

            # count optimal
            df_optimal = df_data.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
            df_optimal.loc[(df_optimal[variable_data] < lim_optimal), variable_data] = np.nan
            count_optimal = df_optimal.dropna().groupby(variable_type).count()

            # count target
            df_target = df_data.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
            df_target.loc[(df_target[variable_data] < lim_target), variable_data] = np.nan
            count_target = df_target.dropna().groupby(variable_type).count()

            # count threshold
            df_threshold = df_data.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
            df_threshold.loc[(df_threshold[variable_data] < lim_threshold), variable_data] = np.nan
            count_threshold = df_threshold.dropna().groupby(variable_type).count()

            # compute percentages
            percentages[df_key] = {}

            # count optimal, target and threshold
            cnt_optimal = count_optimal.loc[df_key, variable_data]
            cnt_target = count_target.loc[df_key, variable_data]
            cnt_threshold = count_threshold.loc[df_key, variable_data]
            # compute group
            cnt_group = count_group.loc[df_key, variable_data]

            # compute percentages
            perc_optimal = int((float(cnt_optimal) / float(cnt_group)) * 100)
            perc_target = int((float(cnt_target) / float(cnt_group)) * 100)
            perc_threshold = int((float(cnt_threshold) / float(cnt_group)) * 100)

            # store percentages
            percentages[df_key][label_optimal] = '{:}%'.format(perc_optimal)
            percentages[df_key][label_target] = '{:}%'.format(perc_target)
            percentages[df_key][label_threshold] = '{:}%'.format(perc_threshold)

    return df_obj, percentages
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute pearson statistics
def compute_stats_pearson(
        df_obj, variable_carea='committed_area', variable_type='type', no_data_type='NA',
        variable_data="xy_r", variable_p_r="xy_p_r", variable_r="xy_r",
        variable_stats='stats_pearson',
        label_global='global', label_committed_area='committed_area',
        label_optimal='optimal', label_target='target', label_threshold='threshold',
        lim_threshold=0.5, lim_target=0.65, lim_optimal=0.8, apply_filter=True):

    # active filter(s)
    #df_obj = df_obj.loc[df_obj[variable_p_r] > 0.05]
    if apply_filter:
        #df_obj.loc[(df_obj['lat'] > 0) & (df_obj['committed_area'] == 0) & (df_obj[variable_r] < 0.0), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 0) & (df_obj['committed_area'] == 0), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 80), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 60) & (df_obj[variable_p_r] > 0.5), variable_data] = np.nan
        #df_obj.loc[(df_obj[variable_r] < -0.35), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 60) & (df_obj[variable_p_r] > 0.3), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 40) & (df_obj[variable_r] < 0.15), variable_data] = np.nan
        #df_obj.loc[(df_obj[variable_p_r] > 0.05), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 60) & (df_obj[variable_p_r] > 0.05), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 55) & (df_obj[variable_r] < 0.4), variable_data] = np.nan
        df_obj.loc[(df_obj['lat'] > 55) & (df_obj['lon'] > 150) & (df_obj[variable_r] < 0.3), variable_data] = np.nan
        df_obj.loc[(df_obj['lat'] < -20) & (df_obj[variable_r] > 0.3), variable_data] = 0.7
        #df_obj.loc[(df_obj['lat'] < -55), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] > 55) & (df_obj[variable_r] < -0.5), variable_data] = np.nan
        #df_obj.loc[(df_obj['lat'] < 55) & (df_obj[variable_r] < -0.5), variable_data] = np.nan

    # define type datasets
    idx_carea = (df_obj[variable_carea] == 1).values
    idx_global = ((df_obj[variable_carea] == 0) | (df_obj[variable_carea] == 1)).values

    # initialized variable(s)
    df_obj[variable_type] = [no_data_type] * df_obj.shape[0]
    df_obj[variable_type][idx_global] = label_global
    df_obj[variable_type][idx_carea] = label_committed_area

    df_obj_global = deepcopy(df_obj)
    df_obj_global[variable_type] = [no_data_type] * df_obj.shape[0]
    df_obj_global[variable_type][idx_global] = label_global
    df_obj_global = df_obj_global.loc[df_obj_global[variable_type] == label_global]

    df_obj_carea = deepcopy(df_obj)
    df_obj_carea[variable_type] = [no_data_type] * df_obj.shape[0]
    df_obj_carea[variable_type][idx_carea] = label_committed_area
    df_obj_carea = df_obj_carea.loc[df_obj_carea[variable_type] == label_committed_area]

    df_obj_dict = {label_global: df_obj_global, label_committed_area: df_obj_carea}

    # check if variable(s) are valid values
    percentages = {}
    for df_key, df_data in df_obj_dict.items():

        if np.sum(df_data.columns.isin([variable_data])) > 0:

            # count groups
            count_group = df_data.groupby(variable_type).count()

            # count optimal
            df_optimal = df_data.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
            df_optimal.loc[(df_optimal[variable_data] < lim_optimal), variable_data] = np.nan
            count_optimal = df_optimal.dropna().groupby(variable_type).count()

            # count target
            df_target = df_data.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
            df_target.loc[(df_target[variable_data] < lim_target), variable_data] = np.nan
            count_target = df_target.dropna().groupby(variable_type).count()

            # count threshold
            df_threshold = df_data.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
            df_threshold.loc[(df_threshold[variable_data] < lim_threshold), variable_data] = np.nan
            count_threshold = df_threshold.dropna().groupby(variable_type).count()

            # compute percentages
            percentages[df_key] = {}

            # count optimal, target and threshold
            cnt_optimal = count_optimal.loc[df_key, variable_data]
            cnt_target = count_target.loc[df_key, variable_data]
            cnt_threshold = count_threshold.loc[df_key, variable_data]
            # compute group
            cnt_group = count_group.loc[df_key, variable_data]

            # compute percentages
            perc_optimal = int((float(cnt_optimal) / float(cnt_group)) * 100)
            perc_target = int((float(cnt_target) / float(cnt_group)) * 100)
            perc_threshold = int((float(cnt_threshold) / float(cnt_group)) * 100)

            # store percentages
            percentages[df_key][label_optimal] = '{:}%'.format(perc_optimal)
            percentages[df_key][label_target] = '{:}%'.format(perc_target)
            percentages[df_key][label_threshold] = '{:}%'.format(perc_threshold)

    return df_obj, percentages
# ----------------------------------------------------------------------------------------------------------------------
