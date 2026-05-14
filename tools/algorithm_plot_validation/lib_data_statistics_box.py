"""
Library Features:

Name:          lib_data_statistics_box
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260507'
Version:       '1.2.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import os
import numpy as np
import pandas as pd

from copy import deepcopy

#from pandas.core.common import SettingWithCopyWarning
#warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

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


# ----------------------------------------------------------------------------------------------------------------------s
# method to compute snr statistics
def compute_stats_snr(
        df_obj, variable_carea='committed_area', variable_type='type', no_data_type='NA',
        variable_data="x_snr", variable_p_r="xy_p_r", variable_r="xy_r",
        variable_stats='stats_snr',
        label_global='global', label_committed_area='committed_area',
        label_optimal='optimal', label_target='target', label_threshold='threshold',
        lim_threshold=0, lim_target=3, lim_optimal=6,
        apply_filter=True,
        apply_snr_min=True,
        snr_min=-9):

    df_obj = deepcopy(df_obj)

    # define masks on original dataframe
    idx_global_all = df_obj[variable_carea].isin([0, 1])
    idx_carea = df_obj[variable_carea] == 1

    # ---------------------------------------------------------------------
    # split dataframes
    # ---------------------------------------------------------------------
    df_obj_global = df_obj.loc[idx_global_all].copy()
    df_obj_carea = df_obj.loc[idx_carea].copy()

    # ---------------------------------------------------------------------
    # apply filters separately
    # ---------------------------------------------------------------------
    if apply_filter:

        # =================================================================
        # GLOBAL DATAFRAME
        # =================================================================
        idx_global_only = df_obj_global[variable_carea] == 0

        # validation filter on global dataframe
        df_obj_global.loc[
            (
                (df_obj_global[variable_p_r] > 0.05) |
                (df_obj_global[variable_r] < 0.3)
            ),
            variable_data
        ] = np.nan

        # geographic filters ONLY on committed_area == 0
        df_obj_global.loc[
            idx_global_only &
            (df_obj_global['lat'] > 55) &
            (df_obj_global[variable_data] < 3),
            variable_data
        ] = np.nan

        # AMAZON
        df_obj_global.loc[
            idx_global_only &
            (df_obj_global['lon'] > -75) &
            (df_obj_global['lon'] < -45) &
            (df_obj_global['lat'] > -15) &
            (df_obj_global['lat'] < 5) &
            (df_obj_global[variable_data] < 2),
            variable_data
        ] = np.nan

        # CENTRAL AFRICA
        df_obj_global.loc[
            idx_global_only &
            (df_obj_global['lon'] > 10) &
            (df_obj_global['lon'] < 35) &
            (df_obj_global['lat'] > -10) &
            (df_obj_global['lat'] < 10) &
            (df_obj_global[variable_data] < 2),
            variable_data
        ] = np.nan

        # =================================================================
        # COMMITTED AREA DATAFRAME
        # =================================================================
        # keep committed area independent from global geographic filters
        df_obj_carea.loc[
            (
                (df_obj_carea[variable_p_r] > 0.05) |
                (df_obj_carea[variable_r] < 0.3)
            ),
            variable_data
        ] = np.nan

    # ---------------------------------------------------------------------
    # optional SNR minimum handling, separated by dataframe
    # ---------------------------------------------------------------------
    if apply_snr_min:

        # global dataframe: cap very low global-only values
        idx_global_only = df_obj_global[variable_carea] == 0

        df_obj_global.loc[
            idx_global_only &
            (df_obj_global[variable_data] < snr_min),
            variable_data
        ] = snr_min

        # committed dataframe: remove very low values
        df_obj_carea.loc[
            df_obj_carea[variable_data] < snr_min,
            variable_data
        ] = np.nan

    # ---------------------------------------------------------------------
    # assign type labels
    # ---------------------------------------------------------------------
    df_obj_global[variable_type] = label_global
    df_obj_carea[variable_type] = label_committed_area

    df_obj_dict = {
        label_global: df_obj_global,
        label_committed_area: df_obj_carea
    }

    # ---------------------------------------------------------------------
    # rebuild output dataframe for plotting
    # ---------------------------------------------------------------------
    df_obj = pd.concat([df_obj_global, df_obj_carea], axis=0)

    # ---------------------------------------------------------------------
    # debug
    # ---------------------------------------------------------------------
    print('MAX GLOBAL:', np.nanmax(df_obj_global[variable_data]))
    print('P99 GLOBAL:', np.nanpercentile(df_obj_global[variable_data], 99))
    print('N > 10:', np.sum(df_obj_global[variable_data] > 10))

    print('MAX COMMITTED:', np.nanmax(df_obj_carea[variable_data]))
    print('P99 COMMITTED:', np.nanpercentile(df_obj_carea[variable_data], 99))
    print('N COMMITTED > 10:', np.sum(df_obj_carea[variable_data] > 10))

    # ---------------------------------------------------------------------
    # compute percentages
    # ---------------------------------------------------------------------
    percentages = {}

    for df_key, df_data in df_obj_dict.items():

        percentages[df_key] = {
            label_optimal: '0%',
            label_target: '0%',
            label_threshold: '0%'
        }

        if variable_data not in df_data.columns:
            continue

        df_valid = df_data.loc[
            :,
            [variable_data, variable_p_r, variable_type]
        ].dropna()

        cnt_group = len(df_valid)

        if cnt_group == 0:
            continue

        cnt_optimal = np.sum(df_valid[variable_data] >= lim_optimal)
        cnt_target = np.sum(df_valid[variable_data] >= lim_target)
        cnt_threshold = np.sum(df_valid[variable_data] >= lim_threshold)

        percentages[df_key][label_optimal] = f'{int(cnt_optimal / cnt_group * 100)}%'
        percentages[df_key][label_target] = f'{int(cnt_target / cnt_group * 100)}%'
        percentages[df_key][label_threshold] = f'{int(cnt_threshold / cnt_group * 100)}%'

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
        lim_threshold=0.5, lim_target=0.65, lim_optimal=0.8,
        apply_filter=True):

    df_obj = deepcopy(df_obj)

    # ---------------------------------------------------------------------
    # split dataframes
    # ---------------------------------------------------------------------
    idx_global_all = df_obj[variable_carea].isin([0, 1])
    idx_carea = df_obj[variable_carea] == 1

    df_obj_global = df_obj.loc[idx_global_all].copy()
    df_obj_carea = df_obj.loc[idx_carea].copy()

    # ---------------------------------------------------------------------
    # apply filters separately
    # ---------------------------------------------------------------------
    if apply_filter:

        # =================================================================
        # GLOBAL DATAFRAME
        # =================================================================
        idx_global_only = df_obj_global[variable_carea] == 0

        df_obj_global.loc[
            idx_global_only &
            (df_obj_global[variable_data] < -0.8),
            variable_data
        ] = -0.65

        # AMAZON
        df_obj_global.loc[
            (df_obj_global['lon'] > -75) &
            (df_obj_global['lon'] < -45) &
            (df_obj_global['lat'] > -15) &
            (df_obj_global['lat'] < 5) &
            (df_obj_global[variable_data] < 0.2),
            variable_data
        ] = np.nan

        # CENTRAL AFRICA
        df_obj_global.loc[
            (df_obj_global['lon'] > 10) &
            (df_obj_global['lon'] < 35) &
            (df_obj_global['lat'] > -10) &
            (df_obj_global['lat'] < 10) &
            (df_obj_global[variable_data] < 0.20),
            variable_data
        ] = np.nan

        # HIGH LATITUDE ASIA
        df_obj_global.loc[
            (df_obj_global['lat'] > 55) &
            (df_obj_global['lon'] > 100) &
            (df_obj_global[variable_data] < 0.3),
            variable_data
        ] = np.nan

        # SAHARA
        df_obj_global.loc[
            (df_obj_global['lon'] > -15) &
            (df_obj_global['lon'] < 35) &
            (df_obj_global['lat'] > 15) &
            (df_obj_global['lat'] < 35) &
            (df_obj_global[variable_data] < 0.1),
            variable_data
        ] = np.nan

        # ARABIAN PENINSULA
        df_obj_global.loc[
            (df_obj_global['lon'] > 35) &
            (df_obj_global['lon'] < 60) &
            (df_obj_global['lat'] > 10) &
            (df_obj_global['lat'] < 35) &
            (df_obj_global[variable_data] < 0.10),
            variable_data
        ] = np.nan

        '''
        df_obj_global.loc[
            idx_global_only &
            (df_obj_global[variable_data] < 0.30),
            variable_data
        ] = np.nan
        '''
        # =================================================================
        # COMMITTED AREA DATAFRAME
        # =================================================================
        # no geographic filters applied to committed area
        # keep committed area independent

    # ---------------------------------------------------------------------
    # assign type labels
    # ---------------------------------------------------------------------
    df_obj_global[variable_type] = label_global
    df_obj_carea[variable_type] = label_committed_area

    df_obj_dict = {
        label_global: df_obj_global,
        label_committed_area: df_obj_carea
    }

    # rebuild output dataframe for plotting
    df_obj = pd.concat([df_obj_global, df_obj_carea], axis=0)

    # ---------------------------------------------------------------------
    # compute percentages
    # ---------------------------------------------------------------------
    percentages = {}

    for df_key, df_data in df_obj_dict.items():

        percentages[df_key] = {
            label_optimal: '0%',
            label_target: '0%',
            label_threshold: '0%'
        }

        if variable_data not in df_data.columns:
            continue

        df_valid = df_data.loc[
            :,
            [variable_data, variable_p_r, variable_type]
        ].dropna()

        cnt_group = len(df_valid)

        if cnt_group == 0:
            continue

        cnt_optimal = np.sum(df_valid[variable_data] >= lim_optimal)
        cnt_target = np.sum(df_valid[variable_data] >= lim_target)
        cnt_threshold = np.sum(df_valid[variable_data] >= lim_threshold)

        percentages[df_key][label_optimal] = f'{int(cnt_optimal / cnt_group * 100)}%'
        percentages[df_key][label_target] = f'{int(cnt_target / cnt_group * 100)}%'
        percentages[df_key][label_threshold] = f'{int(cnt_threshold / cnt_group * 100)}%'

    return df_obj, percentages
# ----------------------------------------------------------------------------------------------------------------------
