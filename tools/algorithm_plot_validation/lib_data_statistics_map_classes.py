"""
Library Features:

Name:          lib_data_statistics_map_classes
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260514'
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

from lib_data_statistics_generic import filter_dataframe_nan

#from pandas.core.common import SettingWithCopyWarning
#warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

logging.getLogger('pandas').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute pearson classes
def compute_pearson_classes(
        df_obj,
        variable_carea='committed_area',
        variable_data='xy_r', variable_classes='classes_r',
        variable_n_obs='obs',
        variable_p_r='xy_p_r',
        variable_r='xy_r',
        variable_lon='lon',
        variable_lat='lat',
        n_sampling=10,
        p_r_thr=0.05,
        r_thr=0.6):

    # ------------------------------------------------------------------------------------------------------------------
    # generic settings
    df_obj = df_obj.copy()

    idx_global = df_obj[variable_carea] == 0
    idx_carea = df_obj[variable_carea] == 1
    idx_weak_r = df_obj[variable_r].abs() < r_thr
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # metric filters
    df_obj.loc[
        df_obj[variable_n_obs] < n_sampling,
        variable_data
    ] = np.nan

    df_obj.loc[
        df_obj[variable_p_r] > p_r_thr,
        variable_data
    ] = np.nan
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # classify Pearson R values for map
    # same classes as compute_stats_pr
    # committed area: 3.5, 2.5, 1.5, 0.5
    # global area:    7.5, 6.5, 5.5, 4.5

    # --- create class variable
    df_obj[variable_classes] = np.nan

    idx_carea = df_obj[variable_carea] == 1
    idx_global = df_obj[variable_carea] == 0
    pr = df_obj[variable_data]

    pr_min, pr_max = np.nanmin(pr), np.nanmax(pr)

    # committed area classes
    df_obj.loc[idx_carea & (pr >= -1) & (pr < 0.5), variable_classes] = 3.5
    df_obj.loc[idx_carea & (pr >= 0.5) & (pr < 0.65), variable_classes] = 2.5
    df_obj.loc[idx_carea & (pr >= 0.65) & (pr < 0.8), variable_classes] = 1.5
    df_obj.loc[idx_carea & (pr >= 0.8) & (pr <= 1), variable_classes] = 0.5

    # non-committed / global classes
    df_obj.loc[idx_global & (pr >= -1) & (pr < 0.5), variable_classes] = 7.5
    df_obj.loc[idx_global & (pr >= 0.5) & (pr < 0.65), variable_classes] = 6.5
    df_obj.loc[idx_global & (pr >= 0.65) & (pr < 0.8), variable_classes] = 5.5
    df_obj.loc[idx_global & (pr >= 0.8) & (pr <= 1), variable_classes] = 4.5
    # ------------------------------------------------------------------------------------------------------------------

    return df_obj
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# compute for snr map
def compute_snr_classes(
        df_obj,
        variable_carea='committed_area',
        variable_data='x_snr', variable_classes='classes_snr',
        variable_n_obs='obs',
        variable_p_r='xy_p_r',
        variable_r='xy_r',
        variable_lon='lon',
        variable_lat='lat',
        n_sampling=10,
        p_r_thr=0.05,
        r_thr=0.6):

    # ------------------------------------------------------------------------------------------------------------------
    # generic settings
    df_obj = df_obj.copy()

    idx_global = df_obj[variable_carea] == 0
    idx_carea = df_obj[variable_carea] == 1
    idx_weak_r = df_obj[variable_r].abs() < r_thr
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # metric filters
    df_obj.loc[
        df_obj[variable_n_obs] < n_sampling,
        variable_data
    ] = np.nan

    df_obj.loc[
        df_obj[variable_p_r] > p_r_thr,
        variable_data
    ] = np.nan
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # lower SNR handling
    df_obj.loc[
        idx_global &
        (df_obj[variable_data] < -14),
        variable_data
    ] = -14

    df_obj.loc[
        idx_carea &
        (df_obj[variable_data] < -14),
        variable_data
    ] = np.nan
    # ------------------------------------------------------------------------------------------------------------------

    '''
    # ------------------------------------------------------------------------------------------------------------------
    # geographic filters, mainly for non-committed/global pixels

    # Russia
    #df_obj.loc[
    #    (df_obj[variable_lat] > 50) &
    #    (df_obj[variable_lon] > 25) &
    #    idx_global &
    #    idx_weak_r,
    #    variable_data
    #] = np.nan

    # Sahara
    df_obj.loc[
        (df_obj[variable_lon] > -20) & (df_obj[variable_lon] < 40) &
        (df_obj[variable_lat] > 17.7) & (df_obj[variable_lat] < 30) &
        idx_global &
        idx_weak_r,
        variable_data
    ] = np.nan

    # Arabia
    df_obj.loc[
        (df_obj[variable_lon] > 40) & (df_obj[variable_lon] < 63) &
        (df_obj[variable_lat] > 11) & (df_obj[variable_lat] < 22) &
        (df_obj[variable_data] < 0),
        variable_data
    ] = np.nan

    # Pakistan
    df_obj.loc[
        (df_obj[variable_lon] > 60) & (df_obj[variable_lon] < 78) &
        (df_obj[variable_lat] > 23) & (df_obj[variable_lat] < 38) &
        (df_obj[variable_data] < 0),
        variable_data
    ] = np.nan

    # Canada
    #df_obj.loc[
    #    (df_obj[variable_lon] > -180) & (df_obj[variable_lon] < -50) &
    #    (df_obj[variable_lat] > 45) & (df_obj[variable_lat] < 90) &
    #    idx_weak_r,
    #    variable_data
    #] = np.nan

    # Norway / Sweden / Finland
    #df_obj.loc[
    #    (df_obj[variable_lon] > 3) & (df_obj[variable_lon] < 36) &
    #    (df_obj[variable_lat] > 60) & (df_obj[variable_lat] < 90) &
    #    idx_weak_r,
    #    variable_data
    #] = np.nan

    # Brazil
    df_obj.loc[
        (df_obj[variable_lon] > -72) & (df_obj[variable_lon] < -42) &
        (df_obj[variable_lat] > -60) & (df_obj[variable_lat] < 6) &
        idx_global &
        idx_weak_r,
        variable_data
    ] = np.nan

    # Australia
    df_obj.loc[
        (df_obj[variable_lon] > 112) & (df_obj[variable_lon] < 154) &
        (df_obj[variable_lat] > -44) & (df_obj[variable_lat] < -10) &
        (df_obj[variable_data] < 0),
        variable_data
    ] = np.nan

    # Central Asia
    df_obj.loc[
        (df_obj[variable_lon] > 50) & (df_obj[variable_lon] < 90) &
        (df_obj[variable_lat] > 35) & (df_obj[variable_lat] < 56) &
        (df_obj[variable_data] < 0),
        variable_data
    ] = np.nan
    # ------------------------------------------------------------------------------------------------------------------
    '''
    # ------------------------------------------------------------------------------------------------------------------
    # classify SNR values for map
    # same classes as compute_stats_snr
    # committed area: 3.5, 2.5, 1.5, 0.5
    # global area:    7.5, 6.5, 5.5, 4.5

    # --- create class variable
    df_obj[variable_classes] = np.nan

    idx_carea = df_obj['committed_area'] == 1
    idx_global = df_obj['committed_area'] == 0
    snr = df_obj[variable_data]

    # committed area classes
    df_obj.loc[idx_carea & (snr < 0), variable_classes] = 3.5
    df_obj.loc[idx_carea & (snr >= 0) & (snr < 3), variable_classes] = 2.5
    df_obj.loc[idx_carea & (snr >= 3) & (snr <= 6), variable_classes] = 1.5
    df_obj.loc[idx_carea & (snr > 6), variable_classes] = 0.5

    # non-committed / global classes
    df_obj.loc[idx_global & (snr < 0), variable_classes] = 7.5
    df_obj.loc[idx_global & (snr >= 0) & (snr < 3), variable_classes] = 6.5
    df_obj.loc[idx_global & (snr >= 3) & (snr <= 6), variable_classes] = 5.5
    df_obj.loc[idx_global & (snr > 6), variable_classes] = 4.5
    # ------------------------------------------------------------------------------------------------------------------

    '''
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pylab as plt

    plt.figure(figsize=(10, 6))

    plt.scatter(
        df_obj['lon'].values,
        df_obj['lat'].values,
        c=df_obj[variable_classes].values,
        s=1,
        cmap='tab20c',
        vmin=0,
        vmax=8,
        rasterized=True
    )

    plt.colorbar(label=variable_classes)

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(variable_classes)

    plt.show(block=True)
    '''

    return df_obj
# ----------------------------------------------------------------------------------------------------------------------
