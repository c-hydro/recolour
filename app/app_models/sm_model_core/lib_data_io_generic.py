"""
Library Features:

Name:          lib_data_io_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260615'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd

from copy import deepcopy

from lib_info_args import logger_name
from lib_utils_time import define_time_frequency

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to combine data over the expected time range
def combine_data_point_by_time(dframe_k1, dframe_k2, dframe_k3,
                               time_tag='time', time_frequency='h', time_reverse=True,
                               time_ref=None, fill_value_missing=-9997,
                               no_data_k1=-9999, no_data_k2=-9999, no_data_k3=-9999,
                               scale_factor_k1=1, scale_factor_k2=1, scale_factor_k3=0.01):

    time_frequency = time_frequency.lower()

    if (dframe_k1 is None) or (dframe_k2 is None) or (dframe_k3 is None):
        log_stream.warning(' ===> Dataframes are not defined; One or more dataframes are defined by None')
        return None

    time_start_k1, time_end_k1 = dframe_k1.index.min(), dframe_k1.index.max()
    time_start_k2, time_end_k2 = dframe_k2.index.min(), dframe_k2.index.max()
    time_start_k3, time_end_k3 = dframe_k3.index.min(), dframe_k3.index.max()

    time_start_common = pd.DatetimeIndex([time_start_k1, time_start_k2, time_start_k3]).min()
    time_end_common = pd.DatetimeIndex([time_end_k1, time_end_k2, time_end_k3]).max()

    if time_ref is not None:
        time_ref = pd.Timestamp(time_ref)
        time_end_common = max(time_end_common, time_ref)

    time_range_common = pd.date_range(time_start_common, time_end_common, freq=time_frequency)

    attrs_common = dframe_k1.attrs

    dframe_k1 = dframe_k1.drop(columns=[time_tag], errors='ignore')
    if 'values_k1' not in dframe_k1.columns:
        log_stream.error(' ===> Dataframe 1 does not have the column "values_k1"')
        raise RuntimeError('Column "values_k1" must be included in the dataframe.')

    dframe_k1.loc[dframe_k1['values_k1'] == no_data_k1, 'values_k1'] = np.nan
    dframe_k1['values_k1'] = dframe_k1['values_k1'].values * scale_factor_k1
    dframe_k1.loc[np.isnan(dframe_k1['values_k1']), 'values_k1'] = no_data_k1

    dframe_k2 = dframe_k2.drop(columns=[time_tag], errors='ignore')
    if 'values_k2' not in dframe_k2.columns:
        log_stream.error(' ===> Dataframe 2 does not have the column "values_k2"')
        raise RuntimeError('Column "values_k2" must be included in the dataframe.')

    dframe_k2.loc[dframe_k2['values_k2'] == no_data_k2, 'values_k2'] = np.nan
    dframe_k2['values_k2'] = dframe_k2['values_k2'].values * scale_factor_k2
    dframe_k2.loc[np.isnan(dframe_k2['values_k2']), 'values_k2'] = no_data_k2

    dframe_k3 = dframe_k3.drop(columns=[time_tag], errors='ignore')
    if 'values_k3' not in dframe_k3.columns:
        log_stream.error(' ===> Dataframe 3 does not have the column "values_k3"')
        raise RuntimeError('Column "values_k3" must be included in the dataframe.')

    dframe_k3.loc[dframe_k3['values_k3'] == no_data_k3, 'values_k3'] = np.nan
    dframe_k3['values_k3'] = dframe_k3['values_k3'].values * scale_factor_k3
    dframe_k3.loc[np.isnan(dframe_k3['values_k3']), 'values_k3'] = no_data_k3

    dframe_common = pd.DataFrame(index=time_range_common)
    dframe_common.index.name = time_tag

    dframe_common = dframe_common.join(dframe_k1)
    dframe_common = dframe_common.join(dframe_k2)
    dframe_common = dframe_common.join(dframe_k3)

    if time_ref is not None:
        dframe_common = dframe_common.loc[dframe_common.index <= time_ref]

        time_range_expected = pd.date_range(
            start=dframe_common.index.min(),
            end=time_ref,
            freq=time_frequency
        )

        dframe_common = dframe_common.reindex(time_range_expected)
        dframe_common.index.name = time_tag

        value_cols = ['values_k1', 'values_k2', 'values_k3']

        missing_rows = dframe_common[value_cols].isna().all(axis=1)
        dframe_common.loc[missing_rows, value_cols] = fill_value_missing

    if time_tag not in dframe_common.columns:
        dframe_common[time_tag] = dframe_common.index

    if time_reverse:
        dframe_common = dframe_common.sort_index(ascending=False)

    dframe_common.attrs = attrs_common

    return dframe_common
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to combine data over the expected time range
def combine_data_point_by_time_OLD(dframe_k1, dframe_k2, dframe_k3,
                               time_tag='time', time_frequency='h', time_reverse=True,
                               no_data_k1=-9999, no_data_k2=-9999, no_data_k3=-9999,
                               scale_factor_k1=1, scale_factor_k2=1, scale_factor_k3=0.01):

    time_frequency = time_frequency.lower()

    # check dataframes
    if (dframe_k1 is None) or (dframe_k2 is None) or (dframe_k3 is None):
        log_stream.warning(' ===> Dataframes are not defined; One or more dataframes are defined by None')
        return None

    # define time common limits
    time_start_k1, time_end_k1 = dframe_k1.index.min(), dframe_k1.index.max()
    time_start_k2, time_end_k2 = dframe_k2.index.min(), dframe_k2.index.max()
    time_start_k3, time_end_k3 = dframe_k3.index.min(), dframe_k3.index.max()

    time_start_common = pd.DatetimeIndex([time_start_k1, time_start_k2, time_start_k3]).min()
    time_end_common = pd.DatetimeIndex([time_end_k1, time_end_k2, time_end_k3]).max()

    time_range_common = pd.date_range(time_start_common, time_end_common, freq=time_frequency)

    # get attributes
    attrs_common = dframe_k1.attrs

    # remove time information from each dataframe and apply scale factor
    dframe_k1 = dframe_k1.drop(columns=[time_tag])

    if 'values_k1' not in dframe_k1.columns:
        log_stream.error(' ===> Dataframe 1 does not have the column "values_k1"')
        raise RuntimeError('Column "values_k1" must be included in the dataframe. '
                           'Check if the variable mapping is correctly defined')

    dframe_k1['values_k1'][dframe_k1['values_k1'] == no_data_k1] = np.nan
    dframe_k1['values_k1'] = dframe_k1['values_k1'].values * scale_factor_k1
    dframe_k1['values_k1'][np.isnan(dframe_k1['values_k1'])] = no_data_k1

    dframe_k2 = dframe_k2.drop(columns=[time_tag])

    if 'values_k2' not in dframe_k2.columns:
        log_stream.error(' ===> Dataframe 2 does not have the column "values_k2"')
        raise RuntimeError('Column "values_k2" must be included in the dataframe. '
                           'Check if the variable mapping is correctly defined')

    dframe_k2['values_k2'][dframe_k2['values_k2'] == no_data_k2] = np.nan
    dframe_k2['values_k2'] = dframe_k2['values_k2'].values * scale_factor_k2
    dframe_k2['values_k2'][np.isnan(dframe_k2['values_k2'])] = no_data_k2

    dframe_k3 = dframe_k3.drop(columns=[time_tag])

    if 'values_k3' not in dframe_k3.columns:
        log_stream.error(' ===> Dataframe 3 does not have the column "values_k3"')
        raise RuntimeError('Column "values_k3" must be included in the dataframe. '
                           'Check if the variable mapping is correctly defined')

    dframe_k3['values_k3'][dframe_k3['values_k3'] == no_data_k3] = np.nan
    dframe_k3['values_k3'] = dframe_k3['values_k3'].values * scale_factor_k3
    dframe_k3['values_k3'][np.isnan(dframe_k3['values_k3'])] = no_data_k3

    # define common dataframe
    dframe_common = pd.DataFrame(index=time_range_common)
    dframe_common.index.name = time_tag
    dframe_common = dframe_common.join(dframe_k1)
    dframe_common = dframe_common.join(dframe_k2)
    dframe_common = dframe_common.join(dframe_k3)
    # add column time to the dataframe
    if time_tag not in dframe_common.columns:
        dframe_common[time_tag] = dframe_common.index

    # time reverse flag
    if time_reverse:
        dframe_common = dframe_common.sort_index(ascending=False)

    dframe_common.attrs = attrs_common

    return dframe_common

# ----------------------------------------------------------------------------------------------------------------------
