"""
Library Features:

Name:          lib_utils_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20241031'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings

import numpy as np
import pandas as pd

from copy import deepcopy

import pytesmo.scaling as scaling
import pytesmo.metrics as metrics

from pytesmo.time_series.filters import exp_filter, boxcar_filter
from scipy import stats

from lib_info_args import logger_name

# logging
warnings.filterwarnings('ignore')
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to join time series
def join_time_series(ts_obj_collections, ts_name_ref='ref',
                     time_start=None, time_end=None, time_range=None, time_freq='D'):

    if time_range is None:
        ts_time_expected = pd.date_range(start=time_start, end=time_end, freq=time_freq)
    else:
        ts_time_expected = deepcopy(time_range)

    # get time-series reference
    ts_reference = ts_obj_collections[ts_name_ref]
    ts_obj_collections.pop(ts_name_ref)
    # get values and times reference
    ts_values_reference = ts_reference.values
    ts_time_reference = ts_reference.index.values

    # initialize dframe expected
    ts_dframe_expected = pd.DataFrame(index=ts_time_expected)

    # initialize dframe collections
    ts_obj_reference = {ts_name_ref: ts_values_reference}
    ts_dframe_ref = pd.DataFrame(data=ts_obj_reference, index=ts_time_reference)
    # join expected with reference
    ts_dframe_expected = ts_dframe_expected.join(ts_dframe_ref)

    # iterate over other time-series
    attrs_collections = {}
    for ts_name_other, ts_series_other in ts_obj_collections.items():

        # check other time-series not empty and not NoneType
        if ts_series_other is not None:

            attrs_series_other = ts_series_other.attrs

            ts_dframe_other = ts_series_other.to_frame()
            ts_dframe_other.rename(columns={'other': ts_name_other}, inplace=True)
            ts_dframe_other.attrs = attrs_series_other

            # join expected with reference
            ts_dframe_expected = ts_dframe_expected.join(ts_dframe_other)
            ts_dframe_expected[ts_name_other].attrs = attrs_series_other

            attrs_collections[ts_name_other] = attrs_series_other

    ts_dframe_expected.attrs = attrs_collections

    return ts_dframe_expected
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to apply time series filter
def apply_time_series_filter(ts_series_data, ts_filter_type='exp', ts_filter_window=3):

    # get ts series
    ts_series_data = ts_series_data.dropna()
    # get julian date of ts series
    ts_series_jd = ts_series_data.index.to_julian_date().values

    # compute filtered ts
    if ts_filter_type == 'exp':
        # exp filter fx
        data_filtered = exp_filter(ts_series_data.values, ts_series_jd, ctime=ts_filter_window)
    elif ts_filter_type == 'boxcar':
        # boxcar filter fx
        data_filtered = boxcar_filter(ts_series_data.values, ts_series_jd, window=ts_filter_window)
    else:
        # filter type not supported
        log_stream.error(' ===> Filter type "' + ts_filter_type + '" not supported')
        raise NotImplementedError('Filter type "' + ts_filter_type + '" not implemented yet')

    ts_data_filtered = pd.Series(data=data_filtered, index=ts_series_data.index)

    return ts_data_filtered
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to scale time series
def apply_time_series_scaling(ts_ref, ts_other, ts_scale_method='cdf_beta_match',
                              ts_min=0.01, ts_max=100):

    # organize time-series reference
    if (ts_min is not None) & (ts_max is not None):
        ts_ref[(ts_ref < ts_min) | (ts_ref > ts_max)] = np.nan
    ts_ref = ts_ref.dropna()

    values_ref = ts_ref.values
    index_ref = ts_ref.index.values
    n_ref_filtered = index_ref.__len__()
    data_ref = {'ref': values_ref}

    dframe_reference = pd.DataFrame(data=data_ref, index=index_ref)

    # organize time-series other
    if (ts_min is not None) & (ts_max is not None):
        ts_other[(ts_other < ts_min) | (ts_other > ts_max)] = np.nan
    ts_other = ts_other.dropna()

    values_other = ts_other.values
    index_other = ts_other.index.values
    n_other_filtered = index_other.__len__()
    data_other = {'other': values_other}
    dframe_other = pd.DataFrame(data=data_other, index=index_other)

    if n_ref_filtered > n_other_filtered:
        log_stream.warning(
            ' ===> Reference time-series (' + str(n_ref_filtered) +
            ') is greater than other time-series (' + str(n_other_filtered) + ')')

        dframe_other = dframe_other.interpolate(method='polynomial', order=4, limit=10, limit_direction='forward')
        n_other_interpolated = dframe_other.index.__len__()

    if n_other_filtered > n_ref_filtered:
        log_stream.warning(
            ' ===> Reference time-series (' + str(n_ref_filtered) +
            ') is less than other time-series (' + str(n_other_filtered) + ')')

        dframe_reference = dframe_reference.interpolate(method='polynomial', order=4, limit=10, limit_direction='forward')
        n_ref_interpolated = dframe_reference.index.__len__()

    # check other time-series
    if not dframe_other.empty:

        # join time-series common
        if dframe_reference.__len__() > dframe_other.__len__():
            dframe_common = dframe_other.join(dframe_reference)
        elif dframe_reference.__len__() < dframe_other.__len__():
            dframe_common = dframe_reference.join(dframe_other)
        else:
            dframe_common = dframe_reference.join(dframe_other)

        # get info common
        dframe_common = dframe_common.dropna()
        fields_common = list(dframe_common.columns)
        idx_common = fields_common.index('ref')

        # scale datasets using a method defined in pytesmo library
        dframe_scaled = scaling.scale(dframe_common, method=ts_scale_method, reference_index=idx_common)
        # return ts scaled
        ts_scaled = dframe_scaled['other']

    else:
        # return ts defined by NoneType
        ts_scaled = None

    return ts_scaled

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to apply time series metrics
def apply_time_series_metrics(ts_ref, ts_other):

    x_df = ts_ref.to_frame()
    if ts_other is not None:
        y_df = ts_other.to_frame()
    else:
        log_stream.error(' ===> Dataframe other is defined by NoneType')
        raise RuntimeError('Check the algorithm. Dataframe must be defined by series or empty series')

    # check dataframe ref
    x_check = deepcopy(x_df)
    x_check = x_check.dropna()
    # check dataframe other
    y_check = deepcopy(y_df)
    y_check = y_check.dropna()

    # check not empty dataframe(s)
    if not x_check.empty and not y_check.empty:

        xy_df = x_df.join(y_df)
        xy_df = xy_df.dropna()

        x_values = xy_df['ref'].values
        y_values = xy_df['other'].values

        with warnings.catch_warnings():
            pearson_r, person_p = stats.pearsonr(x_values, y_values)
            pearson_r, person_p = float('{:.2f}'.format(pearson_r)), float('{:.2e}'.format(person_p))
            spearman_rho, spearman_p = stats.spearmanr(x_values, y_values)
            spearman_rho, spearman_p = float('{:.2f}'.format(spearman_rho)), float('{:.2e}'.format(spearman_p))
            kendall_tau, kendall_p = stats.kendalltau(x_values, y_values)
            kendall_tau, kendall_p = float('{:.2f}'.format(kendall_tau)), float('{:.2e}'.format(kendall_p))
            rmsd = float('{:.2f}'.format(metrics.rmsd(x_values, y_values)))
            bias = float('{:.2f}'.format(metrics.bias(x_values, y_values)))
            nash_sutcliffe = float('{:.2f}'.format(metrics.nash_sutcliffe(x_values, y_values)))

    else:
        if x_check.empty:
            log_stream.warning(' ===> Reference datasets is defined by empty dataframe. Metrics are not computed')
        if y_check.empty:
            log_stream.warning(' ===> Other datasets is defined by empty dataframe. Metrics are not computed')

        pearson_r, person_p = np.nan, np.nan
        spearman_rho, spearman_p =  np.nan, np.nan
        kendall_tau, kendall_p = np.nan, np.nan
        rmsd =np.nan
        bias = np.nan
        nash_sutcliffe = np.nan

    # Calculate correlation coefficients, RMSD, bias, Nash Sutcliffe
    ts_metrics = {'pearson_r': pearson_r, 'pearson_p': person_p,
                  'spearman_rho': spearman_rho, 'spearman_p': spearman_p,
                  'kendall_tau': kendall_tau, 'kendall_p': kendall_p,
                  'rmsd': rmsd,
                  'bias': bias,
                  'nash_sutcliffe': nash_sutcliffe}

    return ts_metrics
# ----------------------------------------------------------------------------------------------------------------------
