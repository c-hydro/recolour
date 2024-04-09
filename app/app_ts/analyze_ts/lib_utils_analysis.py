"""
Library Features:

Name:          lib_utils_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240109'
Version:       '1.0.0'
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
def apply_time_series_scaling(ts_ref, ts_other, ts_scale_method='cdf_beta_match'):

    # organize time-series reference
    ts_ref = ts_ref.dropna()
    ts_values = ts_ref.values
    ts_index = ts_ref.index.values
    ts_n = ts_index.__len__()
    ts_data = {'ref': ts_values}
    ts_dframe_reference = pd.DataFrame(data=ts_data, index=ts_index)
    # organize time-series other
    ts_other = ts_other.dropna()
    ts_values = ts_other.values
    ts_index = ts_other.index.values
    ts_n = ts_index.__len__()
    ts_data = {'other': ts_values}
    ts_dframe_other = pd.DataFrame(data=ts_data, index=ts_index)

    # check other time-series
    if not ts_dframe_other.empty:

        # join time-series common
        if ts_dframe_reference.__len__() > ts_dframe_other.__len__():
            ts_dframe_common = ts_dframe_other.join(ts_dframe_reference)
        elif ts_dframe_reference.__len__() < ts_dframe_other.__len__():
            ts_dframe_common = ts_dframe_reference.join(ts_dframe_other)
        else:
            ts_dframe_common = ts_dframe_reference.join(ts_dframe_other)
        # get info common
        ts_dframe_common = ts_dframe_common.dropna()
        ts_fields_common = list(ts_dframe_common.columns)
        ts_idx_common = ts_fields_common.index('ref')

        # scale datasets using a method defined in pytesmo library
        ts_dframe_scaled = scaling.scale(ts_dframe_common, method=ts_scale_method, reference_index=ts_idx_common)
        # return ts scaled
        ts_scaled = ts_dframe_scaled['other']

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
        log_stream.error()

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

    # Calculate correlation coefficients, RMSD, bias, Nash Sutcliffe
    ts_metrics = {'pearson_r': pearson_r, 'pearson_p': person_p,
                  'spearman_rho': spearman_rho, 'spearman_p': spearman_p,
                  'kendall_tau': kendall_tau, 'kendall_p': kendall_p,
                  'rmsd': rmsd,
                  'bias': bias,
                  'nash_sutcliffe': nash_sutcliffe}

    return ts_metrics
# ----------------------------------------------------------------------------------------------------------------------
