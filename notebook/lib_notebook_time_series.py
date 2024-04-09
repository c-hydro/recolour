"""
Library Features:

Name:          lib_notebook_time_series
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import pickle
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


# ----------------------------------------------------------------------------------------------------------------------
# Method to add year description
def add_dframe_year(dframe_data, column_year='year'):
    dframe_time = dframe_data.index
    grp_year = [pd.Timestamp(t_stamp).year for t_stamp in dframe_time]
    dframe_data[column_year] = grp_year
    return dframe_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to add seasonal description
def add_dframe_seasons(dframe_data, column_season='season', lut_season=None, lut_flag=True):

    # season filter
    if lut_season is None:
        if lut_flag:
            lut_season = {
                1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA',
                7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON', 12: 'DJF'}
        else:
            lut_season = {
                1: 'ALL', 2: 'ALL', 3: 'ALL', 4: 'ALL', 5: 'ALL', 6: 'ALL',
                7: 'ALL', 8: 'ALL', 9: 'ALL', 10: 'ALL', 11: 'ALL', 12: 'ALL'}

    dframe_time = dframe_data.index
    grp_season = [lut_season.get(pd.Timestamp(t_stamp).month) for t_stamp in dframe_time]
    dframe_data[column_season] = grp_season
    return dframe_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to strip time-series
def strip_time_series(ts_dframe, ts_threshold=2):

    # message strip time-series start
    print(' ---> Strip time-series ... ')

    # get time start and time end (according with threshold value)
    time_period_finite = ts_dframe.dropna(thresh=ts_threshold).index.values
    time_start_finite, time_end_finite = time_period_finite[0], time_period_finite[-1]
    # select dframe with start and end time
    ts_dframe = ts_dframe[time_start_finite:time_end_finite]

    # message strip time-series start
    print(' ---> Strip time-series ... DONE')

    return ts_dframe
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to apply time series filter
def apply_time_series_filter(ts_dframe, ts_variable='smap',
                             ts_filter_type='exp', ts_filter_window=3,
                             ts_filter_label='{var_name}_{filter_type}_t{filter_window}'):

    # message filter time-series start
    print(' ---> Filter time-series "' + ts_variable + '" ... ')

    # check ts filter window
    if ts_filter_window is None:
        ts_filter_window = []
    if not isinstance(ts_filter_window, list):
        ts_filter_window = [ts_filter_window]

    # check ts variable
    if ts_variable in list(ts_dframe.columns):

        # get ts series
        ts_series = ts_dframe[ts_variable].dropna()

        # get julian date of ts series
        ts_jd = ts_series.index.to_julian_date().values

        # iterate over filter window(s)
        for step_filter_window in ts_filter_window:

            # string filter window
            string_filter_window = "{:02d}".format(step_filter_window)
            # message filter window start
            print(' ----> Filter window "' + string_filter_window + '" ... ')

            # compute filtered ts
            if ts_filter_type == 'exp':

                # exponential filter fx
                ts_filtered = exp_filter(ts_dframe[ts_variable].values, ts_jd, ctime=step_filter_window)
                # exponential filter label
                ts_label = ts_filter_label.format(
                    var_name=ts_variable, filter_window=string_filter_window, filter_type=ts_filter_type)

            elif ts_filter_type == 'boxcar':

                # boxcar filter fx
                ts_filtered = boxcar_filter(ts_dframe[ts_variable].values, ts_jd, window=step_filter_window)
                # boxcar filter label
                ts_label = (ts_filter_label.format(
                    var_name=ts_variable, filter_window=string_filter_window, filter_type=ts_filter_type))
            else:
                # filter type not supported
                print(' ===> Filter type "' + ts_filter_type + '" not supported')
                raise NotImplementedError('Filter type "' + ts_filter_type + '" not implemented yet')

            # message filter window end
            print(' ----> Filter window "' + string_filter_window + '" ... DONE')

            # update ts dframe
            ts_dframe[ts_label] = ts_filtered

    else:
        # ts variable not available
        print(' ===> Variable "' + ts_variable + '" is not available in the source DataFrame')

    # message filter time-series start
    print(' ---> Filter time-series "' + ts_variable + '" ... DONE')

    return ts_dframe
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to print time series metrics
def print_time_series_metrics(metrics_obj):
    for metrics_datasets, metrics_fields in metrics_obj.items():
        print(' --- Datasets: "' + metrics_datasets + '"')
        for metrics_name, metrics_value in metrics_fields.items():
            print(' ---- ' + metrics_name + ' = ' + metrics_value)
        print(' ')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to dump time series metrics
def dump_time_series_metrics(metrics_obj_src, metrics_point='point', metrics_score_list=None,
                             file_name='metrics_obj.csv', file_sep=';'):

    if metrics_score_list is None:
        metrics_score_list = ['Pearson R', 'Spearman rho', 'Bias', 'Nash Sutcliffe']
    if not isinstance(metrics_score_list, list):
        metrics_score_list = [metrics_score_list]

    if metrics_score_list:
        metrics_obj_collections = {}
        for metrics_type, metrics_fields in metrics_obj_src.items():
            if metrics_score_list:
                metrics_obj_selected = {}
                for metrics_score_name in metrics_score_list:
                    if metrics_score_name in list(metrics_fields.keys()):
                        metrics_score_value = metrics_fields[metrics_score_name]
                        metrics_obj_selected[metrics_score_name] = metrics_score_value
                metrics_obj_collections[metrics_type] = metrics_obj_selected
    else:
        metrics_obj_collections = deepcopy(metrics_obj_src)

    metrics_obj_upd = {"point": metrics_point}
    metrics_obj_upd.update(metrics_obj_collections)

    metrics_obj_dst = pd.DataFrame.from_dict(data=metrics_obj_upd, orient='columns')
    metrics_obj_dst.to_csv(file_name, sep=file_sep)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
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
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
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

    # iterate over time-series
    for ts_other in ts_fields:

        # message scale time-series start
        print(' ---> Scale time-series "' + ts_other + '" ... ')

        # get ref and other time-series
        ts_dframe_ref = ts_dframe_generic[ts_reference]
        ts_dframe_other = ts_dframe_generic[ts_other]

        # analyze reference time-series
        ts_series_ref = ts_dframe_ref.dropna()
        ts_dframe_ref = ts_series_ref.to_frame()
        if ts_dframe_ref.empty:
            print(' ===> Reference time-series "' + ts_reference + '" is empty')
            raise RuntimeError('Reference time-series must be defined in the source DataFrame')

        # analyze other time-series
        ts_series_other = ts_dframe_other.dropna()
        ts_dframe_other = ts_series_other.to_frame()

        # check other time-series
        if not ts_dframe_other.empty:

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

            # message scale time-series end
            print(' ---> Scale time-series "' + ts_other + '" ... DONE')

        else:
            # message scale time-series end (empty)
            print(' ---> Scale time-series "' + ts_other + '" ... SKIPPED. Empty time-series')

    # add control series
    ts_control = np.zeros(shape=[ts_n])
    ts_control[:] = 0
    ts_dframe_control = pd.DataFrame(data={'control':ts_control}, index=ts_index)
    ts_dframe_reference = ts_dframe_reference.join(ts_dframe_control)

    return ts_dframe_reference
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remove time series by variable (or variable root)
def remove_time_series_by_variable(ts_dframe_generic, ts_remove_variables=None):

    # message remove time-series start
    print(' ---> Remove time-series variable ... ')

    # check ts remove variables format
    if ts_remove_variables is None:
        ts_remove_variables = []
    if not isinstance(ts_remove_variables, list):
        ts_remove_variables = [ts_remove_variables]

    # list of dataset variables
    ts_dframe_variables = list(ts_dframe_generic.columns)

    ts_select_variables = []
    for step_remove_variable in ts_remove_variables:
        if step_remove_variable in ts_dframe_variables:
            ts_select_variables.append(step_remove_variable)

    # drop dataframe variable(s)
    ts_dframe_generic = ts_dframe_generic.drop(columns=ts_select_variables)

    # message remove time-series end
    print(' ---> Remove time-series variable ... DONE')

    return ts_dframe_generic
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select time series by variable (or variable root)
def select_time_series_by_variable(ts_dframe_generic, ts_search_variables=None):

    # message select time-series start
    print(' ---> Select time-series variable ... ')

    # check ts search variables format
    if ts_search_variables is None:
        ts_search_variables = []
    if not isinstance(ts_search_variables, list):
        ts_search_variables = [ts_search_variables]

    # list of dataset variables
    ts_dframe_variables = list(ts_dframe_generic.columns)

    ts_select_variables = []
    for step_search_variable in ts_search_variables:
        tmp_select_variables = [step_variable for step_variable in ts_dframe_variables
                                if step_search_variable in step_variable]
        ts_select_variables.extend(tmp_select_variables)

    # iterate over selected variable(s)
    ts_dframe_selected = None
    for step_select_variable in ts_select_variables:

        # message select variable start
        print(' ----> Variable "' + step_select_variable + '" ... ')

        # get values and index
        ts_values = ts_dframe_generic[step_select_variable].values
        ts_index = ts_dframe_generic.index.values

        if ts_dframe_selected is None:
            ts_dframe_selected = pd.DataFrame(index=ts_index, data={step_select_variable: ts_values})
        else:
            ts_dframe_selected[step_select_variable] = ts_values

        # message select variable end
        print(' ----> Variable "' + step_select_variable + '" ... DONE')

    # add seasons and year columns
    ts_dframe_selected = add_dframe_seasons(ts_dframe_selected, column_season='season', lut_season=None, lut_flag=True)
    ts_dframe_selected = add_dframe_year(ts_dframe_selected, column_year='year')

    # message select time-series end
    print(' ---> Select time-series variable ... DONE')

    return ts_dframe_selected
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select time series by season
def select_time_series_by_season(ts_dframe_generic, season_filter='DJF,MAM,JJA,SON'):

    # message select time-series start
    print(' ---> Select time-series seasons ... ')

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

    # message select time-series end
    print(' ---> Select time-series seasons ... DONE')

    return ts_dframe_filtered
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to join time series
def join_time_series(ts_dframe_product, ts_dframe_datasets):

    # message join time-series start
    print(' ---> Join time-series product and datasets ... ')

    if 'season' in list(ts_dframe_datasets.columns):
        ts_dframe_datasets = ts_dframe_datasets.drop(columns=['season'])
    if 'year' in list(ts_dframe_datasets.columns):
        ts_dframe_datasets = ts_dframe_datasets.drop(columns=['year'])
    ts_dframe_common = ts_dframe_product.join(ts_dframe_datasets)

    # message join time-series end
    print(' ---> Join time-series product and datasets ... DONE')

    return ts_dframe_common

# ----------------------------------------------------------------------------------------------------------------------
