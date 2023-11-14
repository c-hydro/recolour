#!/usr/bin/python3
"""
RECOLOUR NOTEBOOK - TIME-SERIES ANALYSIS  - REprocess paCkage for sOiL mOistUre pRoducts
__date__ = '20231023'
__version__ = '1.0.0'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'SM'

General command line:
python3 app_sm_ts_analysis.py -settings_file configuration.json

Version(s):
20231023 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
import os
import datetime

import numpy as np

from argparse import ArgumentParser

from lib_notebook_io_json import (read_time_series_product, read_time_series_datasets,
                                  join_time_series, scale_time_series,
                                  print_time_series_metrics, compute_time_series_metrics,
                                  filter_time_series_by_season)
from lib_notebook_io_generic import read_file_json
from lib_notebook_geo import read_point_data
from lib_notebook_system import fill_tags2string, make_folder

import matplotlib.pylab as plt
import matplotlib.dates as mdates
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Notebook for analyzing time-series reference and k1'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2023-10-23'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main(file_name_settings):

    # -------------------------------------------------------------------------------------
    # set algorithm settings
    data_settings = read_file_json(file_name_settings)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    print(' ============================================================================ ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> START ... ')
    print(' ')

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # get settings object(s)
    obj_parameters = data_settings['parameters']
    obj_template = data_settings['template']
    obj_data = data_settings['data']
    obj_figure = data_settings['figure']
    obj_metrics = data_settings['metrics']
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # get registry info
    folder_name_registry = obj_data['registry']['folder_name']
    file_name_registry = obj_data['registry']['file_name']
    file_path_registry = os.path.join(folder_name_registry, file_name_registry)
    # get registry data
    registry_dframe = read_point_data(file_path_registry)
    registry_name, registry_tag = list(registry_dframe['name']), list(registry_dframe['tag'])

    # select point
    registry_tag_value = obj_parameters['point_name']
    registry_point = registry_dframe.loc[registry_dframe['tag'] == registry_tag_value]
    if registry_point.empty:
        name_list = ','.join(registry_dframe['tag'])
        print('Point available: "' + name_list + '"')
        raise RuntimeError(' ===> Point "' + registry_tag_value + '" is not available in the registry')
    registry_dict = registry_point.to_dict('records')[0]
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # get time-series info
    folder_name_ts = obj_data['time_series']['folder_name']
    file_name_ts = obj_data['time_series']['file_name']
    file_path_ts_template = os.path.join(folder_name_ts, file_name_ts)

    obj_values = {'point_name': registry_tag_value}
    file_path_ts_point = fill_tags2string(file_path_ts_template, obj_template, obj_values)[0]

    # get time-series data datasets
    (ts_dframe_datasets,
     ts_dframe_ecmwf_l1, ts_dframe_ecmwf_l2, ts_dframe_ecmwf_l3,
     ts_dframe_hmc, ts_dframe_smap) = read_time_series_datasets(file_path_ts_point)

    # get time-series data product
    ts_dframe_product, ts_dframe_obs, ts_dframe_tc = read_time_series_product(file_path_ts_point)

    # join time-series
    ts_dframe_common = join_time_series(ts_dframe_datasets, ts_dframe_product)
    # filter time-series by season
    ts_dframe_filtered = filter_time_series_by_season(ts_dframe_common, 'ALL')
    # scale time-series
    ts_dframe_scaled = scale_time_series(ts_dframe_filtered, ts_reference='sm_obs', ts_scale_method='cdf_beta_match')
    # compute time-series metrics
    ts_obj_metrics = compute_time_series_metrics(ts_dframe_scaled, ts_reference='sm_obs')
    # print time-series metrics
    print_time_series_metrics(ts_obj_metrics)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # select finite period (to skip the control series)
    time_period_finite = ts_dframe_scaled.dropna(thresh=2).index.values
    time_start_finite, time_end_finite = time_period_finite[0], time_period_finite[-1]

    ts_dframe_scaled = ts_dframe_scaled[time_start_finite:time_end_finite]
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Plot the datasets
    fig, ax = plt.subplots()

    ts_dframe_scaled['control'].plot(figsize=(10, 5), ax=ax, color='white')
    ts_dframe_scaled['sm_obs'].plot(figsize=(10, 5), ax=ax, marker='o', color='black')
    ts_dframe_scaled['sm_tc'].plot(figsize=(10, 5), ax=ax, color='blue')
    ts_dframe_scaled['sm_ecmwf_layer_0_7'].plot(figsize=(10, 5), ax=ax, color='green')
    ts_dframe_scaled['sm_hmc'].plot(figsize=(10, 5), ax=ax, color='red')
    ts_dframe_scaled['sm_smap'].plot(figsize=(10, 5), ax=ax, color='pink')
    #ts_dframe_scaled.plot(figsize=(10, 5), ax=ax)  # , color='red')
    ax.set_xlabel('time')
    ax.set_ylabel('sm [%]')
    ax.legend(["", "obs", "tc", "ecmwf_l1", "hmc", "smap"])
    ax.grid(True)

    plt.ylim(0, 100)
    plt.title('sm - datasets')
    plt.show()

    # Info
    print(' ==> Plot the time-series figure ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # set figure info
    folder_name_figure = obj_figure['folder_name']
    file_name_figure = obj_figure['file_name']
    file_path_figure_template = os.path.join(folder_name_figure, file_name_figure)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # set metrics info
    folder_name_metrics = obj_metrics['folder_name']
    file_name_metrics = obj_metrics['file_name']
    file_path_metrics_template = os.path.join(folder_name_metrics, file_name_metrics)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    print(' ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    print(' ==> ... END')
    print(' ==> Bye, Bye')
    print(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-time', action="store", dest="alg_time")
    parser_values = parser_handle.parse_args()

    alg_settings, alg_time = 'configuration.json', None
    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    if parser_values.alg_time:
        alg_time = parser_values.alg_time

    return alg_settings, alg_time

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == '__main__':

    folder_name_settings = os.getcwd()
    file_name_settings = 'notebook_recolour_sm_ts_hmc_datasets_tc.json'

    path_name_settings = os.path.join(folder_name_settings, file_name_settings)

    main(path_name_settings)
# ----------------------------------------------------------------------------
