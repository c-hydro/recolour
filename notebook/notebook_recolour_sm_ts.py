#!/usr/bin/python3
"""
RECOLOUR NOTEBOOK - TIME-SERIES ANALYSIS  - REprocess paCkage for sOiL mOistUre pRoducts
__date__ = '2024-02-29'
__version__ = '1.6.0'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'recolour'

General command line:
python3 notebook_recolour_sm_ts.py -settings_file configuration.json

Version(s):
20240229 (1.6.0) --> Update methods to multiyear analysis
20231123 (1.5.0) --> Beta release
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import time
import os
from argparse import ArgumentParser
from copy import deepcopy

from lib_notebook_io_json import read_file_settings, read_time_series_product, read_time_series_datasets

from lib_notebook_time_series import (join_time_series, scale_time_series, strip_time_series,
                                      print_time_series_metrics, dump_time_series_metrics, compute_time_series_metrics,
                                      select_time_series_by_season,
                                      select_time_series_by_variable, remove_time_series_by_variable,
                                      apply_time_series_filter)
from lib_notebook_geo import read_point_data
from lib_notebook_generic import fill_tags2string, make_folder
from lib_notebook_figure import plot_ts_generic, plot_ts_by_year

import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Notebook for analyzing time-series reference and k1'
alg_type = 'Package'
alg_version = '1.6.0'
alg_release = '2024-02-29'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Script Main
def main(file_name_settings):

    # ------------------------------------------------------------------------------------------------------------------
    # set algorithm settings
    data_settings = read_file_settings(file_name_settings)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm
    print(' ============================================================================ ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> START ... ')
    print(' ')

    # Time algorithm information
    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get settings object(s)
    obj_parameters = data_settings['parameters']
    obj_template = data_settings['template']
    obj_data = data_settings['data']
    obj_figure = data_settings['figure']
    obj_metrics = data_settings['metrics']
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set figure info
    folder_name_figure = obj_figure['folder_name']
    file_name_figure = obj_figure['file_name']
    file_path_figure_template = os.path.join(folder_name_figure, file_name_figure)

    # set metrics info
    folder_name_metrics = obj_metrics['folder_name']
    file_name_metrics = obj_metrics['file_name']
    file_path_metrics_template = os.path.join(folder_name_metrics, file_name_metrics)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
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
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
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
    # select time-series by season
    ts_dframe_filtered = select_time_series_by_season(ts_dframe_common, 'ALL')

    # apply time-series filter (to select variable)
    ts_dframe_filtered = apply_time_series_filter(
        ts_dframe_filtered,
        ts_variable='sm_smap', ts_filter_type='exp', ts_filter_window=[2, 3, 4])
    ts_dframe_filtered = apply_time_series_filter(
        ts_dframe_filtered,
        ts_variable='sm_smap', ts_filter_type='boxcar', ts_filter_window=[2, 3, 4])

    # scale time-series
    ts_dframe_scaled = scale_time_series(ts_dframe_filtered, ts_reference='sm_obs', ts_scale_method='cdf_beta_match')
    # compute time-series metrics
    ts_obj_metrics = compute_time_series_metrics(ts_dframe_scaled, ts_reference='sm_obs')
    # print time-series metrics
    print_time_series_metrics(ts_obj_metrics)

    # dump metrics values to file
    obj_values = {'point_name': registry_tag_value, 'season_name': 'ALL'}
    file_path_metrics_def = fill_tags2string(file_path_metrics_template, obj_template, obj_values)[0]
    folder_name_metrics_def, file_name_metrics_def = os.path.split(file_path_metrics_def)
    os.makedirs(folder_name_metrics_def, exist_ok=True)

    dump_time_series_metrics(
        ts_obj_metrics, file_name=file_path_metrics_def,
        metrics_point=registry_tag_value,
        metrics_score_list=['Pearson R', 'Spearman rho', 'Bias', 'Nash Sutcliffe'])

    # hmc datasets
    # select time-series by variable
    ts_dframe_hmc = select_time_series_by_variable(
        ts_dframe_scaled, ts_search_variables=['control', 'sm_obs', 'sm_hmc'])

    # ecmwf datasets
    # select time-series by variable
    ts_dframe_ecmwf = select_time_series_by_variable(
        ts_dframe_scaled, ts_search_variables=['control', 'sm_obs', 'sm_ecmwf'])

    # smap datasets
    # select time-series by variable
    ts_dframe_smap = select_time_series_by_variable(
        ts_dframe_scaled, ts_search_variables=['control', 'sm_obs', 'sm_smap'])
    # remove time-series by variable
    ts_dframe_smap = remove_time_series_by_variable(ts_dframe_smap, ts_remove_variables=['sm_smap'])
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # select finite period (to skip the control series)
    ts_dframe_scaled = strip_time_series(ts_dframe_scaled, ts_threshold=2)
    ts_dframe_hmc = strip_time_series(ts_dframe_hmc, ts_threshold=2)
    ts_dframe_ecmwf = strip_time_series(ts_dframe_ecmwf, ts_threshold=2)
    ts_dframe_smap = strip_time_series(ts_dframe_smap, ts_threshold=2)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # plot yearly ecmwf time-series
    plot_ts_by_year(ts_dframe_ecmwf, tag_ref='sm_obs', tag_kn=['sm_ecmwf_layer_0_7', 'sm_ecmwf_layer_0_28'],
                    tag_year=[2021, 2022, 2023],
                    figure_legend=['obs', 'sm_ecmwf_layer_0_7', 'sm_ecmwf_layer_0_28'],
                    figure_title='Soil Moisture ECMWF', figure_point=registry_tag_value)

    # plot yearly hmc time-series
    plot_ts_by_year(ts_dframe_hmc, tag_ref='sm_obs', tag_kn=['sm_hmc'],
                    tag_year=[2021, 2022, 2023],
                    figure_legend=['obs', 'sm_hmc'],
                    figure_title='Soil Moisture Continuum', figure_point=registry_tag_value)

    # plot yearly smap time-series
    plot_ts_by_year(ts_dframe_smap, tag_ref='sm_obs', tag_kn=['sm_smap_exp_t02', 'sm_smap_exp_t04'],
                    tag_year=[2021, 2022, 2023],
                    figure_legend=['obs', 'sm_smap_exp_t02', 'sm_smap_exp_t04'],
                    figure_title='Soil Moisture SMAP', figure_point=registry_tag_value)

    # plot generic smap time-series
    plot_ts_generic(ts_dframe_smap, tag_ref='sm_obs', tag_kn=['sm_smap_exp_t02', 'sm_smap_exp_t04'],
                    figure_legend=['obs', 'sm_smap_exp_t02', 'sm_smap_exp_t04'],
                    figure_title='Soil Moisture SMAP', figure_point=registry_tag_value)

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set figure info
    folder_name_figure = obj_figure['folder_name']
    file_name_figure = obj_figure['file_name']
    file_path_figure_template = os.path.join(folder_name_figure, file_name_figure)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set metrics info
    folder_name_metrics = obj_metrics['folder_name']
    file_name_metrics = obj_metrics['file_name']
    file_path_metrics_template = os.path.join(folder_name_metrics, file_name_metrics)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    print(' ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    print(' ==> ... END')
    print(' ==> Bye, Bye')
    print(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get script argument(s)
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

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# call script from external library
if __name__ == '__main__':

    folder_name_settings = os.getcwd()
    # case tc italy
    # file_name_settings = 'notebook_recolour_sm_ts_datasets_tc_hmc_italy.json'
    # case obs liguria
    file_name_settings = 'notebook_recolour_sm_ts_datasets_obs_liguria.json'

    path_name_settings = os.path.join(folder_name_settings, file_name_settings)

    main(path_name_settings)
# ----------------------------------------------------------------------------------------------------------------------
