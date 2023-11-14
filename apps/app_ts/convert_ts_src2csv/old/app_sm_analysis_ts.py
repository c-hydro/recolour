#!/usr/bin/python3
"""
ARPAL Analysis Tool - SM TIME-SERIES
__date__ = '20221122'
__version__ = '1.0.0'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org',
        'Francesco Silvestro (francesco.silvestro@cimafoundation.org)'

__library__ = 'SM'

General command line:
python3 app_sm_ts_analysis.py -settings_file configuration.json

Version(s):
20221122 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
import os

from driver_data_src import DriverData as DriverData_Static
from driver_data_io_dynamic_grid import DriverData as DriverData_Dynamic_Grid
from driver_data_io_dynamic_point import DriverData as DriverData_Dynamic_Point

from driver_analysis import DriverAnalysis

from argparse import ArgumentParser

from lib_data_io_json import read_file_json
from lib_utils_time import set_time
from lib_utils_logging import set_logging_file
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_version = '1.0.0'
alg_release = '2022-11-22'
alg_name = 'SM ANALYSIS'
# Algorithm parameter(s)
time_format = '%Y-%m-%d %H:%M'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    # Set algorithm logging
    set_logging_file(
        logger_name=logger_name,
        logger_file=os.path.join(data_settings['log']['folder_name'], data_settings['log']['file_name']))

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    log_stream.info(' ============================================================================ ')
    log_stream.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    log_stream.info(' ==> START ... ')
    log_stream.info(' ')

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Organize time run
    time_reference, time_range = set_time(
        time_ref_args=alg_time,
        time_ref_file=data_settings['time']['time_reference'],
        time_ref_file_start=data_settings['time']['time_start'],
        time_ref_file_end=data_settings['time']['time_end'],
        time_format=time_format,
        time_period=data_settings['time']['time_period'],
        time_frequency=data_settings['time']['time_frequency'],
        time_rounding=data_settings['time']['time_rounding']
    )
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Geographical datasets
    driver_data_static = DriverData_Static(
        src_dict=data_settings['data']['static']['source'],
        anc_dict=data_settings['data']['static']['ancillary'],
        alg_dict=data_settings['algorithm']['ancillary'],
        tmp_dict=data_settings['tmp'],
        template_tags_dict=data_settings['algorithm']['template'],
        flag_data_updating=data_settings['algorithm']['flags']['updating_ancillary_static'])
    data_static_collection = driver_data_static.organize_data()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Iterate over time(s)
    for time_step in time_range:

        # Soil moisture point datasets
        driver_data_dynamic_point = DriverData_Dynamic_Point(
            time_step, time_reference,
            src_dict=data_settings['data']['dynamic']['source'],
            anc_dict=data_settings['data']['dynamic']['ancillary'],
            alg_dict=data_settings['algorithm']['ancillary'],
            time_dict=data_settings['data']['dynamic']['time'],
            geo_dict=data_static_collection,
            tmp_dict=data_settings['tmp'],
            template_tags_dict=data_settings['algorithm']['template'],
            flag_data_updating=data_settings['algorithm']['flags']['updating_ancillary_dynamic_point'])
        driver_data_dynamic_point.organize_data()

        # Soil moisture grid datasets
        driver_data_dynamic_grid = DriverData_Dynamic_Grid(
            time_step, time_reference,
            src_dict=data_settings['data']['dynamic']['source'],
            anc_dict=data_settings['data']['dynamic']['ancillary'],
            alg_dict=data_settings['algorithm']['ancillary'],
            time_dict=data_settings['data']['dynamic']['time'],
            geo_dict=data_static_collection,
            tmp_dict=data_settings['tmp'],
            template_tags_dict=data_settings['algorithm']['template'],
            flag_data_updating=data_settings['algorithm']['flags']['updating_ancillary_dynamic_grid'])
        driver_data_dynamic_grid.organize_data()

        # Soil moisture analysis
        driver_analysis = DriverAnalysis(
            time_step, time_reference,
            anc_dict=data_settings['data']['dynamic']['ancillary'],
            dst_dict=data_settings['data']['dynamic']['destination'],
            alg_dict=data_settings['algorithm']['ancillary'],
            time_dict=data_settings['data']['dynamic']['time'],
            geo_dict=data_static_collection,
            tmp_dict=data_settings['tmp'],
            template_tags_dict=data_settings['algorithm']['template'],
            flag_data_updating=data_settings['algorithm']['flags']['updating_ancillary_analysis'])

        driver_analysis.organize_data()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    log_stream.info(' ')
    log_stream.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    log_stream.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    log_stream.info(' ==> ... END')
    log_stream.info(' ==> Bye, Bye')
    log_stream.info(' ============================================================================ ')
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
    main()
# ----------------------------------------------------------------------------
