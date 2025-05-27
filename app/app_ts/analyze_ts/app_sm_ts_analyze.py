#!/usr/bin/python3
"""
RECOLOUR APP - TIME-SERIES - METRICS ANALYZER  - REprocess paCkage for sOiL mOistUre pRoducts
__date__ = '20250425'
__version__ = '1.3.0'
__author__ = 'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'recolour'

General command line:
python3 app_sm_ts_analyze.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20250425 (1.3.0) --> Fix bugs for the operational release
20240705 (1.2.0) --> Add new features in managing times for data record and nrt time-series; fix bugs in write csv datasets
20240305 (1.1.0) --> Operational release
20240108 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
import os

from driver_data_static import DriverData as DriverData_Static
from driver_data_dynamic import DriverData as DriverData_Dynamic

from argparse import ArgumentParser

from lib_data_io_json import read_file_json
from lib_utils_time import set_time
from lib_utils_logging import set_logging_file
from lib_info_args import logger_name, time_format_algorithm

# Logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for analyzing and computing the time-series metrics'
alg_type = 'Package'
alg_version = '1.3.0'
alg_release = '2025-04-25'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # get algorithm settings
    alg_settings, alg_time = get_args()

    # set algorithm settings
    data_settings = read_file_json(alg_settings)

    # set algorithm logging
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
    time_run, time_range = set_time(
        time_ref_args=alg_time,
        time_ref_file=data_settings['time']['time_reference'],
        time_ref_file_start=data_settings['time']['time_start'],
        time_ref_file_end=data_settings['time']['time_end'],
        time_format=time_format_algorithm,
        time_period=data_settings['time']['time_period'],
        time_frequency=data_settings['time']['time_frequency'],
        time_rounding=data_settings['time']['time_rounding']
    )
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # geographical datasets
    driver_data_static = DriverData_Static(
        src_dict=data_settings['data']['static']['source'],
        dst_dict=data_settings['data']['static']['destination'],
        tmp_dict=data_settings['tmp'],
        template_dict=data_settings['algorithm']['template'],
        flags_dict=data_settings['algorithm']['flags'])
    data_static_collection = driver_data_static.organize_data()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # source datasets
    driver_data_dynamic = DriverData_Dynamic(
        time_run=time_run, time_range=time_range,
        static_obj=data_static_collection,
        source_dict=data_settings['data']['dynamic']['source'],
        ancillary_dict=data_settings['data']['dynamic']['ancillary'],
        destination_dict=data_settings['data']['dynamic']['destination'],
        methods_dict=data_settings['algorithm']['methods'],
        template_dict=data_settings['algorithm']['template'],
        flags_dict=data_settings['algorithm']['flags'],
        tmp_dict=data_settings['tmp'])
    # method to organize data
    data_source_collection = driver_data_dynamic.organize_data()
    # method to analyze data
    data_analyzed_collection = driver_data_dynamic.analyze_data(data_source_collection)
    # method to dump data
    driver_data_dynamic.dump_data(data_analyzed_collection)
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
