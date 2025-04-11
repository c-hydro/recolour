#!/usr/bin/python3

"""
RECOLOUR APPS - CELL2GRID ASCAT - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20240423'
__version__ = '1.1.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_map_cell2grid_ascat.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20240423 (1.1.0) --> Update codes and fix bugs
20240327 (1.0.0) --> First development
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import time

from copy import deepcopy
import argparse

from lib_info_args import logger_name, logger_format, time_format_algorithm

from lib_info_settings import get_data_settings
from lib_utils_time import set_time_info

from drv_map_static import DrvMap as DrvMapStatic
from drv_map_dynamic import DrvMap as DrvMapDynamic

# set logger
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for resampling cells 2 grid ascat'
alg_type = 'Package'
alg_version = '1.1.0'
alg_release = '2024-04-23'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # get file settings
    alg_file_settings, alg_time_settings = get_args()
    # read data settings
    alg_data_settings = get_data_settings(alg_file_settings)
    # set logging
    set_logging(logger_name=logger_name, logger_format=logger_format,
                logger_folder=alg_data_settings['log']['path_log'],
                logger_file=alg_data_settings['log']['file_log'])
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # info algorithm (start)
    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    # time algorithm
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Organize time information
    alg_time_run, alg_time_range, alg_time_chunks = set_time_info(
        time_run_args=alg_time_settings,
        time_run_file=None,
        time_run_file_start=None,time_run_file_end=None,
        time_format=time_format_algorithm,
        time_period=1, time_frequency='H', time_rounding='H', time_reverse=True)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # configure static driver
    drv_data_static = DrvMapStatic(alg_data_settings)
    # organize static datasets
    alg_data_static = drv_data_static.organize_data()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # configure dynamic driver
    drv_data_dynamic = DrvMapDynamic(alg_time_run, alg_data_static, alg_data_settings)
    # organize datasets
    collections_data_dynamic = drv_data_dynamic.organize_data()
    # analyze data
    drv_data_dynamic.analyze_data(collections_data_dynamic)
    # dump datasets
    drv_data_dynamic.dump_data()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(' ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    alg_logger.info(' ==> ... END')
    alg_logger.info(' ==> Bye, Bye')
    alg_logger.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to get script argument(s)
def get_args():
    # parser algorithm arg(s)
    parser_obj = argparse.ArgumentParser()
    parser_obj.add_argument('-settings_file', action="store", dest="settings_file")
    parser_obj.add_argument('-time', action="store", dest="arg_time")
    parser_value = parser_obj.parse_args()
    # set algorithm arg(s)
    settings_file, arg_time = 'configuration.json', None
    if parser_value.settings_file:
        settings_file = parser_value.settings_file
    if parser_value.arg_time:
        arg_time = parser_value.arg_time

    return settings_file, arg_time

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_name='algorithm_logger', logger_folder=None, logger_file='log.txt', logger_format=None):

    if logger_format is None:
        logger_format = deepcopy(logger_format_def)
    if logger_file is None:
        logger_file = deepcopy(logger_file_def)

    if logger_folder is not None:
        logger_path = os.path.join(logger_folder, logger_file)
    else:
        logger_path = deepcopy(logger_file)

    logger_loc = os.path.split(logger_path)
    if logger_loc[0] == '' or logger_loc[0] == "":
        logger_folder_name, logger_file_name = os.path.dirname(os.path.abspath(sys.argv[0])), logger_loc[1]
    else:
        logger_folder_name, logger_file_name = logger_loc[0], logger_loc[1];

    os.makedirs(logger_folder_name, exist_ok=True)

    # define logger path
    logger_path = os.path.join(logger_folder_name, logger_file_name)

    # Remove old logging file
    if os.path.exists(logger_path):
        os.remove(logger_path)

    # Open logger
    logging.getLogger(logger_name)
    logging.root.setLevel(logging.DEBUG)

    # Open logging basic configuration
    logging.basicConfig(level=logging.DEBUG, format=logger_format, filename=logger_path, filemode='w')

    # Set logger handle
    logger_handle_1 = logging.FileHandler(logger_path, 'w')
    logger_handle_2 = logging.StreamHandler()
    # Set logger level
    logger_handle_1.setLevel(logging.DEBUG)
    logger_handle_2.setLevel(logging.DEBUG)
    # Set logger formatter
    logger_formatter = logging.Formatter(logger_format)
    logger_handle_1.setFormatter(logger_formatter)
    logger_handle_2.setFormatter(logger_formatter)

    # Add handle to logging
    logging.getLogger('').addHandler(logger_handle_1)
    logging.getLogger('').addHandler(logger_handle_2)

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# call script from external library
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------------------------
