#!/usr/bin/python3

"""
RECOLOUR TOOLS - ASCAT SWATH2CELL - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20250610'
__version__ = '1.4.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_swath2cell_ascat.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20250610 (1.4.0) --> Fix bugs in managing file list
20240322 (1.3.0) --> Fix bugs in nc file and update code(s) for improving log messages
20240218 (1.2.0) --> Code refactor for the data record and nrt mode
20231128 (1.1.0) --> Code refactor for the recolour package
20230804 (1.0.0) --> First development
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# libraries
import os
import sys
import pandas as pd
import time
import logging

from argparse import ArgumentParser
from copy import deepcopy

from lib_info_args import logger_name
from lib_info_args import logger_file as logger_file_def, logger_format as logger_format_def
from lib_info_args import time_format_algorithm as time_format

from lib_info_settings import get_data_settings
from lib_utils_time import set_time

from drv_fx_wrapper import DrvFxWrapper
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'swath2cell'
alg_type = 'Application'
alg_version = '1.4.0'
alg_release = '2025-06-10'
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# script main
def main_wrapper():

    # -------------------------------------------------------------------------------------
    # method to get script argument(s)
    alg_file_settings, alg_time_settings = get_args()

    # read data settings
    alg_data_settings = get_data_settings(alg_file_settings)

    # set logging
    set_logging(logger_name=logger_name,
                logger_folder=alg_data_settings['log']['folder_name'],
                logger_file=alg_data_settings['log']['file_name'])
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # info algorithm (start)
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ')

    # time algorithm
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Organize time information
    alg_time_run, alg_time_chunks = set_time(
        time_run_args=alg_time_settings,
        time_run_file=alg_data_settings['time']['time_now'],
        time_run_file_start=alg_data_settings['time']['time_start'],
        time_run_file_end=alg_data_settings['time']['time_end'],
        time_format=time_format,
        time_period=alg_data_settings['time']['time_period'],
        time_frequency=alg_data_settings['time']['time_frequency'],
        time_rounding=alg_data_settings['time']['time_rounding'], time_reverse=True)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # iterate over time information
    for time_reference, time_info in alg_time_chunks.items():

        # -------------------------------------------------------------------------------------
        # initialize driver
        drv_fx = DrvFxWrapper(
            alg_settings=alg_data_settings,
            alg_time_reference=time_reference,
            alg_time_start=time_info[0], alg_time_end=time_info[1])
        # organize fx args
        alg_fx_args = drv_fx.organize_fx_args()
        # organize fx classes
        alg_fx_classes, alg_fx_time = drv_fx.organize_fx_classes(alg_fx_args)
        # execute fx
        drv_fx.execute_fx(alg_fx_classes, alg_fx_time)
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
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

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_name='algorithm_logger', logger_folder=None,
                logger_file='log.txt', logger_format=None):

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
        logger_folder_name, logger_file_name = logger_loc[0], logger_loc[1]

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
    main_wrapper()
# -------------------------------------------------------------------------------------




