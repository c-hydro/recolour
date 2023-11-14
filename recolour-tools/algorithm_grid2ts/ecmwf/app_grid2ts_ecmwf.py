#!/usr/bin/python3

"""
RECOLOUR TOOLS - ECMWF GRID2TS - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20230919'
__version__ = '1.5.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
    'Martina Natali (martina01.natali@edu.unife.it)'
__library__ = 'recolour'

General command line:
python app_grid2ts_ecmwf.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20230918 (1.5.0) --> Add "data_record" and "nrt" mode
20230803 (1.0.0) --> First development
"""
# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import time
from argparse import ArgumentParser

from lib_info_args import time_format_datasets as time_format
from lib_utils_time import set_time_info, update_time_info
from lib_info_settings import get_data_settings, parse_data_settings
from lib_utils_ecmwf import create_file_grid

from lib_reshuffle_ecmwf import main as main_runner
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'grid2ts'
alg_type = 'Package'
alg_version = '1.5.0'
alg_release = '2023-09-19'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# script main
def main_wrapper():

    # -------------------------------------------------------------------------------------
    # method to get script argument(s)
    file_settings, time_settings = get_args()
    # set logging
    set_logging()

    # read data settings
    data_settings = get_data_settings(file_settings)
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
    # method to organize time information
    alg_time_now, alg_time_run, alg_time_start, alg_time_end = set_time_info(
        time_run_args=time_settings, time_run_file=data_settings['time']['time_now'],
        time_run_file_start=data_settings['time']['time_start'],
        time_run_file_end=data_settings['time']['time_end'],
        time_format=time_format,
        time_period=data_settings['time']['time_period'],
        time_frequency=data_settings['time']['time_frequency'],
        time_rounding=data_settings['time']['time_rounding'])
    # method to update time information
    data_settings = update_time_info(data_settings, alg_time_run, alg_time_start, alg_time_end)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # parse data settings
    args_settings = parse_data_settings(data_settings)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # create data grid
    create_file_grid(
        grid_path=os.path.join(data_settings['grid']['folder_name'], data_settings['grid']['file_name']),
        data_path=data_settings['data']['path_grid'],
        data_ext='.nc', grid_update=data_settings['flags']['reset_static'],
    )
    # -------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # run procedure to convert grid to time-series datasets
    main_runner(args_settings)
    # ----------------------------------------------------------------------------------------------------------------------

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
def set_logging(logger_file='./log.txt', logger_format=None):
    if logger_format is None:
        logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                        '%(filename)s:[%(lineno)-6s - %(funcName)20s()] %(message)s'

    logger_loc = os.path.split(logger_file)
    if logger_loc[0] == "":
        logger_folder = os.path.dirname(__file__)
    else:
        logger_folder = logger_loc[0]
    os.makedirs(logger_folder, exist_ok=True)

    # Remove old logging file
    if os.path.exists(logger_file):
        os.remove(logger_file)

    # Set level of root debugger
    logging.root.setLevel(logging.DEBUG)

    # Open logging basic configuration
    logging.basicConfig(level=logging.DEBUG, format=logger_format, filename=logger_file, filemode='w')

    # Set logger handle
    logger_handle_1 = logging.FileHandler(logger_file, 'w')
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
