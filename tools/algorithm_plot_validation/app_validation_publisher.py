#!/usr/bin/python3

"""
RECOLOUR TOOLS - PUBLISHER APP - VALIDATION - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20240410'
__version__ = '1.1.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org),
    'Martina Natali (martina01.natali@edu.unife.it)'

__library__ = 'recolour'

General command line:
python app_validation_publisher.py -settings_file configuration.json

Version(s):
20240410 (1.1.0) --> Add codes for boxplot, update codes and fix bugs
20230727 (1.0.0) --> First development
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import time
import sys

import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')

from copy import deepcopy
import argparse

from lib_info_args import logger_name as logger_name_def
from lib_info_args import logger_format as logger_format_def
from lib_info_settings import get_data_settings
from lib_utils_grid import get_grid_cells

from drv_data_publisher import DrvData
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'publisher'
alg_type = 'Package'
alg_version = '1.1.0'
alg_release = '2024-04-10'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# script main
def main():

    # -------------------------------------------------------------------------------------
    # get file settings
    file_settings = get_args()
    # read data settings
    alg_settings = get_data_settings(file_settings)
    # set logging
    set_logging(logger_name=logger_name_def, logger_format=logger_format_def,
                logger_folder=alg_settings['log']['path_log'],
                logger_file=alg_settings['log']['file_log'])
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
    # generate reference cells and gpis
    alg_cell_list, alg_gpis, alg_cell_grid = get_grid_cells(
        cell_start=alg_settings['domain']['cell_start'],
        cell_end=alg_settings['domain']['cell_end'],
        cells_list=alg_settings['domain']['cell_list'],
        path_grid=alg_settings['grid']['folder_name'],
        file_grid=alg_settings['grid']['file_name'])
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # configure driver
    drv_driver_publisher = DrvData(alg_cell_list, alg_cell_grid, alg_settings)
    # organize datasets
    dframe_cells = drv_driver_publisher.organize_data()
    # publish datasets
    drv_driver_publisher.plot_data(dframe_cells)
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
    # parser algorithm arg(s)
    parser_obj = argparse.ArgumentParser()
    parser_obj.add_argument('-settings_file', action="store", dest="settings_file")
    parser_value = parser_obj.parse_args()
    # set algorithm arg(s)
    settings_file = 'configuration.json'
    if parser_value.settings_file:
        settings_file = parser_value.settings_file

    return settings_file

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_name='algorithm_logger', logger_folder=None, logger_file='log.txt', logger_format=None):
    if logger_format is None:
        logger_format = deepcopy(logger_format_def)
    if logger_file is None:
        logger_file = deepcopy(logger_file)

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
    logging.basicConfig(level=logging.DEBUG, format=logger_format, filename=logger_file, filemode='w')

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
