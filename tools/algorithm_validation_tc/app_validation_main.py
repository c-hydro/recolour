#!/usr/bin/python3

"""
RECOLOUR TOOLS - VALIDATION APP - REprocess paCkage for sOiL mOistUre pRoductst

__date__ = '20230719'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)',
    'Martina Natali (martina01.natali@edu.unife.it)'

__library__ = 'recolour'

General command line:
python app_validation_main.py -settings_file configuration.json

Version(s):
20230726 (1.0.1) --> Update to read SMAP data
20230719 (1.0.0) --> First development

Note(s):
https://pytesmo.readthedocs.io/en/latest/examples.html#the-pytesmo-validation-framework
TEST over Italy domain (complete) --> [1357, 1358, 1359, 1393, 1394, 1395, 1429, 1430, 1431]
TEST over Italy domain (partial) --> [1394, 1395]
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import argparse
import os
import sys
import time

from copy import deepcopy

from drv_process_main import DrvProcess

from lib_info_args import logger_name as logger_name_def
from lib_info_args import logger_format as logger_format_def
from lib_info_args import logger_file as logger_file_def

from lib_info_settings import get_data_settings

from lib_utils_generic import get_grid_reference, get_grid_cells
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'validation'
alg_type = 'Package'
alg_version = '1.0.1'
alg_release = '2023-07-19'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# algorithm main
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
    # get reference grid information
    alg_path_grid, alg_file_grid = get_grid_reference(alg_settings['datasets']['source'])

    # generate reference cells and gpis
    alg_cells, alg_gpis = get_grid_cells(
        cell_start=alg_settings['domain']['cell_start'],
        cell_end=alg_settings['domain']['cell_end'],
        cells_list=alg_settings['domain']['cell_list'],
        path_grid=alg_path_grid,
        file_grid=alg_file_grid)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    '''
    # set algorithm for debugging in sequential mode
    alg_cells = [1394, 1395] # italy case
    alg_cells = [31, 32]     # other case
    alg_settings['mode']['mp_mode'] = False
    '''

    # configure driver
    drv_process = DrvProcess(alg_cells, alg_settings)
    # setup process
    drv_process.setup_process()
    # execute process (sequential or parallel execution)
    drv_process.execute_process()
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
        logger_file = deepcopy(logger_file_def)

    if logger_folder is not None:
        logger_path = os.path.join(logger_folder, logger_file)
    else:
        logger_path = deepcopy(logger_file)

    logger_loc = os.path.split(logger_path)
    if logger_loc[0] == '' or logger_loc[0] == "":
        logger_folder_name, logger_file_name = os.path.dirname(os.path.abspath(sys.argv[0])), logger_loc[1]
    else:
        logger_folder_name, logger_file_name = logger_loc[0], logger_loc[1]; print('ea')
        
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
