#!/usr/bin/python3

"""
RECOLOUR TOOLS - SOILGRIDS IMG2CELL - REprocess paCkage for static soil datasets

__date__ = '20260521'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_img2cell_soilgrids.py -settings_file configuration.json

Version(s):
20260521 (1.0.0) --> First development for static SoilGrids datasets
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# libraries
import os
import logging
import time
from argparse import ArgumentParser

from lib_info_settings import get_data_settings, parse_data_settings, get_data_by_tag
from lib_utils_soil_grids import create_file_grid, select_origin_grid

from lib_reshuffle_soil_grids import main as main_runner
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'img2cell_soilgrids'
alg_type = 'Application'
alg_version = '1.0.0'
alg_release = '2026-05-21'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# script main
def main_wrapper():

    # -------------------------------------------------------------------------------------
    # method to get script argument(s)
    file_settings = get_args()

    # read data settings
    data_settings = get_data_settings(file_settings)

    # set logging
    log_settings = get_data_by_tag(
        data_settings,
        data_tag='log',
        data_default={
            'path': os.path.dirname(__file__),
            "file": "log.txt"
        }
    )

    set_logging(logger_file=os.path.join(log_settings['path'],log_settings['file']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # info algorithm (start)
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name +
                 ' (Version: ' + alg_version +
                 ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ')

    # time algorithm
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # parse data settings
    args_settings = parse_data_settings(data_settings)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # define origin grid
    data_ext = select_origin_grid(
        data_path=data_settings['data']['path_grid'],
        data_origin=data_settings['grid'].get("derived_by")
    )

    # create file grid
    create_file_grid(
        grid_path=os.path.join(
            data_settings['grid']['folder_name'],
            data_settings['grid']['file_name']
        ),
        data_path=data_settings['data']['path_grid'],
        data_ext=data_ext,
        grid_update=data_settings['flags']['reset_static']
    )
    # -------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # run procedure to convert grid to cell datasets
    main_runner(args_settings)
    # ----------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name +
                 ' (Version: ' + alg_version +
                 ' Release_Date: ' + alg_release + ')')

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

    parser_handle.add_argument(
        '-settings_file',
        action="store",
        dest="alg_settings"
    )

    parser_values = parser_handle.parse_args()

    alg_settings = 'configuration.json'

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings

    return alg_settings

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_file='./log.txt', logger_format=None):

    if logger_format is None:
        logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                        '%(filename)s:[%(lineno)-6s - %(funcName)20s()] %(message)s'

    logger_loc = os.path.split(logger_file)

    if logger_loc[0] == "":
        logger_folder_name, logger_file_name = (
            os.path.dirname(__file__),
            logger_loc[1]
        )
    else:
        logger_folder_name, logger_file_name = (
            logger_loc[0],
            logger_loc[1]
        )

    os.makedirs(logger_folder_name, exist_ok=True)

    # define logger path
    logger_path = os.path.join(
        logger_folder_name,
        logger_file_name
    )

    # Remove old logging file
    if os.path.exists(logger_path):
        os.remove(logger_path)

    # Set level of root debugger
    logging.root.setLevel(logging.DEBUG)

    # Open logging basic configuration
    logging.basicConfig(
        level=logging.DEBUG,
        format=logger_format,
        filename=logger_path,
        filemode='w'
    )

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
