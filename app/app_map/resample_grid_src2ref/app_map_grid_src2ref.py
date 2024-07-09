#!/usr/bin/python3
"""
RECOLOUR APP - RESAMPLE GRID SOURCE TO REFERENCE  - REprocess paCkage for sOiL mOistUre pRoducts
__date__ = '20240708'
__version__ = '1.0.0'
__author__ = 'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'recolour'

General command line:
python3 app_map_grid_src2ref.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20240708 (1.0.0) --> Beta release
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import sys
import time
import argparse
import os

from copy import deepcopy

from lib_data_utils import apply_method_resample, apply_method_mask
from lib_data_tiff import read_file_tiff, write_file_tiff

from lib_utils_io import get_file_info, set_file_info
from lib_utils_time import set_time_info

from lib_info_settings import get_data_settings
from lib_info_args import logger_name, logger_format, time_format_algorithm

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for resampling grid data from source to reference'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2024-07-09'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# script main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get file settings
    alg_file_settings, alg_time_settings = get_args()
    # read data settings
    alg_data_settings = get_data_settings(alg_file_settings)

    # set logging
    set_logging(
        logger_name=logger_name, logger_format=logger_format,
        logger_folder=alg_data_settings['log']['folder_name'],
        logger_file=alg_data_settings['log']['file_name'])
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (start)
    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    # time algorithm
    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Organize time information
    alg_time_run = set_time_info(
        time_run_args=alg_time_settings,
        time_run_file=alg_data_settings['time']['time_run'],
        time_format=time_format_algorithm,
        time_frequency=alg_data_settings['time']['time_frequency'],
        time_rounding=alg_data_settings['time']['time_rounding'])
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set file reference
    file_path_ref_generic = get_file_info(alg_data_settings['datasets']['static'])
    file_path_ref_defined = set_file_info(
        file_path_ref_generic, template_time=alg_time_run, template_tags=alg_data_settings['template'])

    # set file source
    file_path_src_generic = get_file_info(alg_data_settings['datasets']['dynamic']['source'])
    file_path_src_defined = set_file_info(
        file_path_src_generic, template_time=alg_time_run, template_tags=alg_data_settings['template'])

    # set file destination
    file_path_dst_generic = get_file_info(alg_data_settings['datasets']['dynamic']['destination'])
    file_path_dst_defined = set_file_info(
        file_path_dst_generic, template_time=alg_time_run, template_tags=alg_data_settings['template'])
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # time info start
    alg_logger.info(' ---> Time Run: "' + str(alg_time_run) + '" ... ')

    # get reset flags
    reset_datasets = alg_data_settings['flags']['reset_datasets']
    # check reset flags
    if reset_datasets:
        if os.path.exists(file_path_dst_defined):
            os.remove(file_path_dst_defined)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # check destination file availability
    if not os.path.exists(file_path_dst_defined):

        # ------------------------------------------------------------------------------------------------------------------
        # get file reference
        alg_logger.info(' ----> Get reference file: "' + file_path_ref_defined + '" ... ')
        if os.path.exists(file_path_ref_defined):
            (file_data_ref, file_lons_ref, file_lats_ref,
             file_height_ref, file_width_ref, file_proj_ref, file_geotrans_ref) = read_file_tiff(file_path_ref_defined)
            alg_logger.info(' ----> Get reference file: "' + file_path_ref_defined + '" ... DONE')
        else:
            # message error file not found
            alg_logger.info(' ----> Get reference file: "' + file_path_ref_defined + '" ... FAILED')
            alg_logger.error(' ==> Reference file "' + file_path_ref_defined + '" does not exist')
            raise IOError('File not found')
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # get file source
        alg_logger.info(' ----> Get source file: "' + file_path_src_defined + '" ... ')
        if os.path.exists(file_path_src_defined):
            (file_data_src, file_lons_src, file_lats_src,
             file_height_src, file_width_src, file_proj_src, file_geotrans_src) = read_file_tiff(file_path_src_defined)
            alg_logger.info(' ----> Get source file: "' + file_path_src_defined + '" ... DONE')
        else:
            # message error file not found
            alg_logger.info(' ----> Get source file: "' + file_path_src_defined + '" ... FAILED')
            alg_logger.error(' ==> Source file "' + file_path_src_defined + '" does not exist')
            raise IOError('File not found')
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # resample grid source to reference
        alg_logger.info(' ----> Apply resampling source to reference ... ')
        file_data_tmp = apply_method_resample(
            file_lons_src, file_lats_src, file_data_src, file_lons_ref, file_lats_ref,
            **alg_data_settings['methods']['resample'])
        alg_logger.info(' ----> Apply resampling source to reference ... DONE')

        # mask grid source over reference
        alg_logger.info(' ----> Apply masking source to reference ... ')
        file_data_dst = apply_method_mask(
            file_data_tmp, file_data_ref,
            **alg_data_settings['methods']['mask'])
        alg_logger.info(' ----> Apply masking source to reference ... DONE')
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # save destination file
        alg_logger.info(' ----> Save destination file: "' + file_path_dst_defined + '" ... ')
        write_file_tiff(
            file_path_dst_defined, file_data_dst, file_width_ref, file_height_ref, file_geotrans_ref, file_proj_ref)
        alg_logger.info(' ----> Save destination file: "' + file_path_dst_defined + '" ... DONE')
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # time info end
        alg_logger.info(' ---> Time Run: "' + str(alg_time_run) + '" ... DONE')
        # ------------------------------------------------------------------------------------------------------------------
    else:
        # ------------------------------------------------------------------------------------------------------------------
        # time info end
        alg_logger.info(' ---> Time Run: "' + str(alg_time_run) + '" ... SKIPPED. File already processed')
        # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(' ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    alg_logger.info(' ==> ... END')
    alg_logger.info(' ==> Bye, Bye')
    alg_logger.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_name='algorithm_logger', logger_folder=None, logger_file='log.txt', logger_format=None):

    if logger_format is None:
        logger_format = deepcopy(logger_format)
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

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# call script from external library
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------------------------------
