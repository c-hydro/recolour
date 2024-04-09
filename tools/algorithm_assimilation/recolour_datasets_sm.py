#!/usr/bin/python3

"""
DTE HYDROLOGY - PREPROCESSING TOOL FOR SOIL MOISTURE VARIABLE(S)

__date__ = '20210514'
__version__ = '1.0.0'
__author__ = 'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'DTE'

General command line:
python dte_preprocessing_datasets_sm.py -settings_file configuration.json

Version:
20210312 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
import os
import json

from argparse import ArgumentParser

from lib_utils_system import make_folder
from lib_utils_time import set_time
from lib_info_args import time_format_algorithm

from driver_data_io_static import DriverStatic
from driver_data_io_dynamic import DriverDynamic
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
project_name = 'DTE'
alg_name = 'TOOL FOR SOIL MOISTURE VARIABLES'
alg_type = 'PREPROCESSING'
alg_version = '1.0.0'
alg_release = '2021-03-12'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # Read arguments and settings file
    file_name_args, time_args = get_args()
    file_data_settings = read_file_settings(file_name_args)

    # Set logging
    make_folder(file_data_settings['log']['folder_name'])
    set_logging(logger_file=os.path.join(file_data_settings['log']['folder_name'],
                                         file_data_settings['log']['file_name']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info('[' + project_name + ' ' + alg_type + ' - ' + alg_name + ' (Version ' + alg_version +
                 ' - Release ' + alg_release + ')]')
    logging.info(' ==> START ... ')
    logging.info(' ')

    # Time algorithm information
    alg_time_start = time.time()

    # Organize time run
    time_run, time_range, time_chunks = set_time(
        time_run_args=time_args,
        time_run_file=file_data_settings['time']['time_run'],
        time_run_file_start=file_data_settings['time']['time_start'],
        time_run_file_end=file_data_settings['time']['time_end'],
        time_format=time_format_algorithm,
        time_period=file_data_settings['time']['time_period'],
        time_frequency=file_data_settings['time']['time_frequency'],
        time_rounding=file_data_settings['time']['time_rounding']
    )

    # Organize static datasets
    driver_data_static = DriverStatic(
        src_dict=file_data_settings['data']['static']['source'],
        dst_dict=file_data_settings['data']['static']['destination'],
        parameters_dict=file_data_settings['parameters'],
        flag_cleaning_static=file_data_settings['flags']['cleaning_static_destination'])
    static_data_collections = driver_data_static.organize_static()

    # Organize dynamic datasets
    driver_data_dynamic = DriverDynamic(
        static_data_collection=static_data_collections,
        src_dict=file_data_settings['data']['dynamic']['source'],
        anc_dict=file_data_settings['data']['dynamic']['ancillary'],
        anl_dict=file_data_settings['data']['dynamic']['analysis'],
        dest_dict=file_data_settings['data']['dynamic']['destination'],
        parameters_dict=file_data_settings['parameters'],
        tmp_dict=file_data_settings['tmp'],
        alg_template_tags=file_data_settings['template'],
        flag_active_variables=file_data_settings['flags']['active_variables'],
        flag_cleaning_ancillary_model_collections=file_data_settings['flags']['cleaning_dynamic_ancillary_model_collections'],
        flag_cleaning_ancillary_dset_obj=file_data_settings['flags']['cleaning_dynamic_ancillary_datasets_obj'],
        flag_cleaning_ancillary_dset_ref=file_data_settings['flags']['cleaning_dynamic_ancillary_datasets_reference'],
        flag_cleaning_ancillary_dset_grp=file_data_settings['flags']['cleaning_dynamic_ancillary_datasets_group'],
        flag_cleaning_ancillary_dset_collections=file_data_settings['flags']['cleaning_dynamic_ancillary_datasets_collections'],
        flag_cleaning_analysis=file_data_settings['flags']['cleaning_dynamic_analysis'],
        flag_cleaning_destination=file_data_settings['flags']['cleaning_dynamic_destination'],
    )
    driver_data_dynamic.organize_dynamic_data_obs()
    driver_data_dynamic.organize_dynamic_data_model()

    driver_data_dynamic.analyze_dynamic_data()

    driver_data_dynamic.filter_dynamic_stats()
    driver_data_dynamic.dump_dynamic_data()

    # Info algorithm
    alg_time_elapsed = round(time.time() - alg_time_start, 1)

    logging.info(' ')
    logging.info('[' + project_name + ' ' + alg_type + ' - ' + alg_name + ' (Version ' + alg_version +
                 ' - Release ' + alg_release + ')]')
    logging.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read file settings
def read_file_settings(file_name):
    if os.path.exists(file_name):
        with open(file_name) as file_handle:
            file_data = json.load(file_handle)
    else:
        logging.error(' ===> Error in reading settings file ' + file_name)
        raise IOError('File not found')
    return file_data
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-time', action="store", dest="alg_time")
    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    if parser_values.alg_time:
        alg_time = parser_values.alg_time
    else:
        alg_time = None

    return alg_settings, alg_time

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to set logging information
def set_logging(logger_file='log.txt', logger_format=None):
    if logger_format is None:
        logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                        '%(filename)s:[%(lineno)-6s - %(funcName)20s()] %(message)s'

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
# Call script from external library
if __name__ == '__main__':
    main()
# -------------------------------------------------------------------------------------
