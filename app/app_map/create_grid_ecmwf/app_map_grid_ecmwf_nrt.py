#!/usr/bin/python3

"""
RECOLOUR APP - ECMWF GRID PRODUCT  - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20260429'
__version__ = '2.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_map_grid_ecmwf_nrt.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20260429 (2.0.0) --> Update application for new computing profile methods based on dynamic weighted layers
20230830 (1.0.0) --> First development
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import time

try:
    import matplotlib
    matplotlib.use("Agg")
except ImportError:
    warnings.warn(
        "matplotlib is not installed. Plotting features will be disabled.",
        ImportWarning
    )
import argparse

from lib_utils_time import set_time_info

from lib_info_args import time_format_algorithm as time_format
from lib_info_args import logger_name as logger_name_def
from lib_info_args import logger_format as logger_format_def
from lib_info_settings import get_data_settings
from lib_utils_logging import set_logging

from drv_data_static import DrvData as DrvDataStatic
from drv_data_dynamic import DrvData as DrvDataDynamic
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for computing grid ecmwf product'
alg_type = 'Package'
alg_version = '2.0.0'
alg_release = '2026-04-29'
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
    set_logging(logger_name=logger_name_def, logger_format=logger_format_def,
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
    alg_time_now, alg_time_range, alg_time_chunks = set_time_info(
        time_run_args=alg_time_settings, time_run_file=alg_data_settings['time']['time_now'],
        time_run_file_start=alg_data_settings['time']['time_start'],
        time_run_file_end=alg_data_settings['time']['time_end'],
        time_format=time_format,
        time_period=alg_data_settings['time']['time_period'],
        time_frequency=alg_data_settings['time']['time_frequency'],
        time_rounding=alg_data_settings['time']['time_rounding'], time_reverse=True)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # configure static driver
    drv_data_static = DrvDataStatic(alg_data_settings)
    # organize static datasets
    alg_data_static = drv_data_static.organize_data()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # iterate over time information
    for alg_time_reference, alg_time_datasets in alg_time_chunks.items():

        # -------------------------------------------------------------------------------------
        # configure dynamic driver
        drv_data_dynamic = DrvDataDynamic(
            alg_time_reference, alg_time_datasets, alg_data_static, alg_data_settings)
        # organize datasets
        drv_data_dynamic.organize_data()
        # compute data resampling
        drv_data_dynamic.analyze_data()
        # dump datasets
        drv_data_dynamic.dump_data()
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

    # parser algorithm arg(s)
    parser_obj = argparse.ArgumentParser()
    parser_obj.add_argument('-settings_file', action="store", dest="settings_file")
    parser_obj.add_argument('-time', action="store", dest="settings_time")
    parser_value = parser_obj.parse_args()

    # set algorithm arg(s)
    settings_file, settings_time = 'configuration.json', None
    if parser_value.settings_file:
        settings_file = parser_value.settings_file
    if parser_value.settings_time:
        settings_time = parser_value.settings_time

    return settings_file, settings_time

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# call script from external library
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------------------------
