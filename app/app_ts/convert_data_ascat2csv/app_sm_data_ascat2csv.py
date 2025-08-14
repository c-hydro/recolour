#!/usr/bin/env python3

"""
RECOLOUR APPS - EXTRACT DATA FROM ASCAT POINTS OR GRID TO TIME-SERIES - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20250812'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_map_grid_smr.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20250812 (1.0.0) --> First development
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations
import logging
import os
import time

from lib_utils_generic import parse_args, setup_config, setup_logging
from lib_utils_time import compute_time_window, compute_time_ranges

from drv_geo import DrvGeo
from drv_data_src import DrvData as DrvDataSrc
from drv_data_dst import DrvData as DrvDataDst

from lib_utils_info import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for extracting ASCAT data (converted to grid or points) to CSV'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2025-08-12'
# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# main function
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # parse command-line arguments
    alg_args = parse_args()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # load config from JSON file
    alg_settings = setup_config(alg_args.settings_file)
    # define logging
    path_log = os.path.join(alg_settings["log"]['folder_name'], alg_settings["log"]['file_name'])
    setup_logging(path_log)
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
    # define time window
    alg_time_start, alg_time_end = compute_time_window(
        time_now=alg_args.time_now,
        time_start=alg_settings['time']['time_start'], time_end=alg_settings['time']['time_end'])
    # define time chunks
    alg_time_chunks = compute_time_ranges(alg_time_start, alg_time_end, freq='MS')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # geographical class
    drv_geo = DrvGeo(
        folder_name=alg_settings['geo']['registry']['folder_name'],
        file_name=alg_settings['geo']['registry']['file_name'])
    # geographical data
    obj_data_geo = drv_geo.organize_data()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # iterate over time chunks
    for alg_time_period in alg_time_chunks:

        # class to set data source csv
        drv_data_csv = DrvDataSrc(
            folder_name=alg_settings['source']['data_csv']['folder_name'],
            file_name=alg_settings['source']['data_csv']['file_name'],
            file_variable=alg_settings['source']['data_csv']['variable'],
            file_delimiter=alg_settings['source']['data_csv']['delimiter'],
            file_index=alg_settings['source']['data_csv']['index'],
            file_format='csv',
            time_start=alg_time_period[0], time_end=alg_time_period[1], frequency='D',
            search_radius_km=alg_settings['parameters']['search_radius_km'],
            time_window=alg_settings['source']['data_windows'])

        # method to get data
        obj_data_csv = drv_data_csv.read_data()
        # method to organize data
        obj_data_csv = drv_data_csv.organize_data_ts(
            file_dfs=obj_data_csv, registry_df=obj_data_geo)
        # method to sync data
        obj_data_csv = drv_data_csv.sync_data(obj_data_csv)

        # class to set data source tiff
        drv_data_tiff = DrvDataSrc(
            folder_name=alg_settings['source']['data_tiff']['folder_name'],
            file_name=alg_settings['source']['data_tiff']['file_name'],
            file_variable=alg_settings['source']['data_tiff']['variable'],
            file_delimiter=alg_settings['source']['data_tiff']['delimiter'],
            file_index=alg_settings['source']['data_tiff']['index'],
            file_format='tiff',
            time_start=alg_time_period[0], time_end=alg_time_period[1], frequency='D',
            search_radius_km=alg_settings['parameters']['search_radius_km'],
            time_window=alg_settings['source']['data_windows'])

        # method to get data
        obj_data_tiff = drv_data_tiff.read_data()
        # method to organize data
        obj_data_tiff = drv_data_tiff.organize_data_maps(
            file_das=obj_data_tiff, registry_df=obj_data_geo)
        # method to sync data
        obj_data_tiff = drv_data_tiff.sync_data(obj_data_tiff)

        # class to set data destination
        drv_data_dst = DrvDataDst(
            folder_name=alg_settings['destination']['folder_name'],
            file_name=alg_settings['destination']['file_name'],
            file_variable=alg_settings['destination']['variable'],
            file_delimiter=alg_settings['destination']['delimiter'],
            file_index=alg_settings['destination']['index'],
            file_format='csv',
            time_start=alg_time_period[0], time_end=alg_time_period[1], frequency='MS',
            mapping=alg_settings['destination']['mapping'])

        # method to organize data destination
        obj_data_dst = drv_data_dst.organize_data(obj_data_csv, obj_data_tiff)

        # method to dump data destination
        drv_data_dst.dump_data(obj_data_dst, obj_data_geo)

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
# method to configure main function
if __name__ == "__main__":
    raise SystemExit(main())
# ----------------------------------------------------------------------------------------------------------------------
