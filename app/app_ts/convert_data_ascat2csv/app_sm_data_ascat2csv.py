
#!/usr/bin/env python3

# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations
import logging
import os
import csv, sys

from pathlib import Path
from typing import Dict, Optional

from lib_utils_generic import setup_config, setup_logging
from lib_utils_default import parse_args
from lib_utils_time import compute_time_window, compute_time_ranges

from driver_registry import DriverRegistry
from driver_data_src import DriverData as DriverDataSrc
from driver_data_dst import DriverData as DriverDataDst

from lib_utils_info import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

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
    # define time window
    alg_time_start, alg_time_end = compute_time_window(
        time_now=alg_args.time_now,
        time_start=alg_settings['time']['time_start'], time_end=alg_settings['time']['time_end'])
    # define time chunks
    alg_time_chunks = compute_time_ranges(alg_time_start, alg_time_end, freq='MS')
    # ------------------------------------------------------------------------------------------------------------------

    alg_logger.info("Starting index build")
    alg_logger.debug(f"Time window: {alg_time_start} -> {alg_time_end}")

    # ------------------------------------------------------------------------------------------------------------------
    # registry class
    driver_registry = DriverRegistry(
        folder_name=alg_settings['auxiliary']['registry']['folder_name'],
        file_name=alg_settings['auxiliary']['registry']['file_name'])
    # registry data
    obj_data_registry = driver_registry.read_data()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # iterate over time chunks
    for alg_time_period in alg_time_chunks:

        # class to set data source csv
        driver_data_csv = DriverDataSrc(
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
        obj_data_csv = driver_data_csv.read_data()
        # method to organize data
        obj_data_csv = driver_data_csv.organize_data_ts(file_dfs=obj_data_csv, registry_df=obj_data_registry)
        # method to sync data
        obj_data_csv = driver_data_csv.sync_data(obj_data_csv)

        # class to set data source tiff
        driver_data_tiff = DriverDataSrc(
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
        obj_data_tiff = driver_data_tiff.read_data()
        # method to organize data
        obj_data_tiff = driver_data_tiff.organize_data_maps(
            file_das=obj_data_tiff, registry_df=obj_data_registry)
        # method to sync data
        obj_data_tiff = driver_data_csv.sync_data(obj_data_tiff)

        # class to set data destination
        driver_data_dst = DriverDataDst(
            folder_name=alg_settings['destination']['folder_name'],
            file_name=alg_settings['destination']['file_name'],
            file_variable=alg_settings['destination']['variable'],
            file_delimiter=alg_settings['destination']['delimiter'],
            file_index=alg_settings['destination']['index'],
            file_format='csv',
            time_start=alg_time_period[0], time_end=alg_time_period[1], frequency='MS')

        # method to organize data destination
        obj_data_dst = driver_data_dst.organize_data(obj_data_csv, obj_data_tiff)

        # method to dump data destination
        driver_data_dst.dump_data(obj_data_dst, obj_data_registry)

    # ------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to configure main function
if __name__ == "__main__":
    raise SystemExit(main())
# ----------------------------------------------------------------------------------------------------------------------
