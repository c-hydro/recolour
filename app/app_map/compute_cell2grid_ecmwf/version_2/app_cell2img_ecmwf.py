#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - ECMWF RZEM LAYERS/PROFILE CELL2GRID

General command line:
python app_cell2grid_ecmwf.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import time
import sys
import argparse

from datetime import datetime

from lib_utils_logging import get_logger
from lib_utils_io import read_file_json
from lib_utils_time import parse_reference_time
from lib_data import process, save
from config_info import LOGGER_NAME, ALG_NAME, ALG_RELEASE, ALG_VERSION

# set logger
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()
    # get settings
    settings = read_file_json(args.settings_file)

    # get reference time
    try:
        reference_time, reference_info = parse_reference_time(args.time_run)
    except Exception as exc:
        print(f" ===> ERROR: parsing time: {exc}")
        sys.exit(1)

    # get logger
    get_logger(logger, settings, reference_time=reference_time)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start message
    logger.info(' ============================================================================ ')
    logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger.info(' ==> START ... ')
    logger.info(' ')

    logger.info(f" ---> Settings file:      {args.settings_file}")
    logger.info(f" ---> Source folder:      {settings.get('source', {}).get('folder')}")
    logger.info(f" ---> Source filename:    {settings.get('source', {}).get('filename')}")
    logger.info(f" ---> Destination folder: {settings.get('destination', {}).get('folder')}")
    logger.info(f" ---> Destination file:   {settings.get('destination', {}).get('filename')}")

    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set time
    logger.info(" ---> Set time ... ")
    try:
        logger.info(f" ----> Get reference time: {reference_time.strftime('%Y-%m-%d %H:%M:%S')}")
        if args.time_run:
            logger.info(f" ----> Set by user: {args.time_run} (apply time configuration)")
        else:
            logger.info(" ----> Set by system (aqpply time configuration)")

        # print reference info
        if reference_info is not None:
            logger.info(" ----> Time Info ...")
            for info_key, info_value in reference_info.items():
                if isinstance(info_value, datetime):
                    info_value = info_value.strftime('%Y-%m-%d %H:%M:%S')
                logger.info(f" -----> {info_key}: {info_value}")
            logger.info(" ----> Time Info ... DONE")

        logger.info(" ---> Set time ... DONE")

    except Exception as exc:
        logger.error(f" ===> ERROR: parsing time: {exc}")
        logger.info(" ---> Set time ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # process datasets
    logger.info(' ---> Process datasets ... ')
    try:
        (soil_moisture_map_smooth, soil_moisture_map_interp, time_lag_map,
         profile, stats, time_start, time_end) = process(settings, reference_time)
        logger.info(' ---> Process datasets ... DONE')

    except Exception as exc:
        logger.error(f" ===> ERROR in processing datasets: {exc}")
        logger.info(' ---> Process datasets ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # save datasets
    logger.info(' ---> Save datasets ... ')
    try:
        save(
            soil_moisture_map_smooth, soil_moisture_map_interp,
            time_lag_map,
            profile,
            stats,
            time_start,
            time_end,
            reference_time,
            settings,
            start_time
        )
        logger.info(' ---> Save datasets ... DONE')

    except Exception as exc:
        logger.error(f" ===> ERROR in saving datasets: {exc}")
        logger.info(' ---> Save datasets ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message
    alg_time_elapsed = round(time.time() - start_time, 1)

    logger.info(' ')
    logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    logger.info(' ==> ... END')
    logger.info(' ==> Bye, Bye')
    logger.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# cli
def get_args():

    parser = argparse.ArgumentParser(
        description="Convert SMAP SPL2SMP_E cell files to raster maps"
    )

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        required=True,
        help="Path to JSON settings file"
    )

    parser.add_argument(
        "-time",
        dest="time_run",
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM". Minutes are rounded down to 00.'
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------