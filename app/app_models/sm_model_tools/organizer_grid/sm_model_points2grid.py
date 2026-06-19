#!/usr/bin/env python3
"""
APP - SM MODEL TOOLS - POINTS 2 GRID

__date__ = '20260615'
__version__ = '2.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'sm_model_tools'

General command line:
python sm_model_points2grid.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20260618 (1.0.0) --> Beta release for sm_model libraries
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations

import argparse
import logging
import sys
import time

from datetime import datetime

from lib_utils_logging import get_logger
from lib_utils_io import read_file_settings
from lib_utils_time import parse_reference_time

from lib_data import process, save

from config_info import LOGGER_NAME, ALG_NAME, ALG_RELEASE, ALG_VERSION

# set logger
logger_stream = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = ''
alg_name = 'Soil Moisture Models Tools - Points2Map'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2026-06-18'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main() -> None:

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()
    # get settings
    settings = read_file_settings(args.settings_file)

    # get reference time
    try:
        # get reference objects
        reference_time, reference_info = parse_reference_time(args.time_run)

        # optional hour override
        if args.hour is not None:
            reference_time = reference_time.replace(hour=args.hour,minute=0,second=0,microsecond=0)

            if reference_info is None:
                reference_info = {}

            reference_info["hour_override"] = f"{args.hour:02d}:00"
    except Exception as exc:
        print(f" ===> ERROR: parsing time: {exc}")
        sys.exit(1)

    # get logger
    get_logger(logger_stream, settings, reference_time=reference_time)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start message
    logger_stream.info(' ============================================================================ ')
    logger_stream.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger_stream.info(' ==> START ... ')
    logger_stream.info(' ')
    logger_stream.info(" ")
    logger_stream.info(f" ---> Settings file:      {args.settings_file}")
    logger_stream.info(f" ---> Source folder:      {settings.get('source', {}).get('folder')}")
    logger_stream.info(f" ---> Source filename:    {settings.get('source', {}).get('filename')}")
    logger_stream.info(f" ---> Destination folder: {settings.get('destination', {}).get('folder')}")
    logger_stream.info(f" ---> Destination file:   {settings.get('destination', {}).get('filename')}")

    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set time
    logger_stream.info(" ---> Set time ... ")
    try:
        logger_stream.info(f" ----> Get reference time: {reference_time.strftime('%Y-%m-%d %H:%M:%S')}")
        if args.time_run:
            logger_stream.info(f" ----> Set by user: {args.time_run} (apply time configuration)")
        else:
            logger_stream.info(" ----> Set by system (aqpply time configuration)")

        # print reference info
        if reference_info is not None:
            logger_stream.info(" ----> Time Info ...")
            for info_key, info_value in reference_info.items():
                if isinstance(info_value, datetime):
                    info_value = info_value.strftime('%Y-%m-%d %H:%M:%S')
                logger_stream.info(f" -----> {info_key}: {info_value}")
            logger_stream.info(" ----> Time Info ... DONE")

        logger_stream.info(" ---> Set time ... DONE")

    except Exception as exc:
        logger_stream.error(f" ===> ERROR: parsing time: {exc}")
        logger_stream.info(" ---> Set time ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # process datasets
    logger_stream.info(' ---> Process datasets ... ')
    try:
        grid_data, grid_mask, grid_profile = process(settings, reference_time)
        logger_stream.info(' ---> Process datasets ... DONE')

    except Exception as exc:
        logger_stream.error(f" ===> ERROR in processing datasets: {exc}")
        logger_stream.info(' ---> Process datasets ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # save datasets
    logger_stream.info(' ---> Save datasets ... ')
    try:
        save(grid_data, grid_mask, grid_profile, settings,
             time_reference=reference_time
        )
        logger_stream.info(' ---> Save datasets ... DONE')

    except Exception as exc:
        logger_stream.error(f" ===> ERROR in saving datasets: {exc}")
        logger_stream.info(' ---> Save datasets ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message
    alg_time_elapsed = round(time.time() - start_time, 1)

    logger_stream.info(' ')
    logger_stream.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger_stream.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    logger_stream.info(' ==> ... END')
    logger_stream.info(' ==> Bye, Bye')
    logger_stream.info(' ============================================================================ ')
    # --------------------------------------------------------------------------------------------------------------

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

    parser.add_argument(
        "-hour",
        type=int,
        default=None,
        choices=range(24),
        help="Optional hour override (0-23). If provided, the reference time is forced to HH:00."
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# wrapper algorithm main
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
