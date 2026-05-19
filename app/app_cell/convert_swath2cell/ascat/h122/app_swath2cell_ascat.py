#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SSM SWATH2CELL CONVERT - REprocess paCkage for sOiL mOistUre pRoducts

__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_ascat_swath2cell.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import time
import sys

from datetime import timedelta
from lib_ascat_swath import SwathGridFiles

from lib_utils_base import read_file_json, get_logger, get_args
from lib_utils_obj import map_dict_values
from lib_utils_time import parse_reference_time, compute_time_period
from lib_utils_report import collect_report, save_report
from config_info import LOGGER_NAME, ALG_NAME, ALG_RELEASE, ALG_VERSION

import matplotlib
matplotlib.use("Agg")

# set logger
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get arguments
    args = get_args()
    # get settings
    settings_algorithm = read_file_json(args.settings_file)

    # get reference time
    try:
        time_reference = parse_reference_time(args.time_run)
    except Exception as exc:
        print(f" ===> ERROR: parsing time: {exc}")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # resolve logging folder using reference time
    settings_logging = settings_algorithm.get("logging", {})

    tags_logging = {"time_step": time_reference, "time_start": time_reference,
                    "time_end": time_reference, "time_run": time_reference}

    if settings_logging.get("folder") is not None:
        settings_logging["folder"] = settings_logging.get("folder").format(**tags_logging)

    # update settings
    settings_algorithm["logging"] = settings_logging

    # get logger
    get_logger(settings_algorithm, reference_time=time_reference)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start message
    logger.info(' ============================================================================ ')
    logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger.info(' ==> START ... ')
    logger.info(' ')

    logger.info(f" ---> Settings file:      {args.settings_file}")
    logger.info(f" ---> Source folder:      {settings_algorithm.get('source', {}).get('folder')}")
    logger.info(f" ---> Source filename:    {settings_algorithm.get('source', {}).get('filename')}")
    logger.info(f" ---> Destination folder: {settings_algorithm.get('destination', {}).get('folder')}")
    logger.info(f" ---> Destination file:   {settings_algorithm.get('destination', {}).get('filename')}")
    logger.info(f" ---> Grid name:          {settings_algorithm.get('grid', {}).get('name')}")

    cells_list = settings_algorithm.get('parameters', {}).get('cells_list')
    logger.info(
        f" ---> Cells list:         {', '.join(map(str, cells_list)) if cells_list else 'None'}"
    )

    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # organize settings time
    settings_time = map_dict_values(
        settings_algorithm,
        dict_map={
            'time_start' : ['time', 'time_start'], 'time_end' : ['time', 'time_end'],
            'time_delta' : ['time', 'time_delta'], 'time_frequency' : ['time', 'time_frequency'],
            'floor_start' : ['time', 'floor_start_to_midnight'],
            'floor_end' : ['time', 'floor_end_to_midnight']
        },
        dict_default={
            'time_start': None, 'time_end': None, 'time_delta': '1D', 'time_frequency': 'H',
            'floor_start': True, 'floor_end': False,
        },
        dict_mandatory=None
    )

    # compute time period
    time_period = compute_time_period(time_run=time_reference, **settings_time)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # organize settings source
    settings_swath = map_dict_values(
        settings_algorithm,
        dict_map={'path_in' : ['source', 'folder'], 'filename_in' : ['source', 'filename'],
                  'path_out': ['destination', 'folder'], 'filename_out': ['destination', 'filename'],
                  'product_id': ['parameters', 'product'],
                  'mode': ['parameters', 'mode'],
                  "cells": ['parameters', 'cells_list'],
                  "cell_fn_format": ['parameters', 'cells_format'],
                  'max_nbytes_mb': ['parameters', 'max_nbytes_mb'],
                  "grid_name": ["grid", "name"]
                  },
        dict_default=None,
        dict_mandatory={
            'path_in' : True, 'filename_in' : True, 'product_id' : True,
            'path_out': True, 'filename_out': True, 'max_nbytes_mb': True,
            'cells_list': True, 'mode': True, 'grid_name': True
        },
    )
    settings_swath['dt_run'] = time_reference
    settings_swath['date_range'] = time_period
    settings_swath['print_progress'] = False
    settings_swath['parallel'] = False
    settings_swath['fmt_kwargs'] = {}
    settings_swath['fmt_kwargs']['dt_delta'] = timedelta(hours=1)

    # initialize swath class
    swath_collection = SwathGridFiles.from_product_id(**settings_swath)

    # convert swaths to cell files
    swath_collection.stack_to_cell_files(**settings_swath)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # organize report
    destination_folder_resolved = settings_swath.get("path_out")

    settings_report = {
        "time_period_start": time_period[0].strftime("%Y-%m-%d %H:%M:%S") if len(time_period) > 0 else None,
        "time_period_end": time_period[-1].strftime("%Y-%m-%d %H:%M:%S") if len(time_period) > 0 else None,
        "time_period_steps": len(time_period),
        "cells": settings_swath.get("cells"),
        "mode": settings_swath.get("mode"),
        "product_id": settings_swath.get("product_id")
    }

    # collect report info
    report = collect_report(
        settings=settings_algorithm,
        reference_time=time_reference,
        time_start=time_period[0],
        time_end=time_period[-1],
        destination_folder=destination_folder_resolved,
        stats=settings_report,
        start_time=start_time,
        written_objects=[]
    )

    # save report info
    save_report(settings_algorithm, report, time_reference)
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

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
