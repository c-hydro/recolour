"""
Library Features:

Name:           lib_utils_report
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json
import time

from lib_utils_basic import make_folder
from lib_utils_time import resolve_time_tags
from config_utils import LOGGER_NAME, ALG_NAME, ALG_VERSION, ALG_RELEASE

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to collect report
def collect_report(
        settings, reference_time, time_start, time_end, destination_file,
        stats, start_time):

    elapsed = round(time.time() - start_time, 1)

    report = {
        "algorithm": {
            "name": ALG_NAME,
            "version": ALG_VERSION,
            "release": ALG_RELEASE
        },
        "run": {
            "reference_time": reference_time.strftime("%Y-%m-%d %H:%M:%S"),
            "time_start": time_start.strftime("%Y-%m-%d %H:%M:%S"),
            "time_end": time_end.strftime("%Y-%m-%d %H:%M:%S"),
            "elapsed_seconds": elapsed
        },
        "parameters": settings.get("parameters", {}),
        "source": settings.get("source", {}),
        "destination": {
            **settings.get("destination", {}),
            "resolved_file": destination_file
        },
        "time": settings.get("time", {}),
        "stats": stats
    }

    return report
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to write report
def write_report(report, file_path):
    folder_name = os.path.dirname(file_path)
    if folder_name:
        make_folder(folder_name)

    with open(file_path, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=4)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to save report
def save_report(path, obj):
    write_report(obj, path)
# ----------------------------------------------------------------------------------------------------------------------