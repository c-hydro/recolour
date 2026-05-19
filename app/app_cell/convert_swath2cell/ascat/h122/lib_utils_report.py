"""
Library Features:

Name:          lib_utils_report
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260518'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json
import time

from config_info import ALG_NAME, ALG_VERSION, ALG_RELEASE, LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# report utils
def collect_report(settings, reference_time, time_start, time_end,
                   destination_folder, stats, start_time, written_objects=None):

    if written_objects is None:
        written_objects = []

    elapsed = round(time.time() - start_time, 1)

    return {
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
        "grid": settings.get("grid", {}),
        "source": settings.get("source", {}),
        "destination": {
            **settings.get("destination", {}),
            "resolved_folder": destination_folder
        },
        "time": settings.get("time", {}),
        "stats": stats,
        "written_objects": written_objects
    }


def write_report(report, file_path):
    folder_name = os.path.dirname(file_path)
    if folder_name:
        os.makedirs(folder_name, exist_ok=True)

    with open(file_path, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=4)


def save_report(settings, report, reference_time):
    report_settings = settings.get("report", {})
    if not report_settings.get("enabled", False):
        return

    logger.info(" ----> Save report ... ")

    tags = {
        "time_step": reference_time,
        "time_start": reference_time,
        "time_end": reference_time,
        "time_run": reference_time
    }

    report_folder_template = report_settings.get("folder", "./report")
    report_filename_template = report_settings.get(
        "filename", "swath2cell_h122_report_%Y%m%d_%H%M.json"
    )

    report_folder = report_folder_template.format(**tags)
    report_filename = reference_time.strftime(report_filename_template)
    report_path = os.path.join(report_folder, report_filename)

    write_report(report, report_path)
    logger.info(" ::: Report: " + report_path)
    logger.info(" ----> Save report ... DONE")
# ----------------------------------------------------------------------------------------------------------------------