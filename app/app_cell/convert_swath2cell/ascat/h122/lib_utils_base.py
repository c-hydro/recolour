"""
Library Features:

Name:          lib_utils_base
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260518'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import argparse
import sys
import os
import json

from lib_utils_time import resolve_time_tags
from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read json file
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f'File "{file_name}" not found. Exit')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to configure logger
def get_logger(settings, reference_time=None):

    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    log_folder_template = log_settings.get("folder", "./log")
    log_filename = log_settings.get("filename", "swath2cell_h122.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    if reference_time is None:
        reference_time = datetime.now()

    log_folder = resolve_time_tags(
        log_folder_template,
        {
            "time_step": reference_time,
            "time_start": reference_time,
            "time_end": reference_time,
            "time_run": reference_time
        }
    )

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = reference_time.strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    logger.handlers = []
    logger.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to get arguments
def get_args():
    parser = argparse.ArgumentParser(
        description="Convert H122 swath files to gpi/cell outputs"
    )

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        required=True,
        help="Path to JSON settings file",
    )

    parser.add_argument(
        "-time",
        dest="time_run",
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM". Minutes are rounded down to 00.',
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------
