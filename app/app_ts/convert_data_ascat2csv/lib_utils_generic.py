"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20250813'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import argparse
import logging
import json

from lib_utils_info import logger_name

# set logger obj
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to parse command-line arguments
def parse_args() -> argparse.Namespace:

    p = argparse.ArgumentParser(description="ASCAT daily index builder (JSON-driven).")
    p.add_argument("-settings_file", required=True, help="Path to JSON config.")
    p.add_argument("-time_now", help="YYYY-MM-DD; used only if config lacks time_start/time_end.")
    p.add_argument("-require-exist", action="store_true",
                   help="Emit rows only when TIFF, CSV, and registry exist.")
    p.add_argument("-stdout", action="store_true",
                   help="Write a single combined CSV to stdout, ignoring 'out' pattern.")
    return p.parse_args()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to setup configuration
def setup_config(file_path):
    if os.path.exists(file_path):
        with open(file_path) as file_handle:
            file_cfg = json.load(file_handle)
    else:
        logger_stream.error(' ===> Error in reading settings file "' + file_path + '"')
        raise IOError('File not found')
    return file_cfg
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to setup logging
def setup_logging(log_path):

    # manage log path
    log_folder, log_file = os.path.split(log_path)
    os.makedirs(log_folder, exist_ok=True)

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()  # Avoid duplicates if called multiple times

    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # File handler
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console (prompt) handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

# ----------------------------------------------------------------------------------------------------------------------
