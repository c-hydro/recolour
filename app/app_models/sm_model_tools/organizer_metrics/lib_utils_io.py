"""
Library Features:

Name:           lib_utils_io
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260527'
Version:        '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json

from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read settings file
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f'File "{file_name}" not found. Exit')
# ----------------------------------------------------------------------------------------------------------------------
