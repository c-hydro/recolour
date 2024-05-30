"""
Library Features:

Name:          lib_info_settings
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import json

from lib_info_args import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to get data settings
def get_data_settings(file_name):
    if os.path.exists(file_name):
        with open(file_name) as file_handle:
            data_settings = json.load(file_handle)
    else:
        alg_logger.error(' ===> Error in reading settings file "' + file_name + '"')
        raise IOError('File not found')
    return data_settings
# -------------------------------------------------------------------------------------
