"""
Library Features:

Name:          lib_utils_system
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240209'
Version:       '2.0.8'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to make folder
def make_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
# ----------------------------------------------------------------------------------------------------------------------
