"""
Library Features:

Name:          lib_data_io_ascii
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd

from copy import deepcopy

from lib_utils_generic import fill_tags2string
from lib_utils_obj import map_vars_dframe, sanitize_string, fill_tags_time
from lib_utils_time import replace_time_part
from lib_info_args import logger_name, time_format_algorithm

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to wrap parameters in ascii format
def read_parameters_ascii(file_name, file_column='param_{:}'):

    with open(file_name, 'r') as file:
        file_content = file.readlines()
    file_content = [x.strip() for x in file_content]
    file_content = [float(x) for x in file_content if x]

    file_object = {}
    for i, x in enumerate(file_content):
        file_object[file_column.format(i + 1)] = x

    return file_object
# ----------------------------------------------------------------------------------------------------------------------
