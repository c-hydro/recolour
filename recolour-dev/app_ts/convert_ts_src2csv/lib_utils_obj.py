"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""


# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import unicodedata
import re
from copy import deepcopy

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to sanitize string
def sanitize_string(string_name):

    string_name = string_name.lower()
    string_name = re.sub(r"['.,-]", "", string_name)
    string_name = string_name.replace(' ', '')
    string_name = unicodedata.normalize('NFD', string_name).encode('ascii', 'ignore')
    string_name = string_name.decode("utf-8")

    return string_name
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to map variable(s) in a dataframe object
def map_vars_dframe(var_dframe_in, var_map=None):

    if var_map is not None:
        var_name_excluded = None
        for var_name_in, var_name_out in var_map.items():
            if var_name_in not in list(var_dframe_in.columns):
                if var_name_excluded is None:
                    var_name_excluded = []
                var_name_excluded.append(var_name_in)
                log_stream.warning(' ===> Variable "' + var_name_in + '" not included in the dataframe obj')

        if var_name_excluded is not None:
            for var_name_tmp in var_name_excluded:
                var_map.pop(var_name_tmp)

        var_dframe_out = var_dframe_in.rename(columns=var_map)

    else:
        var_dframe_out = deepcopy(var_dframe_in)

    return var_dframe_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to map variable(s) in a dictionary object
def map_vars_dict(var_data_in, var_map=None):

    if var_map is not None:
        var_data_out = {}
        for var_name_out, var_name_in in var_map.items():
            if var_name_in in list(var_data_in.keys()):
                values_tmp = var_data_in[var_name_in]
                var_data_out[var_name_out] = values_tmp
            else:
                log_stream.warning(' ===> Variable "' + var_name_in + '" is not included in the source obj')
    else:
        var_data_out = deepcopy(var_data_in)

    return var_data_out
# ----------------------------------------------------------------------------------------------------------------------
