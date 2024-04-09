"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20190117'
Version:       '2.0.8'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np

from copy import deepcopy
from random import randint
from netCDF4 import Dataset
from datetime import datetime

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to convert list 2 dictionary
def convert_list_2_dict(var_list, var_fill=None, var_split={'split': False, 'chunk': 2, 'key': 0}):

    var_dict = {}
    for step, var in enumerate(var_list):

        if var_split['split'] is True:

            var_splitted = [var[i:i + var_split['chunk']] for i in range(0, len(var), var_split['chunk'])][0]

            var_keys = deepcopy(var_splitted)
            var_values = deepcopy(var_splitted)
            var_index = int(var_split['key'])
            var_key = var_keys[var_index]
            var_values.pop(var_index)

            var_dict[var_key] = var_values[0]

        else:
            if var_fill is not None:
                var_dict[var_fill] = var
            else:
                var_dict[step] = var

    return var_dict
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Method to add time in a unfilled string (path or filename)
def fill_tags2string(string_raw, tags_format=None, tags_filling=None):

    apply_tags = False
    if string_raw is not None:
        for tag in list(tags_format.keys()):
            if tag in string_raw:
                apply_tags = True
                break

    if apply_tags:
        string_filled = string_raw.format(**tags_format)

        for tag_format_name, tag_format_value in list(tags_format.items()):

            if tag_format_name in list(tags_filling.keys()):
                tag_filling_value = tags_filling[tag_format_name]
                if tag_filling_value is not None:

                    if isinstance(tag_filling_value, datetime):
                        tag_filling_value = tag_filling_value.strftime(tag_format_value)

                    string_filled = string_filled.replace(tag_format_value, tag_filling_value)

        string_filled = string_filled.replace('//', '/')
        return string_filled
    else:
        return string_raw
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to check time-series filename availability
def check_filename(path_ts, path_grid, filename_ts='%04d.nc',
                   filename_grid='TUW_WARP5_grid_info_2_3.nc', var_cell='cell'):

    dset_grid = Dataset(os.path.join(path_grid, filename_grid), 'r')
    cells = np.unique(dset_grid.variables[var_cell][:])

    n = cells.__len__()

    file_available = np.ones(n, dtype=np.bool)
    file_available[:] = False
    for i, cell in enumerate(cells):

        filename_ts_def = filename_ts % cell
        file_ts = os.path.join(path_ts, filename_ts_def)
        if os.path.exists(file_ts):
            file_available[i] = True

    if np.any(file_available == False):
        return False
    else:
        return True
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to create a random string
def random_string(string_root='temporary', string_sep='_', n_rand_min=0, n_rand_max=1000):

    # Generate random number
    n_rand = str(randint(n_rand_min, n_rand_max))
    # Generate time now
    time_now = datetime.now().strftime('%Y%m%d-%H%M%S_%f')

    # Concatenate string(s) with defined separator
    string_rand = string_sep.join([string_root, time_now, n_rand])

    return string_rand
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to update recursively dictionary value
def get_dict_values(d, key, value=[]):

    for k, v in iter(d.items()):

        if isinstance(v, dict):
            if k == key:

                for kk, vv in iter(v.items()):
                    temp = [kk, vv]
                    value.append(temp)

            else:
                vf = get_dict_values(v, key, value)

                if isinstance(vf, list):
                    if vf:
                        vf_end = vf[0]
                    else:
                        vf_end = None

                elif isinstance(vf, np.ndarray):
                    vf_end = vf.tolist()
                else:
                    vf_end = vf

                if vf_end not in value:
                    if vf_end:

                        if isinstance(value, list):
                            value.append(vf_end)
                        elif isinstance(value, str):
                            value = [value, vf_end]

                    else:
                        pass
                else:
                    pass

        else:
            if k == key:

                if isinstance(v, np.ndarray):
                    value = v
                else:
                    value = v
    return value
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to merge dictionaries
def merge_dict(d1, d2):
    dd = {}
    for d in (d1, d2):  # you can list as many input dicts as you want here
        for key, value in iter(d.items()):
            dd[key] = value
    return dd
# ----------------------------------------------------------------------------------------------------------------------
