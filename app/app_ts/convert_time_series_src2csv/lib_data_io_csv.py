"""
Library Features:

Name:          lib_data_io_csv
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd

from lib_utils_obj import map_vars_dframe, sanitize_string
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to wrap registry in csv format
def wrap_registry_csv(file_name, file_fields, file_sep=',', file_decimal='.'):

    # get file fields
    fields_data_raw = pd.read_table(file_name, sep=file_sep, decimal=file_decimal)
    fields_data_raw.columns = fields_data_raw.columns.str.strip()
    # map file fields
    fields_data_map = map_vars_dframe(fields_data_raw, file_fields)

    # get name fields
    var_data_name = fields_data_map['name'].values
    # parser tag
    if 'tag' not in list(fields_data_map.keys()):
        var_data_tag = []
        for string_name in var_data_name:
            string_tag = sanitize_string(string_name)
            var_data_tag.append(string_tag)

        fields_data_map['tag'] = var_data_tag

    # create fields dframe
    fields_dframe = pd.DataFrame(data=fields_data_map)

    if 'valid' in list(fields_dframe.columns):
        fields_dframe = fields_dframe[fields_dframe['valid'] == 1]

    return fields_dframe
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write file csv
def write_file_csv(file_name, file_dframe,
                   dframe_sep=',', dframe_decimal='.', dframe_float_format='%.2f',
                   dframe_index=False, dframe_header=True,
                   dframe_index_label='time', dframe_index_format='%Y-%m-%d %H:%M',
                   dframe_no_data=-9999):

    if np.isfinite(dframe_no_data):
        file_dframe.fillna(dframe_no_data, inplace=True)

    if dframe_index_format is not None:
        file_dframe.index = file_dframe.index.strftime(dframe_index_format)

    file_dframe.to_csv(
        file_name, mode='w',
        index=dframe_index, sep=dframe_sep, decimal=dframe_decimal,
        index_label=dframe_index_label,
        header=dframe_header, float_format=dframe_float_format,  quotechar='"')

# ----------------------------------------------------------------------------------------------------------------------
