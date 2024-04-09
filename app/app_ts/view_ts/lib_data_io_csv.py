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
import pandas as pd
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file csv
def read_file_csv(file_name, dframe_index='time', dframe_date_format='%Y-%m-%d %H:%M',
                  dframe_sep=',', dframe_decimal='.', dframe_float_precision='legacy'):

    file_dframe = pd.read_csv(file_name,
                              decimal=dframe_decimal, sep=dframe_sep, float_precision=dframe_float_precision)

    if dframe_index not in list(file_dframe.columns):
        log_stream.error(' ===> Index column "' + dframe_index +
                         '"  must be available in the source dataframe. Check the source file')
        raise RuntimeError('Including the index column in the source file for skipping this error.')

    if dframe_index == 'time':
        file_dframe[dframe_index] = pd.DatetimeIndex(file_dframe[dframe_index].values).strftime(dframe_date_format)
        file_dframe[dframe_index] = pd.DatetimeIndex(file_dframe[dframe_index])

    file_dframe = file_dframe.reset_index()
    if 'index' in list(file_dframe.columns):
        file_dframe = file_dframe.drop(['index'], axis=1)
    file_dframe = file_dframe.set_index(dframe_index)
    file_dframe.sort_index(inplace=True)

    return file_dframe

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write file csv
def write_file_csv(file_name, file_dframe,
                   dframe_sep=',', dframe_decimal='.', dframe_float_format='%.2f',
                   dframe_index=False, dframe_header=True,
                   dframe_index_label='time'):

    file_dframe.to_csv(
        file_name, mode='w',
        index=dframe_index, sep=dframe_sep, decimal=dframe_decimal,
        index_label=dframe_index_label,
        header=dframe_header, float_format=dframe_float_format,  quotechar='"')

# ----------------------------------------------------------------------------------------------------------------------
