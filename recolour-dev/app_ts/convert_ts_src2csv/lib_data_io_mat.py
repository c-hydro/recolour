"""
Library Features:

Name:          lib_data_io_mat
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'

NOTE:
    matlab file .mat should not have string.java format
    convert string.java to cell;
    >> cell_var = {cell(string_var)}
    save the file:
    >> save file_name.mat -v7
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import scipy.io
import os
import pandas as pd

from datetime import datetime
from copy import deepcopy

from hmc.lib_data_io_pickle import read_obj, write_obj
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file mat
def read_file_mat(file_name, file_mandatory=True, file_fields_excluded=None):

    if file_fields_excluded is None:
        file_fields_excluded = ['__header__', '__version__', '__globals__']

    file_data_raw = None
    if os.path.exists(file_name):
        file_data_raw = scipy.io.loadmat(file_name, struct_as_record=True)
    else:
        if file_mandatory:
            log_stream.error(' ===> File "' + file_name + '" does not exists.')
            raise FileNotFoundError('File is mandatory to run the algorithm')
        else:
            log_stream.warning(' ===> File "' + file_name + '" does not exists. Datasets will be defined by NoneType')

    key_data = deepcopy(list(file_data_raw.keys()))
    if file_data_raw is not None:
        for key_name in key_data:
            if key_name in file_fields_excluded:
                file_data_raw.pop(key_name)

    return file_data_raw
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert variable
def convert_var_mat(var_data_in, var_type='time', string_no_data='-',
                    reset_tmp=True, folder_tmp=None, file_tmp='times.workspace',
                    time_format='%d/%m/%Y %H:%M'):

    path_tmp = os.path.join(folder_tmp, file_tmp)
    if reset_tmp:
        if os.path.exists(path_tmp):
            os.remove(path_tmp)

    if not os.path.exists(path_tmp):
        os.makedirs(folder_tmp, exist_ok=True)

        tmp_values = var_data_in
        tmp_values = tmp_values[:, 0][0][:, 0]

        var_data_out = []
        for tmp_value in tmp_values:

            if tmp_value:
                tmp_value = str(tmp_value[0])
                if var_type == 'time':
                    tmp_time = datetime.strptime(tmp_value, time_format)
                    obj_step = pd.Timestamp(tmp_time)
                else:
                    obj_step = deepcopy(tmp_value)
            else:
                obj_step = string_no_data
            var_data_out.append(obj_step)

        write_obj(path_tmp, var_data_out)

    else:
        var_data_out = read_obj(path_tmp)

    return var_data_out
# ----------------------------------------------------------------------------------------------------------------------
