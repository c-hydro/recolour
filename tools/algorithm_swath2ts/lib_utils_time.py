"""
Library Features:

Name:          lib_utils_time
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import pandas as pd

from copy import deepcopy
from datetime import datetime, timedelta
from datetime import date

from lib_utils_generic import fill_tags2string
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set time information
def set_time(time_run_args=None, time_run_file=None, time_format='%Y-%m-%d %H:$M',
             time_run_file_start=None, time_run_file_end=None,
             time_period=1, time_frequency='24H', time_rounding='D', time_reverse=True):

    logging.info(' ----> Set time info ... ')
    if (time_run_file_start is None) and (time_run_file_end is None):

        logging.info(' -----> Time info defined by "time_run" argument ... ')

        if time_run_args is not None:
            time_run = time_run_args
            logging.info(' ------> Time ' + time_run + ' set by argument')
        elif (time_run_args is None) and (time_run_file is not None):
            time_run = time_run_file
            logging.info(' ------> Time ' + time_run + ' set by user')
        elif (time_run_args is None) and (time_run_file is None):
            time_now = date.today()
            time_run = time_now.strftime(time_format)
            logging.info(' ------> Time ' + time_run + ' set by system')
        else:
            logging.info(' ----> Set time info ... FAILED')
            logging.error(' ===> Argument "time_run" is not correctly set')
            raise IOError('Time type or format is wrong')

        time_tmp = pd.Timestamp(time_run)

        time_day_run = time_tmp.floor(time_rounding)
        time_day_end = time_tmp.floor('H')
        time_day_start = time_tmp.floor(time_rounding)

        time_ramge_end = time_tmp.floor(time_rounding)

        if time_period > 0:
            time_range = pd.date_range(end=time_ramge_end, periods=time_period, freq=time_frequency)
        else:
            logging.warning(' ===> TimePeriod must be greater then 0. TimePeriod is set automatically to 1')
            time_range = pd.date_range(end=time_ramge_end, periods=1, freq=time_frequency)

        if time_reverse:
            time_range = time_range[::-1]

        time_chunks = {time_day_run: [time_day_start, time_day_end]}

        if time_period > 0:
            for time_step in time_range:
                time_tmp_step = deepcopy(time_step)
                time_tmp_end = deepcopy(time_tmp_step) - pd.Timedelta(seconds=1)
                time_tmp_start = pd.date_range(end=time_tmp_step, periods=2, freq=time_frequency)[0]

                time_run_start = deepcopy(time_tmp_start)

                if time_run_start not in list(time_chunks.keys()):
                    time_chunks[time_run_start] = [time_tmp_start, time_tmp_end]

        logging.info(' -----> Time info defined by "time_run" argument ... DONE')

    elif (time_run_file_start is not None) and (time_run_file_end is not None):

        logging.info(' -----> Time info defined by "time_start" and "time_end" arguments ... ')

        time_run = date.today()
        time_run = time_run.strftime(time_format)
        time_run = pd.Timestamp(time_run)

        time_period_run = time_run.floor(time_rounding)

        time_run_file_start = pd.Timestamp(time_run_file_start)
        time_period_start = time_run_file_start.floor(time_rounding)
        time_run_file_end = pd.Timestamp(time_run_file_end)
        time_period_end = time_run_file_end.floor(time_rounding)

        time_chunks = {time_period_run: [time_period_start, time_period_end]}

        logging.info(' -----> Time info defined by "time_start" and "time_end" arguments ... DONE')

    else:
        logging.info(' ----> Set time info ... FAILED')
        logging.error(' ===> Arguments "time_start" and/or "time_end" is/are not correctly set')
        raise IOError('Time type or format is wrong')

    logging.info(' ----> Set time info ... DONE')

    return time_run, time_chunks

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure date string
def parse_time_string(time_string, time_format='%Y-%m-%d %H:%M'):
    time_stamp_obj = pd.Timestamp(time_string)
    time_stamp_str = time_stamp_obj.strftime(time_format)
    return datetime.strptime(time_stamp_str, time_format)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to fill path string
def fill_time_string(path_in, time_info_list=None,
                     time_tags_list=None, time_template_list=None):

    if (time_tags_list is not None) and (time_template_list is not None):

        if not isinstance(time_info_list, list):
            time_info_list = [time_info_list]

        if time_info_list.__len__() == 1:
            time_info_list = time_info_list * time_tags_list.__len__()

        obj_tags_format, obj_tags_value = {}, {}
        for time_tags_step, time_tmpl_step, time_info_step in zip(time_tags_list, time_template_list, time_info_list):
            obj_tags_format[time_tags_step] = time_tmpl_step
            obj_tags_value[time_tags_step] = time_info_step

        path_out = fill_tags2string(path_in, tags_format=obj_tags_format, tags_filling=obj_tags_value)

    else:
        path_out = deepcopy(path_in)

    return path_out

# ----------------------------------------------------------------------------------------------------------------------
