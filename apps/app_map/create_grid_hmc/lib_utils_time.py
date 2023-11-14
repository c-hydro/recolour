"""
Library Features:

Name:          lib_utils_time
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230830'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# Libraries
import logging
import pandas as pd

from copy import deepcopy
from datetime import date

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set time file
def set_time_file(time_now=None, time_start=None, time_end=None,
                  time_frequency='H', time_rounding='D',
                  time_period=None, time_reverse=True, time_direction=None):

    if (time_start is None) and (time_end is None):
        logging.error(' ===> Variable "time_start" and "time_end" are defined by NoneType')
        raise RuntimeError('Variable "time_start" and/or "time_end" must be defined to get correctly time(s)')

    if time_now is not None:
        time_now = time_now.floor(time_rounding)

    if time_start is not None:
        time_start = time_start.round(time_rounding)
    if time_end is not None:
        time_end = time_end.round(time_rounding)

    if time_frequency == 'D':
        if time_start is None:
            time_start = deepcopy(time_end.replace(hour=0, minute=0, second=0))
        if time_end is None:
            time_end = deepcopy(time_start.replace(hour=23, minute=59, second=59))
    elif time_frequency == 'H':
        if time_start is None:
            time_start = deepcopy(time_end.replace(minute=0, second=0))
        if time_end is None:
            time_end = deepcopy(time_start.replace(minute=59, second=59))
    else:
        logging.error(' ===> Variable "time_frequency" "' + str(time_frequency) + '" is not supported')
        raise NotImplemented('Case not implemented yet')

    if time_direction is None:
        time_range = pd.date_range(start=time_start, end=time_end, freq=time_frequency)
    elif time_direction == 'left':
        time_end = time_end.floor(time_rounding)
        time_range = pd.date_range(end=time_end, periods=time_period, freq=time_frequency)
    elif time_direction == 'right':
        time_start = time_start.floor(time_rounding)
        time_range = pd.date_range(start=time_start, periods=time_period, freq=time_frequency)
    else:
        logging.error(' ===> Variable "time_direction" "' + str(time_direction) + '" is not supported')
        raise NotImplemented('Case not implemented yet')

    if time_now is not None:
        if time_now > time_range[-1]:
            tmp_range = pd.date_range(end=time_now, periods=time_period, freq=time_frequency)
            intersection_range = time_range.intersection(tmp_range)
            if not intersection_range.empty:
                if tmp_range[0] > time_range[0]:
                    time_range = time_range.union(tmp_range)
                elif tmp_range[0] <= time_range[0]:
                    pass

    if time_reverse:
        time_range = time_range[::-1]

    return time_range

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to set time info
def set_time_info(time_run_args=None, time_run_file=None, time_format='%Y-%m-%d %H:$M',
                  time_run_file_start=None, time_run_file_end=None,
                  time_period=1, time_frequency='H', time_rounding='H', time_reverse=True):

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
        time_now = deepcopy(time_tmp)
        time_run = time_tmp.floor(time_rounding)

        if time_period > 0:
            time_range = pd.date_range(end=time_run, periods=time_period, freq=time_frequency)
        else:
            logging.warning(' ===> TimePeriod must be greater then 0. TimePeriod is set automatically to 1')
            time_range = pd.DatetimeIndex([time_run], freq=time_frequency)

        logging.info(' -----> Time info defined by "time_run" argument ... DONE')

    elif (time_run_file_start is not None) and (time_run_file_end is not None):

        logging.info(' -----> Time info defined by "time_start" and "time_end" arguments ... ')

        time_run_file_start = pd.Timestamp(time_run_file_start)
        time_run_file_start = time_run_file_start.floor(time_rounding)
        time_run_file_end = pd.Timestamp(time_run_file_end)
        time_run_file_end = time_run_file_end.floor(time_rounding)

        time_now = date.today()
        time_run = time_now.strftime(time_format)
        time_run = pd.Timestamp(time_run)
        time_run = time_run.floor(time_rounding)
        time_range = pd.date_range(start=time_run_file_start, end=time_run_file_end, freq=time_frequency)

        logging.info(' -----> Time info defined by "time_start" and "time_end" arguments ... DONE')

    else:
        logging.info(' ----> Set time info ... FAILED')
        logging.error(' ===> Arguments "time_start" and/or "time_end" is/are not correctly set')
        raise IOError('Time type or format is wrong')

    if time_reverse:
        time_range = time_range[::-1]

    time_chunks = get_time_chunks(time_range, time_reverse=time_reverse)

    logging.info(' ----> Set time info ... DONE')

    return time_now, time_run, time_range, time_chunks

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to get time chunks
def get_time_chunks(time_range, time_period='D', time_reverse=False):

    time_groups = time_range.to_period(time_period)
    tmp_chunks = time_range.groupby(time_groups)

    time_chunks = {}
    for key in sorted(tmp_chunks.keys(), reverse=time_reverse):
        time_chunks[key] = tmp_chunks[key]

    return time_chunks
# ----------------------------------------------------------------------------------------------------------------------
