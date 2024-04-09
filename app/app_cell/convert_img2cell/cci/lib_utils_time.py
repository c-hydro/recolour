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
from datetime import date, timedelta

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to update time info
def update_time_info(data_settings,
                     time_run, time_start, time_end, time_format='%Y-%m-%d %H:%M'):

    # check time_end in data settings
    if 'time_end' in list(data_settings['time']):
        time_tmp = data_settings['time']['time_end']
        if time_tmp is not None:
            logging.warning(' ===> Variable "time_end" is not None and the algorithm will update it on settings')
        if isinstance(time_end, pd.Timestamp):
            time_end = time_end.strftime(time_format)
        data_settings['time']['time_end'] = time_end
    else:
        logging.error(' ===> Variable "time_end" is not defined in the settings')
        raise RuntimeError('Variable must be included in the settings')
    # check time_start in data settings
    if 'time_start' in list(data_settings['time']):
        time_tmp = data_settings['time']['time_start']
        if time_tmp is not None:
            logging.warning(' ===> Variable "time_start" is not None and the algorithm will update it on settings')
        if isinstance(time_start, pd.Timestamp):
            time_start = time_start.strftime(time_format)
        data_settings['time']['time_start'] = time_start
    else:
        logging.error(' ===> Variable "time_start" is not defined in the settings')
        raise RuntimeError('Variable must be included in the settings')

    # add time run to data settings
    if isinstance(time_run, pd.Timestamp):
        time_run = time_run.strftime(time_format)
    data_settings['time']['time_run'] = time_run

    return data_settings
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set time info
def set_time_info(time_run_args=None, time_run_file=None, time_format='%Y-%m-%d %H:$M',
                  time_run_file_start=None, time_run_file_end=None,
                  time_period=1, time_frequency='H', time_rounding='H'):

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
        time_tmp = time_tmp.replace(minute=0, second=0)

        time_now = deepcopy(time_tmp)
        time_day_end = time_tmp.floor(time_rounding)
        if time_frequency == 'D':
            time_day_start = time_now - timedelta(hours=24*time_period)
        if time_frequency == 'H':
            time_day_start = time_now - timedelta(hours=time_period)

        if time_period > 1:
            # time_range_arr = pd.date_range(end=time_day_start, periods=time_period, freq=time_frequency)
            time_range_arr = pd.date_range(end=time_now, periods=time_period, freq=time_frequency)
            time_range_start = time_range_arr[0]
            time_range_end = deepcopy(time_day_end)
        elif time_period == 1:
            # time_range_arr = pd.date_range(end=time_day_start, periods=time_period, freq=time_frequency)
            time_range_arr = pd.date_range(end=time_now, periods=time_period, freq=time_frequency)
            time_range_start = time_range_arr[0]
            time_range_end = deepcopy(time_day_end)
        elif time_period == 0:
            logging.warning(' ===> Variable "time_period" must be greater then 0. It will be set automatically to 1')
            time_range_start = deepcopy(time_day_start)
            time_range_end = deepcopy(time_day_end)
        else:
            logging.error(' ===> Variable "time_period" is not correctly defined.')
            raise RuntimeError('Check your settings file to define the variable')

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
        time_range_arr = pd.date_range(start=time_run_file_start, end=time_run_file_end, freq=time_frequency)

        time_range_start, time_range_end = time_range_arr[0], time_range_arr[-1]

        logging.info(' -----> Time info defined by "time_start" and "time_end" arguments ... DONE')

    else:
        logging.info(' ----> Set time info ... FAILED')
        logging.error(' ===> Arguments "time_start" and/or "time_end" is/are not correctly set')
        raise IOError('Time type or format is wrong')

    # update time range according with minimum time delta
    time_range_start, time_range_end = set_timedelta_min(time_range_start, time_range_end, time_delta_min='24H')

    logging.info(' ----> Set time info ... DONE')

    return time_now, time_run, time_range_start, time_range_end

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set a minimum time delta to compute time start and end
def set_timedelta_min(time_start, time_end, time_delta_min='24H', time_rounding='D', time_approximation='floor'):
    time_delta_upd = None
    time_delta_min = pd.Timedelta(time_delta_min)
    time_delta_elapsed = time_end - time_start
    while time_delta_elapsed < time_delta_min:

        time_start = time_start - pd.Timedelta('1H')
        time_delta_elapsed = time_end - time_start

        if time_delta_upd is None:
            time_delta_upd = True

    if time_delta_upd:
        if time_approximation == 'round':
            time_start = time_start.round(time_rounding)
        elif time_approximation == 'floor':
            time_start = time_start.floor(time_rounding)
        elif time_approximation == 'ceil':
            time_start = time_start.ceil(time_rounding)
        elif time_approximation == 'nothing':
            pass
        else:
            logging.error(' ===> Time approximation "' + time_approximation + '" is not supported')
            raise NotImplemented('Case not implemented yet')

    return time_start, time_end
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
