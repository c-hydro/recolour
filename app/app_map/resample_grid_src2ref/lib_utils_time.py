"""
Library Features:

Name:          lib_utils_time
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240709'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import pandas as pd

from datetime import date

from lib_info_args import logger_name

# logging
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to set time info
def set_time_info(time_run_args=None, time_run_file=None, time_format='%Y-%m-%d %H:$M',
                  time_frequency='H', time_rounding='H'):

    alg_logger.info(' ---> Set time info ... ')

    if (time_run_args is not None) and (time_run_file is None):
        time_run_reference = pd.Timestamp(time_run_args)
        time_run_reference = time_run_reference.floor(time_rounding)
    elif (time_run_args is None) and (time_run_file is not None):
        time_run_reference = pd.Timestamp(time_run_file)
        time_run_reference = time_run_reference.floor(time_rounding)
    elif (time_run_args is None) and (time_run_file is None):
        time_now = date.today()
        time_run_reference = time_now.strftime(time_format)
        time_run_reference = pd.Timestamp(time_run_reference)
        time_run_reference = time_run_reference.floor(time_rounding)
    else:
        alg_logger.error(' ===> Argument "time_reference" is not correctly set')
        raise IOError('Time type or format is wrong')

    alg_logger.info(' ---> Set time info ... DONE')

    return time_run_reference

# ----------------------------------------------------------------------------------------------------------------------
