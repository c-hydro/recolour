"""
Library Features:

Name:          lib_utils_time
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260518'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
from datetime import datetime, timedelta

import pandas as pd

# constants
TIME_FMT_CLI = "%Y-%m-%d %H:%M"
TIME_FMT_FILE = "%Y%m%d%H%M%S"
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to compute time period
def compute_time_period(
        time_run,
        time_start=None, time_end=None,
        time_delta='1D', time_frequency='H',
        floor_start=False, floor_end=False,
        type='datetime', **kwargs):

    # convert delta
    if isinstance(time_delta, str):
        time_delta = pd.Timedelta(time_delta)

    # convert run time
    time_run = pd.Timestamp(time_run)

    # compute end time
    if time_end is None:
        time_end = time_run
    else:
        time_end = pd.Timestamp(time_end)

    # compute start time
    if time_start is None:
        time_start = time_end - time_delta
    else:
        time_start = pd.Timestamp(time_start)

    # optional floor operations
    if floor_start:
        time_start = time_start.floor("D")

    if floor_end:
        time_end = time_end.floor("D")

    # return datetime tuple
    if type == 'datetime':

        time_period = (
            time_start.to_pydatetime(),
            time_end.to_pydatetime()
        )

    # return pandas date range
    elif type == 'range':

        time_period = pd.date_range(
            start=time_start,
            end=time_end,
            freq=time_frequency
        )

    else:
        raise NotImplementedError(
            f'Time period type "{type}" is not supported'
        )

    return time_period
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# template methods
def resolve_time_tags(path_template, time_dict):
    path_resolved = str(path_template)

    for time_key, time_value in time_dict.items():
        if time_value is None:
            continue

        tag_pattern = r"\{" + re.escape(time_key) + r":([^}]+)\}"

        def replace_tag(match):
            time_format = match.group(1)
            return time_value.strftime(time_format)

        path_resolved = re.sub(tag_pattern, replace_tag, path_resolved)

    return path_resolved
# ----------------------------------------------------------------------------------------------------------------------

def normalize_time_to_hour(time_obj):
    return time_obj.replace(minute=0, second=0, microsecond=0)

def normalize_time_to_midnight(time_obj):
    return time_obj.replace(hour=0, minute=0, second=0, microsecond=0)

def parse_time(time_string):
    try:
        return datetime.strptime(time_string, TIME_FMT_CLI)
    except ValueError as exc:
        raise ValueError('Time must have format "YYYY-MM-DD HH:MM"') from exc

def parse_reference_time(time_string=None):
    if time_string is None:
        return normalize_time_to_hour(datetime.now())
    return normalize_time_to_hour(parse_time(time_string))