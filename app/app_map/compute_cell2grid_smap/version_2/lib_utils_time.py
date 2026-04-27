"""
Library Features:

Name:           lib_utils_time
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
from datetime import datetime, timedelta

from config_utils import LOGGER_NAME, TIME_FMT_CLI

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to normalize time to hour
def normalize_time_to_hour(time_obj):
    return time_obj.replace(minute=0, second=0, microsecond=0)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to normalize time to midnight
def normalize_time_to_midnight(time_obj):
    return time_obj.replace(hour=0, minute=0, second=0, microsecond=0)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to parse time string
def parse_time(time_string):
    try:
        return datetime.strptime(time_string, TIME_FMT_CLI)
    except ValueError as exc:
        raise ValueError('Time must have format "YYYY-MM-DD HH:MM"') from exc
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to parse reference time
def parse_reference_time(time_string=None):
    if time_string is None:
        return normalize_time_to_hour(datetime.now())
    return normalize_time_to_hour(parse_time(time_string))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to parse time delta
def parse_time_delta(time_string):
    """
    Supported formats:
    - H, D
    - 1H, 3H, 24H
    - 1D, 2D, 5D
    """
    if time_string is None:
        raise ValueError("Time string cannot be None")

    time_string = str(time_string).strip().upper()

    if time_string == "H":
        time_string = "1H"
    elif time_string == "D":
        time_string = "1D"

    match = re.match(r"^(?P<value>\d+)(?P<unit>[HD])$", time_string)
    if match is None:
        raise ValueError('Time string must have format like "H", "D", "1H", "24H", "2D"')

    time_value = int(match.group("value"))
    time_unit = match.group("unit")

    if time_unit == "H":
        return timedelta(hours=time_value)
    elif time_unit == "D":
        return timedelta(days=time_value)
    else:
        raise ValueError('Unsupported time unit. Use "H" or "D"')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to resolve time tags
def resolve_time_tags(template_string, time_dict):
    template_resolved = str(template_string)

    for time_key, time_value in time_dict.items():
        if time_value is None:
            continue

        tag_pattern = r"\{" + re.escape(time_key) + r":([^}]+)\}"

        def replace_tag(match):
            time_format = match.group(1)
            return time_value.strftime(time_format)

        template_resolved = re.sub(tag_pattern, replace_tag, template_resolved)

    return template_resolved
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to define time window
def resolve_time_window(settings, reference_time):
    time_settings = settings.get("time", {})

    time_start_str = time_settings.get("time_start")
    time_end_str = time_settings.get("time_end")
    time_delta_str = time_settings.get("time_delta", time_settings.get("time_period", "2D"))

    floor_start_to_midnight = bool(time_settings.get("floor_start_to_midnight", False))
    floor_end_to_midnight = bool(time_settings.get("floor_end_to_midnight", False))

    if (time_start_str is not None) and (time_end_str is not None):
        time_start = normalize_time_to_hour(parse_time(time_start_str))
        time_end = normalize_time_to_hour(parse_time(time_end_str))

    elif (time_start_str is None) and (time_end_str is None):
        time_end = normalize_time_to_hour(reference_time)
        time_delta = parse_time_delta(time_delta_str)
        time_start = normalize_time_to_hour(time_end - time_delta)

    else:
        raise RuntimeError(
            'time_start and time_end must be both provided, or neither provided'
        )

    if floor_start_to_midnight:
        time_start = normalize_time_to_midnight(time_start)

    if floor_end_to_midnight:
        time_end = normalize_time_to_midnight(time_end)

    if time_end < time_start:
        raise RuntimeError(f"time_end {time_end} is earlier than time_start {time_start}")

    return time_start, time_end
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to iterate over steps
def iter_time_steps(settings, time_start, time_end):
    time_settings = settings.get("time", {})
    time_frequency_str = time_settings.get("time_frequency", "H")
    time_frequency_delta = parse_time_delta(time_frequency_str)

    current_time = normalize_time_to_hour(time_start)
    final_time = normalize_time_to_hour(time_end)

    while current_time <= final_time:
        yield current_time
        current_time += time_frequency_delta
# ----------------------------------------------------------------------------------------------------------------------