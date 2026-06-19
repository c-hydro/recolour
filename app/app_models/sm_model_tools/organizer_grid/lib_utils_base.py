"""
Library Features:

Name:           lib_utils_base
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260619'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
from datetime import datetime
from pathlib import Path

from config_info import LOGGER_NAME

# logger
logger_stream = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to format file_path
def format_file_path(file_pattern: str, time_now: datetime, point_tag: str) -> str:

    point = str(point_tag).replace("point_", "")

    return file_pattern.format(
        yyyy=time_now.strftime("%Y"),
        mm=time_now.strftime("%m"),
        dd=time_now.strftime("%d"),
        HH=time_now.strftime("%H"),
        hh=time_now.strftime("%H"),
        yyyymmdd=time_now.strftime("%Y%m%d"),
        yyyymmddhh=time_now.strftime("%Y%m%d%H"),
        yyyymmddhhmm=time_now.strftime("%Y%m%d%H%M"),
        point=point,
        point_tag=point_tag,
        time=time_now,
        time_now=time_now,
    )
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to extract tag from filename
def extract_tag_from_filename(path: str, pattern: str) -> str:
    name = Path(path).name
    match = re.search(pattern, name)
    if not match:
        raise RuntimeError(f"Cannot extract point tag from file name: {name}; regex={pattern}")
    return match.group("tag") if "tag" in match.groupdict() else match.group(1)
# ----------------------------------------------------------------------------------------------------------------------
