"""
Library Features:

Name:           lib_utils_basic
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np

from pathlib import Path
from datetime import datetime
from typing import Any, Dict, Optional
from glob import glob

from lib_utils_time import resolve_time_tags
from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to make folder
def build_path(
    dataset_cfg: Dict[str, Any],
    *,
    time: Optional[datetime] = None,
    time_start: Optional[datetime] = None,
    time_end: Optional[datetime] = None,
    time_run: Optional[datetime] = None,
) -> Path:
    values = {
        key: value
        for key, value in {
            "time": time,
            "time_start": time_start,
            "time_end": time_end,
            "time_run": time_run,
        }.items()
        if value is not None
    }

    folder = fill_time_tags(dataset_cfg["folder"], **values)
    filename = fill_time_tags(dataset_cfg["filename"], **values)
    return Path(folder) / filename
# ----------------------------------------------------------------------------------------------------------------------
