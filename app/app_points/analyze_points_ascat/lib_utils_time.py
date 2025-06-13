"""
Library Features:

Name:          lib_utils_time
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20250612'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
from datetime import datetime, timedelta
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to get time from datetime string
def get_time(dt_in, format='%Y-%m-%d %H:%M:%S'):
    dt_out = datetime.fromisoformat(dt_in)
    return dt_out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to round datetime to nearest 00:00 or 12:00
def round_time(dt: datetime) -> datetime:
    """
    Round datetime to nearest 00:00 or 12:00 on the same date.
    """
    date = dt.date()
    dt0 = datetime.combine(date, datetime.min.time())
    dt12 = dt0 + timedelta(hours=12)
    if abs((dt - dt0).total_seconds()) <= abs((dt - dt12).total_seconds()):
        return dt0
    return dt12
# ----------------------------------------------------------------------------------------------------------------------
