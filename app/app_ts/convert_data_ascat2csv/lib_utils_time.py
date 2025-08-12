# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import pandas as pd

from datetime import datetime, timedelta
from typing import Optional, Tuple
from pandas.tseries.frequencies import to_offset

from lib_utils_info import logger_name

# set logger
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to get feasible dates within a range based on frequency
def get_feasible_dates(start_ts, end_ts, frequency):
    """
    Generate valid dates between start_ts and end_ts according to frequency.
    """
    # Generate all valid dates at the given frequency
    valid_dates = pd.date_range(start=start_ts, end=end_ts, freq=frequency)

    # If there’s no valid date, return None
    if len(valid_dates) == 0:
        return None

    return valid_dates
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to get the latest feasible date within a range
def get_latest_feasible_date(start_ts, end_ts, frequency):
    """
    Get the last feasible date within the given range for saving data.
    """
    valid_dates = get_feasible_dates(start_ts, end_ts, frequency)
    return valid_dates[-1] if valid_dates is not None else None
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to define a time reference based on a timestamp and frequency
def define_time_reference(time_stamp, time_frequency='D', time_format=None):
    """
    Truncate a pandas Timestamp to the specified reference frequency
    and output it in the given format (if suitable).

    Parameters
    ----------
    time_stamp : pandas.Timestamp
        The input timestamp.
    time_frequency : str, optional
        Truncation level:
        - 'D' : Day (default)
        - 'M' : Month
        - 'Y' : Year
    time_format : str, optional
        Desired output format string (strftime syntax).
        If not suitable for the given frequency, a warning is issued
        and a default format will be used.

    Returns
    -------
    str
        Truncated date string according to the specified frequency and format.
    """

    if time_frequency == 'MS':
        log_stream.warning(" ===> 'MS' frequency is not supported. Use 'M' for month start.")
        time_frequency = 'M'

    # Define default formats for each frequency
    default_formats = {
        'D': "%Y-%m-%d",
        'M': "%Y-%m",
        'Y': "%Y"
    }

    # Check format suitability
    if time_format is None:
        time_format = default_formats[time_frequency]
    else:
        # Simple suitability check: ensure format matches frequency granularity
        if time_frequency == 'Y' and ('%m' in time_format or '%d' in time_format):
            log_stream.warning(f" ===> Provided format '{time_format}' is too detailed for yearly frequency. "
                  f"Defaulting to '{default_formats['Y']}'")
            time_format = default_formats['Y']
        elif time_frequency == 'M' and '%d' in time_format:
            log_stream.warning(f" ===> Provided format '{time_format}' is too detailed for monthly frequency. "
                  f"Defaulting to '{default_formats['M']}'")
            time_format = default_formats['M']

    # Truncate timestamp according to frequency
    if time_frequency == 'D':
        truncated = time_stamp.normalize()
    elif time_frequency == 'M':
        truncated = pd.Timestamp(year=time_stamp.year, month=time_stamp.month, day=1)
    elif time_frequency == 'Y':
        truncated = pd.Timestamp(year=time_stamp.year, month=1, day=1)
    else:
        log_stream.error(f" ===> Invalid time frequency '{time_frequency}'. ")
        raise ValueError("Invalid time_frequency. Use 'D', 'M', or 'Y'.")

    return truncated.strftime(time_format)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to validate frequency for time ranges
def _validate_frequency(freq: str, start: pd.Timestamp, end: pd.Timestamp,
                        require_multiple: bool = True,
                        strict_align: bool = False) -> None:
    """
    Validate that `freq` is a usable pandas offset for [start, end].

    - Checks freq is parsable.
    - Checks it advances time (i.e., next boundary > start).
    - Optionally requires start to align to freq anchor (strict_align).
    - Optionally requires that the period contains >=2 chunks (require_multiple).
    """
    if end < start:
        raise ValueError("time_end must be >= time_start")

    # 1) Parsable / valid frequency
    try:
        off = to_offset(freq)
    except Exception as e:
        raise ValueError(f"Invalid frequency '{freq}': {e}") from e

    # 2) It must advance time
    #    Build two boundaries from 'start' to ensure the second is strictly greater.
    bounds_test = pd.date_range(start=start, periods=2, freq=off)
    if len(bounds_test) < 2 or not (bounds_test[1] > bounds_test[0]):
        raise ValueError(f"Frequency '{freq}' does not advance time from start={start}.")

    # 3) Alignment check (optional)
    if strict_align:
        # If the first boundary generated by freq isn't exactly 'start',
        # then start isn't aligned to the anchor (e.g., MS, W-MON, etc.)
        aligned_first = pd.date_range(start=start, periods=1, freq=off)[0]
        if aligned_first != start:
            raise ValueError(
                f"Start {start} is not aligned to freq '{freq}'. "
                f"Set strict_align=False to allow partial first chunk."
            )

    # 4) Require multiple chunks (optional)
    if require_multiple:
        # Construct boundaries across the span; if < 2 intervals, reject.
        # (We include start and end explicitly to be robust to anchors.)
        bounds = pd.date_range(start=start, end=end, freq=off)
        if len(bounds) == 0 or bounds[0] != start:
            bounds = pd.DatetimeIndex([start]).append(bounds)
        if bounds[-1] != end:
            bounds = bounds.append(pd.DatetimeIndex([end]))
        if len(bounds) < 3:
            raise ValueError(
                f"Frequency '{freq}' is too coarse for [{start}, {end}] "
                f"to form multiple chunks. Try a finer freq."
            )
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute time ranges based on start, end, and frequency
def compute_time_ranges(
    time_start, time_end, freq='D',
    closed='right',
    require_multiple: bool = True,
    strict_align: bool = False):
    """
    Create contiguous ranges over [time_start, time_end] using `freq`,
    after validating that the frequency is suitable for chunking.

    Raises ValueError if validation fails (i.e., unsuitable frequency).
    """
    start = pd.to_datetime(time_start)
    end   = pd.to_datetime(time_end)

    _validate_frequency(freq, start, end,
                        require_multiple=require_multiple,
                        strict_align=strict_align)

    # Build boundaries (anchor-aware), then force inclusion of start/end.
    off = to_offset(freq)
    bounds = pd.date_range(start=start, end=end, freq=off)
    if len(bounds) == 0 or bounds[0] != start:
        bounds = pd.DatetimeIndex([start]).append(bounds)
    if bounds[-1] != end:
        bounds = bounds.append(pd.DatetimeIndex([end]))

    # Emit intervals with chosen closure
    ns = pd.Timedelta('1ns')
    ranges = []
    for left, right in zip(bounds[:-1], bounds[1:]):
        if closed == 'right':          # [left, right)
            r_start, r_end = left, min(right - ns, end)
        elif closed == 'left':         # (left, right]
            r_start, r_end = min(left + ns, right), right
        elif closed == 'both':         # [left, right] (shrink interior segments by 1ns)
            r_start, r_end = left, right if right == end else right - ns
        elif closed == 'neither':      # (left, right)
            r_start = min(left + ns, right)
            r_end   = right if right == end else right - ns
        else:
            raise ValueError("closed must be one of {'right','left','both','neither'}")

        if r_start > r_end:
            # Pathological edge case; fallback to hard bounds
            r_start, r_end = left, right
        ranges.append((r_start, r_end))

    return ranges
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute a time window based on given parameters
def compute_time_window(
    time_now: Optional[datetime] = None,
    time_start: Optional[datetime] = None,
    time_end: Optional[datetime] = None,
    *,
    period: int = 1,
    frequency: str = 'D',
    enforce_span: bool = False
) -> Optional[Tuple[datetime, datetime]]:
    """
    Compute a time window with time_now as the primary reference.
    If time_now is given, it becomes the main anchor point for window computation.
    """

    # Frequency → timedelta
    freq_map = {
        'S': timedelta(seconds=period),
        'M': timedelta(minutes=period),
        'H': timedelta(hours=period),
        'D': timedelta(days=period),
        'W': timedelta(weeks=period)
    }
    delta = freq_map.get(frequency.upper())
    if delta is None:
        raise ValueError(f"Unsupported frequency: {frequency}")

    # ---- Priority logic ----
    # 1) No start/end → window anchored to time_now
    if time_start is None and time_end is None:
        start = time_now - delta
        end = time_now

    # 2) Only start given
    elif time_start is not None and time_end is None:
        if time_now < time_start:
            # time_now before start → shift window back from time_now
            end = time_now
            start = time_now - delta
        else:
            # Anchor to time_now if within range
            end = min(time_now, time_start + delta)
            start = time_start

    # 3) Only end given
    elif time_start is None and time_end is not None:
        if time_now < time_end:
            # Anchor end at time_now
            end = time_now
            start = time_now - delta
        else:
            # time_now after end → normal backward window
            end = time_end
            start = time_end - delta

    # 4) Both start/end given
    else:

        if time_now is not None:
            # Clamp end to time_now if time_now < time_end
            if time_now < time_start:
                # time_now before range → shift window back from time_now
                end = time_now
                start = time_now - delta
            elif time_now < time_end:
                # Use time_now as end
                end = time_now
                start = time_start
            else:
                # time_now after range
                start = time_start
                end = time_end
        else:
            # No time_now, use given start and end
            start = time_start
            end = time_end

    # ---- Enforce fixed span if requested ----
    if enforce_span:
        end = min(end, time_now)  # never go past now
        start = end - delta

    # Final validation
    if start > end:
        return None

    return start, end
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to generate time points for a given date based on labels
def compute_time_by_labels(base_date: pd.Timestamp, labels, ref_point='end'):
    """
    Generate pd.Timestamp values for a single date based on labels.

    Parameters
    ----------
    base_date : pd.Timestamp
        The date to generate times for.
    labels : list of str
        Time ranges in format 'startHour_endHour'.
    ref_point : str
        'start', 'mid_point', or 'end'.

    Returns
    -------
    list of pd.Timestamp
    """
    result = []
    for label in labels:
        start_hour, end_hour = map(int, label.split('_'))

        if ref_point == 'start':
            offset_hours = start_hour
        elif ref_point == 'end':
            offset_hours = end_hour
        elif ref_point == 'mid_point':
            offset_hours = (start_hour + end_hour) / 2
        else:
            raise ValueError("ref_point must be 'start', 'mid_point', or 'end'")

        result.append(base_date.normalize() + pd.Timedelta(hours=offset_hours))

    return result
# ----------------------------------------------------------------------------------------------------------------------
