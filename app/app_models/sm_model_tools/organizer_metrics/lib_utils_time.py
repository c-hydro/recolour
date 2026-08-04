"""
Library Features:

Name:           lib_utils_time
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260713'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
import pandas as pd

from typing import Dict, List, Optional

from config_info import LOGGER_NAME, TIME_FMT_CLI

logger = logging.getLogger(LOGGER_NAME)

# set seasons tags and reference
SEASON_MONTHS = {"DJF": [12, 1, 2],"MAM": [3, 4, 5],"JJA": [6, 7, 8],"SON": [9, 10, 11],}
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# remove empty seasons
def remove_seasons_empty(season_periods: Dict[str, pd.DatetimeIndex],) -> Dict[str, pd.DatetimeIndex]:

    season_periods_filtered = {}
    for season_name, season_period in season_periods.items():

        if season_period is None or len(season_period) == 0:
            logger.warning(f" ===> Season {season_name}: no time steps available. Season removed.")
            continue
        season_periods_filtered[season_name] = season_period

    return season_periods_filtered
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create seasons period
def create_seasons_period(
        time_period_common: pd.DatetimeIndex,
        time_period_reference: pd.DatetimeIndex,
        time_period_other: pd.DatetimeIndex,
        seasons: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Split three aligned time periods into seasonal triplets.

    Parameters
    ----------
    time_period_common : pd.DatetimeIndex
        Common time period.

    time_period_reference : pd.DatetimeIndex
        Reference time period.

    time_period_other : pd.DatetimeIndex
        Other time period.

    seasons : list[str], optional
        Seasons to extract.
        Supported values are:

            ALL
            DJF
            MAM
            JJA
            SON

        If None, only ALL is returned.

    Returns
    -------
    dict

    Example
    -------
    {
        "ALL": DataFrame(
            common,
            reference,
            other
        ),
        "DJF": DataFrame(...),
        "JJA": DataFrame(...)
    }
    """

    # ---------------------------------------------------------------------
    # Check inputs
    for name, period in (
            ("time_period_common", time_period_common),
            ("time_period_reference", time_period_reference),
            ("time_period_other", time_period_other),
    ):

        if not isinstance(period, pd.DatetimeIndex):
            raise TypeError(
                f"'{name}' must be a pandas.DatetimeIndex."
            )

    if len(time_period_common) != len(time_period_reference):
        raise RuntimeError(
            "Common and reference periods have different lengths."
        )

    if len(time_period_common) != len(time_period_other):
        raise RuntimeError(
            "Common and other periods have different lengths."
        )

    if seasons is None or len(seasons) == 0:
        seasons = ["ALL"]

    seasons = [season.upper() for season in seasons]

    # ---------------------------------------------------------------------
    # Build aligned table
    time_table = pd.DataFrame(
        {
            "common": time_period_common,
            "reference": time_period_reference,
            "other": time_period_other,
        }
    )

    # ---------------------------------------------------------------------
    # Split by season
    output = {}

    for season in seasons:

        if season == "ALL":
            output["ALL"] = time_table.copy()
            continue

        if season not in SEASON_MONTHS:
            raise ValueError(
                f"Unknown season '{season}'. "
                "Supported seasons are: ALL, DJF, MAM, JJA, SON."
            )

        months = SEASON_MONTHS[season]

        output[season] = (
            time_table[
                time_table["common"].dt.month.isin(months)
            ]
            .reset_index(drop=True)
        )

    return output
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create shifted time period
def _create_shifted_time_period(
        time_period_common: pd.DatetimeIndex,
        time_cfg: dict,
        tag_prefix: str,
) -> pd.DatetimeIndex:
    """
    Create a derived time period from the common time period.

    Month and day values are interpreted as relative offsets:
        ref_month = -1    -> previous month
        ref_day = -1      -> previous day
        other_day = 1     -> next day

    Hour and minute values replace the corresponding clock components:
        ref_hour = 23
        ref_minutes = 0
    """

    month_offset_raw = time_cfg.get(f"{tag_prefix}_month")
    day_offset_raw = time_cfg.get(f"{tag_prefix}_day")
    hour_raw = time_cfg.get(f"{tag_prefix}_hour")
    minute_raw = time_cfg.get(f"{tag_prefix}_minutes")

    month_offset = (
        int(month_offset_raw)
        if month_offset_raw is not None
        else 0
    )

    day_offset = (
        int(day_offset_raw)
        if day_offset_raw is not None
        else 0
    )

    hour = (
        int(hour_raw)
        if hour_raw is not None
        else None
    )

    minute = (
        int(minute_raw)
        if minute_raw is not None
        else None
    )

    if hour is not None and not 0 <= hour <= 23:
        raise RuntimeError(
            f"Invalid '{tag_prefix}_hour': {hour}. "
            "Expected a value between 0 and 23."
        )

    if minute is not None and not 0 <= minute <= 59:
        raise RuntimeError(
            f"Invalid '{tag_prefix}_minutes': {minute}. "
            "Expected a value between 0 and 59."
        )

    time_period_derived = []

    for time_step_common in time_period_common:

        # Apply relative month/day offsets
        time_step_derived = (
            time_step_common
            + pd.DateOffset(
                months=month_offset,
                days=day_offset
            )
        )

        # Replace hour/minute only when configured
        replace_components = {
            "second": 0,
            "microsecond": 0,
        }

        if hour is not None:
            replace_components["hour"] = hour

        if minute is not None:
            replace_components["minute"] = minute

        try:
            time_step_derived = time_step_derived.replace(
                **replace_components
            )
        except ValueError as exc:
            raise RuntimeError(
                f"Unable to create '{tag_prefix}' timestamp from "
                f"'{time_step_common}'. "
                f"Offsets: month={month_offset}, day={day_offset}; "
                f"clock: hour={hour}, minute={minute}."
            ) from exc

        time_period_derived.append(time_step_derived)

    return pd.DatetimeIndex(time_period_derived)


# method to create time periods
def create_time_period(
        time_cfg: dict, time_start_cli: Optional[str] = None,time_end_cli: Optional[str] = None,) -> (
        tuple)[pd.DatetimeIndex,pd.DatetimeIndex,pd.DatetimeIndex,dict]:
    """
    Create three time periods:

    1. time_period_common
       Created directly from time_start, time_end and frequency.

    2. time_period_reference
       Created from time_period_common using the ref_* settings.

    3. time_period_other
       Created from time_period_common using the other_* settings.

    The month/day tags are relative offsets, while hour/minute tags
    replace the corresponding timestamp components.

    Examples
    --------
    ref_day = -1, ref_hour = 23:
        2018-01-02 00:00 -> 2018-01-01 23:00

    other_day = 0, other_hour = 00:
        2018-01-02 00:00 -> 2018-01-02 00:00
    """

    time_format = time_cfg.get("format","%Y-%m-%d %H:%M")
    time_frequency = time_cfg.get("frequency","D")

    # CLI values have priority over JSON values
    time_start_raw = (
        time_start_cli
        if time_start_cli is not None
        else time_cfg.get("time_start")
    )

    time_end_raw = (
        time_end_cli
        if time_end_cli is not None
        else time_cfg.get("time_end")
    )

    if time_start_raw is None:
        raise RuntimeError(
            "Time start is not defined. "
            "Provide it using the CLI or set "
            "'time.time_start' in the JSON."
        )

    if time_end_raw is None:
        raise RuntimeError(
            "Time end is not defined. "
            "Provide it using the CLI or set "
            "'time.time_end' in the JSON."
        )

    try:
        time_start = pd.to_datetime(
            time_start_raw,
            format=time_format
        )
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Invalid time_start '{time_start_raw}'. "
            f"Expected format: '{time_format}'."
        ) from exc

    try:
        time_end = pd.to_datetime(
            time_end_raw,
            format=time_format
        )
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Invalid time_end '{time_end_raw}'. "
            f"Expected format: '{time_format}'."
        ) from exc

    if time_start > time_end:
        raise RuntimeError(
            f"time_start '{time_start}' must be earlier than "
            f"or equal to time_end '{time_end}'."
        )

    # -------------------------------------------------------------------------
    # 1. Common period
    try:
        time_period_common = pd.date_range(
            start=time_start,
            end=time_end,
            freq=time_frequency,
            inclusive="both"
        )
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Invalid time frequency '{time_frequency}'. "
            "Use values such as 'h', '3h', 'D', or '2D'."
        ) from exc

    if time_period_common.empty:
        raise RuntimeError(
            "The generated common time period is empty. "
            "Check the time bounds and frequency."
        )

    # -------------------------------------------------------------------------
    # 2. Reference period
    time_period_reference = _create_shifted_time_period(
        time_period_common=time_period_common,
        time_cfg=time_cfg,
        tag_prefix="ref"
    )

    # -------------------------------------------------------------------------
    # 3. Other period
    time_period_other = _create_shifted_time_period(
        time_period_common=time_period_common,
        time_cfg=time_cfg,
        tag_prefix="other"
    )

    # -------------------------------------------------------------------------
    # Collect information
    time_info = {
        "common": {
            "start": time_period_common[0].strftime(time_format),
            "end": time_period_common[-1].strftime(time_format),
            "steps": len(time_period_common)
        },
        "reference": {
            "start": time_period_reference[0].strftime(time_format),
            "end": time_period_reference[-1].strftime(time_format),
            "steps": len(time_period_reference),
            "month_offset": time_cfg.get("ref_month"),
            "day_offset": time_cfg.get("ref_day"),
            "hour": time_cfg.get("ref_hour"),
            "minutes": time_cfg.get("ref_minutes"),
        },
        "other": {
            "start": time_period_other[0].strftime(time_format),
            "end": time_period_other[-1].strftime(time_format),
            "steps": len(time_period_other),
            "month_offset": time_cfg.get("other_month"),
            "day_offset": time_cfg.get("other_day"),
            "hour": time_cfg.get("other_hour"),
            "minutes": time_cfg.get("other_minutes"),
        },
        "frequency": time_frequency,
        "format": time_format
    }

    return (
        time_period_common,
        time_period_reference,
        time_period_other,
        time_info
    )
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
