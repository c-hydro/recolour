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
        reference_time_period: pd.DatetimeIndex,
        seasons: Optional[List[str]] = None,
) -> Dict[str, pd.DatetimeIndex]:
    """
    Split a DatetimeIndex into seasonal DatetimeIndexes.

    Parameters
    ----------
    reference_time_period : pd.DatetimeIndex
        Input time period.

    seasons : list[str], optional
        Seasons to extract. Supported values are:
            ALL, DJF, MAM, JJA, SON

        If None or empty, only ALL is returned.

    Returns
    -------
    dict
        Example:
        {
            "ALL": DatetimeIndex(...),
            "DJF": DatetimeIndex(...),
            "JJA": DatetimeIndex(...)
        }
    """

    if not isinstance(reference_time_period, pd.DatetimeIndex):
        raise TypeError("'reference_time_period' must be a pandas.DatetimeIndex.")

    if seasons is None or len(seasons) == 0:
        seasons = ["ALL"]

    seasons = [season.upper() for season in seasons]

    output = {}

    for season in seasons:

        if season == "ALL":
            output["ALL"] = reference_time_period
            continue

        if season not in SEASON_MONTHS:
            raise ValueError(
                f"Unknown season '{season}'. "
                f"Supported seasons are: ALL, DJF, MAM, JJA, SON."
            )

        months = SEASON_MONTHS[season]
        output[season] = reference_time_period[
            reference_time_period.month.isin(months)
        ]

    return output
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create time period
def create_time_period(
        time_cfg: dict,
        time_start_cli: Optional[str] = None,
        time_end_cli: Optional[str] = None,
) -> tuple[pd.DatetimeIndex, dict]:
    """
    Create a DatetimeIndex using CLI time bounds when provided,
    otherwise using the bounds defined in the JSON configuration.

    Reference keys override datetime components when defined:
        ref_month
        ref_day
        ref_hour
        ref_minutes

    Supported frequency examples:
        h, 1h, 3h, 6h
        D, 1D, 2D
    """

    time_format = time_cfg.get("format", "%Y-%m-%d %H:%M")
    time_frequency = time_cfg.get("frequency", "D")

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
            "Provide it using the CLI or set 'time.time_start' in the JSON."
        )

    if time_end_raw is None:
        raise RuntimeError(
            "Time end is not defined. "
            "Provide it using the CLI or set 'time.time_end' in the JSON."
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

    # Reference values override the corresponding datetime components
    reference_components = {}

    if time_cfg.get("ref_month") is not None:
        reference_components["month"] = int(time_cfg["ref_month"])

    if time_cfg.get("ref_day") is not None:
        reference_components["day"] = int(time_cfg["ref_day"])

    if time_cfg.get("ref_hour") is not None:
        reference_components["hour"] = int(time_cfg["ref_hour"])

    if time_cfg.get("ref_minutes") is not None:
        reference_components["minute"] = int(time_cfg["ref_minutes"])

    if reference_components:
        reference_components["second"] = 0
        reference_components["microsecond"] = 0

        try:
            time_start = time_start.replace(**reference_components)
            time_end = time_end.replace(**reference_components)
        except ValueError as exc:
            raise RuntimeError(
                "Invalid reference time components: "
                f"{reference_components}. "
                "Check ref_month, ref_day, ref_hour and ref_minutes."
            ) from exc

    if time_start > time_end:
        raise RuntimeError(
            f"time_start '{time_start}' must be earlier than or equal to "
            f"time_end '{time_end}'."
        )

    try:
        time_period = pd.date_range(
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

    if time_period.empty:
        raise RuntimeError(
            "The generated time period is empty. "
            "Check the time bounds and frequency."
        )

    time_info = {
        "start": time_start.strftime(time_format),
        "end": time_end.strftime(time_format),
        "first": time_period[0].strftime(time_format),
        "last": time_period[-1].strftime(time_format),
        "frequency": time_frequency,
        "steps": len(time_period),
        "format": time_format,
        "ref_month": time_cfg.get("ref_month"),
        "ref_day": time_cfg.get("ref_day"),
        "ref_hour": time_cfg.get("ref_hour"),
        "ref_minutes": time_cfg.get("ref_minutes"),
    }

    return time_period, time_info
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
