# ----------------------------------------------------------------------------------------------------------------------
# libraries
import numpy as np
import pandas as pd
import xarray as xr
from functools import wraps
from typing import Iterable, Optional, List

from lib_utils_decoretors import simplify_list, iterate_dict
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to adapt a DataFrame to a specified time range
@iterate_dict
def adapt_dataframe_to_range(
    df: pd.DataFrame,
    time_start,
    time_end,
    freq: str = "D",
    tz: (str, None) = None,
    inclusive: str = "both",
):
    """
    Align a DataFrame to a generated time index from time_start to time_end
    at the given frequency, filling missing timestamps with NaN.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with a DateTimeIndex (or index convertible to one).
    time_start : str | pd.Timestamp | datetime-like
        Start of the target range (e.g., "2025-01-01 00:00").
    time_end : str | pd.Timestamp | datetime-like
        End of the target range (e.g., "2025-01-10 23:00").
    freq : str, default "D"
        Pandas frequency string (e.g., "D", "H", "15min", "MS").
    tz : str | None, default None
        Time zone name (e.g., "Europe/Rome"). If provided, the result index
        will be in this timezone. The input df index will be localized
        (if naive) or converted (if already tz-aware).
    inclusive : {"both","left","right","neither"}, default "both"
        Passed to `pd.date_range` to control boundary inclusion.

    Returns
    -------
    pd.DataFrame
        Reindexed DataFrame over the generated range, with missing rows as NaN.
    """
    out = df.copy()

    # Ensure DateTimeIndex
    if not isinstance(out.index, pd.DatetimeIndex):
        out.index = pd.to_datetime(out.index, format="%Y-%m-%d", errors="coerce")
    # Drop rows with unparseable timestamps to avoid reindex errors
    out = out[~out.index.isna()]

    # Handle timezone alignment
    if tz is not None:
        if out.index.tz is None:
            out.index = out.index.tz_localize(tz)
        else:
            out.index = out.index.tz_convert(tz)

    # Build expected index
    expected_index = pd.date_range(
        start=pd.to_datetime(time_start),
        end=pd.to_datetime(time_end),
        freq=freq,
        inclusive=inclusive,
        tz=tz,
    )

    adapted = out.reindex(expected_index)

    return adapted
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compose paths based on a template and time
class PreserveMissing(dict):
    """A dict that returns {key} if the key is missing."""
    def __missing__(self, key):
        return "{" + key + "}"

@simplify_list
def compose_paths(
    path_template: str,
    path_time: (str, pd.Timestamp) = None,
    path_labels: List[str] = None,
    zero_pad: bool = True
) -> List[str]:
    """
    Compose a list of formatted paths from a template, preserving unfilled tags.

    Parameters
    ----------
    path_template : str
        Template containing placeholders: {yyyy}, {mm}, {dd}, {label}, etc.
    path_time : str | pd.Timestamp | None
        Date for formatting. If None, leaves date placeholders intact.
    path_labels : list of str, optional
        List of labels for {label}. Defaults to [""].
    zero_pad : bool, optional
        Whether to zero-pad month/day.

    Returns
    -------
    list of str
        Formatted paths with unfilled tags preserved.
    """
    if path_labels is None:
        path_labels = [""]

    # Start with empty dictionary but preserve placeholders
    values = PreserveMissing()

    if path_time is not None:
        if not isinstance(path_time, pd.Timestamp):
            path_time = pd.Timestamp(path_time)
        values["yyyy"] = path_time.year
        values["mm"] = f"{path_time.month:02d}" if zero_pad else str(path_time.month)
        values["dd"] = f"{path_time.day:02d}" if zero_pad else str(path_time.day)

    path_list = []
    for label in path_labels:
        if label != "":
            values["label"] = label
        else:
            values.pop("label", None)  # so {label} stays unfilled
        path_list.append(path_template.format_map(values))

    return path_list
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def _normalize_to_list(obj) -> List[pd.DataFrame]:
    """Turn input into a clean list of DataFrames, dropping None."""
    if obj is None:
        return []
    if isinstance(obj, pd.DataFrame):
        return [obj]
    if isinstance(obj, Iterable) and not isinstance(obj, (str, bytes)):
        cleaned = [x for x in obj if x is not None]
        # type check
        for i, x in enumerate(cleaned):
            if not isinstance(x, pd.DataFrame):
                raise TypeError(f"Item at index {i} is not a pandas DataFrame (got {type(x)}).")
        return cleaned
    raise TypeError("Input must be a DataFrame, a list/tuple of DataFrames, or None.")

def accept_df_or_list(func):
    """Decorator to allow a function to accept a single DF, a list of DFs, or None."""
    @wraps(func)
    def wrapper(obj, *args, **kwargs):
        dfs = _normalize_to_list(obj)
        return func(dfs, *args, **kwargs)
    return wrapper

@accept_df_or_list
def merge_by_rows(
    dfs: List[pd.DataFrame],
    ignore_index: bool = True,
    sort: bool = False,
    return_empty_df: bool = False,
    fill_value: dict = None,
    drop_rows_all_na: bool = False,
):
    """
    Row-bind DataFrames (like rbind).
    - Accepts DF | [DF, ...] | None.
    - Ignores None inside lists.
    - If nothing remains, returns None (or an empty DF if return_empty_df=True).
    - Mismatched columns are aligned like pandas.concat.
    """
    if not dfs:
        return pd.DataFrame() if return_empty_df else None

    out = pd.concat(dfs, ignore_index=ignore_index, sort=sort)

    if drop_rows_all_na:
        out = out.dropna(how='all')

    if fill_value is not None:
        # fill_value can be a scalar or a dict per-column
        out = out.fillna(fill_value)

    return out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to merge DataArrays by time
def merge_by_time(df_1, df_2, time_dim: str = "time") -> Optional[xr.DataArray]:
    if df_1 is None and df_2 is None:
        return None
    elif df_1 is None:
        return df_2
    elif df_2 is None:
        return df_1

    # Merge when both are present
    merged = xr.concat([df_1, df_2], dim=time_dim)
    merged = merged.groupby(time_dim).first()  # Remove duplicates
    merged = merged.sortby(time_dim)
    return merged
# ----------------------------------------------------------------------------------------------------------------------


