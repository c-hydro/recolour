"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20250813'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations
import logging
import pandas as pd
import xarray as xr
from typing import Optional, List, Dict

from lib_utils_decoretors import simplify_list, iterate_items

from lib_utils_info import logger_name

# set logger stream
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to filter and rename DataFrame columns
def filter_dataframe(
    df: pd.DataFrame,
    rename_map: Optional[Dict[str, str]] = None,
    filter_map: Optional[Dict[str, bool]] = None,
    *,
    keep_order: str = "filter"  # "filter" -> order by filter_map; "dataframe" -> original df order
) -> pd.DataFrame:
    """
    Filter columns using a boolean map and (optionally) rename them.

    - Only columns marked True in `filter_map` AND present in `df` are kept.
    - Missing columns in `df` are silently skipped (no errors).
    - If `rename_map` is provided, only kept columns are renamed.
    - Column order:
        * "filter": follows the order of keys in `filter_map` that are True and exist in df.
        * "dataframe": keeps the original df column order among the kept columns.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    rename_map : dict[str, str] | None
        Mapping from original column names -> new names (optional).
    filter_map : dict[str, bool] | None
        Mapping from original column names -> whether to keep (True) or drop (False).
        If None, all columns are considered kept.
    keep_order : {"filter", "dataframe"}
        How to order the resulting columns.

    Returns
    -------
    pd.DataFrame
        Filtered (and possibly renamed) DataFrame.
    """
    if filter_map is None:
        cols_to_keep = list(df.columns)
    else:
        # keep True entries that exist in df
        requested = [col for col, keep in filter_map.items() if keep]
        present = set(df.columns)
        cols_to_keep = [col for col in requested if col in present]

        if keep_order == "dataframe":
            cols_to_keep = [c for c in df.columns if c in cols_to_keep]

    out = df.loc[:, cols_to_keep].copy()

    if rename_map:
        # rename only the kept columns that appear in the map
        rename_subset = {k: v for k, v in rename_map.items() if k in out.columns}
        out.rename(columns=rename_subset, inplace=True)

    return out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to adapt a DataFrame to a specified time range
@iterate_items(strict_zip=True, dict_key_source='first')
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
# method to merge DataFrames by rows
def merge_by_rows(
    df1: Optional[pd.DataFrame],
    df2: Optional[pd.DataFrame],
    ignore_index: bool = True,
    sort: bool = False,
    return_empty_df: bool = False,
    fill_value: Union[dict, int, float, str, None] = None,
    drop_rows_all_na: bool = False,
) -> Optional[pd.DataFrame]:

    # Filter out None
    dfs = [df for df in (df1, df2) if df is not None]

    if not dfs:
        return pd.DataFrame() if return_empty_df else None
    if len(dfs) == 1:
        out = dfs[0].copy()
    else:
        out = pd.concat(dfs, ignore_index=ignore_index, sort=sort)

    if drop_rows_all_na:
        out = out.dropna(how='all')

    if fill_value is not None:
        out = out.fillna(fill_value)

    return out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to merge DataArrays by time
def merge_by_time(df1, df2, time_dim: str = "time") -> Optional[xr.DataArray]:
    if df1 is None and df2 is None:
        return None
    elif df1 is None:
        return df2
    elif df2 is None:
        return df1

    # Merge when both are present
    merged = xr.concat([df1, df2], dim=time_dim)
    merged = merged.groupby(time_dim).first()  # Remove duplicates
    merged = merged.sortby(time_dim)
    return merged
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to merge DataFrames by data
@iterate_items(strict_zip=True, dict_key_source='first')
def merge_by_data(df_1, df_2, no_data=-9999) -> pd.DataFrame:

    if df_1 is not None and df_2 is not None:
        if not isinstance(df_1, pd.DataFrame) or not isinstance(df_2, pd.DataFrame):
            log_stream.error("Both inputs must be pandas DataFrames.")
            raise NotImplementedError("Case not implemented: both inputs must be pandas DataFrames.")

    # check if the DataFrame is empty
    if df_1.empty and not df_2.empty:
        log_stream.warning(f" ===> Empty DataFrame 1")
        return df_2
    elif df_2.empty and not df_1.empty:
        log_stream.warning(f" ===> Empty DataFrame 2")
        return df_1
    elif df_1.empty and df_2.empty:
        log_stream.warning(f" ===> Both DataFrames are empty")
        return pd.DataFrame()

    # merge by time index (outer join to keep all dates)
    df_merged = pd.merge(df_1, df_2, left_index=True, right_index=True, how='outer')
    # fill NaN values (e.g., with 0)
    merged_df = df_merged.fillna(no_data)

    return merged_df
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to check if all dictionaries have the same keys
def check_dict_keys(dicts):
    if not dicts:
        log_stream.warning(" ===> No dictionaries provided.")
        return

    # Get the keys from the first dictionary
    reference_keys = set(dicts[0].keys())

    check_keys = True
    for i, d in enumerate(dicts[1:], start=1):
        current_keys = set(d.keys())
        if current_keys != reference_keys:
            missing = reference_keys - current_keys
            extra = current_keys - reference_keys

            log_stream.warning(f" ===> Dictionary {i} key mismatch:")
            if missing:
                log_stream.warning(f" ===> Missing keys: {missing}")
            if extra:
                log_stream.warning(f" ===> Extra keys: {extra}")
            check_keys = False

    return check_keys
# ----------------------------------------------------------------------------------------------------------------------
