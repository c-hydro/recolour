# ----------------------------------------------------------------------------------------------------------------------
# libraries
import numpy as np
import pandas as pd
import xarray as xr
from sklearn.neighbors import BallTree
from typing import Dict, Literal, Union, Optional

from lib_utils_decoretors import iterate_dict
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to resample values by period (daily, monthly, etc.)
@iterate_dict
def aggregate_values_ts_by_frequency(df, time_col='timestamp', value_cols=None, freq='D', skipna=True):
    """
    Sorts by time and computes mean values and counts for a given period.
    Works whether time is an index or a column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    time_col : str
        Name of the datetime column (if not already index).
    value_cols : str, list, or None
        Columns to average and count. If None, all numeric columns are used.
    period : str
        'hourly', 'daily', 'monthly' (case-insensitive).
    skipna : bool
        Whether to ignore NaNs in the mean calculation.

    Returns
    -------
    pd.DataFrame
        DataFrame with mean and count values per period, including the time column.
    """

    # get frequency
    freq = freq.lower()
    # copy DataFrame to avoid modifying the original
    df = df.copy()

    # Ensure datetime index
    if not isinstance(df.index, pd.DatetimeIndex):
        if time_col not in df.columns:
            raise ValueError(f"'{time_col}' not found in DataFrame columns.")
        df[time_col] = pd.to_datetime(df[time_col])
        df = df.sort_values(time_col).set_index(time_col)
    else:
        if not pd.api.types.is_datetime64_any_dtype(df.index):
            df.index = pd.to_datetime(df.index)
        df = df.sort_index()

    # Normalize value_cols to a list
    if isinstance(value_cols, str):
        value_cols = [value_cols]
    elif value_cols is None:
        value_cols = df.select_dtypes(include=np.number).columns.tolist()
    if not value_cols:
        raise ValueError("No numeric columns to aggregate.")

    # Build aggregation dictionary
    agg_funcs = {}
    for col in value_cols:
        if skipna:
            agg_funcs[col] = [
                ('mean', lambda s: s.mean()),
                ('count', lambda s: s.count())
            ]
        else:
            agg_funcs[col] = [
                ('mean', lambda s: s.mean(skipna=False)),
                ('count', lambda s: s.notna().sum())
            ]

    # Perform resample + aggregation
    resampled = df[value_cols].resample(freq).agg(agg_funcs)

    # Flatten MultiIndex columns
    resampled.columns = [f"{col}_{stat}" for col, stat in resampled.columns]
    resampled = resampled.reset_index()

    # Restore time column name
    resampled.rename(columns={resampled.columns[0]: time_col}, inplace=True)

    # Set as index
    resampled = resampled.set_index(time_col)
    # (Optional) sort by time
    resampled = resampled.sort_index()

    return resampled
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to search values by point within a specified radius
@iterate_dict
def search_values_points(
    target_df: pd.DataFrame,
    source_df: pd.DataFrame,
    radius_m: float = 12_500,
    target_lat_col: str = "lat",
    target_lon_col: str = "lon",
    source_lat_col: str = "lat",
    source_lon_col: str = "lon",      # column in target_df to use as the dict key; if None uses target_df.index
    include_distance: bool = True,
    time_col: str = "time",
    distance_col: str = "distance_m",
) -> pd.DataFrame:
    """
    For each row in target_df, collect all rows from source_df within `radius_m` meters.

    Returns a dict mapping target key -> DataFrame of neighbors.

    Each DataFrame:
      - Contains the selected rows from source_df (preserving original columns).
      - Optionally adds distance in meters.
      - Stores the **full original target row** in `.attrs["target_row"]` for metadata.
    """
    # constants
    R = 6_371_000.0  # Earth radius (m)
    r_rad = radius_m / R

    # coords in radians (BallTree + haversine expects [lat, lon] in radians)
    src_coords = np.deg2rad(source_df[[source_lat_col, source_lon_col]].to_numpy())
    tgt_coords = np.deg2rad(target_df[[target_lat_col, target_lon_col]].to_numpy())

    # tree + query
    tree = BallTree(src_coords, metric="haversine")
    idx_lists, dist_lists = tree.query_radius(
        tgt_coords, r=r_rad, return_distance=True, sort_results=True
    )

    # build a tidy result
    rows = []
    for t_i, (src_idx, dists) in enumerate(zip(idx_lists, dist_lists)):
        if src_idx.size == 0:
            continue  # no neighbors for this target
        block = source_df.iloc[src_idx].copy()
        block[time_col] = source_df.index[src_idx]

        if include_distance:
            block[distance_col] = dists * R  # meters
        rows.append(block)

    df_neighbors = pd.concat(rows, ignore_index=True) if rows else None

    # store original target row(s) as attrs if you really want (optional; can be large)
    if df_neighbors is not None:
        df_neighbors.attrs = target_df.to_dict()

    # Set as index
    df_neighbors = df_neighbors.set_index(time_col)
    # (Optional) sort by time
    df_neighbors = df_neighbors.sort_index()

    return df_neighbors
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to search values in a DataArray map within a specified radius
@iterate_dict
def search_values_maps(
    target_df: pd.DataFrame,
    source_da: xr.DataArray,              # dims: (time, lat, lon) by default
    radius_m: float = 12_500,
    target_lat_col: str = "lat",
    target_lon_col: str = "lon",
    time_name: str = "time",
    lat_name: str = "lat",
    lon_name: str = "lon",
    include_distance: bool = True,
    distance_col: str = "distance_m",
    value_col: Optional[str] = None,      # column name for values; defaults to source_da.name or "value"
    return_mode: Literal["stacked", "dict"] = "stacked",
    target_id_col: Optional[str] = None,  # id column in target_df; defaults to index
) -> Union[pd.DataFrame, Dict[object, pd.DataFrame]]:
    """
    For each row in target_df, collect all grid cells from source_da within `radius_m` meters
    around (lat, lon) and include values for **all times**.

    Returns
    -------
    'stacked' : single tidy DataFrame indexed by time with columns:
                [lat, lon, value, target_id, (distance_m)]
    'dict'    : {target_id: DataFrame indexed by time}
    """
    # Checks
    if any(d not in source_da.dims for d in (time_name, lat_name, lon_name)) or source_da.ndim != 3:
        raise ValueError(f"source_da must be 3-D with dims ({time_name}, {lat_name}, {lon_name})")
    if target_df.empty:
        return {} if return_mode == "dict" else pd.DataFrame()

    # Names / constants
    if value_col is None:
        value_col = source_da.name or "value"
    R = 6_371_000.0
    r_rad = float(radius_m) / R

    # IDs
    target_ids = target_df.index.to_numpy() if target_id_col is None else target_df[target_id_col].to_numpy()

    # Coordinates (regular 1-D lat/lon expected)
    times = pd.DatetimeIndex(source_da[time_name].values)
    lats = np.asarray(source_da[lat_name].values)
    lons = np.asarray(source_da[lon_name].values)
    if lats.ndim != 1 or lons.ndim != 1:
        raise ValueError("Expect 1-D lat and lon coords (regular grid).")

    lat2d, lon2d = np.meshgrid(lats, lons, indexing="ij")  # (nlat, nlon)

    # Normalize target longitudes to grid convention
    lon_min, lon_max = float(lons.min()), float(lons.max())
    grid_is_0360 = (lon_min >= 0.0) and (lon_max <= 360.0)
    def _normalize_lon(lon_val: float) -> float:
        if grid_is_0360 and lon_val < 0:
            return lon_val % 360.0
        if (not grid_is_0360) and lon_val > 180:
            return ((lon_val + 180) % 360) - 180
        return lon_val

    # Spatial index on grid cell centers (radians, order: [lat, lon])
    grid_coords_rad = np.deg2rad(np.c_[lat2d.ravel(), lon2d.ravel()])
    tree = BallTree(grid_coords_rad, metric="haversine")

    if return_mode == "dict":
        out: Dict[object, pd.DataFrame] = {}

    stacked_blocks = []
    nlat, nlon = lat2d.shape
    T = times.size

    for t_id, (_, row) in zip(target_ids, target_df.iterrows()):
        tlat = float(row[target_lat_col])
        tlon = _normalize_lon(float(row[target_lon_col]))

        # Query neighbors once per target
        tgt_rad = np.deg2rad([[tlat, tlon]])
        idx_lists, dist_lists = tree.query_radius(tgt_rad, r=r_rad, return_distance=True, sort_results=True)
        src_idx = idx_lists[0]
        dists = dist_lists[0]  # radians

        if src_idx.size == 0:
            if return_mode == "dict":
                out[t_id] = pd.DataFrame().set_index(pd.Index([], name=time_name))
            continue

        # Map flat -> (ilat, ilon)
        ilat, ilon = np.unravel_index(src_idx, (nlat, nlon))

        # Slice all times at those (lat, lon) points: result dims (time, points)
        # Use xarray vectorized indexing for efficiency
        da_sub = source_da.isel({lat_name: xr.DataArray(ilat), lon_name: xr.DataArray(ilon)})  # (time, points)
        vals = np.asarray(da_sub.values)                    # shape (T, P)
        P = vals.shape[1]

        # Prepare tidy columns
        time_col_vals = np.tile(times.values, P)            # T repeated for each point
        lat_col_vals  = np.repeat(lat2d[ilat, ilon], T)     # each point repeated across all times
        lon_col_vals  = np.repeat(lon2d[ilat, ilon], T)
        val_col_vals  = vals.reshape(T * P, order="F")      # column-major to align with tiling
        tgt_ids_vals  = np.repeat(t_id, T * P)
        if include_distance:
            dist_m_vals = np.repeat(dists * R, T)           # meters

        block = pd.DataFrame({
            time_name: pd.to_datetime(time_col_vals),
            lat_name:  lat_col_vals,
            lon_name:  lon_col_vals,
            value_col: val_col_vals,
            "target_id": tgt_ids_vals,
        })
        if include_distance:
            block[distance_col] = dist_m_vals

        block = block.set_index(time_name).sort_index()

        if return_mode == "dict":
            out[t_id] = block
        else:
            stacked_blocks.append(block)

    if return_mode == "dict":
        return out

    if not stacked_blocks:
        cols = [lat_name, lon_name, value_col, "target_id"] + ([distance_col] if include_distance else [])
        return pd.DataFrame(columns=cols).set_index(pd.Index([], name=time_name))

    df = pd.concat(stacked_blocks, axis=0, ignore_index=False)
    df = df.sort_index()

    # remove useless columns
    df = df.drop(columns=['target_id'], errors='ignore')

    return df

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to aggregate values maps using frequency
@iterate_dict
def aggregate_values_maps_by_frequency(df, time_col='time', value_col='ssm_filtered', freq="D", min_frac=0.75):
    """
    Resample to `freq` and compute mean/max/min only when the fraction of finite
    values in the period is >= min_frac. Also return counts and coverage.

    Parameters
    ----------
    df : DataFrame
        Input data for a single station or series.
    time_col : str
        Timestamp column name.
    value_col : str, default 'ssm_filtered'
        Column with numeric values to aggregate.
    freq : str, default "D"
        Resample frequency (e.g., "D" for daily).
    min_frac : float, default 0.75
        Minimum fraction of finite values required for a valid aggregate.

    Returns
    -------
    DataFrame
        Index: resample period
        Columns:
            {value_col}_mean, {value_col}_max, {value_col}_min  (masked if below threshold)
            n_total  (rows in period)
            n_used   (finite values used)
            pct_used (100 * n_used / n_total)
    """
    # Handle empty input
    if df is None or len(df) == 0:
        cols = [f"{value_col}_mean", f"{value_col}_max", f"{value_col}_min",
                "n_total", "n_used", "pct_used"]
        return pd.DataFrame(columns=cols)

    # Normalize types and index
    g = (
        df.copy()
        .reset_index(drop=False)
        .assign(**{
            time_col: lambda d: pd.to_datetime(d[time_col]),
            value_col: lambda d: pd.to_numeric(d[value_col], errors='coerce')
        })
        .set_index(time_col)
        .sort_index()
    )

    # Counts per period
    n_total = g[value_col].resample(freq).size().rename("n_total")
    n_used = g[value_col].apply(np.isfinite).resample(freq).sum().astype("Int64").rename("n_used")

    # Coverage fraction/percent
    frac = n_used / n_total.replace(0, np.nan)
    pct_used = (frac * 100.0).rename("pct_used")

    # Aggregates
    agg = g[value_col].resample(freq).agg(["mean", "max", "min"])
    agg = agg.rename(columns={
        "mean": f"{value_col}_mean",
        "max":  f"{value_col}_max",
        "min":  f"{value_col}_min",
    })

    # Mask aggregates if below threshold, but keep diagnostics
    agg = agg.where(frac >= min_frac)

    # Combine
    out = pd.concat([agg, n_total, n_used, pct_used], axis=1)

    return out
# ----------------------------------------------------------------------------------------------------------------------