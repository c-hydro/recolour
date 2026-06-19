"""
Library Features:

Name:           lib_utils_io
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json
import rasterio
import csv

from pathlib import Path
import pandas as pd
import numpy as np

from datetime import datetime, timedelta
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from dataclasses import dataclass

from config_info import LOGGER_NAME, VALUE_NODATA_DEFAULT
from lib_utils_base import format_file_path
from lib_utils_time import format_time_for_csv

# logger stream
logger_stream = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# point dataclass
@dataclass
class PointValue:
    tag: str
    lon: float
    lat: float
    value: float
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to collect data
def collect_data(
    path_raw: str,
    registry: Dict[str, Dict[str, Any]],
    time_now: datetime, params: Dict[str, Any],
) -> (List[PointValue], List, List):

    points: List[PointValue] = []
    times: List[datetime] = []
    lags: List[int] = []

    # iterate over registry
    for tag in registry.keys():

        # info data start
        logger_stream.info(f' -----> Get data for point {tag} ... ')

        # define path
        path_def = format_file_path(file_pattern=path_raw, time_now=time_now, point_tag=tag)

        # check if file exists
        if not os.path.exists(path_def):
            logger_stream.warning(f" ===> File not found: {path_def}")

            # info data end (no file)
            logger_stream.info(f' -----> Get data for point {tag} ... SKIPPED.')
            continue

        # read data from file
        point, time_selection, time_lag = read_file_data(path_def, tag, registry, params, time_now)

        # collect point, lags and times
        if point is not None:
            points.append(point)
            lags.append(time_lag)
            times.append(time_selection)

            # info data end
            logger_stream.info(f' -----> Get data for point {tag} ... DONE')

        else:
            # info data end (no data)
            logger_stream.info(f' -----> Get data for point {tag} ... SKIPPED. NO DATA')

    # check valid point(s)
    if not points:
        logger_stream.error(' ===> Point data are not defined')
        raise RuntimeError("No valid point values found")

    return points, times, lags
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read file data
def read_file_data(
    path: str,
    tag: str,
    registry: Dict[str, Dict[str, Any]],
    params: Dict[str, Any],
    time_ref: datetime,
) -> Optional[Tuple[PointValue, datetime, float]]:

    # get params
    time_format = params.get("time_format")
    step_hours = params.get("step_hours")
    max_previous_steps = params.get("max_previous_steps")
    nodata_values = params.get("nodata_values")
    valid_min = params.get("valid_min")
    valid_max = params.get("valid_max")
    delimiter = params.get("delimiter", ",")
    time_col = params.get("time_col", "time")
    value_col = params.get("value_col", "theta_simulated")
    lon_col = params.get("registry_lon_col", "longitude")
    lat_col = params.get("registry_lat_col", "latitude")

    # check path in registry
    if tag not in registry:
        logger_stream.warning(f" ===> Point {tag}: not found in registry, skipped")
        return None

    # read data
    df = pd.read_csv(path, sep=delimiter)

    # check expected cols
    if time_col not in df.columns:
        logger_stream.error(f' ===> Column {time_col} is mandatory. Not found.')
        raise KeyError(f"Column '{time_col}' not found in {path}. Available: {list(df.columns)}")
    if value_col not in df.columns:
        logger_stream.error(f' ===> Column {value_col} is mandatory. Not found.')
        raise KeyError(f"Column '{value_col}' not found in {path}. Available: {list(df.columns)}")

    # iterate to select the time
    selected_time, selected_value = None, None
    for step in range(max_previous_steps + 1):

        candidate_time = time_ref - timedelta(hours=step * step_hours)
        target_time = format_time_for_csv(candidate_time, time_format)

        row = df.loc[df[time_col].astype(str) == target_time]

        if row.empty:
            continue

        value = float(row.iloc[0][value_col])

        if any(np.isclose(value, float(v)) for v in nodata_values):
            continue
        if valid_min is not None and value < float(valid_min):
            continue
        if valid_max is not None and value > float(valid_max):
            continue
        if not np.isfinite(value):
            continue

        selected_time = candidate_time
        selected_value = value
        break

    # check selected value
    if selected_value is None:
        logger_stream.warning(f" ===> Point {tag}: no valid value found within {max_previous_steps} previous steps")
        return None

    # compute lag hours
    lag_hours = (time_ref - selected_time).total_seconds() / 3600.0

    # define point obj
    rec = registry[tag]
    point = PointValue(tag=tag,lon=float(rec[lon_col]),lat=float(rec[lat_col]),value=selected_value,)

    return point, selected_time, lag_hours
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read settings file
def read_file_settings(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        logger_stream.error(f' ===> File {file_name} not found')
        raise FileNotFoundError(f'File "{file_name}" not found. Exit')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read file registry
def read_file_registry(path: str, params: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:

    # get params
    rows: Dict[str, Dict[str, Any]] = {}
    tag_col = params.get("tag_col", "tag")
    lon_col = params.get("lon_col", "longitude")
    lat_col = params.get("lat_col", "latitude")
    valid_col = params.get("valid_col", "valid")
    only_valid = params.get("use_only_valid", True)
    encoding = params.get("encoding", "utf-8-sig")
    delimiter = params['delimiter']

    if delimiter == "auto":
        first = Path(path).read_text(encoding=cfg.get("encoding", "utf-8-sig")).splitlines()[0]
        delimiter = max([",", ";", "\t"], key=lambda d: first.count(d))

    # open csv
    with open(path, newline="", encoding=encoding) as fp:
        reader = csv.DictReader(fp, delimiter=delimiter)
        for row in reader:
            clean = {str(k).strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k is not None}
            if not clean:
                continue
            if only_valid and valid_col in clean:
                try:
                    if int(float(clean[valid_col])) != 1:
                        continue
                except Exception:
                    continue
            tag = str(clean[tag_col]).strip()
            clean[lon_col] = float(clean[lon_col])
            clean[lat_col] = float(clean[lat_col])
            rows[tag] = clean

    # check points
    if not rows:
        logger_stream.error(f' ===> Regitry points are not defined')
        raise RuntimeError(f"No valid registry points found in {path}")

    return rows
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to write geo tiff file
def write_file_geotiff(
        path: str, array: np.ndarray, ref_profile: Dict[str, Any], nodata: float,
        compress: str = "deflate") -> None:

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    profile = ref_profile.copy()
    profile.update(
        driver="GTiff",
        count=1,
        dtype="float32",
        nodata=float(nodata),
        compress=compress,
        tiled=False,
    )
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype(np.float32), 1)
# ----------------------------------------------------------------------------------------------------------------------
