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

from typing import Any, Dict

import numpy as np
import pandas as pd
from netCDF4 import Dataset, chartostring, num2date

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
    format_data: str = 'csv'
) -> (List[PointValue], List, List):

    points: List[PointValue] = []
    times: List[datetime] = []
    lags: List[int] = []

    # read data from file
    if format_data == 'csv':

        # iterate over registry
        for tag in registry.keys():

            # info data start
            logger_stream.info(f' -----> Get point from file {tag} ... ')

            # define path
            path_def = format_file_path(file_pattern=path_raw, time_now=time_now, point_tag=tag)

            # check if file exists
            if not os.path.exists(path_def):
                logger_stream.warning(f" ===> File not found: {path_def}")

                # info data end (no file)
                logger_stream.info(f' -----> Get point from file {tag} ... SKIPPED.')
                continue

            point, time_selection, time_lag = read_file_point(path_def, tag, registry, params, time_now)

            # collect point, lags and times
            if point is not None:
                points.append(point)
                lags.append(time_lag)
                times.append(time_selection)

                # info data end
                logger_stream.info(f' -----> Get point from file {tag} ... DONE')

            else:
                # info data end (no data)
                logger_stream.info(f' -----> Get point from file {tag} ... SKIPPED. NO DATA')

    elif format_data == 'netcdf':

        # info data start
        logger_stream.info(f' -----> Read data collections ... ')

        # define path
        path_def = format_file_path(file_pattern=path_raw, time_now=time_now, point_tag='collections')

        # check if file exists
        if not os.path.exists(path_def):
            logger_stream.warning(f" ===> File not found: {path_def}")

        collections = read_file_collections(path_def, params)

        # info data start
        logger_stream.info(f' -----> Read data collections ... DONE')

        points, lags, times = [], [], []
        for tag in registry.keys():

            # info data start
            logger_stream.info(f' -----> Get point from collections {tag} ... ')

            point, time_selection, time_lag = read_collection_point(
                collections=collections,
                tag=tag,
                registry=registry,
                params=params,
                time_ref=time_now
            )

            if point is not None:
                points.append(point)
                lags.append(time_lag)
                times.append(time_selection)

                # info data end
                logger_stream.info(f' -----> Get point from collections {tag} ... DONE')

            else:

                # info data end (no data)
                logger_stream.info(f' -----> Get point from collections {tag} ... SKIPPED. NO DATA')

    else:
        raise ValueError(f"Format data are not supported: {format_data} ")

    # check valid point(s)
    if not points:
        logger_stream.error(' ===> Point data are not defined')
        raise RuntimeError("No valid point values found")

    return points, times, lags
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
def read_collection_point(
    collections: Dict[str, Any],
    tag: str,
    registry: Dict[str, Dict[str, Any]],
    params: Dict[str, Any],
    time_ref: datetime,
) -> Tuple[
    Optional[PointValue],
    Optional[datetime],
    Optional[float]
]:

    # get params
    step_hours = params.get("step_hours")
    max_previous_steps = params.get("max_previous_steps")
    nodata_values = params.get("nodata_values", [])
    valid_min = params.get("valid_min")
    valid_max = params.get("valid_max")

    time_col = params.get("time_col", "time")
    value_col = params.get("value_col", "theta_simulated")
    point_col = params.get("point_name_var", "point_name")

    lon_col = params.get(
        "registry_lon_col",
        "longitude"
    )

    lat_col = params.get(
        "registry_lat_col",
        "latitude"
    )

    # check registry
    if tag not in registry:
        logger_stream.warning(
            f" ===> Point {tag}: not found in registry, skipped"
        )
        return None, None, None

    # check collection variables
    required_variables = [
        time_col,
        point_col,
        value_col
    ]

    missing_variables = [
        variable_name
        for variable_name in required_variables
        if variable_name not in collections
    ]

    if missing_variables:
        raise KeyError(
            f"Missing collection variables: {missing_variables}. "
            f"Available variables: {list(collections.keys())}"
        )

    # get collection data
    time_values = pd.DatetimeIndex(
        pd.to_datetime(
            collections[time_col]
        )
    )

    point_names = np.asarray(
        collections[point_col]
    ).astype(str)

    value_values = np.asarray(
        collections[value_col],
        dtype=np.float64
    )

    # values must be organized as [time, point]
    if value_values.ndim != 2:
        raise RuntimeError(
            f'Variable "{value_col}" must be two-dimensional '
            f"with shape (time, point). "
            f"Current shape: {value_values.shape}"
        )

    if value_values.shape[0] != time_values.size:
        raise RuntimeError(
            f'Time dimension mismatch for "{value_col}": '
            f"{value_values.shape[0]} data steps versus "
            f"{time_values.size} time steps"
        )

    if value_values.shape[1] != point_names.size:
        raise RuntimeError(
            f'Point dimension mismatch for "{value_col}": '
            f"{value_values.shape[1]} data points versus "
            f"{point_names.size} point names"
        )

    # search point
    point_indexes = np.where(
        point_names == str(tag)
    )[0]

    if point_indexes.size == 0:
        logger_stream.warning(
            f' ===> Point "{tag}" not found in collection'
        )
        return None, None, None

    if point_indexes.size > 1:
        logger_stream.warning(
            f' ===> Point "{tag}" occurs {point_indexes.size} times '
            f"in collection; using the first occurrence"
        )

    point_index = int(
        point_indexes[0]
    )

    point_series = value_values[
        :,
        point_index
    ]

    # iterate backward to select time
    selected_time = None
    selected_value = None

    for step in range(max_previous_steps + 1):

        candidate_time = (
            pd.Timestamp(time_ref) -
            pd.Timedelta(
                hours=step * step_hours
            )
        )

        time_indexes = np.where(
            time_values == candidate_time
        )[0]

        if time_indexes.size == 0:
            continue

        value = float(
            point_series[
                int(time_indexes[0])
            ]
        )

        if not np.isfinite(value):
            continue

        if any(
                np.isclose(
                    value,
                    float(nodata_value)
                )
                for nodata_value in nodata_values
        ):
            continue

        if (
                valid_min is not None and
                value < float(valid_min)
        ):
            continue

        if (
                valid_max is not None and
                value > float(valid_max)
        ):
            continue

        selected_time = candidate_time.to_pydatetime()
        selected_value = value
        break

    # no valid value found
    if selected_value is None:
        logger_stream.warning(
            f" ===> Point {tag}: no valid value found within "
            f"{max_previous_steps} previous steps"
        )
        return None, None, None

    # compute lag
    lag_hours = (
        pd.Timestamp(time_ref) -
        pd.Timestamp(selected_time)
    ).total_seconds() / 3600.0

    # create point
    registry_record = registry[tag]

    point = PointValue(
        tag=tag,
        lon=float(
            registry_record[lon_col]
        ),
        lat=float(
            registry_record[lat_col]
        ),
        value=selected_value,
    )

    return point, selected_time, lag_hours

# helper to read file collections
def read_file_collections(
        path: str,
        params: Dict[str, Any]
) -> Dict[str, Any]:

    """
    Read the collection variables from a NetCDF file.

    Variables read:
        - time
        - point_name
        - theta_simulated

    Returns
    -------
    Dict[str, Any]

        {
            "time": pandas.DatetimeIndex,
            "point_name": numpy.ndarray,
            "theta_simulated": numpy.ndarray
        }

    The theta_simulated array is always organized as:

        theta_simulated[time, point]
    """

    # get variable names
    time_var_name = params.get(
        "time_col",
        "time"
    )

    point_var_name = params.get(
        "point_name_var",
        "point_name"
    )

    value_var_name = params.get(
        "value_col",
        "theta_simulated"
    )

    point_dim_name = params.get(
        "point_dim_name",
        "point"
    )

    time_dim_name = params.get(
        "time_dim_name",
        "time"
    )

    # open NetCDF file
    with Dataset(path, mode="r") as file_handle:

        # ---------------------------------------------------------------------
        # Check variables

        required_variables = [
            time_var_name,
            point_var_name,
            value_var_name
        ]

        missing_variables = [
            variable_name
            for variable_name in required_variables
            if variable_name not in file_handle.variables
        ]

        if missing_variables:
            raise KeyError(
                f'Missing variables in NetCDF file "{path}": '
                f"{missing_variables}. Available variables are: "
                f"{list(file_handle.variables.keys())}"
            )

        # ---------------------------------------------------------------------
        # Read and convert time

        time_variable = file_handle.variables[
            time_var_name
        ]

        if not hasattr(time_variable, "units"):
            raise RuntimeError(
                f'Time variable "{time_var_name}" does not define '
                f'the mandatory "units" attribute'
            )

        time_values_raw = time_variable[:]

        if np.ma.isMaskedArray(time_values_raw):
            time_values_raw = time_values_raw.filled(
                np.nan
            )

        time_values_decoded = num2date(
            time_values_raw,
            units=time_variable.units,
            calendar=getattr(
                time_variable,
                "calendar",
                "standard"
            ),
            only_use_python_datetimes=True,
            only_use_cftime_datetimes=False
        )

        time_values = pd.DatetimeIndex(
            pd.to_datetime(
                np.asarray(time_values_decoded)
            )
        )

        # ---------------------------------------------------------------------
        # Read point names

        point_variable = file_handle.variables[
            point_var_name
        ]

        point_values_raw = point_variable[:]

        if np.ma.isMaskedArray(point_values_raw):

            if point_values_raw.dtype.kind in [
                    "S", "U", "O"]:

                point_values_raw = point_values_raw.filled(
                    ""
                )

            else:
                point_values_raw = point_values_raw.filled(
                    np.nan
                )



        # Classic NetCDF character array:
        #
        # point_name(point, string_length)
        if (
                point_values_raw.ndim == 2 and
                point_values_raw.dtype.kind in ["S", "U"]):

            point_values = chartostring(
                point_values_raw
            )

        # NetCDF VLEN string array:
        #
        # point_name(point)
        elif (
                point_values_raw.ndim == 1 and
                point_values_raw.dtype.kind in ["S", "U", "O"]):

            point_values = np.asarray([
                value.decode(
                    "utf-8",
                    errors="ignore"
                )
                if isinstance(value, bytes)
                else str(value)
                for value in point_values_raw
            ])

        # Numeric point identifiers
        elif (
                point_values_raw.ndim == 1 and
                point_values_raw.dtype.kind in ["i", "u", "f"]):

            point_values = np.asarray([
                str(int(value))
                if (
                    np.isfinite(value) and
                    float(value).is_integer()
                )
                else str(value)
                for value in point_values_raw
            ])

        else:
            raise RuntimeError(
                f'Unsupported point variable "{point_var_name}": '
                f"shape={point_values_raw.shape}, "
                f"dtype={point_values_raw.dtype}"
            )

        point_values = np.asarray([
            str(value).strip().strip("\x00")
            for value in point_values
        ])

        # ---------------------------------------------------------------------
        # Read theta_simulated

        value_variable = file_handle.variables[
            value_var_name
        ]

        value_dimensions = list(
            value_variable.dimensions
        )

        if time_dim_name not in value_dimensions:
            raise RuntimeError(
                f'Variable "{value_var_name}" does not use time '
                f'dimension "{time_dim_name}". '
                f"Dimensions are {value_dimensions}"
            )

        if point_dim_name not in value_dimensions:
            raise RuntimeError(
                f'Variable "{value_var_name}" does not use point '
                f'dimension "{point_dim_name}". '
                f"Dimensions are {value_dimensions}"
            )

        time_axis = value_dimensions.index(
            time_dim_name
        )

        point_axis = value_dimensions.index(
            point_dim_name
        )

        value_values = value_variable[:]

        logger_stream.info(
            " -----> Raw variable: shape=%s dtype=%s min=%s max=%s",
            value_values.shape,
            value_values.dtype,
            np.nanmin(np.asarray(value_values, dtype=np.float64)),
            np.nanmax(np.asarray(value_values, dtype=np.float64)),
        )

        if np.ma.isMaskedArray(value_values):
            value_values = value_values.filled(
                np.nan
            )

        logger_stream.info(
            " -----> After fill: min=%s max=%s",
            np.nanmin(value_values),
            np.nanmax(value_values),
        )

        value_values = np.asarray(
            value_values,
            dtype=np.float64
        )

        logger_stream.info(
            " -----> After float conversion: min=%s max=%s",
            np.nanmin(value_values),
            np.nanmax(value_values),
        )

        # Organize variable as [time, point]
        value_values = np.moveaxis(
            value_values,
            [time_axis, point_axis],
            [0, 1]
        )

        logger_stream.info(
            " -----> After moveaxis: shape=%s min=%s max=%s",
            value_values.shape,
            np.nanmin(value_values),
            np.nanmax(value_values),
        )

        # Allow only singleton additional dimensions
        if value_values.ndim > 2:

            extra_shape = value_values.shape[2:]

            if any(size != 1 for size in extra_shape):
                raise RuntimeError(
                    f'Variable "{value_var_name}" has unsupported '
                    f"additional dimensions. Shape after reordering: "
                    f"{value_values.shape}"
                )

            value_values = value_values.reshape(
                value_values.shape[0],
                value_values.shape[1]
            )

        if value_values.ndim != 2:
            raise RuntimeError(
                f'Variable "{value_var_name}" cannot be organized as '
                f"(time, point). Shape is {value_values.shape}"
            )

        # ---------------------------------------------------------------------
        # Check dimensions

        if value_values.shape[0] != time_values.size:
            raise RuntimeError(
                f'Time size mismatch: variable "{time_var_name}" has '
                f"{time_values.size} elements, variable "
                f'"{value_var_name}" has {value_values.shape[0]} steps'
            )

        if value_values.shape[1] != point_values.size:
            raise RuntimeError(
                f'Point size mismatch: variable "{point_var_name}" has '
                f"{point_values.size} elements, variable "
                f'"{value_var_name}" has {value_values.shape[1]} points'
            )

        # ---------------------------------------------------------------------
        # Convert explicit no-data values to NaN

        nodata_values = []

        if hasattr(value_variable, "_FillValue"):
            fill_value = float(value_variable._FillValue)

            logger_stream.info(
                " -----> _FillValue = %s", fill_value
            )

            before = np.count_nonzero(np.isnan(value_values))

            value_values[
                np.isclose(value_values, fill_value)
            ] = np.nan

            after = np.count_nonzero(np.isnan(value_values))

            logger_stream.info(
                " -----> NaNs after _FillValue replacement: %d -> %d",
                before,
                after
            )

            logger_stream.info(
                " -----> Range after _FillValue: min=%s max=%s",
                np.nanmin(value_values),
                np.nanmax(value_values),
            )

        if hasattr(value_variable, "missing_value"):
            nodata_values.extend(
                np.atleast_1d(
                    value_variable.missing_value
                ).tolist()
            )

        for nodata_value in nodata_values:

            try:
                nodata_value = float(
                    nodata_value
                )
            except (TypeError, ValueError):
                continue

            value_values[
                np.isclose(
                    value_values,
                    nodata_value,
                    equal_nan=False
                )
            ] = np.nan

    # return collection variables
    collections = {
        time_var_name: time_values,
        point_var_name: point_values,
        value_var_name: value_values
    }

    logger_stream.info(
        " -----> Final array:"
    )
    logger_stream.info(
        "        shape      : %s", value_values.shape
    )
    logger_stream.info(
        "        finite     : %d", np.isfinite(value_values).sum()
    )
    logger_stream.info(
        "        nan        : %d", np.isnan(value_values).sum()
    )

    if np.isfinite(value_values).any():
        logger_stream.info(
            "        min/max    : %s / %s",
            np.nanmin(value_values),
            np.nanmax(value_values),
        )
    else:
        logger_stream.warning(
            " ===> Variable '%s' contains only NaN values",
            value_var_name
        )

    return collections
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read file point
def read_file_point(
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
