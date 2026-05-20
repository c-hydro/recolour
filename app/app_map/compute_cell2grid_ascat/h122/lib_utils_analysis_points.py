"""
Library Features:

Name:           lib_utils_analysis_points
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""


# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd
import xarray as xr

from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to check 1d arrays
def check_1d_locations_var(dataset: xr.Dataset, var_name: (str, list) = None):

    # check variable name
    if var_name is None:
        raise KeyError('Variable name must be defined. Found NoneType.')

    # allow single string or list/tuple of candidate names
    if isinstance(var_name, str):
        var_names = [var_name]
    else:
        var_names = list(var_name)

    # find first existing variable/coord
    selected_name = None
    for name in var_names:
        if name in dataset.variables or name in dataset.coords:
            selected_name = name
            break
    if selected_name is None:
        raise KeyError(f"None of the variables {var_names} found")

    # get data array
    data_array = dataset[selected_name]

    # check dimensions
    valid_dims = [("locations",), ("obs",)]
    if data_array.dims not in valid_dims:
        raise ValueError(
            f"Variable '{selected_name}' must have dims "
            f"{valid_dims}, got {data_array.dims}"
        )

    return selected_name
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# point collection and dedup
def collect_points_to_dataframe(settings, file_list):
    source_settings = settings.get("source", {})
    value_var = source_settings.get("value_col", "surface_soil_moisture")

    frames = []

    for file_path in file_list:
        logger.info(f" ----> Read source file: {file_path}")

        logger.info(f" ::: GET POINTS ... ")
        with xr.open_dataset(file_path) as dataset:
            var_name_geo_x = check_1d_locations_var(dataset, ["lon", "longitude", "x", "X"])
            logger.info(f" :::: Variable GeoX: {var_name_geo_x}")
            var_name_geo_y = check_1d_locations_var(dataset, ["lat", "latitude", "y", "Y"])
            logger.info(f" :::: Variable GeoY: {var_name_geo_y}")
            var_name_loc = check_1d_locations_var(dataset, ["gpi", "location_id"])
            logger.info(f" :::: Variable GPI/Locations: {var_name_loc}")
            var_name_time = check_1d_locations_var(dataset, "time")
            logger.info(f" :::: Variable Time: {var_name_time}")
            var_name_data = check_1d_locations_var(dataset, value_var)
            logger.info(f" :::: Variable Data: {var_name_data}")

            lon = np.asarray(dataset[var_name_geo_x].values, dtype=np.float64)
            lat = np.asarray(dataset[var_name_geo_y].values, dtype=np.float64)
            gpi = np.asarray(dataset[var_name_loc].values, dtype=np.int64)
            obs_time = pd.to_datetime(dataset[var_name_time].values)
            obs_value = np.asarray(dataset[var_name_data].values, dtype=np.float64)

            if not (lon.shape == lat.shape == gpi.shape == obs_time.shape == obs_value.shape):
                raise ValueError(
                    f'{file_path}: shape mismatch among lon/lat/gpi/time/{value_var}: '
                    f'{lon.shape}, {lat.shape}, {gpi.shape}, {obs_time.shape}, {obs_value.shape}'
                )

            valid_mask = (
                np.isfinite(lon) &
                np.isfinite(lat) &
                np.isfinite(obs_value)
            )

            if valid_mask.any():
                frame = pd.DataFrame(
                    {
                        "gpi": gpi[valid_mask],
                        "lon": lon[valid_mask],
                        "lat": lat[valid_mask],
                        "time": obs_time[valid_mask],
                        value_var: obs_value[valid_mask],
                        "source_file": os.path.basename(file_path)
                    }
                )
                frames.append(frame)

        logger.info(f" ::: GET POINTS ... DONE")

    # manage merge points
    logger.info(f" ::: MERGE POINTS ... ")
    if not frames:
        merge_frames = pd.DataFrame(columns=["gpi", "lon", "lat", "time", value_var, "source_file"])
    else:
        merge_frames = pd.concat(frames, ignore_index=True)
    logger.info(f" ::: MERGE POINTS ... DONE")

    return merge_frames

# remove duplicates from points
def deduplicate_latest_points(df, value_var, reference_time=None, reference_flag=False):

    info = {
        "initial_rows": 0,
        "reference_active": bool(reference_flag),
        "reference_time": reference_time,
        "filtered_by_reference_time": 0,
        "duplicates_by_value": 0,
        "duplicates_by_gpi": 0,
        "final_rows": 0,
    }

    if df.empty:
        logger.warning(" ===> DataFrame is empty")
        return df.copy(), info

    info["initial_rows"] = len(df)

    logger.info(" ::: DEDUPLICATE POINTS ... ")
    logger.info(f" :::: Initial rows: {info['initial_rows']}")

    df = df.copy()

    # filter by reference time if active
    if reference_flag:
        if reference_time is None:
            logger.warning(" ===> Reference flag is active but reference_time is None")
        else:
            rows_before_filter = len(df)

            df = df.loc[df["time"] <= reference_time].copy()

            info["filtered_by_reference_time"] = rows_before_filter - len(df)

            logger.info(f" :::: Reference time active: {reference_time}")
            logger.info(
                f" :::: Rows filtered by reference time: "
                f"{info['filtered_by_reference_time']}"
            )

    if df.empty:
        logger.warning(" ===> DataFrame is empty after reference time filtering")
        info["final_rows"] = 0
        return df.copy(), info

    # sort by gpi and time
    df = df.sort_values(["gpi", "time"]).copy()

    # duplicates based on ["gpi", value_var]
    dup_mask_value = df.duplicated(subset=["gpi", value_var], keep="last")
    dup_points_value = df.loc[dup_mask_value, ["gpi", "time", value_var]]

    info["duplicates_by_value"] = len(dup_points_value)
    logger.info(f" :::: Duplicates by value: {info['duplicates_by_value']}")

    # remove duplicates by value, keeping latest
    df = df.drop_duplicates(subset=["gpi", value_var], keep="last")

    # duplicates based on ["gpi"]
    dup_mask_gpi = df.duplicated(subset=["gpi"], keep="last")
    dup_points_gpi = df.loc[dup_mask_gpi, ["gpi", "time", value_var]]

    info["duplicates_by_gpi"] = len(dup_points_gpi)
    logger.info(f" :::: Duplicates by gpi/locations: {info['duplicates_by_gpi']}")

    # keep latest point per gpi
    df = (
        df.sort_values(["gpi", "time"])
          .drop_duplicates(subset=["gpi"], keep="last")
    )

    # final sorting
    df = df.sort_values("gpi").reset_index(drop=True)

    info["final_rows"] = len(df)
    logger.info(f" :::: Final rows: {info['final_rows']}")

    logger.info(" ::: DEDUPLICATE POINTS ... DONE")

    return df, info
# ----------------------------------------------------------------------------------------------------------------------