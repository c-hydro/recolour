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

from config_utils import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to check 1d arrays
def check_1d_locations_var(dataset, var_name):
    if var_name not in dataset.variables and var_name not in dataset.coords:
        raise KeyError(f"Variable '{var_name}' not found")

    data_array = dataset[var_name]
    if data_array.dims != ("locations",):
        raise ValueError(f"Variable '{var_name}' must have dims ('locations',), got {data_array.dims}")
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# point collection and dedup
def collect_points_to_dataframe(settings, file_list):
    source_settings = settings.get("source", {})
    value_var = source_settings.get("value_col", "surface_soil_moisture")

    frames = []

    for file_path in file_list:
        logger.info(f" ----> Read source file: {file_path}")

        with xr.open_dataset(file_path) as dataset:
            check_1d_locations_var(dataset, "lon")
            check_1d_locations_var(dataset, "lat")
            check_1d_locations_var(dataset, "gpi")
            check_1d_locations_var(dataset, "time")
            check_1d_locations_var(dataset, value_var)

            lon = np.asarray(dataset["lon"].values, dtype=np.float64)
            lat = np.asarray(dataset["lat"].values, dtype=np.float64)
            gpi = np.asarray(dataset["gpi"].values, dtype=np.int64)
            obs_time = pd.to_datetime(dataset["time"].values)
            obs_value = np.asarray(dataset[value_var].values, dtype=np.float64)

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

    if not frames:
        return pd.DataFrame(columns=["gpi", "lon", "lat", "time", value_var, "source_file"])

    return pd.concat(frames, ignore_index=True)


def deduplicate_latest_points(df, value_var):
    if df.empty:
        return df.copy()

    df = df.sort_values(["gpi", "time"]).copy()
    df = df.drop_duplicates(subset=["gpi", value_var], keep="last")
    df = df.sort_values(["gpi", "time"]).drop_duplicates(subset=["gpi"], keep="last")
    df = df.sort_values("gpi").reset_index(drop=True)

    return df
# ----------------------------------------------------------------------------------------------------------------------