# libraries
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import logging
from lib_utils_time import datetime_to_jd

# method to decode NetCDF files and extract soil moisture data
def decode_netcdf(filepath):

    try:
        ds = xr.open_dataset(filepath)

        # Extract core variables
        lat = ds["latitude"].values.astype(np.float32)
        lon = ds["longitude"].values.astype(np.float32)
        lat = np.clip(lat, -90.0, 90.0)
        lon = np.clip(lon, -180.0, 180.0)

        v_ssm = ds["surface_soil_moisture"].values
        v_proc = ds["processing_flag"].values if "processing_flag" in ds else np.full_like(v_ssm, np.nan)
        v_corr = ds["correction_flag"].values if "correction_flag" in ds else np.full_like(v_ssm, np.nan)
        v_qual = ds["surface_flag"].values if "surface_flag" in ds else np.full_like(v_ssm, np.nan)

        # Time handling
        v_time_index = pd.to_datetime(ds["time"].values)
        v_dates = np.array([[t.year, t.month, t.day, t.hour, t.minute, t.second] for t in v_time_index])

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Discarding nonzero nanoseconds")
            v_jd = np.array([datetime_to_jd(t.to_pydatetime()) for t in v_time_index])

        #v_jd = np.array([datetime_to_jd(t.to_pydatetime()) for t in v_time_index])

        # Build output in same structure as decode_bufr
        workspace = {
            "longitude": lon,
            "latitude": lat,
            "ssm": v_ssm,
            "flag_processing": v_proc,
            "flag_corrections": v_corr,
            "quality": v_qual,
            "jd": v_jd,
            "time": v_time_index,
            "date": v_dates
        }

        return workspace

    except Exception as e:
        logging.error(f"Error decoding NetCDF: {e}")
        return {
            "longitude": np.array([]),
            "latitude": np.array([]),
            "ssm": np.array([]),
            "flag_processing": np.array([]),
            "flag_corrections": np.array([]),
            "quality": np.array([]),
            "jd": np.array([]),
            "time": pd.to_datetime([]),
            "date": np.array([])
        }
