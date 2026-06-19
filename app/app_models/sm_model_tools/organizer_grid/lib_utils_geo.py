"""
Library Features:

Name:           lib_utils_geo
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import rasterio
import math
import numpy as np

from rasterio.transform import xy as rio_xy
from rasterio.warp import transform as rio_transform

from typing import Optional, Sequence, Any, Tuple

from config_info import LOGGER_NAME, EARTH_RADIUS_M

# logger
logger_stream = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to radius
def radius_to_chord(radius_m: float) -> float:
    return 2.0 * math.sin(float(radius_m) / (2.0 * EARTH_RADIUS_M))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to transform lonlat 2 xyz
def lonlat_to_unit_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    return np.column_stack([x.ravel(), y.ravel(), z.ravel()])
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to compute radius
def great_circle_distance_m(lon0: np.ndarray, lat0: np.ndarray, lon1: np.ndarray, lat1: np.ndarray) -> np.ndarray:
    lon0r = np.deg2rad(lon0)
    lat0r = np.deg2rad(lat0)
    lon1r = np.deg2rad(lon1)
    lat1r = np.deg2rad(lat1)
    dlon = lon1r - lon0r
    dlat = lat1r - lat0r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat0r) * np.cos(lat1r) * np.sin(dlon / 2.0) ** 2
    return EARTH_RADIUS_M * 2.0 * np.arcsin(np.sqrt(a))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to compute radius
def get_grid_lonlat(src: Any) -> Tuple[np.ndarray, np.ndarray]:
    rows = np.arange(src.height)
    cols = np.arange(src.width)
    cc, rr = np.meshgrid(cols, rows)
    xs, ys = rio_xy(src.transform, rr, cc, offset="center")
    xs = np.asarray(xs, dtype=np.float64)
    ys = np.asarray(ys, dtype=np.float64)
    if src.crs and str(src.crs).upper() not in ("EPSG:4326", "OGC:CRS84"):
        lons, lats = rio_transform(src.crs, "EPSG:4326", xs.ravel().tolist(), ys.ravel().tolist())
        return np.asarray(lons).reshape(xs.shape), np.asarray(lats).reshape(ys.shape)
    return xs, ys
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read file grid
def read_file_grid(path: str, params: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:

    band = int(params.get("band", 1))
    no_data = params.get("no_data", None)

    if not os.path.exists(path):
        logger_stream.error(f"File grid {path} does not exist")
        raise RuntimeError("File grid is mandatory")

    with rasterio.open(path) as src:

        dem = src.read(band).astype(np.float64)
        mask = np.isfinite(dem)

        if src.nodata is not None:
            mask &= ~np.isclose(dem, float(src.nodata), equal_nan=False)

        if no_data is not None:
            if not isinstance(no_data, (list, tuple, np.ndarray)):
                no_data = [no_data]

            for value in no_data:
                mask &= ~np.isclose(dem, float(value), equal_nan=False)

        lon, lat = get_grid_lonlat(src)

        profile = src.profile

    return mask, lon, lat, profile
# ----------------------------------------------------------------------------------------------------------------------
