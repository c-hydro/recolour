"""
Library Features:

Name:          lib_utils_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260618'
Version:       '1.5.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations

import numpy as np

from pyresample import geometry, kd_tree
from typing import Any, Dict, Sequence

from lib_utils_io import PointValue
from config_info import LOGGER_NAME, VALUE_NODATA_DEFAULT
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to interpolate points to grid
def interpolate_points2grid(
    points: Sequence[PointValue],
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
    dem_valid_mask: np.ndarray,
    cfg: Dict[str, Any],
) -> np.ndarray:

    nodata = float(cfg.get("nodata", VALUE_NODATA_DEFAULT))
    roi_m = float(cfg.get("roi_m", cfg.get("radius_m", 10000)))
    method = str(cfg.get("method", "idw")).lower()
    power = float(cfg.get("idw_power", 2.0))
    min_points = int(cfg.get("min_points", 1))
    max_points = cfg.get("max_points", None)
    max_points = int(max_points) if max_points is not None else None

    if lon_grid.shape != dem_valid_mask.shape:
        lon_grid = lon_grid.reshape(dem_valid_mask.shape)

    if lat_grid.shape != dem_valid_mask.shape:
        lat_grid = lat_grid.reshape(dem_valid_mask.shape)

    out = np.full(dem_valid_mask.shape, nodata, dtype=np.float32)

    if not points:
        return out

    point_lon = np.asarray([p.lon for p in points], dtype=np.float64)
    point_lat = np.asarray([p.lat for p in points], dtype=np.float64)
    point_val = np.asarray([p.value for p in points], dtype=np.float32)

    valid_points = (
        np.isfinite(point_lon)
        & np.isfinite(point_lat)
        & np.isfinite(point_val)
    )

    if not np.any(valid_points):
        return out

    point_lon = point_lon[valid_points]
    point_lat = point_lat[valid_points]
    point_val = point_val[valid_points]

    source_def = geometry.SwathDefinition(
        lons=point_lon,
        lats=point_lat,
    )

    target_def = geometry.SwathDefinition(
        lons=lon_grid,
        lats=lat_grid,
    )

    neighbours = max_points if max_points is not None else len(point_val)
    neighbours = max(neighbours, min_points)

    if method == "nearest":

        data = kd_tree.resample_nearest(
            source_def,
            point_val,
            target_def,
            radius_of_influence=roi_m,
            fill_value=nodata,
        )

    elif method == "mean":

        def weight_mean(dist):
            return np.ones_like(dist, dtype=np.float64)

        data = kd_tree.resample_custom(
            source_def,
            point_val,
            target_def,
            radius_of_influence=roi_m,
            neighbours=neighbours,
            weight_funcs=weight_mean,
            fill_value=nodata,
        )

    elif method == "idw":

        def weight_idw(dist):
            return 1.0 / np.maximum(dist, 1.0) ** power

        data = kd_tree.resample_custom(
            source_def,
            point_val,
            target_def,
            radius_of_influence=roi_m,
            neighbours=neighbours,
            weight_funcs=weight_idw,
            fill_value=nodata,
        )

    else:
        raise ValueError(f"Unsupported interpolation method: {method}")

    data = np.asarray(data, dtype=np.float32)

    valid_target = (
        dem_valid_mask
        & np.isfinite(lon_grid)
        & np.isfinite(lat_grid)
        & np.isfinite(data)
        & ~np.isclose(data, nodata)
    )

    out[valid_target] = data[valid_target]

    return out
# ----------------------------------------------------------------------------------------------------------------------
