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

import logging
import numpy as np

from pyresample import geometry, kd_tree
from typing import Any, Dict, Sequence

from astropy.convolution import (
    convolve,
    Box2DKernel,
    Gaussian2DKernel,
    Tophat2DKernel,
    TrapezoidDisk2DKernel,
)

from lib_utils_io import PointValue
from config_info import LOGGER_NAME, VALUE_NODATA_DEFAULT

# logger stream
logger_stream = logging.getLogger(LOGGER_NAME)
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

# ----------------------------------------------------------------------------------------------------------------------
# method to smooth grid using Astropy
def smooth_grid(
        data: np.ndarray,
        valid_mask: np.ndarray,
        cfg: Dict[str, Any],
) -> np.ndarray:
    """
    Smooth interpolated gridded data using Astropy convolution kernels.

    Invalid cells are handled using normalized convolution:

        smoothed = convolution(data * weights) / convolution(weights)

    This avoids Astropy warnings caused by contiguous NaN regions larger
    than the smoothing kernel.

    Supported methods:
        - gaussian
        - box
        - tophat
        - trapezoid
    """

    active = bool(cfg.get("active", False))
    method = str(cfg.get("method", "gaussian")).lower()
    nodata = float(cfg.get("nodata", VALUE_NODATA_DEFAULT))
    iterations = int(cfg.get("iterations", 1))

    preserve_original_nodata = bool(
        cfg.get("preserve_original_nodata", True)
    )

    preserve_range = bool(
        cfg.get("preserve_range", True)
    )

    original_weight = float(
        cfg.get("original_weight", 0.35)
    )

    # Minimum convolved support required to consider an output cell valid.
    # A small positive value is usually sufficient.
    minimum_weight = float(
        cfg.get("minimum_weight", 1.0e-6)
    )

    # -------------------------------------------------------------------------
    # Check parameters

    data = np.asarray(
        data,
        dtype=np.float64
    )

    valid_mask = np.asarray(
        valid_mask,
        dtype=bool
    )

    if data.ndim != 2:
        raise ValueError(
            f"Input data must be two-dimensional. "
            f"Received: {data.ndim}D"
        )

    if data.shape != valid_mask.shape:
        raise ValueError(
            f"Data shape {data.shape} differs from valid-mask shape "
            f"{valid_mask.shape}"
        )

    if not 0.0 <= original_weight <= 1.0:
        raise ValueError(
            f'"original_weight" must be between 0 and 1. '
            f"Received: {original_weight}"
        )

    if not active or iterations < 1:
        logger_stream.warning(
            " ===> Smoothing is not activated or iterations "
            "are less than 1"
        )

        return data.astype(
            np.float32
        )

    # -------------------------------------------------------------------------
    # Define original valid cells

    original_valid = (
        valid_mask
        & np.isfinite(data)
        & ~np.isclose(data, nodata)
    )

    output = np.full(
        data.shape,
        nodata,
        dtype=np.float32
    )

    if not np.any(original_valid):
        logger_stream.warning(
            " ===> All datasets are invalid. Skip smoothing"
        )

        return output

    # -------------------------------------------------------------------------
    # Create kernel

    kernel = get_smoothing_kernel(
        cfg
    )

    data_min = float(
        np.nanmin(data[original_valid])
    )

    data_max = float(
        np.nanmax(data[original_valid])
    )

    logger_stream.info(
        " ----> Smooth grid using method=%s iterations=%d "
        "original_weight=%.3f",
        method,
        iterations,
        original_weight
    )

    logger_stream.info(
        " -----> Input valid cells: %d / %d; min=%.6f max=%.6f",
        np.count_nonzero(original_valid),
        original_valid.size,
        data_min,
        data_max
    )

    # Keep original values as NaN outside the valid input footprint.
    original_data = np.where(
        original_valid,
        data,
        np.nan
    )

    data_work = original_data.copy()

    # -------------------------------------------------------------------------
    # Iterative normalized convolution

    for iteration in range(iterations):

        current_valid = (
            valid_mask
            & np.isfinite(data_work)
        )

        if preserve_original_nodata:
            current_valid &= original_valid

        if not np.any(current_valid):
            logger_stream.warning(
                " ===> No valid cells available at smoothing "
                "iteration %d",
                iteration + 1
            )
            break

        # Replace invalid cells by zero for the numerator.
        data_numerator = np.where(
            current_valid,
            data_work,
            0.0
        )

        # Binary support weights.
        data_weights = current_valid.astype(
            np.float64
        )

        # Convolve values.
        convolved_values = convolve(
            data_numerator,
            kernel,
            boundary="extend",
            nan_treatment="fill",
            fill_value=0.0,
            normalize_kernel=False,
            preserve_nan=False
        )

        # Convolve support weights using the same kernel.
        convolved_weights = convolve(
            data_weights,
            kernel,
            boundary="extend",
            nan_treatment="fill",
            fill_value=0.0,
            normalize_kernel=False,
            preserve_nan=False
        )

        # Compute normalized convolution only where the kernel has support.
        data_smoothed = np.full(
            data.shape,
            np.nan,
            dtype=np.float64
        )

        supported_cells = (
            valid_mask
            & np.isfinite(convolved_weights)
            & (convolved_weights > minimum_weight)
        )

        data_smoothed[supported_cells] = (
            convolved_values[supported_cells]
            / convolved_weights[supported_cells]
        )

        # Blend the original values and the smoothed values.
        blend_valid = (
            supported_cells
            & np.isfinite(data_smoothed)
        )

        data_next = np.full(
            data.shape,
            np.nan,
            dtype=np.float64
        )

        # Where original data exist, blend original and smoothed values.
        original_blend = (
            blend_valid
            & original_valid
        )

        data_next[original_blend] = (
            original_weight
            * original_data[original_blend]
            + (1.0 - original_weight)
            * data_smoothed[original_blend]
        )

        # When extension into original nodata areas is allowed, use only
        # the smoothed value because no original value exists there.
        if not preserve_original_nodata:

            extended_cells = (
                blend_valid
                & ~original_valid
            )

            data_next[extended_cells] = data_smoothed[
                extended_cells
            ]

        # Never extend values outside the DEM domain.
        data_next[~valid_mask] = np.nan

        if preserve_original_nodata:
            data_next[~original_valid] = np.nan

        data_work = data_next

        finite_iteration = np.isfinite(
            data_work
        )

        if np.any(finite_iteration):
            logger_stream.info(
                " -----> Iteration %d/%d: valid=%d min=%.6f max=%.6f",
                iteration + 1,
                iterations,
                np.count_nonzero(finite_iteration),
                float(np.nanmin(data_work)),
                float(np.nanmax(data_work))
            )
        else:
            logger_stream.warning(
                " ===> Iteration %d/%d produced no finite values",
                iteration + 1,
                iterations
            )

    # -------------------------------------------------------------------------
    # Preserve original range

    if preserve_range:
        data_work = np.clip(
            data_work,
            data_min,
            data_max
        )

    # -------------------------------------------------------------------------
    # Define output

    output_valid = (
        valid_mask
        & np.isfinite(data_work)
    )

    if preserve_original_nodata:
        output_valid &= original_valid

    output[output_valid] = data_work[
        output_valid
    ].astype(np.float32)

    logger_stream.info(
        " -----> Smoothed output valid cells: %d / %d",
        np.count_nonzero(output_valid),
        output_valid.size
    )

    if np.any(output_valid):
        logger_stream.info(
            " -----> Smoothed output range: %.6f / %.6f",
            float(np.nanmin(output[output_valid])),
            float(np.nanmax(output[output_valid]))
        )

    return output
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to smooth grid using Astropy
def smooth_grid_OLD(
        data: np.ndarray,
        valid_mask: np.ndarray,
        cfg: Dict[str, Any],
) -> np.ndarray:
    """
    Smooth interpolated gridded data using Astropy convolution kernels.

    Supported methods:
        - gaussian
        - box
        - tophat
        - trapezoid

    Parameters
    ----------
    data : np.ndarray
        Interpolated two-dimensional grid.

    valid_mask : np.ndarray
        Boolean mask defining valid DEM cells.

    cfg : Dict[str, Any]
        Configuration dictionary. Smoothing settings are read from
        cfg["smoothing"].

    Returns
    -------
    np.ndarray
        Smoothed grid with the same shape as the input.
    """

    active = bool(cfg.get("active", False))
    method = str(cfg.get("method", "gaussian")).lower()
    nodata = float(cfg.get("nodata", VALUE_NODATA_DEFAULT))
    iterations = int(cfg.get("iterations", 1))
    preserve_original_nodata = bool(cfg.get("preserve_original_nodata", True))
    preserve_range = bool(cfg.get("preserve_range", True))

    # Weight assigned to the original interpolated values:
    # 0.0 = completely smoothed
    # 0.35 -> balanced
    # 0.50 -> centers and original values preserved more
    # 0.70 -> limited smoothing
    # 1.0 = completely original
    original_weight = float(cfg.get("original_weight", 0.35))

    data = np.asarray(data, dtype=np.float64)
    valid_mask = np.asarray(valid_mask, dtype=bool)
    if data.ndim != 2:
        raise ValueError(
            f"Input data must be two-dimensional. Received: {data.ndim}D"
        )

    if data.shape != valid_mask.shape:
        raise ValueError(
            f"Data shape {data.shape} differs from valid-mask shape "
            f"{valid_mask.shape}"
        )

    if not active or iterations < 1:
        logger_stream.warning(f' ===> Smoothing is not activated or iterations are less than 1')
        return data.astype(np.float32)

    original_valid = (
            valid_mask
            & np.isfinite(data)
            & ~np.isclose(data, nodata)
    )

    output = np.full(
        data.shape,
        nodata,
        dtype=np.float32,
    )

    if not np.any(original_valid):
        logger_stream.warning(f' ===> All datasets are not valid. Skip smoothing')
        return output

    kernel = get_smoothing_kernel(cfg)

    data_min = float(np.nanmin(data[original_valid]))
    data_max = float(np.nanmax(data[original_valid]))

    # Astropy uses NaN values to identify cells excluded from convolution.
    original_data = np.where(
        original_valid,
        data,
        np.nan,
    )

    data_work = original_data.copy()

    for _ in range(iterations):

        data_smoothed = convolve(
            data_work,
            kernel,
            boundary="extend",
            nan_treatment="interpolate",
            normalize_kernel=True,
            preserve_nan=False,
        )

        data_work = (
                original_weight * original_data
                + (1.0 - original_weight) * data_smoothed
        )

        # Prevent convolution from extending values outside the DEM domain.
        data_work[~valid_mask] = np.nan

        if preserve_original_nodata:
            data_work[~original_valid] = np.nan

    if preserve_range:
        data_work = np.clip(
            data_work,
            data_min,
            data_max,
        )

    output_valid = (
            valid_mask
            & np.isfinite(data_work)
    )

    if preserve_original_nodata:
        output_valid &= original_valid

    output[output_valid] = data_work[output_valid].astype(np.float32)

    return output

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get a configuration value, using a fallback for missing or None values
def get_cfg_value(cfg, key, default):
    value = cfg.get(key, default)

    if value is None:
        value = default

    return value
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create Gaussian smoothing kernel
def create_gaussian_kernel(cfg):

    sigma = float(get_cfg_value(cfg, "sigma", 1.5))

    sigma_x = float(get_cfg_value(cfg, "sigma_x", sigma))
    sigma_y = float(get_cfg_value(cfg, "sigma_y", sigma))
    theta_deg = float(get_cfg_value(cfg, "theta_deg", 0.0))

    if sigma_x <= 0 or sigma_y <= 0:
        raise ValueError(
            f"Gaussian sigma must be greater than zero: "
            f"sigma_x={sigma_x}, sigma_y={sigma_y}"
        )

    return Gaussian2DKernel(
        x_stddev=sigma_x,
        y_stddev=sigma_y,
        theta=np.deg2rad(theta_deg),
    )
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create box smoothing kernel
def create_box_kernel(cfg):

    width = int(get_cfg_value(cfg, "width", 3))

    if width < 1:
        raise ValueError(f"Box width must be at least 1: width={width}")

    # Use an odd width to have a central pixel.
    if width % 2 == 0:
        width += 1

    return Box2DKernel(width=width)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create top-hat smoothing kernel
def create_tophat_kernel(cfg):

    radius = float(get_cfg_value(cfg, "radius", 2.0))

    if radius <= 0:
        raise ValueError(
            f"Top-hat radius must be greater than zero: radius={radius}"
        )

    return Tophat2DKernel(radius=radius)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create trapezoid smoothing kernel
def create_trapezoid_kernel(cfg):

    radius = float(get_cfg_value(cfg, "radius", 2.0))
    slope = float(get_cfg_value(cfg, "slope", 1.0))

    if radius <= 0:
        raise ValueError(
            f"Trapezoid radius must be greater than zero: radius={radius}"
        )

    if slope <= 0:
        raise ValueError(
            f"Trapezoid slope must be greater than zero: slope={slope}"
        )

    return TrapezoidDisk2DKernel(
        radius=radius,
        slope=slope,
    )
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# smoothing-kernel registry
SMOOTHING_KERNELS = {
    "gaussian": create_gaussian_kernel,
    "box": create_box_kernel,
    "tophat": create_tophat_kernel,
    "trapezoid": create_trapezoid_kernel,
}
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get smoothing kernel
def get_smoothing_kernel(cfg):

    method = str(get_cfg_value(cfg, "method", "gaussian")).lower().strip()

    if method not in SMOOTHING_KERNELS:
        raise ValueError(
            f"Unknown smoothing method '{method}'. "
            f"Available methods: {list(SMOOTHING_KERNELS.keys())}"
        )

    try:
        return SMOOTHING_KERNELS[method](cfg)

    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid configuration for smoothing method '{method}': "
            f"{cfg}. Error: {exc}"
        ) from exc
# ----------------------------------------------------------------------------------------------------------------------
