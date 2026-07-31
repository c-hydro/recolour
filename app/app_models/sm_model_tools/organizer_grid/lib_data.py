"""
Library Features:

Name:          lib_data
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260723'
Version:       '1.2.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations

import logging
import os
import string
from pathlib import Path
from typing import Any
from datetime import timedelta

import numpy as np

from lib_utils_geo import read_file_grid
from lib_utils_io import (
    read_file_registry,
    collect_data,
    write_file_geotiff,
)
from lib_utils_analysis import interpolate_points2grid, smooth_grid
from config_info import LOGGER_NAME, VALUE_NODATA_DEFAULT

# logger stream
logger_stream = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class used to preserve unresolved template fields
class UnresolvedField:
    """
    Preserve a field that cannot yet be resolved.

    For example:

        {point}
        {dataset}
        {variable:03d}

    This allows time fields to be formatted first while keeping other fields
    available for a subsequent formatting operation.
    """

    def __init__(self, field_name: str):
        self.field_name = field_name

    def __format__(self, format_spec: str) -> str:
        if format_spec:
            return f"{{{self.field_name}:{format_spec}}}"

        return f"{{{self.field_name}}}"

    def __str__(self) -> str:
        return f"{{{self.field_name}}}"

    def __repr__(self) -> str:
        return f"{{{self.field_name}}}"


# ----------------------------------------------------------------------------------------------------------------------
# formatter used to resolve only available template fields
class PartialFormatter(string.Formatter):
    """
    Python string formatter that preserves unknown fields.

    Example:

        template = (
            "/data/{time:%Y/%m/%d}/"
            "file_{time:%Y%m%d%H}_point_{point}.csv"
        )

        result = formatter.format(
            template,
            time=time_reference,
        )

    Result:

        /data/2026/06/18/file_2026061804_point_{point}.csv
    """

    def get_value(
            self,
            key: Any,
            args: tuple,
            kwargs: dict,
    ) -> Any:

        if isinstance(key, str):
            if key in kwargs:
                return kwargs[key]

            return UnresolvedField(key)

        return super().get_value(key, args, kwargs)


# formatter instance
partial_formatter = PartialFormatter()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to format a path template
def format_path_template(
        path_template: str,
        time_reference,
        preserve_unknown: bool = False,
        **kwargs,
) -> str:
    """
    Format a path using a generic datetime placeholder.

    Supported examples:

        {time:%Y}
        {time:%m}
        {time:%d}
        {time:%H}
        {time:%M}
        {time:%Y%m}
        {time:%Y%m%d}
        {time:%Y%m%d%H}
        {time:%Y%m%d%H%M}
        {time:%Y/%m/%d}

    Additional placeholders can be passed using kwargs:

        {var_name}
        {point}
        {dataset}

    Parameters
    ----------
    path_template : str
        Path containing Python format placeholders.
    time_reference : datetime-like
        Reference datetime used for the ``time`` placeholder.
    preserve_unknown : bool
        Preserve unknown placeholders when True. This is useful for source
        paths where placeholders such as ``{point}`` are resolved later.
    **kwargs
        Additional values used to format the path.

    Returns
    -------
    str
        Formatted path.
    """

    if not isinstance(path_template, str):
        raise TypeError(
            f"Path template must be a string. "
            f"Received: {type(path_template).__name__}"
        )

    format_values = {
        "time": time_reference,
        "time_now": time_reference,
        **kwargs,
    }

    try:
        if preserve_unknown:
            return partial_formatter.format(
                path_template,
                **format_values,
            )

        return path_template.format(**format_values)

    except (KeyError, ValueError, TypeError, AttributeError) as exc:
        raise RuntimeError(
            f"Unable to format path template '{path_template}': {exc}"
        ) from exc


# ----------------------------------------------------------------------------------------------------------------------
# helper to execute data process
def process(settings, time_reference):

    # ------------------------------------------------------------------------------------------------------------------
    # start message - script
    logger_stream.info(" ----> Process execution ... ")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get registry info
    folder_name_reg = settings["data_static"]["registry"]["folder"]
    file_name_reg = settings["data_static"]["registry"]["filename"]
    path_name_reg = os.path.join(folder_name_reg, file_name_reg)

    params_reg = {
        "delimiter": settings["data_static"]["registry"].get(
            "delimiter", ","
        ),
        "tag_col": settings["data_static"]["registry"].get(
            "tag_col", "tag"
        ),
        "lon_col": settings["data_static"]["registry"].get(
            "lon_col", "longitude"
        ),
        "lat_col": settings["data_static"]["registry"].get(
            "lat_col", "latitude"
        ),
        "valid_col": settings["data_static"]["registry"].get(
            "valid_col", "valid"
        ),
        "use_only_valid": bool(
            settings["data_static"]["registry"].get(
                "use_only_valid", True
            )
        ),
    }

    # get grid info
    folder_name_grid = settings["data_static"]["grid"]["folder"]
    file_name_grid = settings["data_static"]["grid"]["filename"]
    path_name_grid = os.path.join(folder_name_grid, file_name_grid)

    params_grid = {
        "band": settings["data_static"]["grid"].get(
            "band", 1
        ),
        "no_data": settings["data_static"]["grid"].get(
            "no_data",
            settings["data_static"]["grid"].get(
                "nodata_values", [-9999]
            )[0],
        ),
    }

    # get source info
    folder_name_src = settings["data_dynamic"]["source"]["folder"]
    file_name_src = settings["data_dynamic"]["source"]["filename"]
    path_name_src = os.path.join(folder_name_src, file_name_src)

    file_format_src = settings["data_dynamic"]["source"].get("format", "csv")

    # Format all source time fields while preserving fields resolved later,
    # such as {point or collections}.
    path_name_src = format_path_template(path_template=path_name_src, time_reference=time_reference, preserve_unknown=True,)

    params_src = {
        "time_format": settings["data_dynamic"]["source"].get(
            "time_format", "%Y-%m-%d %H:%M"
        ),
        "step_hours": float(
            settings["data_dynamic"]["source"].get(
                "step_hours", 1
            )
        ),
        "max_previous_steps": int(
            settings["data_dynamic"]["source"].get(
                "max_previous_steps", 10
            )
        ),
        "nodata_values": settings["data_dynamic"]["source"].get(
            "nodata_values",
            [VALUE_NODATA_DEFAULT, -9998],
        ),
        "valid_min": settings["data_dynamic"]["source"].get(
            "valid_min", 0
        ),
        "valid_max": settings["data_dynamic"]["source"].get(
            "valid_max", 1
        ),
        "delimiter": settings["data_dynamic"]["source"].get(
            "delimiter", ","
        ),
        "time_col": settings["data_dynamic"]["source"].get(
            "time_col", "time"
        ),
        "value_col": settings["data_dynamic"]["source"].get(
            "value_col", "theta_simulated"
        ),
        "selection_method": settings["data_dynamic"]["source"].get(
            "selection_method", "last_mod"
        ),
        "selection_columns": settings["data_dynamic"]["source"].get(
            "selection_columns",
            {
                "theta_simulated": "theta_simulated",
                "theta_observed": "theta_observed",
            },
        ),
        "strict_nodata": settings["data_dynamic"]["source"].get(
            "strict_nodata", VALUE_NODATA_DEFAULT
        ),
        "registry_lon_col": settings["data_dynamic"]["source"].get(
            "registry_lon_col", "longitude"
        ),
        "registry_lat_col": settings["data_dynamic"]["source"].get(
            "registry_lat_col", "latitude"
        ),
    }

    # Validate type, extension and filename template
    file_ext = Path(path_name_src).suffix.lower()
    if file_format_src == "csv":
        if file_ext != ".csv":
            raise ValueError(f"CSV output type requires a '.csv' destination file, found '{file_ext}'.")
    elif file_format_src == "netcdf":
        if file_ext != ".nc":
            raise ValueError(f"NetCDF output type requires a '.nc' destination file, found '{file_ext}'.")
    else:
        raise ValueError(f"Unrecognized file format '{file_format_src}'.")

    # get destination info for logging
    folder_name_dst = settings["data_dynamic"]["destination"]["folder"]
    file_name_dst = settings["data_dynamic"]["destination"]["filename"]
    path_name_dst = os.path.join(folder_name_dst, file_name_dst)

    name_dst = settings["data_dynamic"]["destination"].get("name", "soil_moisture")
    path_name_dst = format_path_template(path_template=path_name_dst,time_reference=time_reference, var_name=name_dst,)

    # get interpolation info
    params_interp = {
        "no_data": settings["interpolation"].get(
            "no_data", VALUE_NODATA_DEFAULT
        ),
        "method": settings["interpolation"].get(
            "method", "idw"
        ),
        "roi_m": settings["interpolation"].get(
            "roi_m", 10000
        ),
        "min_points": settings["interpolation"].get(
            "min_points", 1
        ),
        "max_points": settings["interpolation"].get(
            "max_points", 8
        ),
        "idw_power": settings["interpolation"].get(
            "idw_power", 8
        ),
    }

    # get smooth info
    params_smooth = {
        "no_data": settings["smoothing"].get(
            "no_data", VALUE_NODATA_DEFAULT
        ),
        "active": settings["smoothing"].get(
            "active", True
        ),
        "method": settings["smoothing"].get(
            "method", "gaussian"
        ),
        "iterations": settings["smoothing"].get(
            "iterations", 1
        ),
        "preserve_original_nodata": settings["smoothing"].get(
            "preserve_original_nodata", True
        ),
        "preserve_range": settings["smoothing"].get(
            "preserve_range", True
        ),
        "original_weight": settings["smoothing"].get(
            "original_weight", 0.35
        ),

        # Gaussian
        "sigma": settings["smoothing"].get(
            "sigma", 1.5
        ),
        "sigma_x": settings["smoothing"].get(
            "sigma_x", None
        ),
        "sigma_y": settings["smoothing"].get(
            "sigma_y", None
        ),
        "theta_deg": settings["smoothing"].get(
            "theta_deg", 0.0
        ),

        # Circular kernels
        "radius": settings["smoothing"].get(
            "radius", 2.0
        ),
        "slope": settings["smoothing"].get(
            "slope", 1.0
        ),

        # Box kernel
        "width": settings["smoothing"].get(
            "width", 3
        ),

        # Kernel generation
        "oversampling_factor": settings["smoothing"].get(
            "oversampling_factor", 10
        ),
    }
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # log paths
    logger_stream.info(f" ----> Source template:   {path_name_src}")
    logger_stream.info(f" ----> Destination path: {path_name_dst}")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # read file registry
    logger_stream.info(f" ----> Read registry {path_name_reg} ... ")
    registry_obj = read_file_registry(path_name_reg,params_reg,)
    logger_stream.info(f" ----> Read registry {path_name_reg} ... DONE")

    # read file grid
    logger_stream.info(f" ----> Read grid {path_name_grid} ... ")
    grid_mask, grid_lon, grid_lat, grid_profile = read_file_grid(path_name_grid,params_grid,)
    logger_stream.info( f" ----> Read grid {path_name_grid} ... DONE")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # collect data points
    logger_stream.info(" ----> Collect data points ... ")

    points_data, points_times, points_lags = collect_data(
        path_name_src,
        registry_obj, time_reference, params_src, format_data=file_format_src
    )

    # check valid points for interpolation
    if len(points_data) == 0:
        logger_stream.warning(
            " ====> No valid point values available for interpolation"
        )
        logger_stream.info(
            " ----> Collect data points ... FAILED"
        )

        return None

    logger_stream.info(
        f" ----> Valid points: {len(points_data)}"
    )
    logger_stream.info(
        " ----> Collect data points ... DONE"
    )

    # interpolate data points
    logger_stream.info(
        " ----> Interpolate data points to grid ... "
    )

    grid_data = interpolate_points2grid(
        points_data,
        grid_lon,
        grid_lat,
        grid_mask,
        params_interp,
    )

    logger_stream.info(
        " ----> Interpolate data points to grid ... DONE"
    )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # smooth grid
    logger_stream.info(" ----> Smooth grid ... ")

    grid_data = smooth_grid(
        data=grid_data,
        valid_mask=grid_mask,
        cfg=params_smooth,
    )

    logger_stream.info(" ----> Smooth grid ... DONE")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message - script
    logger_stream.info(" ----> Process execution ... DONE")

    return grid_data, grid_mask, grid_profile
    # ------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to save data
def save(
        grid_data,
        grid_mask,
        grid_profile,
        settings,
        time_reference,
):

    # ------------------------------------------------------------------------------------------------------------------
    # start message - script
    logger_stream.info(" ----> Dumping execution ... ")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get destination info
    folder_name_dst = settings["data_dynamic"]["destination"]["folder"]
    file_name_dst = settings["data_dynamic"]["destination"]["filename"]
    path_name_dst = os.path.join(folder_name_dst, file_name_dst)

    name_dst = settings["data_dynamic"]["destination"].get(
        "name", "soil_moisture"
    )
    compress_dst = settings["data_dynamic"]["destination"].get(
        "compress", "deflate"
    )
    no_data_dst = settings["data_dynamic"]["destination"].get(
        "no_data", np.nan
    )
    
    # optional time shift
    shift_cfg = settings["data_dynamic"]["destination"].get(
        "time_shift", {}
    )
    # ------------------------------------------------------------------------------------------------------------------
	
    # ------------------------------------------------------------------------------------------------------------------
    # compute time shifted (if needed by the algorithm)
    time_destination = time_reference + timedelta(
        days=shift_cfg.get("days", 0),
        hours=shift_cfg.get("hours", 0),
        minutes=shift_cfg.get("minutes", 0),
        seconds=shift_cfg.get("seconds", 0),
    )
    # ------------------------------------------------------------------------------------------------------------------
	
    # ------------------------------------------------------------------------------------------------------------------
    # check no data
    if no_data_dst is None:
        no_data_dst = np.nan

    # define destination path
    path_name_dst = format_path_template(
        path_template=path_name_dst,
        time_reference=time_destination,
        var_name=name_dst,
    )

    # create destination folder
    folder_name_dst = os.path.dirname(path_name_dst)

    if folder_name_dst:
        os.makedirs(
            folder_name_dst,
            exist_ok=True,
        )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # dump data
    logger_stream.info(
        f" -----> Reference time: {time_reference}"
    )
    logger_stream.info(
        f" -----> File path:      {path_name_dst}"
    )

    write_file_geotiff(
        path_name_dst,
        grid_data,
        grid_profile,
        no_data_dst,
        compress=compress_dst,
    )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message - script
    logger_stream.info(" ----> Dumping execution ... DONE")
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
