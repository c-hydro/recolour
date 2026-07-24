
"""
Library Features:

Name:          lib_io_ascii
Author(s):     Fabio Delogu
Date:          '20260715'
Version:       '1.1.0'

Purpose:
    Validate, load and organize geographical datasets configured in JSON.

Supported types:
    - ascii_grid
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import rasterio

from typing import Any, Dict

from lib_utils_geo import create_coords
from lib_io_base import process_values
from config_info import LOGGER_NAME, TIME_FMT_CLI

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read a ascii grid file
def read_file_ascii_grid(
        file_path, file_band: int = 1, file_crs: str = 'epsg:4326',
        file_scale_factor: float = 1, file_offset: float = 0,
        file_valid_min: float = None, file_valid_max: float = None,
        file_nodata: float = None, no_data_values: list = [],) -> Dict[str, Any]:

    # check format of no_data_values
    if no_data_values is not None:
        if not isinstance(no_data_values, list):
            no_data_values = [no_data_values]

    # try to read file
    try:
        with rasterio.open(file_path) as file_handle:

            if file_band > file_handle.count:
                raise ValueError(
                    f"Band {file_band} was requested but the file contains only {file_handle.count} band(s).")

            values_raw = file_handle.read(file_band)

            transform = file_handle.transform
            width = file_handle.width
            height = file_handle.height

            crs = file_handle.crs
            nodata = file_handle.nodata

            profile = file_handle.profile.copy()
            tags = file_handle.tags(file_band)

    except Exception as exc:
        raise RuntimeError(f"Unable to read file'{file_path}'.") from exc

    # check crs
    if crs is None:
        crs = file_crs
        logger.warning(f" ===> Geo file {file_path} does not define a CRS. Configured CRS {file_crs} will be used.",)
    else:
        if crs != file_crs:
            logger.warning(
                f" ===> CRS mismatch; file CRS is {crs}, configured CRS is {file_crs}. File CRS will be used.")

    # process data
    values = process_values(
        values=values_raw, dtype=np.float64,
        scale_factor=file_scale_factor, offset=file_offset,
        valid_min=file_valid_min, valid_max=file_valid_max,
        nodata_values=no_data_values, file_nodata=file_nodata,
    )

    # create coordinates
    longitude, latitude = create_coords(transform=transform, width=width, height=height)

    # create finite mask and apply to values
    finite_mask = np.isfinite(values)
    finite_values = values[finite_mask]

    # recompute min, max and mean
    if finite_values.size > 0:
        value_min = float(np.nanmin(finite_values))
        value_max = float(np.nanmax(finite_values))
        value_mean = float(np.nanmean(finite_values))
    else:
        value_min, valid_max, valid_mean = np.nan, np.nan, np.nan
        logger.warning(f"Dataset contains no valid values.")

    # create datasets object
    geo_dataset = {

        "values": values,
        "longitude": longitude,
        "latitude": latitude,

        "transform": transform,
        "crs": crs,
        "width": width,
        "height": height,

        "band": file_band,
        "nodata": nodata,
        "file_path": file_path,

        "profile": profile,
        "tags": tags,

        "valid_count": int(np.count_nonzero(finite_mask)),
        "invalid_count": int(np.count_nonzero(~finite_mask)),

        "value_min": value_min,
        "value_max": value_max,
        "value_mean": value_mean
    }

    return geo_dataset
# ----------------------------------------------------------------------------------------------------------------------
