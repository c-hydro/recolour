"""
Library Features:

Name:          lib_io_tiff
Author(s):     Fabio Delogu
Date:          '20260716'
Version:       '1.0.0'

Purpose:
    Validate, load and organize soil-moisture GeoTIFF datasets.

Supported types:
    - geotiff
    - tiff

Notes:
    The TIFF internal compression, such as DEFLATE, LZW or ZSTD, is handled
    transparently by rasterio. No external decompression step is required.
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from datetime import datetime
from typing import Any, Dict, Optional, Sequence, Union

import numpy as np
import rasterio

from rasterio.crs import CRS

from lib_utils_geo import create_coords
from lib_io_base import process_values
from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# constants
TIFF_DEFAULT_CRS = "EPSG:4326"
TIFF_DEFAULT_NODATA = -9999.0
TIFF_DEFAULT_VARIABLE = "soil_moisture"
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to normalize scalar or sequence values to a list
def _normalize_list(
        values: Optional[Union[Any, Sequence[Any]]],
) -> list:
    """
    Normalize scalar or sequence values to a Python list.
    """

    if values is None:
        return []

    if isinstance(values, list):
        return values

    if isinstance(values, (tuple, set, np.ndarray)):
        return list(values)

    return [values]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to parse time from TIFF filename
def parse_time_from_filename(
        file_path: str,
        time_pattern: str = "%Y%m%d%H",
) -> Optional[datetime]:
    """
    Parse the reference time from filenames such as:

        soil_moisture_2026071123.tif

    Parameters
    ----------
    file_path : str
        TIFF file path.

    time_pattern : str
        Datetime pattern used by the filename timestamp.

    Returns
    -------
    datetime or None
        Parsed reference time.
    """

    file_name = os.path.basename(
        file_path
    )

    file_stem = os.path.splitext(
        file_name
    )[0]

    # The timestamp is expected after the final underscore.
    time_string = file_stem.rsplit(
        "_",
        maxsplit=1,
    )[-1]

    # Remove optional filename suffixes such as "(1)".
    if "(" in time_string:
        time_string = time_string.split(
            "(",
            maxsplit=1,
        )[0]

    try:
        return datetime.strptime(
            time_string,
            time_pattern,
        )

    except ValueError:

        logger.warning(
            " ===> Unable to parse reference time from TIFF filename '%s' "
            "using format '%s'.",
            file_name,
            time_pattern,
        )

        return None
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to collect TIFF nodata values
def get_tiff_nodata_values(
        file_nodata_metadata: Optional[float],
        file_nodata: Optional[float] = None,
        no_data_values: Optional[Union[Any, Sequence[Any]]] = None,
) -> list:
    """
    Collect nodata values from TIFF metadata and configuration.
    """

    nodata_values = _normalize_list(
        no_data_values
    )

    if file_nodata_metadata is not None:
        nodata_values.append(
            file_nodata_metadata
        )

    if file_nodata is not None:
        nodata_values.append(
            file_nodata
        )

    nodata_unique = []

    for nodata_value in nodata_values:

        try:
            nodata_value = float(
                nodata_value
            )
        except (TypeError, ValueError):
            continue

        already_available = False

        for existing_value in nodata_unique:

            if np.isnan(nodata_value) and np.isnan(existing_value):
                already_available = True
                break

            if np.isclose(
                    nodata_value,
                    existing_value,
                    equal_nan=True,
            ):
                already_available = True
                break

        if not already_available:
            nodata_unique.append(
                nodata_value
            )

    return nodata_unique
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to validate the TIFF grid orientation
def validate_grid_orientation(
        transform,
) -> None:
    """
    Validate the expected north-up raster orientation.
    """

    if not np.isclose(transform.b, 0.0):
        logger.warning(
            " ===> TIFF transform contains x rotation: %s.",
            transform.b,
        )

    if not np.isclose(transform.d, 0.0):
        logger.warning(
            " ===> TIFF transform contains y rotation: %s.",
            transform.d,
        )

    if transform.a <= 0:
        logger.warning(
            " ===> TIFF x resolution is not positive: %s.",
            transform.a,
        )

    if transform.e >= 0:
        logger.warning(
            " ===> TIFF y resolution is not negative: %s. "
            "A north-up raster normally has a negative y resolution.",
            transform.e,
        )
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read soil-moisture TIFF file
def read_file_tiff_cnr(
        file_path: str,
        file_band: int = 1,
        file_variable: str = TIFF_DEFAULT_VARIABLE,
        file_crs: str = TIFF_DEFAULT_CRS,
        file_scale_factor: float = 1.0,
        file_offset: float = 0.0,
        file_valid_min: Optional[float] = None,
        file_valid_max: Optional[float] = None,
        file_nodata: Optional[float] = TIFF_DEFAULT_NODATA,
        no_data_values: Optional[Union[Any, Sequence[Any]]] = None,
        parse_time: bool = True,
        time_pattern: str = "%Y%m%d%H",
) -> Dict[str, Any]:
    """
    Read and organize a soil-moisture TIFF dataset.

    Parameters
    ----------
    file_path : str
        Input TIFF file.

    file_band : int
        TIFF band to read. The uploaded soil-moisture file has one band.

    file_variable : str
        Logical variable name assigned to the output dataset.

    file_crs : str
        Fallback CRS. The uploaded TIFF does not define its CRS, but the
        coordinate range and grid match EPSG:4326.

    file_scale_factor : float
        Scale factor applied to raw values.

    file_offset : float
        Offset applied after scaling.

    file_valid_min, file_valid_max : float, optional
        Valid data range after scaling.

    file_nodata : float, optional
        Configured nodata value. The uploaded TIFF uses -9999 in the data,
        although its TIFF nodata metadata contains NaN.

    no_data_values : scalar or sequence, optional
        Additional nodata values.

    parse_time : bool
        Parse the reference time from the filename.

    time_pattern : str
        Datetime pattern used in the filename.

    Returns
    -------
    dict
        Organized geographical dataset.
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(
            f"TIFF file '{file_path}' was not found."
        )

    if not os.path.isfile(file_path):
        raise RuntimeError(
            f"TIFF path '{file_path}' is not a regular file."
        )

    try:

        with rasterio.open(file_path) as file_handle:

            if file_band < 1:
                raise ValueError(
                    f"Band index must be greater than zero. "
                    f"Requested band: {file_band}."
                )

            if file_band > file_handle.count:
                raise ValueError(
                    f"Band {file_band} was requested, but TIFF file "
                    f"'{file_path}' contains only "
                    f"{file_handle.count} band(s)."
                )

            # Read raw values without applying the rasterio mask because
            # the file stores NaN in its nodata metadata but uses -9999
            # inside the raster values.
            values_raw = file_handle.read(
                file_band
            )

            transform = file_handle.transform
            width = file_handle.width
            height = file_handle.height

            file_crs_metadata = file_handle.crs
            file_nodata_metadata = file_handle.nodata

            bounds = file_handle.bounds
            resolution = file_handle.res

            profile = file_handle.profile.copy()
            global_tags = file_handle.tags()
            band_tags = file_handle.tags(
                file_band
            )

            band_description = file_handle.descriptions[
                file_band - 1
            ]

            band_dtype = file_handle.dtypes[
                file_band - 1
            ]

            file_driver = file_handle.driver

            file_compression = profile.get(
                "compress"
            )

            block_shapes = list(
                file_handle.block_shapes
            )

    except Exception as exc:
        raise RuntimeError(
            f"Unable to read TIFF file '{file_path}'."
        ) from exc

    # validate grid orientation
    validate_grid_orientation(
        transform=transform
    )

    # define CRS
    configured_crs = CRS.from_user_input(
        file_crs
    )

    if file_crs_metadata is None:

        crs = configured_crs

        logger.warning(
            " ===> TIFF file '%s' does not define a CRS. "
            "Configured CRS '%s' will be used.",
            file_path,
            configured_crs,
        )

    else:

        crs = file_crs_metadata

        if crs != configured_crs:
            logger.warning(
                " ===> CRS mismatch; TIFF CRS is '%s', configured CRS is "
                "'%s'. TIFF CRS will be used.",
                crs,
                configured_crs,
            )

    # collect nodata values
    nodata_values = get_tiff_nodata_values(
        file_nodata_metadata=file_nodata_metadata,
        file_nodata=file_nodata,
        no_data_values=no_data_values,
    )

    # process values
    values = process_values(
        values=values_raw,
        dtype=np.float64,
        scale_factor=file_scale_factor,
        offset=file_offset,
        valid_min=file_valid_min,
        valid_max=file_valid_max,
        nodata_values=nodata_values,
        file_nodata=file_nodata,
    )

    # create longitude and latitude coordinates
    longitude, latitude = create_coords(
        transform=transform,
        width=width,
        height=height,
    )

    if longitude.shape != values.shape:
        raise ValueError(
            f"Longitude shape {longitude.shape} does not match "
            f"TIFF values shape {values.shape}."
        )

    if latitude.shape != values.shape:
        raise ValueError(
            f"Latitude shape {latitude.shape} does not match "
            f"TIFF values shape {values.shape}."
        )

    # create finite-data mask
    finite_mask = np.isfinite(
        values
    )

    finite_values = values[
        finite_mask
    ]

    # compute statistics
    if finite_values.size > 0:

        value_min = float(
            np.nanmin(finite_values)
        )

        value_max = float(
            np.nanmax(finite_values)
        )

        value_mean = float(
            np.nanmean(finite_values)
        )

        value_median = float(
            np.nanmedian(finite_values)
        )

        value_std = float(
            np.nanstd(finite_values)
        )

    else:

        value_min = np.nan
        value_max = np.nan
        value_mean = np.nan
        value_median = np.nan
        value_std = np.nan

        logger.warning(
            " ===> TIFF dataset '%s' contains no valid values.",
            file_path,
        )

    # parse time from filename
    if parse_time:

        reference_time = parse_time_from_filename(
            file_path=file_path,
            time_pattern=time_pattern,
        )

    else:
        reference_time = None

    var_dataset = {

        "values": values,
        "longitude": longitude,
        "latitude": latitude,

        "transform": transform,
        "crs": crs,

        "width": width,
        "height": height,
        "shape": (
            height,
            width,
        ),

        "bounds": {
            "left": float(bounds.left),
            "bottom": float(bounds.bottom),
            "right": float(bounds.right),
            "top": float(bounds.top),
        },

        "resolution": {
            "x": float(resolution[0]),
            "y": float(resolution[1]),
        },

        "variable": file_variable,
        "band": file_band,
        "band_count": profile.get("count"),
        "band_description": band_description,
        "dtype": band_dtype,

        "reference_time": reference_time,

        "nodata": file_nodata,
        "file_nodata": file_nodata_metadata,
        "nodata_values": nodata_values,

        "scale_factor": file_scale_factor,
        "offset": file_offset,
        "valid_min": file_valid_min,
        "valid_max": file_valid_max,

        "file_path": file_path,
        "driver": file_driver,
        "compression": file_compression,
        "block_shapes": block_shapes,

        "profile": profile,
        "tags": global_tags,
        "band_tags": band_tags,

        "valid_count": int(
            np.count_nonzero(finite_mask)
        ),

        "invalid_count": int(
            np.count_nonzero(~finite_mask)
        ),

        "value_min": value_min,
        "value_max": value_max,
        "value_mean": value_mean,
        "value_median": value_median,
        "value_std": value_std,
    }

    return var_dataset
# ----------------------------------------------------------------------------------------------------------------------
