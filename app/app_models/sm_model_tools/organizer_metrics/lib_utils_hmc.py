"""
Library Features:

Name:           lib_utils_hmc
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260715'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging

from datetime import datetime
from typing import Any, Dict, Optional, Tuple

import numpy as np
import xarray as xr

from affine import Affine

from lib_utils_geo import compute_volume_max
from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)

# constants
HMC_LONGITUDE_NAME = "Longitude"
HMC_LATITUDE_NAME = "Latitude"

HMC_DIM_Y = "south_north"
HMC_DIM_X = "west_east"

HMC_DEFAULT_CRS = "EPSG:4326"

HMC_TIME_ATTRIBUTE = "time_coverage_end"
HMC_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"

HMC_DEFAULT_NODATA = -8.999999815811072e15
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute HMC soil moisture using VTot and Vmax derived from CN
def compute_hmc_soil_moisture(
        values_vtot,
        values_cn,
        scale_factor=1.0,
        clip_values=True,
        nodata_values=(-9999, -9998),
):
    """
    Compute HMC relative soil moisture using:

        SM = VTot / Vmax

    where:

        Vmax = 25400 / CN - 254

    Parameters
    ----------
    values_vtot : array-like
        Total water volume.

    values_cn : array-like
        Curve Number grid.

    scale_factor : float
        Output scale:
            1.0   -> soil moisture in the range 0-1
            100.0 -> soil moisture in the range 0-100

    clip_values : bool
        Clip relative soil moisture to the physical range 0-1 before
        applying the scale factor.

    nodata_values : tuple
        Values treated as missing in both input datasets.

    Returns
    -------
    np.ndarray
        Soil-moisture grid. Invalid cells are returned as NaN.
    """

    # convert input values
    values_vtot = np.asarray(values_vtot, dtype=np.float64,)
    values_cn = np.asarray(values_cn,dtype=np.float64,)

    # check grid shapes
    if values_vtot.shape != values_cn.shape:
        raise ValueError(f"VTot shape {values_vtot.shape} differs from Curve Number shape {values_cn.shape}.")

    # compute maximum volume
    values_vmax = compute_volume_max(values_cn=values_cn,nodata_values=nodata_values,)

    # initialize output
    values_sm = np.full(values_vtot.shape,np.nan,dtype=np.float64,)

    # identify valid cells
    valid_mask = (np.isfinite(values_vtot)& np.isfinite(values_vmax) & (values_vmax > 0.0))

    # exclude configured nodata values from VTot
    for nodata_value in nodata_values:
        valid_mask &= ~np.isclose(values_vtot,nodata_value,)

    # compute relative soil moisture
    values_sm[valid_mask] = (values_vtot[valid_mask]/ values_vmax[valid_mask])

    # constrain values to their physical range
    if clip_values:
        values_sm[valid_mask] = np.clip(values_sm[valid_mask],0.0,1.0,)

    # apply output scaling
    values_sm[valid_mask] *= float(scale_factor)

    vtot_finite = np.isfinite(values_vtot)
    cn_finite = np.isfinite(values_cn)
    vmax_finite = np.isfinite(values_vmax)
    sm_finite = np.isfinite(values_sm)

    vtot_valid = np.isfinite(values_vtot)
    for nodata_value in nodata_values:
        vtot_valid &= ~np.isclose(values_vtot, nodata_value)

    logger.info(" ::: SM statistics")
    logger.info("  CN valid    : %d", np.count_nonzero(cn_finite))
    logger.info("  VTOT valid  : %d", np.count_nonzero(vtot_valid))
    logger.info("  VMAX valid  : %d", np.count_nonzero(vmax_finite))
    logger.info("  MASK valid  : %d", np.count_nonzero(valid_mask))
    logger.info("  SM valid    : %d", np.count_nonzero(sm_finite))


    logger.info("  CN    : finite=%d/%d  min=%.3f  max=%.3f",
        np.count_nonzero(cn_finite),
        values_cn.size,
        np.nanmin(values_cn[cn_finite]) if np.any(cn_finite) else np.nan,
        np.nanmax(values_cn[cn_finite]) if np.any(cn_finite) else np.nan,
    )

    logger.info("  VTot  : finite=%d/%d  min=%.3f  max=%.3f",
        np.count_nonzero(vtot_finite),
        values_vtot.size,
        np.nanmin(values_vtot[vtot_finite]) if np.any(vtot_finite) else np.nan,
        np.nanmax(values_vtot[vtot_finite]) if np.any(vtot_finite) else np.nan,
    )

    logger.info("  VMax  : finite=%d/%d  min=%.3f  max=%.3f",
        np.count_nonzero(vmax_finite),
        values_vmax.size,
        np.nanmin(values_vmax[vmax_finite]) if np.any(vmax_finite) else np.nan,
        np.nanmax(values_vmax[vmax_finite]) if np.any(vmax_finite) else np.nan,
    )

    logger.info("  SM    : finite=%d/%d  min=%.3f  max=%.3f",
        np.count_nonzero(sm_finite),
        values_sm.size,
        np.nanmin(values_sm[sm_finite]) if np.any(sm_finite) else np.nan,
        np.nanmax(values_sm[sm_finite]) if np.any(sm_finite) else np.nan,
    )

    logger.info("  Pixels removed from CN mask      : %d",
        np.count_nonzero(cn_finite & ~sm_finite),
    )

    logger.info("  Pixels valid in both CN and SM   : %d",
        np.count_nonzero(cn_finite & sm_finite),
    )

    logger.info("  Valid computation mask           : %d",
        np.count_nonzero(valid_mask),
    )

    logger.info(
        "  Lost because CN invalid      : %d",
        np.count_nonzero(vtot_valid & ~np.isfinite(values_cn)),
    )

    logger.info("  Lost because Vmax invalid    : %d",
        np.count_nonzero(vtot_valid & ~np.isfinite(values_vmax)),
    )

    logger.info("  Lost because Vmax <= 0       : %d",
        np.count_nonzero(
            vtot_valid &
            np.isfinite(values_vmax) &
            (values_vmax <= 0)
        ),
    )

    logger.info("  Final valid SM              : %d",
        np.count_nonzero(valid_mask),
    )

    return values_sm
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to parse HMC reference time
def parse_hmc_reference_time(global_attributes: Dict[str, Any],) -> Optional[datetime]:
    """
    Parse the HMC time_coverage_end global attribute.
    """

    time_raw = global_attributes.get(HMC_TIME_ATTRIBUTE)
    if time_raw is None:
        return None

    time_raw = convert_attribute(time_raw)
    try:
        return datetime.strptime(str(time_raw), HMC_TIME_FORMAT,)
    except ValueError:
        logger.warning(f" ===> Unable to parse HMC reference time {time_raw} using format {HMC_TIME_FORMAT}.")
        return None
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create affine transform from HMC global attributes
def create_hmc_transform(global_attributes: Dict[str, Any], longitude: np.ndarray, latitude: np.ndarray,) -> Affine:
    """
    Create the affine transform using HMC global attributes.

    The file defines:
        xllcorner
        yllcorner
        xcellsize
        ycellsize
        nrows
        ncols

    HMC arrays are ordered from north to south. Therefore, the raster
    transform uses a negative vertical resolution.
    """

    xllcorner = global_attributes.get("xllcorner")
    yllcorner = global_attributes.get("yllcorner")
    xcellsize = global_attributes.get("xcellsize")
    ycellsize = global_attributes.get("ycellsize")
    nrows = global_attributes.get("nrows")

    required_values = [xllcorner, yllcorner, xcellsize, ycellsize, nrows,]

    if all(value is not None for value in required_values):

        xllcorner = float(convert_attribute(xllcorner))
        yllcorner = float(convert_attribute(yllcorner))
        xcellsize = float(convert_attribute(xcellsize))
        ycellsize = float(convert_attribute(ycellsize))
        nrows = int(convert_attribute(nrows))

        # yllcorner is the lower-left boundary. Rasterio requires
        # the upper-left boundary.
        upper_left_y = (yllcorner + nrows * ycellsize)

        return Affine(xcellsize,0.0,xllcorner,0.0, -ycellsize, upper_left_y,)

    logger.warning(" ===> HMC grid attributes are incomplete. Transform will be estimated from Longitude and Latitude.")

    longitude_step = float(np.nanmedian(np.diff(longitude[0, :])))
    latitude_step = float(np.nanmedian(np.diff(latitude[:, 0])))

    longitude_origin = (float(longitude[0, 0])- longitude_step / 2.0)
    latitude_origin = (float(latitude[0, 0])- latitude_step / 2.0)

    return Affine(longitude_step, 0.0, longitude_origin,0.0, latitude_step, latitude_origin,)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to retrieve file nodata values
def get_hmc_nodata_values(data_array: xr.DataArray, file_nodata: Optional[float] = None,
                          no_data_values: Any = None,) -> list:
    """
    Collect nodata values from configuration, attributes and encoding.
    """

    nodata_values = _normalize_list(no_data_values)
    for attribute_name in ["missing_value", "_FillValue",]:

        attribute_value = data_array.attrs.get(attribute_name)

        if attribute_value is not None:
            nodata_values.extend(_normalize_list(attribute_value))

        encoding_value = data_array.encoding.get(attribute_name)

        if encoding_value is not None:
            nodata_values.extend(_normalize_list(encoding_value))

    if file_nodata is not None:
        nodata_values.append(file_nodata)

    # The HMC missing value stored in this file is approximately -9e15.
    # Add it as a fallback even if a variable does not expose the attribute.
    nodata_values.append(HMC_DEFAULT_NODATA)

    nodata_unique = []
    for nodata_value in nodata_values:

        nodata_value = convert_attribute(nodata_value)

        try:
            nodata_value = float(nodata_value)
        except (TypeError, ValueError):
            continue

        already_available = any(np.isclose(nodata_value, existing_value, equal_nan=True,)
                                for existing_value in nodata_unique)

        if not already_available:
            nodata_unique.append(nodata_value)

    return nodata_unique
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select a layer from a three-dimensional HMC variable
def select_hmc_layer(
        data_array: xr.DataArray,
        file_layer: Optional[int] = None,
) -> Tuple[xr.DataArray, Optional[str], Optional[int]]:
    """
    Select a layer from a 3D HMC variable.

    Examples:
        T24       -> dimensions day1_steps, south_north, west_east
        T_1Days   -> dimensions day1_steps, south_north, west_east
        T_5Days   -> dimensions day5_steps, south_north, west_east
        Tmk       -> dimensions tmarked_steps, south_north, west_east
    """

    spatial_dimensions = {HMC_DIM_Y, HMC_DIM_X}

    additional_dimensions = [dimension_name
                             for dimension_name in data_array.dims if dimension_name not in spatial_dimensions]

    if not additional_dimensions:
        return data_array, None, None

    if len(additional_dimensions) > 1:
        raise ValueError(
            f"Variable '{data_array.name}' has more than one non-spatial dimension: {additional_dimensions}.")

    layer_dimension = additional_dimensions[0]
    layer_count = data_array.sizes[layer_dimension]

    if file_layer is None:
        if layer_count == 1:
            file_layer = 0
        else:
            raise ValueError(
                f"Variable '{data_array.name}' contains {layer_count} layers "
                f"along dimension '{layer_dimension}'. "
                f"Parameter 'file_layer' must be provided."
            )

    if not -layer_count <= file_layer < layer_count:
        raise IndexError(
            f"Layer index {file_layer} is outside dimension "
            f"'{layer_dimension}' having size {layer_count}."
        )

    data_array = data_array.isel({layer_dimension: file_layer})

    return data_array, layer_dimension, file_layer
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to normalize an object to a list
def _normalize_list(values: Any) -> list:
    """
    Normalize a scalar, tuple, set or array to a Python list.
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
# method to convert attributes to standard Python values
def convert_attribute(attribute_value: Any) -> Any:
    """
    Convert NumPy and byte attributes to ordinary Python objects.
    """

    if isinstance(attribute_value, bytes):
        return attribute_value.decode("utf-8", errors="replace")

    if isinstance(attribute_value, np.ndarray):

        if attribute_value.size == 1:
            return convert_attribute(attribute_value.reshape(-1)[0])

        return [convert_attribute(value) for value in attribute_value.tolist()]

    if isinstance(attribute_value, np.generic):
        return attribute_value.item()

    return attribute_value
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to clean attributes dictionary
def clean_attributes(attributes: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert all attributes to standard Python values.
    """
    return {key: convert_attribute(value) for key, value in attributes.items()}
# ----------------------------------------------------------------------------------------------------------------------
