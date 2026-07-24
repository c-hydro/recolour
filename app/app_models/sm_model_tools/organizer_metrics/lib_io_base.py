
"""
Library Features:

Name:          lib_geo
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

from copy import deepcopy
from typing import Any, Dict, Mapping, Optional, Tuple
from affine import Affine

from config_info import LOGGER_NAME, TIME_FMT_CLI

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to process values
def process_values(
        values: np.ndarray,
        dtype: type[np.generic] = np.float64,
        scale_factor: float = 1.0,
        offset: float = 0.0,
        valid_min: Optional[float] = None,
        valid_max: Optional[float] = None,
        nodata_values: Optional[list[float]] = None,
        file_nodata: Optional[float] = None,) -> np.ndarray:

    # cast values to dtype
    values_out = np.asarray(values, dtype=dtype,).copy()

    # apply not finite values to create invalid mask
    invalid_mask = ~np.isfinite(values_out)
    for nodata_value in nodata_values:
        invalid_mask |= np.isclose(values_out, nodata_value, equal_nan=False,)

    # apply no data values to create invalid mask
    if file_nodata is not None:
        try:
            if np.isfinite(file_nodata):
                invalid_mask |= np.isclose( values_out,file_nodata, equal_nan=False,)
        except TypeError:
            logger.warning(f" ===> File nodata value {file_nodata} is not numeric.")

    # use min and mask to create invalid mask
    if valid_min is not None:
        invalid_mask |= values_out < valid_min
    if valid_max is not None:
        invalid_mask |= values_out > valid_max

    # apply scale factor and offset
    values_out = (values_out * scale_factor + offset)
    # apply invalid mask
    values_out[invalid_mask] = np.nan

    return values_out
# ----------------------------------------------------------------------------------------------------------------------
