"""
Library Features:

Name:           lib_utils_geo
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260715'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np

from typing import Any, Dict, Mapping, Optional, Tuple
from affine import Affine

from config_info import LOGGER_NAME, TIME_FMT_CLI

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to check grids
def check_grids(
        reference_data: Mapping[str, Any],
        dataset_data: Mapping[str, Any],
        transform_tolerance: float = 1e-8, raise_error: bool = True) -> bool:

    # initialize errors
    errors = []

    # get reference and dataset labels
    reference_label = f"{reference_data['metadata']['group']}_{reference_data['metadata']['key']}"
    dataset_label = f"{dataset_data['metadata']['group']}_{dataset_data['metadata']['key']}"

    # check shape
    if reference_data["values"].shape != dataset_data["values"].shape:
        errors.append(f"shape differs: {reference_data['values'].shape} != {dataset_data['values'].shape}")
    # check crs
    if reference_data["crs"] != dataset_data["crs"]:
        errors.append(f"CRS differs: {reference_data['crs']} != {dataset_data['crs']}")

    # check transform reference
    reference_transform_obj = reference_data["transform"]
    # get transform dataset
    dataset_transform_obj = dataset_data["transform"]


    # check transform objects
    if reference_transform_obj is None or dataset_transform_obj is None:
        if reference_transform_obj is not dataset_transform_obj:
            errors.append("transform is missing from one dataset")
    else:

        # compute reference info
        reference_transform = np.asarray(reference_transform_obj[:6], dtype=np.float64, )
        reference_resolution = min(abs(reference_transform_obj.a), abs(reference_transform_obj.e),)
        # compute dataset info
        dataset_transform = np.asarray(dataset_transform_obj[:6], dtype=np.float64)
        dataset_resolution = min(abs(dataset_transform_obj.a), abs(dataset_transform_obj.e),)

        # compute min resolution
        minimum_resolution = min(reference_resolution, dataset_resolution,)

        # At least the configured absolute tolerance, or 0.1% of a pixel.
        effective_tolerance = max(transform_tolerance, minimum_resolution * 1e-3,)

        # check for errors
        if not np.allclose(reference_transform, dataset_transform, atol=effective_tolerance,rtol=0.0,):

            # compute checks
            transform_difference = (dataset_transform - reference_transform)
            max_difference = float(np.max(np.abs(transform_difference)))
            pixel_difference = (max_difference / minimum_resolution)

            # store errors
            errors.append(
                "transform differs: "
                f"reference={reference_transform.tolist()}, "
                f"dataset={dataset_transform.tolist()}, "
                f"difference={transform_difference.tolist()}, "
                f"max_difference={max_difference:.16e}, "
                f"pixel_difference={pixel_difference:.16e}, "
                f"tolerance={effective_tolerance:.16e}"
            )

    # exit if grids are compatibles
    if not errors:
        return True
    else:
        # return errors
        error_message = (
                f" ===> Grid '{dataset_label}' is not compatible with reference grid '{reference_label}': "
                + "; ".join(errors))

        # exit with different code
        if raise_error:
            logger.error(error_message)
            raise ValueError(error_message)
        else:
            logger.warning(error_message)

        return False

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create coords
def create_coords(transform: Affine, width: int, height: int,) -> Tuple[np.ndarray, np.ndarray]:

    rows, cols = np.meshgrid(np.arange(height), np.arange(width), indexing="ij",)

    longitude = (transform.c + (cols + 0.5) * transform.a + (rows + 0.5) * transform.b)
    latitude = (transform.f + (cols + 0.5) * transform.d + (rows + 0.5) * transform.e)

    return longitude, latitude
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to compute max volume from curve number
def compute_volume_max(
        values_cn,
        nodata_values=(-9999, -9998),
):
    """
    Compute maximum water volume from Curve Number.

    Formula:
        Vmax = 25400 / CN - 254

    Parameters
    ----------
    values_cn : array-like
        Curve Number values.

    nodata_values : tuple
        Values treated as missing data.

    Returns
    -------
    np.ndarray
        Maximum volume values. Invalid cells are returned as NaN.
    """

    # convert input values
    values_cn = np.asarray(
        values_cn,
        dtype=np.float64,
    )

    # initialize output
    values_vmax = np.full(
        values_cn.shape,
        np.nan,
        dtype=np.float64,
    )

    # identify valid cells
    valid_mask = np.isfinite(
        values_cn
    )

    # exclude configured nodata values
    for nodata_value in nodata_values:
        valid_mask &= ~np.isclose(
            values_cn,
            nodata_value,
        )

    # CN must be greater than zero and lower than 100.
    #
    # CN == 100 produces Vmax == 0 and cannot be used
    # as a denominator for soil-moisture computation.
    valid_mask &= (
        (values_cn > 0.0)
        & (values_cn < 100.0)
    )

    # compute maximum volume only for valid cells
    values_vmax[valid_mask] = (
        25400.0 / values_cn[valid_mask]
        - 254.0
    )

    return values_vmax
# ----------------------------------------------------------------------------------------------------------------------
