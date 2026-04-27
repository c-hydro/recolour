"""
Library Features:

Name:           lib_utils_io
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json
import rasterio
import numpy as np

from lib_utils_basic import make_folder
from config_utils import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read settings file
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f'File "{file_name}" not found. Exit')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to read target grid
def load_target_grid(mask_file, mask_band=1):

    with rasterio.open(mask_file) as dataset:
        mask = dataset.read(mask_band).astype(np.float32)
        profile = dataset.profile.copy()
        transform = dataset.transform
        height = dataset.height
        width = dataset.width

        cols = np.arange(width, dtype=np.float64)
        rows = np.arange(height, dtype=np.float64)

        x_coords = transform.c + (cols + 0.5) * transform.a
        y_coords = transform.f + (rows + 0.5) * transform.e

        grid_lons, grid_lats = np.meshgrid(x_coords, y_coords)

        logger.info(f" ----> Mask file: {mask_file}")
        logger.info(f" ----> Mask CRS: {dataset.crs}")
        logger.info(f" ----> Mask shape: {mask.shape}")

    domain_mask = np.isfinite(mask) & (mask == 1)

    logger.info(f" ----> Valid domain pixels: {int(domain_mask.sum())}")
    logger.info(f" ----> Invalid domain pixels: {int((~domain_mask).sum())}")

    return (
        grid_lons.astype(np.float64),
        grid_lats.astype(np.float64),
        domain_mask,
        profile
    )
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to write output map
def write_output_map(
        destination_file,
        soil_moisture_map,
        time_lag_map,
        profile,
        value_var,
        fill_value):

    out_profile = profile.copy()
    out_profile.update(
        count=2,
        dtype="float32",
        compress="deflate",
        nodata=np.nan if np.isnan(fill_value) else float(fill_value),
    )

    folder_name = os.path.dirname(destination_file)
    if folder_name:
        make_folder(folder_name)

    with rasterio.open(destination_file, "w", **out_profile) as dst:
        dst.write(soil_moisture_map.astype(np.float32), 1)
        dst.write(time_lag_map.astype(np.float32), 2)
        dst.set_band_description(1, value_var)
        dst.set_band_description(2, "time_lag_hours")
# ----------------------------------------------------------------------------------------------------------------------
