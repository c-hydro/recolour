# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import rasterio
import xarray as xr
import pandas as pd

from lib_utils_decoretors import iterate_items

from lib_utils_info import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read file tif
@iterate_items(iter_types=(list, tuple), strict_zip=True)
def read_tiff(file_name: str,
              file_time: pd.Timestamp = None, variable: str = 'ssm_filtered') -> (xr.DataArray, None):
    if os.path.exists(file_name):
        with rasterio.open(file_name) as src:
            descriptions = list(src.descriptions)

            # Try matching by description first
            if variable in descriptions:
                band_index = descriptions.index(variable) + 1
            else:
                try:
                    band_index = int(variable)
                except ValueError:
                    raise ValueError(f"Variable '{variable}' not found in band descriptions.")

            band_data = src.read(band_index)
            transform = src.transform

            coords = {
                "latitude": [transform[5] + i * transform[4] for i in range(src.height)],
                "longitude": [transform[2] + i * transform[0] for i in range(src.width)],
            }

            # Add time coordinate if provided
            if file_time is not None:
                coords["time"] = [file_time]
            else:
                alg_logger.warning(f" ===> Time not defined for file: {file_name}")

            file_da = xr.DataArray(
                band_data if file_time is None else band_data[None, :, :],  # expand dims if time is given
                dims=("latitude", "longitude") if file_time is None else ("time", "latitude", "longitude"),
                coords=coords,
                attrs={
                    "crs": src.crs.to_string(),
                    "transform": transform,
                    "band": descriptions[band_index - 1] if descriptions[band_index - 1] else f"Band {band_index}"
                }
            )

    else:
        alg_logger.warning(f" ===> File data tiff not found: {file_name}")
        file_da = None

    return file_da
# ----------------------------------------------------------------------------------------------------------------------
