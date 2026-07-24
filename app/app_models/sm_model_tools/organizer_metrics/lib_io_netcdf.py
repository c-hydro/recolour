"""
Library Features:

Name:          lib_io_netcdf
Author(s):     Fabio Delogu
Date:          '20260716'
Version:       '1.0.0'

Purpose:
    Decompress, validate, load and organize HMC state-grid NetCDF datasets.

Supported formats:
    - NetCDF
    - gzip-compressed NetCDF (.nc.gz)

Supported HMC variables:
    Two-dimensional variables:
        - AgeS
        - AlbedoS
        - DFE
        - HydroLevel
        - LST
        - RhoS
        - Routing
        - SWE
        - VRet
        - VTot
        - WS
        - WTLevel

    Three-dimensional variables:
        - T24
        - T_1Days
        - T_5Days
        - Tmk
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import gzip
import logging
import os
import shutil
import tempfile
import numpy as np
import xarray as xr

from rasterio.crs import CRS
from typing import Any, Dict, Optional, Tuple

from lib_io_base import process_values
from lib_utils_hmc import (select_hmc_layer, parse_hmc_reference_time,
                           create_hmc_transform, get_hmc_nodata_values, clean_attributes)
from lib_utils_zip import unzip_file_gzip
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
# method to prepare compressed or uncompressed input
def prepare_file_netcdf(
        file_path: str,
        file_compression: Optional[bool] = None,
        compression_folder: Optional[str] = None,
        compression_temporary: bool = True,
        compression_overwrite: bool = False,
) -> Tuple[str, Optional[str]]:
    """
    Prepare a NetCDF input file.

    Compression can be detected automatically from the '.gz' suffix.

    Returns
    -------
    tuple
        Reading path and optional temporary folder.
    """

    if file_compression is None:
        file_compression = file_path.lower().endswith(".gz")

    if not file_compression:
        return file_path, None

    temporary_folder = None

    if compression_temporary:
        temporary_folder = tempfile.mkdtemp(prefix="hmc_file_")
        reading_path = unzip_file_gzip(file_path=file_path, output_folder=temporary_folder, overwrite=True,)
    else:
        reading_path = unzip_file_gzip(file_path=file_path, output_folder=compression_folder,
                                       overwrite=compression_overwrite)

    return reading_path, temporary_folder
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read an HMC NetCDF state-grid
def read_file_netcdf_hmc(
        file_path: str,
        file_variable: str = "VTot",
        file_layer: Optional[int] = None,
        file_crs: str = HMC_DEFAULT_CRS,
        file_scale_factor: Optional[float] = None,
        file_offset: Optional[float] = None,
        file_valid_min: Optional[float] = None,
        file_valid_max: Optional[float] = None,
        file_nodata: Optional[float] = None,
        no_data_values: Any = None,
        file_compression: Optional[bool] = None,
        compression_folder: Optional[str] = None,
        compression_temporary: bool = True,
        compression_overwrite: bool = False,
        use_file_scale_factor: bool = True,
        engine: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Read an HMC state-grid NetCDF dataset.

    Parameters
    ----------
    file_path : str
        Input .nc or .nc.gz file.

    file_variable : str
        Variable to read. Default is VTot.

    file_layer : int, optional
        Layer index for 3D variables such as T24, T_1Days, T_5Days and Tmk.

    file_crs : str
        Configured CRS. HMC state-grid files use EPSG:4326.

    file_scale_factor : float, optional
        Configured scale factor. When omitted, the variable scale_factor
        attribute is used if use_file_scale_factor=True.

    file_offset : float, optional
        Configured offset. When omitted, the add_offset variable attribute
        is used if available.

    file_valid_min, file_valid_max : float, optional
        Valid data limits applied after scaling.

    file_nodata : float, optional
        Additional configured nodata value.

    no_data_values : scalar or list, optional
        Additional nodata values.

    file_compression : bool, optional
        Compression flag. When None, it is detected from the .gz suffix.

    compression_temporary : bool
        Decompress to a temporary directory and remove it after reading.

    use_file_scale_factor : bool
        Read scale_factor and add_offset directly from NetCDF attributes.

    engine : str, optional
        Xarray engine, for example netcdf4 or h5netcdf.

    Returns
    -------
    dict
        Organized HMC geographical dataset.
    """

    # check file availability
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"HMC NetCDF file '{file_path}' was not found.")

    reading_path = file_path
    temporary_folder = None
    dataset = None

    try:

        # prepare netcdf file (if compressed or not)
        reading_path, temporary_folder = prepare_file_netcdf(
            file_path=file_path, file_compression=file_compression,
            compression_folder=compression_folder, compression_temporary=compression_temporary,
            compression_overwrite=compression_overwrite,)

        # define kwargs to open dataset
        open_kwargs = {"decode_cf": False,"mask_and_scale": False,}

        if engine is not None:
            open_kwargs["engine"] = engine

        # open dataset
        dataset = xr.open_dataset(reading_path, **open_kwargs,)

        # check if variable available in the dataset
        if file_variable not in dataset.data_vars:
            raise KeyError(
                f"Variable '{file_variable}' is not available in HMC file "
                f"'{file_path}'. Available data variables are: "
                f"{list(dataset.data_vars)}."
            )

        # check longitude and latitude
        if HMC_LONGITUDE_NAME not in dataset.variables:
            raise KeyError(f"HMC longitude variable '{HMC_LONGITUDE_NAME}' is missing.")
        if HMC_LATITUDE_NAME not in dataset.variables:
            raise KeyError(f"HMC latitude variable '{HMC_LATITUDE_NAME}' is missing.")

        # get variable values
        data_array = dataset[file_variable]

        # manage variable values
        data_array, layer_dimension, layer_index = select_hmc_layer(data_array=data_array, file_layer=file_layer,)
        data_array = data_array.squeeze(drop=True)

        if data_array.ndim != 2:
            raise ValueError(
                f"HMC variable '{file_variable}' must be 2D after layer "
                f"selection. Current dimensions are {dict(data_array.sizes)}."
            )

        if HMC_DIM_Y not in data_array.dims:
            raise ValueError(
                f"Expected HMC dimension '{HMC_DIM_Y}' is not available in variable '{file_variable}'.")
        if HMC_DIM_X not in data_array.dims:
            raise ValueError(
                f"Expected HMC dimension '{HMC_DIM_X}' is not available in variable '{file_variable}'.")

        # Ensure the output order is always south_north, west_east.
        data_array = data_array.transpose(HMC_DIM_Y, HMC_DIM_X,)

        values_raw = np.asarray(data_array.values)
        values_raw = np.flipud(values_raw)

        longitude = np.asarray(dataset[HMC_LONGITUDE_NAME].values, dtype=np.float64,)
        latitude = np.asarray(dataset[HMC_LATITUDE_NAME].values, dtype=np.float64,)

        if longitude.shape != values_raw.shape:
            raise ValueError(
                f"Longitude shape {longitude.shape} does not match variable shape {values_raw.shape}.")

        if latitude.shape != values_raw.shape:
            raise ValueError(
                f"Latitude shape {latitude.shape} does not match variable shape {values_raw.shape}.")

        variable_attributes = clean_attributes(dict(data_array.attrs))
        global_attributes = clean_attributes(dict(dataset.attrs))
        encoding = clean_attributes(dict(data_array.encoding))

        nodata_values = get_hmc_nodata_values(
            data_array=data_array,
            file_nodata=file_nodata,
            no_data_values=no_data_values,
        )

        # scale factor
        if file_scale_factor is None:
            if use_file_scale_factor:
                file_scale_factor = variable_attributes.get("scale_factor", 1.0,)
            else:
                file_scale_factor = 1.0
        file_scale_factor = float(file_scale_factor)

        # offset
        if file_offset is None:
            if use_file_scale_factor:
                file_offset = variable_attributes.get("add_offset",0.0 )
            else:
                file_offset = 0.0
        file_offset = float(file_offset)

        # reference time
        reference_time = parse_hmc_reference_time(global_attributes=global_attributes)
        # transform
        transform = create_hmc_transform(global_attributes=global_attributes, longitude=longitude, latitude=latitude,)

    except Exception as exc:
        raise RuntimeError(f"Unable to read HMC NetCDF file '{file_path}'.") from exc

    finally:

        if dataset is not None:
            dataset.close()
        if temporary_folder is not None:
            shutil.rmtree(temporary_folder, ignore_errors=True,)

    # process HMC values
    values = process_values(
        values=values_raw,dtype=np.float64,
        scale_factor=file_scale_factor, offset=file_offset,
        valid_min=file_valid_min, valid_max=file_valid_max,
        nodata_values=nodata_values, file_nodata=file_nodata,)

    finite_mask = np.isfinite(values)
    finite_values = values[finite_mask]

    if finite_values.size > 0:
        value_min = float(np.nanmin(finite_values))
        value_max = float(np.nanmax(finite_values))
        value_mean = float(np.nanmean(finite_values))
    else:
        value_min, value_max, value_mean = np.nan, np.nan, np.nan
        logger.warning(f" ===> HMC variable {file_variable} contains no valid values.")

    height, width = values.shape
    configured_crs = CRS.from_user_input(file_crs)

    var_dataset = {

        "values": values,
        "longitude": longitude,
        "latitude": latitude,

        "transform": transform,
        "crs": configured_crs,
        "width": width,
        "height": height,

        "variable": file_variable,
        "layer_dimension": layer_dimension,
        "layer_index": layer_index,

        "dimensions": {HMC_DIM_Y: height, HMC_DIM_X: width,},

        "reference_time": reference_time,
        "reference_time_raw": global_attributes.get(HMC_TIME_ATTRIBUTE),

        "scale_factor": file_scale_factor,
        "offset": file_offset,

        "nodata": file_nodata,
        "nodata_values": nodata_values,

        "file_path": file_path,
        "compression": (file_path.lower().endswith(".gz") if file_compression is None else file_compression),

        "attributes": variable_attributes,
        "global_attributes": global_attributes,
        "encoding": encoding,

        "valid_count": int(np.count_nonzero(finite_mask)),
        "invalid_count": int(np.count_nonzero(~finite_mask)),

        "value_min": value_min,
        "value_max": value_max,
        "value_mean": value_mean,
    }

    return var_dataset
# ----------------------------------------------------------------------------------------------------------------------
