
"""
Library Features:

Name:          lib_utils_datasets
Author(s):     Fabio Delogu
Date:          '20260715'
Version:       '1.0.0'

Purpose:
    Validate datasets

Supported types:
    - ascii_grid
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import string
import numpy as np
import pandas as pd

from rasterio.crs import CRS
from copy import deepcopy
from typing import Any, Dict, Mapping, Optional

from lib_io_ascii import read_file_ascii_grid
from lib_io_netcdf import read_file_netcdf_hmc
from lib_io_tiff import read_file_tiff_cnr
from lib_utils_hmc import compute_hmc_soil_moisture

from config_info import LOGGER_NAME, TIME_FMT_CLI

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# constants
SUPPORTED_DATASET_TYPES = {"ascii_grid", "netcdf_hmc_state", "tiff", "tif"}
DEFAULT_TYPE = "ascii_grid"
DEFAULT_BAND = 1
DEFAULT_CRS = "EPSG:4326"
DEFAULT_SCALE_FACTOR = 1.0
DEFAULT_OFFSET = 0.0
DEFAULT_NODATA_VALUES = [-9999.0]
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to load one dataset according to its type
def read_datasets(
        dataset_cfg: Mapping[str, Any],
        group_name: str = 'NA', dataset_key: str = 'NA', debug: bool = False,
        **kwargs) -> Dict[str, Any]:

    # get information
    dataset_type = dataset_cfg["type"]
    dataset_label = f"{group_name}_{dataset_key}"
    file_variable = dataset_cfg["variable"]
    file_path = dataset_cfg["file_path"]
    file_band = dataset_cfg["band"]
    file_crs = dataset_cfg["crs"]
    file_compression = dataset_cfg["compression"]
    file_scale_factor = dataset_cfg["scale_factor"]
    file_offset = dataset_cfg["offset"]
    file_valid_min = dataset_cfg["valid_min"]
    file_valid_max = dataset_cfg["valid_max"]
    file_nodata_values = dataset_cfg["nodata_values"]
    file_scaling = dataset_cfg["scaling"]

    # get group geo in set in kwargs
    group_geo = kwargs.get('group_geo', {})

    # set extra metadata
    extra_cfg = {'group': group_name, 'key': dataset_key, 'label': dataset_label}

    # info start
    logger.info(f" --------> Read datasets '{dataset_label}' from {file_path} ...")

    # check datasets type
    if dataset_type == "ascii_grid":

        # get variable from ascii grid file
        datasets_obj = read_file_ascii_grid(
            file_path=file_path, file_band=file_band, file_crs=file_crs,
            file_valid_min=file_valid_min, file_valid_max=file_valid_max,
            file_scale_factor=file_scale_factor, file_offset=file_offset,
            no_data_values=file_nodata_values,)

        # add variable metadata
        datasets_obj = add_metadata(
            geo_dataset=datasets_obj, cfg_dataset=dataset_cfg, cfg_extra=extra_cfg)

    # check datasets type
    elif dataset_type == "netcdf_hmc_state":

        # get variable from netcdf file
        datasets_obj = read_file_netcdf_hmc(
            file_path=file_path, file_variable=file_variable,
            file_compression=file_compression,
            file_valid_min=file_valid_min, file_valid_max=file_valid_max)

        # check variable name
        if file_variable == 'VTot':

            # get vtot values
            values_vtot = datasets_obj['values']
            # get cn values
            if 'cn' in group_geo.keys():
                values_cn = group_geo['cn']['values']
            else:
                logger.error(' ===> Curve number values must be defined to scale VTot')
                raise RuntimeError("Curve number must be defined. Exit")

            values_sm = compute_hmc_soil_moisture(values_vtot, values_cn)

            datasets_obj['values'] = values_sm
            datasets_obj['variable'] = 'soil_moisture'

            # debug
            if debug:
                import matplotlib.pyplot as plt
                values_var = datasets_obj["values"]

                plt.figure(figsize=(8, 6))
                plt.imshow(values_vtot, cmap="viridis", origin="upper")
                plt.colorbar(label="vtot")
                plt.tight_layout()

                plt.figure(figsize=(8, 6))
                plt.imshow(values_cn, cmap="viridis", origin="upper")
                plt.colorbar(label="cn")
                plt.tight_layout()

                plt.figure(figsize=(8, 6))
                plt.imshow(values_var, cmap="viridis", origin="upper")
                plt.colorbar(label="sm")
                plt.tight_layout()
                plt.show(block=True)

        # add variable metadata
        datasets_obj = add_metadata(
            geo_dataset=datasets_obj, cfg_dataset=dataset_cfg, cfg_extra=extra_cfg)

    # check datasets type
    elif dataset_type == "tiff" or dataset_type == "tif":

        # get variable from tiff file
        datasets_obj = read_file_tiff_cnr(
            file_path=file_path, file_variable=file_variable,
            file_valid_min=file_valid_min, file_valid_max=file_valid_max)

        # add variable metadata
        datasets_obj = add_metadata(
            geo_dataset=datasets_obj, cfg_dataset=dataset_cfg, cfg_extra=extra_cfg)

        # debug
        if debug:
            import matplotlib.pyplot as plt
            values_var = datasets_obj["values"]

            plt.figure(figsize=(8, 6))
            plt.imshow(values_var, cmap="viridis", origin="upper")
            plt.colorbar(label="sm")
            plt.tight_layout()
            plt.show(block=True)

    else:
        # error dataset type is not expected
        raise NotImplementedError(f"Loader for geo type '{dataset_type}' is not implemented.")

    # info end
    logger.info(f" --------> Read datasets '{dataset_label}' from {file_path} ... DONE")

    return datasets_obj
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to add metadata to datasets
def add_metadata(
        geo_dataset: Dict[str, Any],
        cfg_dataset: Mapping[str, Any],
        cfg_extra: Mapping[str, Any]) -> Dict[str, Any]:

    if cfg_extra is None:
        cfg_extra = {}

    metadata = {
        "scale_factor": cfg_dataset["scale_factor"],
        "offset": cfg_dataset["offset"],
        "valid_min": cfg_dataset["valid_min"],
        "valid_max": cfg_dataset["valid_max"],
        "nodata_values": cfg_dataset["nodata_values"],
    }

    metadata = {**metadata, **cfg_extra}

    geo_dataset["metadata"] = metadata

    return geo_dataset
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to check whether a template contains a time placeholder
def has_time_placeholder(template: str) -> bool:
    """
    Check whether a string contains a ``{time:...}`` or ``{time}``
    placeholder.

    Examples
    --------
    /data/hmc/{time:%Y/%m/%d}
    hmc.state-grid.{time:%Y%m%d%H%M}.nc.gz
    """

    if not isinstance(template, str):
        return False

    formatter = string.Formatter()

    try:
        for _, field_name, _, _ in formatter.parse(template):

            if field_name == "time":
                return True

    except ValueError as exc:
        raise ValueError(
            f"Invalid string template: '{template}'."
        ) from exc

    return False
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resolve a time-dependent template
def resolve_time_template(
        template: str,
        time_step: Any,
        field_name: str,
) -> str:
    """
    Resolve a template containing a ``time`` placeholder.
    """

    if time_step is None:
        raise ValueError(
            f"A time value is required to resolve '{field_name}': "
            f"'{template}'."
        )

    try:
        time_obj = pd.Timestamp(time_step)

    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Unable to parse time '{time_step}' while resolving "
            f"'{field_name}'."
        ) from exc

    if pd.isna(time_obj):
        raise ValueError(
            f"Time is NaT while resolving '{field_name}'."
        )

    try:
        return template.format(
            time=time_obj.to_pydatetime(),
        )

    except (KeyError, IndexError, ValueError) as exc:
        raise ValueError(
            f"Unable to resolve template '{field_name}': "
            f"'{template}' using time '{time_obj}'."
        ) from exc
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to validate one dynamic dataset
def validate_datasets(dataset_cfg: Mapping[str, Any],
                      group_name: str, dataset_key: str,
                      time_step_common: Optional[Any] = None, time_step_datasets: Optional[Any] = None,
                      check_file: bool = True,) -> Dict[str, Any]:

    # create dataset label
    dataset_label = f"{group_name}.{dataset_key}"

    # mandatory expected fields
    mandatory_fields = ["folder","filename",]

    # check format of dataset configuration
    if not isinstance(dataset_cfg, Mapping):
        raise TypeError(f"Dataset configuration '{dataset_label}' must be a dictionary.")

    # get and format configuration
    cfg = deepcopy(dict(dataset_cfg))

    # check mandatory fields
    for field_name in mandatory_fields:
        if field_name not in cfg:
            raise KeyError(f"Field '{field_name}' is missing from dataset configuration '{dataset_label}'.")

    # check dataset name
    dataset_name = cfg.get("name",dataset_key,)
    if not isinstance(dataset_name, str) or not dataset_name.strip():
        raise ValueError(f"Field 'name' in '{dataset_label}' must be a non-empty string.")

    # check dataset type
    dataset_type = str(cfg.get("type", DEFAULT_TYPE)).strip().lower()

    if dataset_type not in SUPPORTED_DATASET_TYPES:
        raise ValueError(
            f"Unsupported dataset type '{dataset_type}' in '{dataset_label}'. Supported types are: "
            f"{sorted(SUPPORTED_DATASET_TYPES)}."
        )

    # check compression
    compression = cfg.get("compression",False,)
    if not isinstance(compression, bool):
        raise TypeError(f"Field 'compression' in '{dataset_label}' must be boolean.")

    # get raw folder and filename templates
    folder_template = os.path.expandvars(os.path.expanduser(str(cfg["folder"])))
    if not folder_template.strip():
        raise ValueError(f"Field 'folder' in '{dataset_label}' cannot be empty.")
    # detect dynamic placeholders
    folder_is_dynamic = has_time_placeholder(folder_template)

    filename_template = str(cfg["filename"])
    if not filename_template.strip():
        raise ValueError(f"Field 'filename' in '{dataset_label}' cannot be empty.")
    # detect dynamic placeholders
    filename_is_dynamic = has_time_placeholder(filename_template)

    # tag dynamic or not
    is_dynamic = (folder_is_dynamic or filename_is_dynamic)

    # resolve paths only when a time step is available
    if is_dynamic and time_step_common is None:
        folder_name, file_name, file_path = None, None, None
    else:

        if folder_is_dynamic:
            folder_name = resolve_time_template(
                template=folder_template, time_step=time_step_datasets, field_name=f"{dataset_label}.folder",)
        else:
            folder_name = folder_template

        if filename_is_dynamic:
            file_name = resolve_time_template(
                template=filename_template, time_step=time_step_datasets, field_name=f"{dataset_label}.filename",)
        else:
            file_name = filename_template

        file_path = os.path.join(folder_name, file_name,)

        # check file system only when requested
        if check_file:
            if not os.path.isdir(folder_name):
                raise FileNotFoundError(f"Dataset folder does not exist for '{dataset_label}': '{folder_name}'.")
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"Dataset file does not exist for '{dataset_label}': '{file_path}'.")

    # check dataset-specific fields
    band, variable, longitude_name, latitude_name, engine = None, None, None, None, None

    # define info by type
    if dataset_type == "tiff" or dataset_type == 'tif':
        band = cfg.get("band",DEFAULT_BAND,)
        if not isinstance(band, int) or band < 1:
            raise ValueError(f"Field 'band' in '{dataset_label}' must be an integer >= 1.")

    elif dataset_type == "ascii_grid":
        band = cfg.get("band",DEFAULT_BAND,)
        if not isinstance(band, int) or band < 1:
            raise ValueError(f"Field 'band' in '{dataset_label}' must be an integer >= 1.")

    elif dataset_type == "netcdf_hmc_state":

        variable = cfg.get("variable")
        if not isinstance(variable, str) or not variable.strip():
            raise ValueError(f"Field 'variable' in '{dataset_label}' must be a non-empty string for NetCDF datasets.")

        longitude_name = cfg.get("longitude","longitude",)
        latitude_name = cfg.get("latitude","latitude",)

        engine = cfg.get("engine",None,)

        if not isinstance(longitude_name, str) or not longitude_name.strip():
            raise ValueError(f"Field 'longitude' in '{dataset_label}' must be a non-empty string.")

        if not isinstance(latitude_name, str) or not latitude_name.strip():
            raise ValueError(f"Field 'latitude' in '{dataset_label}' must be a non-empty string.")

        if engine is not None and not isinstance(engine, str):
            raise TypeError(f"Field 'engine' in '{dataset_label}' must be a string or null.")

    # check CRS override
    crs_raw = cfg.get("crs",cfg.get("crs", DEFAULT_CRS),)
    try:
        crs = CRS.from_user_input(crs_raw)
    except Exception as exc:
        raise ValueError(f"Invalid CRS '{crs_raw}' in '{dataset_label}'.") from exc

    # check scale factor
    scale_factor = cfg.get("scale_factor",DEFAULT_SCALE_FACTOR,)
    _check_numeric_value(value=scale_factor,field_name="scale_factor",dataset_label=dataset_label,)

    # check offset
    offset = cfg.get("offset",DEFAULT_OFFSET,)
    _check_numeric_value(value=offset, field_name="offset", dataset_label=dataset_label,)

    # check valid min and max
    valid_min = cfg.get("valid_min",None,)
    valid_max = cfg.get("valid_max",None,)

    if valid_min is not None:
        _check_numeric_value(value=valid_min, field_name="valid_min", dataset_label=dataset_label,)
    if valid_max is not None:
        _check_numeric_value(value=valid_max,field_name="valid_max",dataset_label=dataset_label,)

    if valid_min is not None and valid_max is not None and float(valid_min) > float(valid_max):
        raise ValueError(f"Field 'valid_min' is greater than 'valid_max' in '{dataset_label}'.")

    # check nodata values
    nodata_values = cfg.get("nodata_values",DEFAULT_NODATA_VALUES,)
    if nodata_values is None: nodata_values = []

    if not isinstance(nodata_values,(list, tuple),):
        raise TypeError(f"Field 'nodata_values' in '{dataset_label}' must be a list or tuple.")

    nodata_values_normalized = []
    for nodata_value in nodata_values:
        _check_numeric_value(value=nodata_value, field_name="nodata_values", dataset_label=dataset_label,)
        nodata_values_normalized.append(float(nodata_value))

    # check scaling
    scaling = cfg.get("scaling",{},)

    # organize configuration
    cfg.update({
        "name": dataset_name.strip(),
        "type": dataset_type,
        "compression": compression,

        # Preserve original templates for every time step.
        "folder_template": folder_template,
        "filename_template": filename_template,

        # These are None during initial validation of dynamic datasets.
        "folder": folder_name,
        "filename": file_name,
        "file_path": file_path,

        "is_dynamic": is_dynamic,
        "folder_is_dynamic": folder_is_dynamic,
        "filename_is_dynamic": filename_is_dynamic,

        "time": (None if time_step_common is None else pd.Timestamp(time_step_datasets)),

        "band": band,
        "variable": (variable.strip() if variable is not None else None),
        "longitude": (longitude_name.strip() if longitude_name is not None else None),
        "latitude": (latitude_name.strip() if latitude_name is not None else None),
        "engine": (engine.strip() if isinstance(engine, str) else None),

        "crs": crs,

        "scale_factor": float(scale_factor),
        "offset": float(offset),

        "valid_min": (None if valid_min is None else float(valid_min)),
        "valid_max": (None if valid_max is None else float(valid_max)),

        "nodata_values": nodata_values_normalized,

        "scaling": scaling
    })

    return cfg
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to check a numeric configuration field
def _check_numeric_value(
        value: Any, field_name: str,dataset_label: str, ) -> None:

    if not isinstance(
            value,
            (
                int,
                float,
                np.integer,
                np.floating,
            ),
    ):
        raise TypeError(f"Field '{field_name}' in '{dataset_label}' "f"must be numeric.")
# ----------------------------------------------------------------------------------------------------------------------
