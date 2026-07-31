"""
Library Features:

Name:          lib_data_io_csv
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import json
import logging
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from netCDF4 import Dataset, date2num, num2date

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to check datasets
def check_datasets_nc(
        file_name,
        current_point,
        first_point_expected=None,
        last_point_expected=None,
        point_variable="point_name"):
    """
    Check whether point names exist in the NetCDF point_name variable.
    """

    if not os.path.exists(file_name):
        return (
            False,
            None if first_point_expected is None else False,
            None if last_point_expected is None else False,
        )

    with Dataset(file_name, mode="r") as file_handle:

        if point_variable not in file_handle.variables:
            raise RuntimeError(
                f'Point variable "{point_variable}" not found in "{file_name}"'
            )

        point_values = file_handle.variables[point_variable][:]

        if np.ma.isMaskedArray(point_values):
            point_values = point_values.filled("")

        point_names = []

        for value in np.asarray(point_values).reshape(-1):

            if isinstance(value, bytes):
                value = value.decode("utf-8")

            point_names.append(str(value).strip())

        current_exists = str(current_point) in point_names

        first_exists = (
            None
            if first_point_expected is None
            else str(first_point_expected) in point_names
        )

        last_exists = (
            None
            if last_point_expected is None
            else str(last_point_expected) in point_names
        )

    return current_exists, first_exists, last_exists
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to select datasets by point
def select_datasets_by_point(
        datasets_obj,
        point_name,
        file_fields=None,
        point_name_field="point_name",
        time_name="time",
        restore_source_names=True,
        **kwargs):
    """
    Select one point from the object returned by read_datasets_nc.

    Parameters
    ----------
    datasets_obj : dict
        Object returned by ``read_datasets_nc``.

    point_name : str
        Name of the point to select.

    file_fields : dict, optional
        Mapping between DataFrame source names and NetCDF names.

        Example:

            {
                "time": "time",
                "values_k1": "rain",
                "values_k2": "air_temperature",
                "values_k3": "soil_moisture"
            }

    point_name_field : str
        Name of the point variable containing point names.

    restore_source_names : bool
        If True, convert NetCDF variable names back to DataFrame names:

            rain -> values_k1

        If False, keep NetCDF names:

            rain -> rain

    Returns
    -------
    point_dframe : pandas.DataFrame
        DataFrame containing all time-series variables for the selected point.

        Point information is stored in:

            point_dframe.attrs
    """

    # ---------------------------------------------------------------------------------
    # Check input object
    if not isinstance(datasets_obj, dict):
        raise TypeError(
            '"datasets_obj" must be the object returned by '
            "read_datasets_nc"
        )

    required_keys = {
        "time",
        "point",
        "data"
    }

    missing_keys = required_keys.difference(
        datasets_obj.keys()
    )

    if missing_keys:
        raise KeyError(
            "Invalid datasets object. Missing keys: "
            f"{sorted(missing_keys)}"
        )

    point_obj = datasets_obj["point"]
    data_obj = datasets_obj["data"]
    time_index = pd.DatetimeIndex(
        datasets_obj["time"]
    )

    # ---------------------------------------------------------------------------------
    # Check point name
    if point_name is None:
        raise ValueError(
            '"point_name" is mandatory'
        )

    point_name = str(point_name).strip()

    if not point_name:
        raise ValueError(
            '"point_name" cannot be empty'
        )

    if point_name_field not in point_obj:
        raise KeyError(
            f'Point variable "{point_name_field}" is not available. '
            f"Available point variables: {list(point_obj.keys())}"
        )

    # ---------------------------------------------------------------------------------
    # Find point row
    saved_point_names = np.asarray(
        point_obj[point_name_field],
        dtype=object
    )

    point_matches = np.where(
        np.asarray(
            [
                str(saved_name).strip() == point_name
                for saved_name in saved_point_names
            ],
            dtype=bool
        )
    )[0]

    if point_matches.size == 0:
        raise KeyError(
            f'Point "{point_name}" was not found. '
            f"Available points: "
            f"{[str(value) for value in saved_point_names]}"
        )

    if point_matches.size > 1:
        log_stream.warning(
            ' ===> Point name "%s" occurs %d times. '
            "The first occurrence will be selected",
            point_name,
            point_matches.size
        )

    point_index = int(
        point_matches[0]
    )

    # ---------------------------------------------------------------------------------
    # Define output variable names
    if file_fields is None:

        variable_mapping = {
            variable_name: variable_name
            for variable_name in data_obj
        }

    else:

        if not isinstance(file_fields, dict):
            raise TypeError(
                '"file_fields" must be a dictionary or None'
            )

        variable_mapping = {}

        for source_name, destination_name in file_fields.items():

            source_name = str(source_name).strip()
            destination_name = str(destination_name).strip()

            if (
                    source_name == time_name
                    or destination_name == time_name):
                continue

            if destination_name not in data_obj:
                continue

            output_name = (
                source_name
                if restore_source_names
                else destination_name
            )

            variable_mapping[output_name] = destination_name

        if not variable_mapping:
            raise RuntimeError(
                "No configured variables are available in the "
                f"NetCDF object. Available variables: "
                f"{list(data_obj.keys())}"
            )

    # ---------------------------------------------------------------------------------
    # Extract point data
    point_data = {}

    for output_name, variable_name in variable_mapping.items():

        variable_values = np.asarray(
            data_obj[variable_name]
        )

        if variable_values.ndim != 2:
            raise RuntimeError(
                f'Data variable "{variable_name}" must have shape '
                f"(point, time). Found: {variable_values.shape}"
            )

        if point_index >= variable_values.shape[0]:
            raise IndexError(
                f'Point index {point_index} is outside variable '
                f'"{variable_name}" with shape {variable_values.shape}'
            )

        point_data[output_name] = variable_values[
            point_index,
            :
        ].copy()

    # ---------------------------------------------------------------------------------
    # Create DataFrame
    point_dframe = pd.DataFrame(
        data=point_data,
        index=time_index
    )

    point_dframe.index.name = time_name

    # ---------------------------------------------------------------------------------
    # Add all point information to DataFrame attrs
    point_attributes = {}

    for variable_name, variable_values in point_obj.items():

        variable_values = np.asarray(
            variable_values
        )

        if point_index >= variable_values.size:
            continue

        variable_value = variable_values[
            point_index
        ]

        if isinstance(variable_value, np.generic):
            variable_value = variable_value.item()

        point_attributes[variable_name] = variable_value

    point_attributes["point_index"] = point_index
    point_attributes["file_name"] = datasets_obj.get(
        "file_name"
    )

    point_dframe.attrs = point_attributes

    log_stream.info(
        ' ::: Select point "%s" selected at row %d with '
        "%d time step(s) and %d variable(s)",
        point_name,
        point_index,
        len(point_dframe.index),
        len(point_dframe.columns)
    )

    return point_dframe
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read datasets
def read_datasets_nc(
        file_name,
        time_name="time",
        point_dim_name="point",
        point_name_field="point_name",
        file_fill_value=-9999.0,
        **kwargs):
    """
    Read the complete NetCDF point collection.

    Expected NetCDF structure
    -------------------------
    dimensions:
        point = N
        time = T

    variables:
        time(time)

        point_id(point)
        point_name(point)
        longitude(point)
        latitude(point)

        rain(point, time)
        air_temperature(point, time)
        soil_moisture(point, time)

    Returns
    -------
    datasets_obj : dict
        {
            "time": pandas.DatetimeIndex,
            "point": {
                "point_id": numpy.ndarray,
                "point_name": numpy.ndarray,
                "longitude": numpy.ndarray,
                "latitude": numpy.ndarray,
                ...
            },
            "data": {
                "rain": numpy.ndarray,             # shape (point, time)
                "air_temperature": numpy.ndarray,  # shape (point, time)
                "soil_moisture": numpy.ndarray     # shape (point, time)
            },
            "attributes": {...},
            "file_name": "..."
        }
    """

    # ---------------------------------------------------------------------------------
    # Check source file
    file_path = Path(file_name)

    if not file_path.exists():
        raise FileNotFoundError(
            f'NetCDF file "{file_path}" does not exist'
        )

    # ---------------------------------------------------------------------------------
    # Decode a scalar string value
    def _decode_string(value):

        if isinstance(value, bytes):
            return value.decode(
                "utf-8",
                errors="replace"
            ).strip()

        if isinstance(value, np.bytes_):
            return value.tobytes().decode(
                "utf-8",
                errors="replace"
            ).strip()

        if value is None:
            return ""

        return str(value).strip()

    # ---------------------------------------------------------------------------------
    # Read a NetCDF variable
    def _read_variable(variable_obj):

        values = variable_obj[:]

        # Convert masked arrays
        if np.ma.isMaskedArray(values):

            if values.dtype.kind in {"i", "u", "f"}:
                values = values.filled(np.nan)
            else:
                values = values.filled("")

        values = np.asarray(values)

        # Decode string variables
        if values.dtype.kind in {"O", "S", "U"}:

            output_values = np.empty(
                values.shape,
                dtype=object
            )

            for value_index in np.ndindex(values.shape):
                output_values[value_index] = _decode_string(
                    values[value_index]
                )

            return output_values

        # Convert numeric variables
        if values.dtype.kind in {"i", "u", "f"}:

            values = values.astype(
                np.float64,
                copy=True
            )

            invalid_mask = ~np.isfinite(values)

            variable_fill_value = getattr(
                variable_obj,
                "_FillValue",
                None
            )

            if variable_fill_value is not None:
                try:
                    invalid_mask |= np.isclose(
                        values,
                        float(variable_fill_value),
                        rtol=0.0,
                        atol=0.0
                    )
                except (TypeError, ValueError):
                    pass

            if file_fill_value is not None:
                try:
                    invalid_mask |= np.isclose(
                        values,
                        float(file_fill_value),
                        rtol=0.0,
                        atol=0.0
                    )
                except (TypeError, ValueError):
                    pass

            values[invalid_mask] = np.nan

        return values

    # ---------------------------------------------------------------------------------
    # Read file
    try:

        with Dataset(
                str(file_path),
                mode="r") as file_handle:

            # -------------------------------------------------------------------------
            # Check mandatory dimensions and variables
            if point_dim_name not in file_handle.dimensions:
                raise RuntimeError(
                    f'Point dimension "{point_dim_name}" is missing '
                    f'from "{file_path}"'
                )

            if time_name not in file_handle.variables:
                raise RuntimeError(
                    f'Time variable "{time_name}" is missing '
                    f'from "{file_path}"'
                )

            if point_name_field not in file_handle.variables:
                raise RuntimeError(
                    f'Point-name variable "{point_name_field}" is missing '
                    f'from "{file_path}"'
                )

            number_of_points = len(
                file_handle.dimensions[point_dim_name]
            )

            # -------------------------------------------------------------------------
            # Read time
            time_variable = file_handle.variables[
                time_name
            ]

            time_units = getattr(
                time_variable,
                "units",
                None
            )

            time_calendar = getattr(
                time_variable,
                "calendar",
                "standard"
            )

            if time_units is None:
                raise RuntimeError(
                    f'Time variable "{time_name}" does not have '
                    'the mandatory "units" attribute'
                )

            time_values = num2date(
                time_variable[:],
                units=time_units,
                calendar=time_calendar,
                only_use_cftime_datetimes=False,
                only_use_python_datetimes=True
            )

            time_index = pd.DatetimeIndex(
                pd.to_datetime(
                    np.asarray(time_values)
                )
            )

            number_of_time_steps = len(
                time_index
            )

            # -------------------------------------------------------------------------
            # Read variables
            point_obj = {}
            data_obj = {}

            for variable_name, variable_obj in file_handle.variables.items():

                if variable_name == time_name:
                    continue

                variable_dimensions = variable_obj.dimensions

                # Point information: variable(point)
                if variable_dimensions == (point_dim_name,):

                    variable_values = _read_variable(
                        variable_obj
                    )

                    if variable_values.shape != (number_of_points,):
                        raise RuntimeError(
                            f'Point variable "{variable_name}" has shape '
                            f"{variable_values.shape}; expected "
                            f"({number_of_points},)"
                        )

                    point_obj[variable_name] = variable_values

                # Time-series data: variable(point, time)
                elif variable_dimensions == (
                        point_dim_name,
                        time_name):

                    variable_values = _read_variable(
                        variable_obj
                    )

                    expected_shape = (
                        number_of_points,
                        number_of_time_steps
                    )

                    if variable_values.shape != expected_shape:
                        raise RuntimeError(
                            f'Data variable "{variable_name}" has shape '
                            f"{variable_values.shape}; expected "
                            f"{expected_shape}"
                        )

                    data_obj[variable_name] = variable_values

            # -------------------------------------------------------------------------
            # Read global attributes
            file_attributes = {
                attribute_name: file_handle.getncattr(
                    attribute_name
                )
                for attribute_name in file_handle.ncattrs()
            }

    except OSError as exc:
        raise RuntimeError(
            f'Unable to read NetCDF file "{file_path}"'
        ) from exc

    # ---------------------------------------------------------------------------------
    # Create output object
    datasets_obj = {
        "file_name": str(file_path),
        "time": time_index,
        "point": point_obj,
        "data": data_obj,
        "attributes": file_attributes
    }

    log_stream.info(
        ' ::: Read file "%s" read: %d point(s), '
        "%d time step(s), %d data variable(s)",
        file_path,
        number_of_points,
        number_of_time_steps,
        len(data_obj)
    )

    return datasets_obj


# -------------------------------------------------------------------------------------
def _get_netcdf_field_pairs(
        file_dframe, file_fields,
        time_name="time"):
    """
    Organize source and destination variable names.

    Example
    -------
    file_fields = {
        "time": "time",
        "values_k1": "rain",
        "values_k2": "air_temperature",
        "values_k3": "soil_moisture"
    }
    """

    if file_fields is None:
        file_fields = {
            variable_name: variable_name
            for variable_name in file_dframe.columns
        }
    elif not isinstance(file_fields, Mapping):
        raise TypeError(
            '"file_fields" must be a dictionary, mapping or None'
        )

    field_pairs = []

    for source_name, destination_name in file_fields.items():

        source_name = str(source_name).strip()
        destination_name = str(destination_name).strip()

        if not source_name:
            raise ValueError(
                '"file_fields" contains an empty source field name'
            )

        if not destination_name:
            raise ValueError(
                f'No destination field defined for "{source_name}"'
            )

        if (
                source_name == time_name
                or destination_name == time_name):
            continue

        field_pairs.append(
            (
                source_name,
                destination_name
            )
        )

    if not field_pairs:
        raise RuntimeError(
            'No data variables are defined in "file_fields"'
        )

    source_names = [
        source_name
        for source_name, destination_name in field_pairs
    ]

    destination_names = [
        destination_name
        for source_name, destination_name in field_pairs
    ]

    if len(source_names) != len(set(source_names)):
        raise RuntimeError(
            f"Duplicated source variable names: {source_names}"
        )

    if len(destination_names) != len(set(destination_names)):
        raise RuntimeError(
            f"Duplicated destination variable names: "
            f"{destination_names}"
        )

    return field_pairs


# -------------------------------------------------------------------------------------
def _normalize_required_point_metadata(file_dframe):
    """
    Validate and normalize mandatory point metadata.

    Mandatory DataFrame attributes
    ------------------------------
    point_id
    point_name
    longitude
    latitude
    """

    mandatory_fields = [
        "point_id",
        "point_name",
        "longitude",
        "latitude"
    ]

    missing_fields = [
        field_name
        for field_name in mandatory_fields
        if field_name not in file_dframe.attrs
    ]

    if missing_fields:
        raise KeyError(
            "Missing mandatory DataFrame attributes: "
            f"{missing_fields}. "
            f"Available attributes: {list(file_dframe.attrs.keys())}"
        )

    point_metadata = dict(
        file_dframe.attrs
    )

    # ---------------------------------------------------------------------------------
    # Point identifier
    point_id = point_metadata.get(
        "point_id"
    )

    if point_id is None:
        raise ValueError(
            'Mandatory attribute "point_id" cannot be None'
        )

    point_id = str(
        point_id
    ).strip()

    if not point_id:
        raise ValueError(
            'Mandatory attribute "point_id" cannot be empty'
        )

    point_metadata["point_id"] = point_id

    # ---------------------------------------------------------------------------------
    # Point name
    point_name = point_metadata.get(
        "point_name"
    )

    if point_name is None:
        raise ValueError(
            'Mandatory attribute "point_name" cannot be None'
        )

    point_name = str(
        point_name
    ).strip()

    if not point_name:
        raise ValueError(
            'Mandatory attribute "point_name" cannot be empty'
        )

    point_metadata["point_name"] = point_name

    # ---------------------------------------------------------------------------------
    # Longitude and latitude
    for coordinate_name in [
            "longitude",
            "latitude"]:

        coordinate_value = point_metadata.get(
            coordinate_name
        )

        try:
            coordinate_value = float(
                coordinate_value
            )

        except (TypeError, ValueError) as exc:
            raise TypeError(
                f'Mandatory attribute "{coordinate_name}" must be numeric. '
                f'Received: {coordinate_value!r}'
            ) from exc

        if not np.isfinite(
                coordinate_value):

            raise ValueError(
                f'Mandatory attribute "{coordinate_name}" must be finite. '
                f'Received: {coordinate_value!r}'
            )

        point_metadata[
            coordinate_name
        ] = coordinate_value

    if not -180.0 <= point_metadata["longitude"] <= 180.0:
        raise ValueError(
            '"longitude" must be between -180 and 180 degrees. '
            f'Received: {point_metadata["longitude"]}'
        )

    if not -90.0 <= point_metadata["latitude"] <= 90.0:
        raise ValueError(
            '"latitude" must be between -90 and 90 degrees. '
            f'Received: {point_metadata["latitude"]}'
        )

    return point_metadata


# -------------------------------------------------------------------------------------
def _prepare_datasets_dataframe(
        file_dframe,
        field_pairs,
        time_name="time"):
    """
    Validate and organize the input DataFrame.
    """

    if not isinstance(
            file_dframe,
            pd.DataFrame):

        raise TypeError(
            '"file_dframe" must be a pandas.DataFrame'
        )

    if file_dframe.empty:
        raise ValueError(
            '"file_dframe" cannot be empty'
        )

    organized_dframe = file_dframe.copy(
        deep=True
    )

    file_attrs = dict(
        file_dframe.attrs
    )

    # ---------------------------------------------------------------------------------
    # Organize time
    if time_name in organized_dframe.columns:

        time_values = organized_dframe.pop(
            time_name
        )

        time_index = pd.DatetimeIndex(
            pd.to_datetime(
                time_values,
                errors="coerce"
            )
        )

    else:

        time_index = pd.DatetimeIndex(
            pd.to_datetime(
                organized_dframe.index,
                errors="coerce"
            )
        )

    valid_time_mask = ~time_index.isna()

    if not valid_time_mask.all():

        n_invalid = int(
            (~valid_time_mask).sum()
        )

        log_stream.warning(
            " ===> Removing %d row(s) with invalid timestamps",
            n_invalid
        )

        valid_positions = np.where(
            valid_time_mask
        )[0]

        organized_dframe = organized_dframe.iloc[
            valid_positions
        ].copy()

        time_index = time_index[
            valid_time_mask
        ]

    if time_index.empty:
        raise RuntimeError(
            "No valid timestamps are available"
        )

    organized_dframe.index = time_index
    organized_dframe.index.name = time_name

    # ---------------------------------------------------------------------------------
    # Remove duplicated timestamps
    duplicated_mask = organized_dframe.index.duplicated(
        keep="last"
    )

    if duplicated_mask.any():

        n_duplicated = int(
            duplicated_mask.sum()
        )

        log_stream.warning(
            " ===> Removing %d duplicated timestamp(s)",
            n_duplicated
        )

        organized_dframe = organized_dframe.loc[
            ~duplicated_mask
        ].copy()

    organized_dframe = organized_dframe.sort_index()

    # ---------------------------------------------------------------------------------
    # Identify source or final variable names
    source_names = [
        source_name
        for source_name, destination_name in field_pairs
    ]

    destination_names = [
        destination_name
        for source_name, destination_name in field_pairs
    ]

    available_columns = [
        str(column_name)
        for column_name in organized_dframe.columns
    ]

    source_names_available = all(
        variable_name in available_columns
        for variable_name in source_names
    )

    destination_names_available = all(
        variable_name in available_columns
        for variable_name in destination_names
    )

    if source_names_available:
        selected_names = source_names

    elif destination_names_available:
        selected_names = destination_names

    else:

        missing_source_names = [
            variable_name
            for variable_name in source_names
            if variable_name not in available_columns
        ]

        missing_destination_names = [
            variable_name
            for variable_name in destination_names
            if variable_name not in available_columns
        ]

        raise KeyError(
            "Unable to identify input DataFrame variables. "
            f"Available columns: {available_columns}. "
            f"Missing source variables: {missing_source_names}. "
            f"Missing destination variables: "
            f"{missing_destination_names}"
        )

    organized_dframe = organized_dframe[
        selected_names
    ].copy()

    organized_dframe.columns = destination_names

    # ---------------------------------------------------------------------------------
    # Convert values to numeric
    for variable_name in destination_names:

        organized_dframe[
            variable_name
        ] = pd.to_numeric(
            organized_dframe[variable_name],
            errors="coerce"
        )

    organized_dframe.attrs = file_attrs

    return organized_dframe


# -------------------------------------------------------------------------------------
def _write_global_attribute(
        file_handle,
        attribute_name,
        attribute_value):
    """
    Write a generic DataFrame attribute as a NetCDF global attribute.
    """

    if attribute_value is None:
        return

    if isinstance(
            attribute_value,
            np.generic):

        attribute_value = attribute_value.item()

    if isinstance(
            attribute_value,
            pd.Timestamp):

        attribute_value = attribute_value.isoformat()

    if isinstance(
            attribute_value,
            (str, int, float, bool)):

        file_handle.setncattr(
            str(attribute_name),
            attribute_value
        )

        return

    if isinstance(
            attribute_value,
            np.ndarray):

        attribute_value = attribute_value.tolist()

    try:

        attribute_value = json.dumps(
            attribute_value,
            default=str
        )

    except (TypeError, ValueError):

        attribute_value = str(
            attribute_value
        )

    file_handle.setncattr(
        str(attribute_name),
        attribute_value
    )
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to write datasets
def write_datasets_nc(
        file_name,
        file_dframe,
        file_fields,
        point_id,
        point_name,
        longitude,
        latitude,
        point_attrs=None,
        time_name="time",
        point_dim_name="point",
        time_units="hours since 1970-01-01 00:00:00",
        time_calendar="standard",
        file_dtype="f4",
        file_fill_value=-9999.0,
        file_compression=True,
        file_compression_level=4,
        file_shuffle=True,
        file_rebuild_legacy=True,
        **kwargs):
    """
    Write or append one point time series to a NetCDF collection.

    The first call creates the NetCDF file. Subsequent calls append points
    along the unlimited ``point`` dimension.

    NetCDF structure
    ----------------
    dimensions:

        point = UNLIMITED
        time = N

    variables:

        time(time)

        point_id(point)
        point_name(point)
        longitude(point)
        latitude(point)

        optional_attribute(point)

        rain(point, time)
        air_temperature(point, time)
        soil_moisture(point, time)

    Parameters
    ----------
    file_name : str
        Destination NetCDF filename.

    file_dframe : pandas.DataFrame
        DataFrame containing one point time series.

    file_fields : dict
        Mapping between DataFrame fields and NetCDF variable names.

        Example:

            {
                "time": "time",
                "values_k1": "rain",
                "values_k2": "air_temperature",
                "values_k3": "soil_moisture"
            }

    point_id : str or int
        Point identifier.

    point_name : str
        Point name.

    longitude : float
        Point longitude.

    latitude : float
        Point latitude.

    point_attrs : dict, optional
        Additional point metadata. Each attribute is stored as a variable
        using the unlimited point dimension.

    file_rebuild_legacy : bool
        If True, remove an existing file created with a fixed point dimension
        and recreate it using the unlimited point dimension.
    """

    # ---------------------------------------------------------------------------------
    # Check DataFrame
    if not isinstance(file_dframe, pd.DataFrame):
        raise TypeError(
            '"file_dframe" must be a pandas.DataFrame'
        )

    if file_dframe.empty:
        raise ValueError(
            '"file_dframe" cannot be empty'
        )

    # ---------------------------------------------------------------------------------
    # Check mandatory point metadata
    if point_id is None:
        raise ValueError(
            '"point_id" is mandatory'
        )

    point_id = str(point_id).strip()

    if not point_id:
        raise ValueError(
            '"point_id" cannot be empty'
        )

    if point_name is None:
        raise ValueError(
            '"point_name" is mandatory'
        )

    point_name = str(point_name).strip()

    if not point_name:
        raise ValueError(
            '"point_name" cannot be empty'
        )

    try:
        longitude = float(longitude)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            f'"longitude" must be numeric. Received: {longitude!r}'
        ) from exc

    try:
        latitude = float(latitude)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            f'"latitude" must be numeric. Received: {latitude!r}'
        ) from exc

    if not np.isfinite(longitude):
        raise ValueError(
            f'"longitude" must be finite. Received: {longitude}'
        )

    if not np.isfinite(latitude):
        raise ValueError(
            f'"latitude" must be finite. Received: {latitude}'
        )

    if not -180.0 <= longitude <= 180.0:
        raise ValueError(
            f'"longitude" must be between -180 and 180. '
            f"Received: {longitude}"
        )

    if not -90.0 <= latitude <= 90.0:
        raise ValueError(
            f'"latitude" must be between -90 and 90. '
            f"Received: {latitude}"
        )

    # ---------------------------------------------------------------------------------
    # Organize additional point attributes
    if point_attrs is None:
        point_attrs = {}

    if not isinstance(point_attrs, dict):
        raise TypeError(
            '"point_attrs" must be a dictionary or None'
        )

    reserved_names = {
        point_dim_name,
        time_name,
        "point_id",
        "point_name",
        "longitude",
        "latitude"
    }

    invalid_attribute_names = [
        attr_name
        for attr_name in point_attrs
        if attr_name in reserved_names
    ]

    if invalid_attribute_names:
        raise ValueError(
            "Reserved variable names used in point_attrs: "
            f"{invalid_attribute_names}"
        )

    point_metadata = {
        "point_id": point_id,
        "point_name": point_name,
        "longitude": longitude,
        "latitude": latitude,
        **point_attrs
    }

    # ---------------------------------------------------------------------------------
    # Organize field mapping
    field_pairs = _get_netcdf_field_pairs(
        file_dframe=file_dframe,
        file_fields=file_fields,
        time_name=time_name
    )

    destination_names = [
        destination_name
        for _, destination_name in field_pairs
    ]

    # Remove variables that conflict with point metadata
    conflicting_names = [
        variable_name
        for variable_name in destination_names
        if variable_name in reserved_names
           or variable_name in point_attrs
    ]

    if conflicting_names:
        log_stream.warning(
            ' ===> The following variables are reserved as point metadata '
            'and will not be saved as time-series variables: %s',
            conflicting_names
        )

        destination_names = [
            variable_name
            for variable_name in destination_names
            if variable_name not in conflicting_names
        ]
    # ---------------------------------------------------------------------------------
    # Organize DataFrame
    organized_dframe = _prepare_datasets_dataframe(
        file_dframe=file_dframe,
        field_pairs=field_pairs,
        time_name=time_name
    )

    time_index = pd.DatetimeIndex(
        organized_dframe.index
    )

    if time_index.empty:
        raise RuntimeError(
            "No valid timestamps are available"
        )

    number_of_time_steps = len(time_index)

    # ---------------------------------------------------------------------------------
    # Check NetCDF options
    file_compression_level = int(
        file_compression_level
    )

    if not 0 <= file_compression_level <= 9:
        raise ValueError(
            '"file_compression_level" must be between 0 and 9'
        )

    file_fill_value = float(
        file_fill_value
    )

    output_dtype = np.dtype(
        file_dtype
    )

    if not np.issubdtype(output_dtype, np.number):
        raise TypeError(
            f'Unsupported NetCDF data type: "{file_dtype}"'
        )

    # ---------------------------------------------------------------------------------
    # Encode time
    encoded_time = np.asarray(
        date2num(
            time_index.to_pydatetime(),
            units=time_units,
            calendar=time_calendar
        ),
        dtype=np.float64
    )

    # ---------------------------------------------------------------------------------
    # Prepare data arrays
    data_collection = {}

    for variable_name in destination_names:

        data_values = organized_dframe[
            variable_name
        ].to_numpy(
            dtype=output_dtype,
            copy=True
        )

        data_values[
            ~np.isfinite(data_values)
        ] = file_fill_value

        data_collection[variable_name] = data_values

    # ---------------------------------------------------------------------------------
    # Prepare destination path
    file_path = Path(file_name)

    file_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    # ---------------------------------------------------------------------------------
    # Check old fixed-dimension NetCDF layout
    if file_path.exists():

        try:

            with Dataset(
                    str(file_path),
                    mode="r") as file_handle:

                if point_dim_name not in file_handle.dimensions:
                    raise RuntimeError(
                        f'Point dimension "{point_dim_name}" is missing '
                        f'from "{file_path}"'
                    )

                point_dimension = file_handle.dimensions[
                    point_dim_name
                ]

                if not point_dimension.isunlimited():

                    if not file_rebuild_legacy:
                        raise RuntimeError(
                            f'Point dimension "{point_dim_name}" is fixed '
                            f'in "{file_path}" with size '
                            f"{len(point_dimension)}. The file cannot be "
                            "extended"
                        )

                    rebuild_legacy_file = True

                else:
                    rebuild_legacy_file = False

        except OSError as exc:
            raise RuntimeError(
                f'Unable to inspect NetCDF file "{file_path}"'
            ) from exc

        if rebuild_legacy_file:

            log_stream.warning(
                ' ===> NetCDF file "%s" uses a fixed "%s" dimension. '
                "The legacy file will be removed and recreated",
                file_path,
                point_dim_name
            )

            file_path.unlink()

    file_exists = file_path.exists()

    # ---------------------------------------------------------------------------------
    # Define metadata type
    def _get_metadata_type(metadata_value):

        if isinstance(metadata_value, (bool, np.bool_)):
            return "integer"

        if isinstance(metadata_value, (int, np.integer)):
            return "integer"

        if isinstance(metadata_value, (float, np.floating)):
            return "float"

        return "string"

    # ---------------------------------------------------------------------------------
    # Create metadata variable
    def _create_metadata_variable(
            file_handle,
            metadata_name,
            metadata_value):

        metadata_type = _get_metadata_type(
            metadata_value
        )

        if metadata_type == "integer":

            metadata_variable = file_handle.createVariable(
                metadata_name,
                "i8",
                (point_dim_name,),
                fill_value=-9999
            )

        elif metadata_type == "float":

            metadata_variable = file_handle.createVariable(
                metadata_name,
                "f8",
                (point_dim_name,),
                fill_value=np.nan
            )

        else:

            metadata_variable = file_handle.createVariable(
                metadata_name,
                str,
                (point_dim_name,)
            )

        return metadata_variable

    # ---------------------------------------------------------------------------------
    # Write metadata value
    def _write_metadata_value(
            metadata_variable,
            point_index,
            metadata_value):

        variable_dtype = metadata_variable.dtype

        if variable_dtype == str or variable_dtype == object:

            if metadata_value is None:
                metadata_variable[point_index] = ""
            else:
                metadata_variable[point_index] = str(metadata_value)

            return

        if np.issubdtype(variable_dtype, np.integer):

            if metadata_value is None:
                metadata_variable[point_index] = -9999
            else:
                try:
                    metadata_variable[point_index] = int(metadata_value)
                except (TypeError, ValueError, OverflowError):
                    metadata_variable[point_index] = -9999

            return

        if metadata_value is None:

            metadata_variable[point_index] = np.nan

        else:

            try:
                metadata_variable[point_index] = float(metadata_value)
            except (TypeError, ValueError, OverflowError):
                metadata_variable[point_index] = np.nan

    # ---------------------------------------------------------------------------------
    # Create new NetCDF file
    if not file_exists:

        with Dataset(
                str(file_path),
                mode="w",
                format="NETCDF4") as file_handle:

            # -------------------------------------------------------------------------
            # Dimensions
            file_handle.createDimension(
                point_dim_name,
                None
            )

            file_handle.createDimension(
                time_name,
                number_of_time_steps
            )

            # -------------------------------------------------------------------------
            # Time variable
            time_variable = file_handle.createVariable(
                time_name,
                "f8",
                (time_name,)
            )

            time_variable[:] = encoded_time

            time_variable.units = time_units
            time_variable.calendar = time_calendar
            time_variable.standard_name = "time"
            time_variable.long_name = "reference time"
            time_variable.axis = "T"

            # -------------------------------------------------------------------------
            # Point metadata variables
            for metadata_name, metadata_value in point_metadata.items():

                metadata_variable = _create_metadata_variable(
                    file_handle=file_handle,
                    metadata_name=metadata_name,
                    metadata_value=metadata_value
                )

                _write_metadata_value(
                    metadata_variable=metadata_variable,
                    point_index=0,
                    metadata_value=metadata_value
                )

            point_id_variable = file_handle.variables[
                "point_id"
            ]

            point_id_variable.long_name = "point identifier"
            point_id_variable.cf_role = "timeseries_id"

            point_name_variable = file_handle.variables[
                "point_name"
            ]

            point_name_variable.long_name = "point name"

            longitude_variable = file_handle.variables[
                "longitude"
            ]

            longitude_variable.standard_name = "longitude"
            longitude_variable.long_name = "point longitude"
            longitude_variable.units = "degrees_east"

            latitude_variable = file_handle.variables[
                "latitude"
            ]

            latitude_variable.standard_name = "latitude"
            latitude_variable.long_name = "point latitude"
            latitude_variable.units = "degrees_north"

            # -------------------------------------------------------------------------
            # Time-series variables
            for variable_name, data_values in data_collection.items():

                data_variable = file_handle.createVariable(
                    variable_name,
                    file_dtype,
                    (
                        point_dim_name,
                        time_name
                    ),
                    fill_value=file_fill_value,
                    zlib=bool(file_compression),
                    complevel=file_compression_level,
                    shuffle=bool(file_shuffle)
                )

                data_variable[0, :] = data_values

                data_variable.long_name = variable_name
                data_variable.coordinates = (
                    f"{time_name} longitude latitude"
                )

            # -------------------------------------------------------------------------
            # Global attributes
            file_handle.title = "Point time-series collection"
            file_handle.Conventions = "CF-1.8"
            file_handle.featureType = "timeSeries"

            file_handle.number_of_points = 1
            file_handle.number_of_time_steps = number_of_time_steps
            file_handle.number_of_variables = len(destination_names)

            file_handle.time_coverage_start = (
                time_index[0].strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
            )

            file_handle.time_coverage_end = (
                time_index[-1].strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
            )

        point_index = 0

    # ---------------------------------------------------------------------------------
    # Append point to existing NetCDF file
    else:

        with Dataset(
                str(file_path),
                mode="a") as file_handle:

            # -------------------------------------------------------------------------
            # Check point dimension
            if point_dim_name not in file_handle.dimensions:
                raise RuntimeError(
                    f'Point dimension "{point_dim_name}" is missing '
                    f'from "{file_path}"'
                )

            point_dimension = file_handle.dimensions[
                point_dim_name
            ]

            if not point_dimension.isunlimited():
                raise RuntimeError(
                    f'Point dimension "{point_dim_name}" is not unlimited '
                    f'in "{file_path}"'
                )

            # -------------------------------------------------------------------------
            # Check time variable
            if time_name not in file_handle.variables:
                raise RuntimeError(
                    f'Time variable "{time_name}" is missing '
                    f'from "{file_path}"'
                )

            saved_time_variable = file_handle.variables[
                time_name
            ]

            saved_time_values = np.asarray(
                saved_time_variable[:],
                dtype=np.float64
            )

            saved_time_units = getattr(
                saved_time_variable,
                "units",
                time_units
            )

            saved_time_calendar = getattr(
                saved_time_variable,
                "calendar",
                time_calendar
            )

            # -------------------------------------------------------------------------
            # Decode the time axis stored in the NetCDF file
            saved_time_decoded = num2date(
                saved_time_values,
                units=saved_time_units,
                calendar=saved_time_calendar,
                only_use_cftime_datetimes=False,
                only_use_python_datetimes=True
            )

            saved_time_index = pd.DatetimeIndex(
                pd.to_datetime(
                    np.asarray(saved_time_decoded)
                )
            )

            incoming_time_index = pd.DatetimeIndex(
                organized_dframe.index
            )

            # Remove timezone information, if present
            if saved_time_index.tz is not None:
                saved_time_index = saved_time_index.tz_convert(
                    None
                )

            if incoming_time_index.tz is not None:
                incoming_time_index = incoming_time_index.tz_convert(
                    None
                )

            # -------------------------------------------------------------------------
            # Find missing and extra timestamps
            missing_time_index = saved_time_index.difference(
                incoming_time_index
            )

            extra_time_index = incoming_time_index.difference(
                saved_time_index
            )

            time_axis_matches = incoming_time_index.equals(
                saved_time_index
            )

            if not time_axis_matches:

                log_stream.warning(
                    ' ===> Time axis adjusted for point "%s": '
                    "stored file has %d step(s), incoming series has %d step(s)",
                    point_name,
                    saved_time_index.size,
                    incoming_time_index.size
                )

                if missing_time_index.size > 0:
                    missing_percentage = (
                            missing_time_index.size
                            / saved_time_index.size
                            * 100.0
                    )

                    log_stream.warning(
                        ' ===> Point "%s": %d missing time step(s) '
                        "(%.2f%%) will be filled with the NetCDF fill value",
                        point_name,
                        missing_time_index.size,
                        missing_percentage
                    )

                    log_stream.warning(
                        ' ===> Point "%s": missing time range is '
                        "%s --> %s",
                        point_name,
                        missing_time_index.min().strftime(
                            "%Y-%m-%d %H:%M:%S"
                        ),
                        missing_time_index.max().strftime(
                            "%Y-%m-%d %H:%M:%S"
                        )
                    )

                if extra_time_index.size > 0:
                    extra_percentage = (
                            extra_time_index.size
                            / incoming_time_index.size
                            * 100.0
                    )

                    log_stream.warning(
                        ' ===> Point "%s": %d extra time step(s) '
                        "(%.2f%%) will be removed",
                        point_name,
                        extra_time_index.size,
                        extra_percentage
                    )

                    log_stream.warning(
                        ' ===> Point "%s": removed time range is '
                        "%s --> %s",
                        point_name,
                        extra_time_index.min().strftime(
                            "%Y-%m-%d %H:%M:%S"
                        ),
                        extra_time_index.max().strftime(
                            "%Y-%m-%d %H:%M:%S"
                        )
                    )

                if (
                        missing_time_index.size == 0
                        and extra_time_index.size == 0):
                    log_stream.warning(
                        ' ===> Point "%s": timestamps are the same, '
                        "but their order differs; data will be reordered",
                        point_name
                    )

                # Reorder values, remove extra times and add missing times
                organized_dframe = organized_dframe.reindex(
                    saved_time_index
                )

            # Use the saved time axis as the final reference
            encoded_time = saved_time_values
            time_index = saved_time_index

            # -------------------------------------------------------------------------
            # Rebuild data arrays after time-axis rearrangement
            data_collection = {}

            for variable_name in destination_names:

                if variable_name not in organized_dframe.columns:
                    raise RuntimeError(
                        f'Data variable "{variable_name}" is missing '
                        f'for point "{point_name}". Available columns are '
                        f"{list(organized_dframe.columns)}"
                    )

                data_values = pd.to_numeric(
                    organized_dframe[variable_name],
                    errors="coerce"
                ).to_numpy(
                    dtype=output_dtype,
                    copy=True
                )

                invalid_mask = ~np.isfinite(
                    data_values
                )

                if np.any(invalid_mask):
                    log_stream.warning(
                        ' ===> Point "%s", variable "%s": '
                        "%d missing or invalid value(s) will be replaced "
                        "with fill value %s",
                        point_name,
                        variable_name,
                        int(np.count_nonzero(invalid_mask)),
                        str(file_fill_value)
                    )

                    data_values[
                        invalid_mask
                    ] = file_fill_value

                data_collection[
                    variable_name
                ] = data_values

            # -------------------------------------------------------------------------
            # Check data variables
            missing_variables = [
                variable_name
                for variable_name in destination_names
                if variable_name not in file_handle.variables
            ]

            if missing_variables:
                raise RuntimeError(
                    "Missing data variables in existing NetCDF file: "
                    f"{missing_variables}"
                )

            expected_dimensions = (
                point_dim_name,
                time_name
            )

            for variable_name in destination_names:

                data_variable = file_handle.variables[
                    variable_name
                ]

                if data_variable.dimensions != expected_dimensions:
                    raise RuntimeError(
                        f'Variable "{variable_name}" has dimensions '
                        f"{data_variable.dimensions}; expected "
                        f"{expected_dimensions}"
                    )

                if data_collection[variable_name].size != saved_time_values.size:
                    raise RuntimeError(
                        f'Variable "{variable_name}" for point '
                        f'"{point_name}" has '
                        f"{data_collection[variable_name].size} values after "
                        f"reindexing; expected {saved_time_values.size}"
                    )

            # -------------------------------------------------------------------------
            # Define next point row
            point_index = len(
                point_dimension
            )

            # -------------------------------------------------------------------------
            # Create metadata variables that do not exist yet
            for metadata_name, metadata_value in point_metadata.items():

                if metadata_name in file_handle.variables:
                    continue

                metadata_variable = _create_metadata_variable(
                    file_handle=file_handle,
                    metadata_name=metadata_name,
                    metadata_value=metadata_value
                )

                for previous_point_index in range(point_index):
                    _write_metadata_value(
                        metadata_variable=metadata_variable,
                        point_index=previous_point_index,
                        metadata_value=None
                    )

            # -------------------------------------------------------------------------
            # Write metadata for the new point
            metadata_variable_names = [
                variable_name
                for variable_name, variable_obj
                in file_handle.variables.items()
                if variable_obj.dimensions == (point_dim_name,)
            ]

            for metadata_name in metadata_variable_names:
                metadata_variable = file_handle.variables[
                    metadata_name
                ]

                metadata_value = point_metadata.get(
                    metadata_name
                )

                _write_metadata_value(
                    metadata_variable=metadata_variable,
                    point_index=point_index,
                    metadata_value=metadata_value
                )

            # -------------------------------------------------------------------------
            # Write time-series row
            for variable_name, data_values in data_collection.items():
                file_handle.variables[
                    variable_name
                ][point_index, :] = data_values

            # -------------------------------------------------------------------------
            # Update global information
            file_handle.number_of_points = point_index + 1

    log_stream.info(
        ' ::: Save point "%s" [%s] saved at row %d '
        "with %d time step(s) and %d data variable(s)",
        point_name,
        point_id,
        point_index,
        number_of_time_steps,
        len(destination_names)
    )
# ----------------------------------------------------------------------------------------------------------------------
