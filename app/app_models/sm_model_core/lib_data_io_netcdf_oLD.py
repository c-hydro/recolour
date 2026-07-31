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
import logging

from pathlib import Path
from typing import Mapping, Any, Dict

import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num, num2date
from lib_info_args import logger_name, time_format_algorithm

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helpers
# method to convert a NetCDF string value to a standard Python string
def _convert_nc_string(value):
    if np.ma.is_masked(value):
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)

# method to check whether a value is numeric
def _is_numeric_value(value):
    if value is None or pd.isna(value):
        return True
    if isinstance(value, (bool, np.bool_)):
        return False
    return isinstance(value,(int,float,complex,np.integer,np.floating,np.complexfloating))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to write auxiliary point metrics in NetCDF format
# method to write auxiliary metrics organized by point
def write_auxiliary_nc(
        file_name,
        file_dframe,
        point_tag,
        point_name=None,
        point_lon=None,
        point_lat=None,
        file_fields=None,
        file_dtype="f4",
        file_no_data=-9999.0,
        compression=True,
        compression_level=4,
        overwrite_point=True,
        **kwargs):
    """
    Save auxiliary metrics using one NetCDF variable for each point.

    NetCDF structure
    ----------------
    dimensions:
        auxiliary
        point

    auxiliary information:
        auxiliary_index(auxiliary)
        auxiliary_name(auxiliary)

    point metadata:
        point_id(point)
        point_name(point)
        longitude(point)
        latitude(point)

    point variables:
        point_1(auxiliary)
        point_2(auxiliary)
        point_3(auxiliary)

    Each point variable contains an array with shape:

        (n_auxiliary_variables,)

    Important
    ---------
    The DataFrame columns are expected to already use the final auxiliary
    variable names.

    Example
    -------
    file_fields = {
        "metric_1": "r",
        "metric_2": "bias",
        "metric_3": "ubrmsd"
    }

    The DataFrame must already contain:

        r, bias, ubrmsd

    The resulting NetCDF structure is:

        auxiliary_name[:] = ["r", "bias", "ubrmsd"]

        point_1[0] -> r
        point_1[1] -> bias
        point_1[2] -> ubrmsd
    """

    # -------------------------------------------------------------------------------------
    # Validate output file
    if not isinstance(file_name, (str, os.PathLike)):
        raise TypeError(
            'The argument "file_name" must be a string or path-like object'
        )

    file_name = os.fspath(file_name)

    if not file_name:
        raise ValueError(
            'The argument "file_name" must be defined'
        )

    # -------------------------------------------------------------------------------------
    # Validate dataframe
    if not isinstance(file_dframe, pd.DataFrame):
        raise TypeError(
            'The object "file_dframe" must be a pandas DataFrame'
        )

    if file_dframe.empty:
        raise RuntimeError(
            f'No auxiliary metrics available for point "{point_tag}"'
        )

    # -------------------------------------------------------------------------------------
    # Validate point tag
    if point_tag is None:
        raise RuntimeError(
            'The argument "point_tag" must be defined'
        )

    point_tag = str(point_tag).strip()

    if not point_tag:
        raise ValueError(
            'The argument "point_tag" must be a non-empty value'
        )

    if point_name is None:
        point_name = point_tag
    else:
        point_name = str(point_name)

    # -------------------------------------------------------------------------------------
    # Validate compression
    if not isinstance(compression_level, int):
        raise TypeError(
            'The argument "compression_level" must be an integer'
        )

    if compression_level < 0 or compression_level > 9:
        raise ValueError(
            'The argument "compression_level" must be between 0 and 9'
        )

    # -------------------------------------------------------------------------------------
    # Validate no-data value
    try:
        file_no_data = float(file_no_data)

    except (TypeError, ValueError) as exc:

        raise ValueError(
            'The argument "file_no_data" must be numeric'
        ) from exc

    # -------------------------------------------------------------------------------------
    # Work on a copy
    metrics_dframe = file_dframe.copy()

    # Only one row is expected
    if metrics_dframe.shape[0] > 1:

        log_stream.warning(
            ' ===> Auxiliary dataframe for point "%s" contains %d rows. '
            'Only the first row will be saved',
            point_tag,
            metrics_dframe.shape[0]
        )

    # -------------------------------------------------------------------------------------
    # Define final auxiliary variable names
    #
    # The dataframe already contains the final names. The values from
    # file_fields are used only to define the names and their order.
    if file_fields is None:

        auxiliary_names = [
            str(column_name)
            for column_name in metrics_dframe.columns
        ]

    else:

        if not isinstance(file_fields, Mapping):
            raise TypeError(
                'The argument "file_fields" must be a mapping'
            )

        auxiliary_names = [
            str(destination_name)
            for destination_name in file_fields.values()
        ]

    if not auxiliary_names:
        raise RuntimeError(
            f'No auxiliary variables configured for point "{point_tag}"'
        )

    # -------------------------------------------------------------------------------------
    # Remove reserved variables
    reserved_variable_names = {
        "point_id",
        "point_name",
        "longitude",
        "latitude",
        "auxiliary_index",
        "auxiliary_name"
    }

    invalid_variable_names = [
        variable_name
        for variable_name in auxiliary_names
        if variable_name in reserved_variable_names
    ]

    if invalid_variable_names:

        log_stream.warning(
            ' ===> Removing reserved auxiliary variable(s) '
            'for point "%s": %s',
            point_tag,
            invalid_variable_names
        )

        auxiliary_names = [
            variable_name
            for variable_name in auxiliary_names
            if variable_name not in reserved_variable_names
        ]

    if not auxiliary_names:
        raise RuntimeError(
            f'No auxiliary variables available for point "{point_tag}" '
            f'after removing reserved fields'
        )

    # -------------------------------------------------------------------------------------
    # Check duplicated auxiliary variable names
    if len(auxiliary_names) != len(set(auxiliary_names)):
        raise ValueError(
            f'Duplicated auxiliary variable names detected: '
            f'{auxiliary_names}'
        )

    # -------------------------------------------------------------------------------------
    # Check variables in dataframe
    missing_variables = [
        variable_name
        for variable_name in auxiliary_names
        if variable_name not in metrics_dframe.columns
    ]

    if missing_variables:
        raise KeyError(
            f'Missing final auxiliary variables for point "{point_tag}": '
            f'{missing_variables}. Available fields are '
            f'{list(metrics_dframe.columns)}'
        )

    # Select only configured auxiliary variables
    metrics_dframe = metrics_dframe[
        auxiliary_names
    ].copy()

    # -------------------------------------------------------------------------------------
    # Use first row and convert values to numeric
    metrics_series = metrics_dframe.iloc[0].copy()

    for auxiliary_name in auxiliary_names:

        metrics_series[auxiliary_name] = pd.to_numeric(
            metrics_series[auxiliary_name],
            errors="coerce"
        )

    auxiliary_values = metrics_series[
        auxiliary_names
    ].to_numpy(
        dtype=np.dtype(file_dtype),
        copy=True
    )

    auxiliary_values[
        ~np.isfinite(auxiliary_values)
    ] = file_no_data

    expected_shape = (
        len(auxiliary_names),
    )

    if auxiliary_values.shape != expected_shape:
        raise RuntimeError(
            f'Unexpected auxiliary array shape for point "{point_tag}": '
            f'{auxiliary_values.shape}; expected {expected_shape}'
        )

    # -------------------------------------------------------------------------------------
    # Organize coordinates
    if point_lon is None:
        point_lon = np.nan
    else:

        try:
            point_lon = float(point_lon)

        except (TypeError, ValueError):
            point_lon = np.nan

        if not np.isfinite(point_lon):
            point_lon = np.nan

    if point_lat is None:

        point_lat = np.nan

    else:

        try:
            point_lat = float(point_lat)

        except (TypeError, ValueError):
            point_lat = np.nan

        if not np.isfinite(point_lat):
            point_lat = np.nan

    # -------------------------------------------------------------------------------------
    # Create output folder
    folder_name = os.path.dirname(
        file_name
    )

    if folder_name:
        os.makedirs(
            folder_name,
            exist_ok=True
        )

    file_exists = os.path.exists(
        file_name
    )

    auxiliary_indices = np.arange(
        len(auxiliary_names),
        dtype=np.int32
    )

    # -------------------------------------------------------------------------------------
    # Create a new NetCDF file
    if not file_exists:

        with Dataset(
                file_name,
                mode="w",
                format="NETCDF4") as file_handle:

            # -------------------------------------------------------------------------
            # Dimensions
            file_handle.createDimension(
                "auxiliary",
                len(auxiliary_names)
            )

            file_handle.createDimension(
                "point",
                None
            )

            # -------------------------------------------------------------------------
            # Auxiliary index
            auxiliary_index_var = file_handle.createVariable(
                "auxiliary_index",
                "i4",
                ("auxiliary",)
            )

            auxiliary_index_var[:] = auxiliary_indices

            auxiliary_index_var.long_name = (
                "zero-based auxiliary metric position "
                "in each point array"
            )

            # -------------------------------------------------------------------------
            # Auxiliary variable name
            auxiliary_name_var = file_handle.createVariable(
                "auxiliary_name",
                str,
                ("auxiliary",)
            )

            auxiliary_name_var[:] = np.asarray(
                auxiliary_names,
                dtype=object
            )

            auxiliary_name_var.long_name = (
                "auxiliary metric name associated with each array position"
            )

            # -------------------------------------------------------------------------
            # Point identifier
            point_id_var = file_handle.createVariable(
                "point_id",
                str,
                ("point",)
            )

            point_id_var.long_name = "point identifier"
            point_id_var.cf_role = "timeseries_id"

            # -------------------------------------------------------------------------
            # Point descriptive name
            point_name_var = file_handle.createVariable(
                "point_name",
                str,
                ("point",)
            )

            point_name_var.long_name = "point descriptive name"

            # -------------------------------------------------------------------------
            # Longitude
            longitude_var = file_handle.createVariable(
                "longitude",
                "f8",
                ("point",),
                fill_value=np.nan
            )

            longitude_var.units = "degrees_east"
            longitude_var.standard_name = "longitude"
            longitude_var.long_name = "point longitude"
            longitude_var.axis = "X"

            # -------------------------------------------------------------------------
            # Latitude
            latitude_var = file_handle.createVariable(
                "latitude",
                "f8",
                ("point",),
                fill_value=np.nan
            )

            latitude_var.units = "degrees_north"
            latitude_var.standard_name = "latitude"
            latitude_var.long_name = "point latitude"
            latitude_var.axis = "Y"

            # -------------------------------------------------------------------------
            # Global attributes
            file_handle.title = "Point auxiliary metrics collection"
            file_handle.featureType = "point"
            file_handle.Conventions = "CF-1.8"

            file_handle.number_of_auxiliary_variables = len(
                auxiliary_names
            )

    # -------------------------------------------------------------------------------------
    # Append or update point
    with Dataset(
            file_name,
            mode="a") as file_handle:

        # ---------------------------------------------------------------------------------
        # Validate required dimensions
        required_dimensions = [
            "auxiliary",
            "point"
        ]

        for dimension_name in required_dimensions:

            if dimension_name not in file_handle.dimensions:
                raise RuntimeError(
                    f'Dimension "{dimension_name}" is missing '
                    f'from NetCDF file "{file_name}"'
                )

        # ---------------------------------------------------------------------------------
        # Validate required variables
        required_variables = [
            "auxiliary_index",
            "auxiliary_name",
            "point_id",
            "point_name",
            "longitude",
            "latitude"
        ]

        for required_variable in required_variables:

            if required_variable not in file_handle.variables:
                raise RuntimeError(
                    f'Variable "{required_variable}" is missing '
                    f'from NetCDF file "{file_name}"'
                )

        # ---------------------------------------------------------------------------------
        # Validate auxiliary dimension
        file_auxiliary_size = len(
            file_handle.dimensions["auxiliary"]
        )

        if file_auxiliary_size != len(auxiliary_names):
            raise RuntimeError(
                f'Auxiliary dimension mismatch for point "{point_tag}": '
                f'file has {file_auxiliary_size} variable(s), '
                f'point has {len(auxiliary_names)} variable(s)'
            )

        # ---------------------------------------------------------------------------------
        # Validate stored auxiliary indexes
        stored_auxiliary_indices = np.asarray(
            file_handle.variables["auxiliary_index"][:],
            dtype=np.int32
        )

        if not np.array_equal(
                stored_auxiliary_indices,
                auxiliary_indices):

            raise RuntimeError(
                f'Auxiliary-index mismatch for point "{point_tag}": '
                f'stored={stored_auxiliary_indices.tolist()}, '
                f'current={auxiliary_indices.tolist()}'
            )

        # ---------------------------------------------------------------------------------
        # Validate stored auxiliary names
        stored_auxiliary_names = [
            _convert_nc_string(value)
            for value in file_handle.variables["auxiliary_name"][:]
        ]

        if stored_auxiliary_names != auxiliary_names:
            raise RuntimeError(
                f'Auxiliary-name mismatch for point "{point_tag}": '
                f'stored={stored_auxiliary_names}, '
                f'current={auxiliary_names}'
            )

        # ---------------------------------------------------------------------------------
        # Validate metadata dimensions
        metadata_dimensions = {
            "point_id": ("point",),
            "point_name": ("point",),
            "longitude": ("point",),
            "latitude": ("point",)
        }

        for variable_name, expected_dimensions in metadata_dimensions.items():

            current_dimensions = file_handle.variables[
                variable_name
            ].dimensions

            if current_dimensions != expected_dimensions:
                raise RuntimeError(
                    f'Invalid dimensions for variable "{variable_name}" '
                    f'in "{file_name}": {current_dimensions}; '
                    f'expected {expected_dimensions}'
                )

        # ---------------------------------------------------------------------------------
        # Find the point index
        existing_point_ids = [
            _convert_nc_string(value)
            for value in file_handle.variables["point_id"][:]
        ]

        if point_tag in existing_point_ids:

            point_idx = existing_point_ids.index(
                point_tag
            )

            if not overwrite_point:
                raise RuntimeError(
                    f'Point "{point_tag}" already exists '
                    f'in NetCDF file "{file_name}"'
                )

            log_stream.warning(
                ' ===> Point "%s" already exists at index %d. '
                'Existing auxiliary values will be overwritten',
                point_tag,
                point_idx
            )

        else:

            point_idx = len(
                existing_point_ids
            )

        # ---------------------------------------------------------------------------------
        # Write point metadata
        file_handle.variables[
            "point_id"
        ][point_idx] = point_tag

        file_handle.variables[
            "point_name"
        ][point_idx] = point_name

        file_handle.variables[
            "longitude"
        ][point_idx] = point_lon

        file_handle.variables[
            "latitude"
        ][point_idx] = point_lat

        # ---------------------------------------------------------------------------------
        # Create or validate point auxiliary variable
        if point_tag not in file_handle.variables:

            point_var = file_handle.createVariable(
                point_tag,
                file_dtype,
                ("auxiliary",),
                fill_value=file_no_data,
                zlib=bool(compression),
                complevel=(
                    compression_level
                    if compression
                    else 0
                ),
                shuffle=bool(compression)
            )

            point_var.long_name = (
                f"auxiliary metrics for point {point_name}"
            )

            point_var.point_id = point_tag
            point_var.point_name = point_name
            point_var.point_index = point_idx
            point_var.longitude = point_lon
            point_var.latitude = point_lat

            point_var.coordinates = (
                "auxiliary_index auxiliary_name"
            )

        else:

            point_var = file_handle.variables[
                point_tag
            ]

            expected_dimensions = (
                "auxiliary",
            )

            if point_var.dimensions != expected_dimensions:
                raise RuntimeError(
                    f'Invalid dimensions for point variable "{point_tag}": '
                    f'{point_var.dimensions}; expected '
                    f'{expected_dimensions}'
                )

            if point_var.shape != expected_shape:
                raise RuntimeError(
                    f'Invalid shape for point variable "{point_tag}": '
                    f'{point_var.shape}; expected {expected_shape}'
                )

            point_var.point_id = point_tag
            point_var.point_name = point_name
            point_var.point_index = point_idx
            point_var.longitude = point_lon
            point_var.latitude = point_lat

        # ---------------------------------------------------------------------------------
        # Write complete auxiliary array
        point_var[:] = auxiliary_values

        # ---------------------------------------------------------------------------------
        # Update global attributes
        file_handle.number_of_points = len(
            file_handle.dimensions["point"]
        )

        file_handle.number_of_auxiliary_variables = len(
            file_handle.dimensions["auxiliary"]
        )

    log_stream.info(
        ' ---> Auxiliary metrics for point "%s" saved in NetCDF file "%s" '
        'at point index %d with array shape (%d,)',
        point_tag,
        file_name,
        point_idx,
        auxiliary_values.shape[0]
    )

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read point results from a NetCDF collection
def read_results_nc(
        file_name,
        file_fields=None,
        registry_fields=None,
        time_reference=None,
        time_start=None,
        time_end=None,
        time_rounding="H",
        time_frequency="H",
        ascending_index=False,
        sort_index=True,
        time_name="time",
        variable_dimension_name="variable",
        variable_index_name="variable_index",
        variable_name_name="variable_name",
        point_dimension_name="point",
        point_id_name="point_id",
        point_name_name="point_name",
        longitude_name="longitude",
        latitude_name="latitude",
        masked_to_nan=True,
        **kwargs) -> Dict[str, pd.DataFrame]:
    """
    Read all point time series from a NetCDF collection.

    The NetCDF file is opened only once.

    Expected NetCDF structure
    -------------------------
    dimensions:
        time
        variable
        point

    coordinate and metadata variables:
        time(time)
        variable_index(variable)
        variable_name(variable)

        point_id(point)
        point_name(point)
        longitude(point)
        latitude(point)

    point variables:
        point_1(time, variable)
        point_2(time, variable)
        point_3(time, variable)

    Returned object
    ---------------
    {
        "point_1": DataFrame,
        "point_2": DataFrame,
        ...
    }

    Each DataFrame contains:

                         rain    airt      sm
        time
        2026-01-01 00:00   0.0    12.4    45.2
        2026-01-01 01:00   1.2    11.8    44.9

    DataFrame attributes contain:

        point_id
        point_name
        longitude
        latitude
        point_index
        file_name
        variable_names
    """

    # -------------------------------------------------------------------------------------
    # Initialize optional arguments
    if file_fields is None:
        file_fields = {}

    if registry_fields is None:
        registry_fields = {}

    if not isinstance(file_fields, Mapping):
        raise TypeError(
            'The argument "file_fields" must be a mapping'
        )

    if not isinstance(registry_fields, Mapping):
        raise TypeError(
            'The argument "registry_fields" must be a mapping'
        )

    # -------------------------------------------------------------------------------------
    # Validate file
    file_path_obj = Path(
        file_name
    )

    if not file_path_obj.exists():
        raise FileNotFoundError(
            f'NetCDF file not found: "{file_path_obj}"'
        )

    if not file_path_obj.is_file():
        raise FileNotFoundError(
            f'NetCDF path is not a file: "{file_path_obj}"'
        )

    # -------------------------------------------------------------------------------------
    # Initialize output object
    points_obj: Dict[str, pd.DataFrame] = {}

    expected_dimensions = (
        time_name,
        variable_dimension_name
    )

    # -------------------------------------------------------------------------------------
    # Open NetCDF once
    try:

        with Dataset(
                str(file_path_obj),
                mode="r") as file_nc:

            # ---------------------------------------------------------------------------------
            # Validate required dimensions
            required_dimensions = [
                time_name,
                variable_dimension_name,
                point_dimension_name
            ]

            missing_dimensions = [
                dimension_name
                for dimension_name in required_dimensions
                if dimension_name not in file_nc.dimensions
            ]

            if missing_dimensions:
                raise RuntimeError(
                    f'Missing dimensions in "{file_path_obj}": '
                    f'{missing_dimensions}. Available dimensions are '
                    f'{list(file_nc.dimensions.keys())}'
                )

            number_of_time_steps = len(
                file_nc.dimensions[time_name]
            )

            number_of_variables = len(
                file_nc.dimensions[variable_dimension_name]
            )

            number_of_points = len(
                file_nc.dimensions[point_dimension_name]
            )

            if number_of_time_steps == 0:
                raise RuntimeError(
                    f'Dimension "{time_name}" is empty '
                    f'in "{file_path_obj}"'
                )

            if number_of_variables == 0:
                raise RuntimeError(
                    f'Dimension "{variable_dimension_name}" is empty '
                    f'in "{file_path_obj}"'
                )

            if number_of_points == 0:
                log_stream.warning(
                    ' ===> Point dimension "%s" is empty in "%s"',
                    point_dimension_name,
                    file_path_obj
                )

            # ---------------------------------------------------------------------------------
            # Validate required variables
            required_variables = [
                time_name,
                variable_index_name,
                variable_name_name,
                point_id_name,
                point_name_name,
                longitude_name,
                latitude_name
            ]

            missing_variables = [
                variable_name
                for variable_name in required_variables
                if variable_name not in file_nc.variables
            ]

            if missing_variables:
                raise RuntimeError(
                    f'Missing variables in "{file_path_obj}": '
                    f'{missing_variables}. Available variables are '
                    f'{list(file_nc.variables.keys())}'
                )

            # ---------------------------------------------------------------------------------
            # Read time coordinate
            time_variable = file_nc.variables[
                time_name
            ]

            if not hasattr(time_variable, "units"):
                raise RuntimeError(
                    f'Time variable "{time_name}" does not have '
                    f'a "units" attribute in "{file_path_obj}"'
                )

            time_calendar = getattr(
                time_variable,
                "calendar",
                "standard"
            )

            time_values_raw = time_variable[:]

            try:

                time_values = num2date(
                    time_values_raw,
                    units=time_variable.units,
                    calendar=time_calendar,
                    only_use_cftime_datetimes=False,
                    only_use_python_datetimes=True
                )

                time_values = pd.to_datetime(
                    np.asarray(time_values)
                )

            except Exception:

                time_values_cftime = num2date(
                    time_values_raw,
                    units=time_variable.units,
                    calendar=time_calendar,
                    only_use_cftime_datetimes=True
                )

                time_values = pd.to_datetime([
                    value.strftime("%Y-%m-%d %H:%M:%S")
                    for value in time_values_cftime
                ])

            time_index = pd.DatetimeIndex(
                time_values,
                name="time"
            )

            if len(time_index) != number_of_time_steps:
                raise RuntimeError(
                    f'Time variable "{time_name}" contains '
                    f'{len(time_index)} values, but dimension '
                    f'"{time_name}" has length {number_of_time_steps}'
                )

            if time_index.has_duplicates:

                duplicated_times = time_index[
                    time_index.duplicated(
                        keep=False
                    )
                ].unique()

                raise RuntimeError(
                    f'Duplicated timestamps found in "{file_path_obj}": '
                    f'{duplicated_times[:10].tolist()}'
                )

            # ---------------------------------------------------------------------------------
            # Read variable indices
            variable_indices = np.asarray(
                file_nc.variables[variable_index_name][:],
                dtype=np.int32
            ).reshape(-1)

            if variable_indices.size != number_of_variables:
                raise RuntimeError(
                    f'Variable "{variable_index_name}" contains '
                    f'{variable_indices.size} values, but dimension '
                    f'"{variable_dimension_name}" has length '
                    f'{number_of_variables}'
                )

            if np.any(variable_indices < 0):
                raise RuntimeError(
                    f'Negative variable indices found in '
                    f'"{variable_index_name}": '
                    f'{variable_indices.tolist()}'
                )

            if np.any(variable_indices >= number_of_variables):
                raise RuntimeError(
                    f'Variable indices exceed dimension length: '
                    f'{variable_indices.tolist()}. Maximum allowed '
                    f'index is {number_of_variables - 1}'
                )

            if len(np.unique(variable_indices)) != len(variable_indices):
                raise RuntimeError(
                    f'Duplicated variable indices found: '
                    f'{variable_indices.tolist()}'
                )

            # ---------------------------------------------------------------------------------
            # Read final variable names
            variable_names = _decode_netcdf_strings(
                file_nc.variables[variable_name_name][:]
            )

            variable_names = [
                str(variable_name).strip()
                for variable_name in variable_names
            ]

            if len(variable_names) != number_of_variables:
                raise RuntimeError(
                    f'Variable "{variable_name_name}" contains '
                    f'{len(variable_names)} values, but dimension '
                    f'"{variable_dimension_name}" has length '
                    f'{number_of_variables}'
                )

            empty_variable_names = [
                variable_position
                for variable_position, variable_name
                in enumerate(variable_names)
                if not variable_name
            ]

            if empty_variable_names:
                raise RuntimeError(
                    f'Empty variable names found at positions: '
                    f'{empty_variable_names}'
                )

            if len(set(variable_names)) != len(variable_names):
                raise RuntimeError(
                    f'Duplicated variable names found: '
                    f'{variable_names}'
                )

            # ---------------------------------------------------------------------------------
            # Optionally map NetCDF final names back to application field names
            #
            # Example:
            #
            # file_fields = {
            #     "values_1": "rain",
            #     "values_2": "airt",
            #     "values_3": "sm"
            # }
            #
            # NetCDF:
            #     rain, airt, sm
            #
            # Returned DataFrame:
            #     values_1, values_2, values_3
            variable_output_names = []

            for stored_variable_name in variable_names:

                output_variable_name = stored_variable_name

                for field_name, final_variable_name in file_fields.items():

                    if field_name == time_name:
                        continue

                    if str(final_variable_name) == stored_variable_name:
                        output_variable_name = str(field_name)
                        break

                variable_output_names.append(
                    output_variable_name
                )

            if len(set(variable_output_names)) != len(variable_output_names):
                raise RuntimeError(
                    f'Duplicated output names after applying '
                    f'"file_fields": {variable_output_names}'
                )

            # ---------------------------------------------------------------------------------
            # Read point metadata
            point_ids = _decode_netcdf_strings(
                file_nc.variables[point_id_name][:]
            )

            point_names = _decode_netcdf_strings(
                file_nc.variables[point_name_name][:]
            )

            point_longitudes = np.asarray(
                file_nc.variables[longitude_name][:],
                dtype=np.float64
            ).reshape(-1)

            point_latitudes = np.asarray(
                file_nc.variables[latitude_name][:],
                dtype=np.float64
            ).reshape(-1)

            # ---------------------------------------------------------------------------------
            # Validate point metadata lengths
            point_metadata_lengths = {
                point_id_name: len(point_ids),
                point_name_name: len(point_names),
                longitude_name: len(point_longitudes),
                latitude_name: len(point_latitudes)
            }

            for metadata_name, metadata_length in point_metadata_lengths.items():

                if metadata_length != number_of_points:
                    raise RuntimeError(
                        f'Variable "{metadata_name}" contains '
                        f'{metadata_length} values, but dimension '
                        f'"{point_dimension_name}" has length '
                        f'{number_of_points}'
                    )

            point_ids = [
                str(point_id).strip()
                for point_id in point_ids
            ]

            point_names = [
                str(point_name).strip()
                for point_name in point_names
            ]

            if len(set(point_ids)) != len(point_ids):
                raise RuntimeError(
                    f'Duplicated point IDs found in "{file_path_obj}": '
                    f'{point_ids}'
                )

            # ---------------------------------------------------------------------------------
            # Define optional requested time range
            if time_start is not None or time_end is not None:

                if time_start is None:
                    time_start_expected = time_index.min()
                else:
                    time_start_expected = pd.Timestamp(
                        time_start
                    ).floor(
                        str(time_rounding).lower()
                    )

                if time_end is None:
                    time_end_expected = time_index.max()
                else:
                    time_end_expected = pd.Timestamp(
                        time_end
                    ).floor(
                        str(time_rounding).lower()
                    )

                if time_start_expected > time_end_expected:
                    raise RuntimeError(
                        f'Invalid time range: start '
                        f'"{time_start_expected}" is after end '
                        f'"{time_end_expected}"'
                    )

                time_range = pd.date_range(
                    start=time_start_expected,
                    end=time_end_expected,
                    freq=str(time_frequency).lower(),
                    name="time"
                )

            else:

                time_range = None

            # ---------------------------------------------------------------------------------
            # Metadata and coordinate variables to exclude
            excluded_variables = {
                time_name,
                variable_index_name,
                variable_name_name,
                point_id_name,
                point_name_name,
                longitude_name,
                latitude_name
            }

            expected_shape = (
                number_of_time_steps,
                number_of_variables
            )

            # ---------------------------------------------------------------------------------
            # Read all point variables
            for netcdf_variable_name, variable_obj in file_nc.variables.items():

                # Skip NetCDF coordinates and metadata
                if netcdf_variable_name in excluded_variables:
                    continue

                # A point variable must have dimensions:
                #
                # point_1(time, variable)
                if variable_obj.dimensions != expected_dimensions:
                    continue

                # ---------------------------------------------------------------------------------
                # Determine point metadata index
                point_index_attribute = getattr(
                    variable_obj,
                    "point_index",
                    None
                )

                point_id_attribute = str(
                    getattr(
                        variable_obj,
                        "point_id",
                        netcdf_variable_name
                    )
                ).strip()

                point_metadata_index = None

                if point_index_attribute is not None:

                    try:
                        point_metadata_index = int(
                            point_index_attribute
                        )
                    except (TypeError, ValueError):
                        point_metadata_index = None

                    if (
                            point_metadata_index is not None
                            and (
                                point_metadata_index < 0
                                or point_metadata_index >= number_of_points
                            )):

                        log_stream.warning(
                            ' ===> Invalid point_index=%s for variable "%s" '
                            'in "%s". Point ID matching will be used',
                            str(point_index_attribute),
                            netcdf_variable_name,
                            file_path_obj
                        )

                        point_metadata_index = None

                # Fallback using point ID
                if point_metadata_index is None:

                    if point_id_attribute in point_ids:

                        point_metadata_index = point_ids.index(
                            point_id_attribute
                        )

                    elif netcdf_variable_name in point_ids:

                        point_metadata_index = point_ids.index(
                            netcdf_variable_name
                        )

                # ---------------------------------------------------------------------------------
                # Read point matrix
                point_values = variable_obj[:]

                if np.ma.isMaskedArray(point_values):

                    if masked_to_nan:
                        point_values = point_values.filled(
                            np.nan
                        )
                    else:
                        point_values = point_values.data

                point_values = np.array(
                    point_values,
                    dtype=np.float32,
                    copy=True
                )

                if point_values.shape != expected_shape:
                    raise RuntimeError(
                        f'Unexpected shape for point variable '
                        f'"{netcdf_variable_name}": '
                        f'{point_values.shape}; expected '
                        f'{expected_shape}'
                    )

                # ---------------------------------------------------------------------------------
                # Convert explicit fill values to NaN
                if masked_to_nan:

                    fill_values = []

                    if hasattr(variable_obj, "_FillValue"):

                        fill_values.append(
                            getattr(
                                variable_obj,
                                "_FillValue"
                            )
                        )

                    if hasattr(variable_obj, "missing_value"):

                        missing_values = np.asarray(
                            getattr(
                                variable_obj,
                                "missing_value"
                            )
                        ).reshape(-1)

                        fill_values.extend(
                            missing_values.tolist()
                        )

                    for fill_value in fill_values:

                        try:
                            fill_value_float = float(
                                fill_value
                            )
                        except (TypeError, ValueError):
                            continue

                        if np.isnan(fill_value_float):
                            continue

                        point_values[
                            np.isclose(
                                point_values,
                                fill_value_float,
                                rtol=0.0,
                                atol=1.0e-6,
                                equal_nan=False
                            )
                        ] = np.nan

                    point_values[
                        ~np.isfinite(point_values)
                    ] = np.nan

                # ---------------------------------------------------------------------------------
                # Organize point DataFrame
                point_data = {}

                for variable_index, output_variable_name in zip(
                        variable_indices,
                        variable_output_names):

                    point_data[
                        output_variable_name
                    ] = point_values[
                        :,
                        int(variable_index)
                    ].copy()

                point_dframe = pd.DataFrame(
                    data=point_data,
                    index=time_index.copy()
                )

                point_dframe.index.name = "time"

                # ---------------------------------------------------------------------------------
                # Reindex requested time range
                if time_range is not None:

                    point_dframe = point_dframe.reindex(
                        time_range
                    )

                # ---------------------------------------------------------------------------------
                # Sort time
                if sort_index:

                    point_dframe = point_dframe.sort_index(
                        ascending=ascending_index
                    )

                # ---------------------------------------------------------------------------------
                # Organize point metadata
                point_attrs = dict(
                    registry_fields
                )

                point_attrs["time_reference"] = time_reference
                point_attrs["file_name"] = str(file_path_obj)
                point_attrs["netcdf_variable"] = netcdf_variable_name
                point_attrs["variable_names"] = list(variable_names)
                point_attrs["output_variable_names"] = list(
                    variable_output_names
                )

                if point_metadata_index is not None:

                    point_attrs["point_index"] = point_metadata_index
                    point_attrs["point_id"] = point_ids[
                        point_metadata_index
                    ]

                    point_attrs["point_name"] = point_names[
                        point_metadata_index
                    ]

                    point_attrs["longitude"] = float(
                        point_longitudes[
                            point_metadata_index
                        ]
                    )

                    point_attrs["latitude"] = float(
                        point_latitudes[
                            point_metadata_index
                        ]
                    )

                else:

                    log_stream.warning(
                        ' ===> Metadata index was not found for point '
                        'variable "%s" in "%s"',
                        netcdf_variable_name,
                        file_path_obj
                    )

                    point_attrs["point_index"] = None
                    point_attrs["point_id"] = point_id_attribute
                    point_attrs["point_name"] = str(
                        getattr(
                            variable_obj,
                            "point_name",
                            netcdf_variable_name
                        )
                    )

                    point_attrs["longitude"] = float(
                        getattr(
                            variable_obj,
                            "longitude",
                            np.nan
                        )
                    )

                    point_attrs["latitude"] = float(
                        getattr(
                            variable_obj,
                            "latitude",
                            np.nan
                        )
                    )

                # Add point variable attributes
                for attribute_name in variable_obj.ncattrs():

                    try:

                        point_attrs[
                            attribute_name
                        ] = variable_obj.getncattr(
                            attribute_name
                        )

                    except Exception:

                        pass

                point_dframe.attrs = point_attrs

                # ---------------------------------------------------------------------------------
                # Use stored point ID as dictionary key
                point_key = point_attrs.get(
                    "point_id",
                    netcdf_variable_name
                )

                point_key = str(
                    point_key
                ).strip()

                if not point_key:
                    point_key = netcdf_variable_name

                if point_key in points_obj:
                    raise RuntimeError(
                        f'Duplicated output point key "{point_key}" '
                        f'in "{file_path_obj}"'
                    )

                points_obj[
                    point_key
                ] = point_dframe

    except OSError as exc:

        log_stream.error(
            ' ===> Error reading NetCDF file "%s": %s',
            file_path_obj,
            exc
        )

        raise RuntimeError(
            f'Unable to read NetCDF file "{file_path_obj}"'
        ) from exc

    # -------------------------------------------------------------------------------------
    # Ensure points were found
    if not points_obj:
        raise RuntimeError(
            f'No point variables with dimensions '
            f'{expected_dimensions} were found in '
            f'"{file_path_obj}". Available variables are not compatible '
            f'with the expected point-oriented format'
        )

    return points_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write point datasets in a NetCDF collection
def write_results_nc(
        file_name,
        file_dframe,
        point_tag,
        point_name=None,
        point_lon=None,
        point_lat=None,
        file_fields=None,
        time_fields=None,
        time_index_label="time",
        time_units="hours since 1970-01-01 00:00:00",
        time_calendar="standard",
        file_no_data=-9999.0,
        file_dtype="f4",
        compression=True,
        compression_level=4,
        ascending_index=True,
        sort_index=True,
        mismatch_threshold_percentage=50.0,
        overwrite_point=True,
        **kwargs):
    """
    Save one point time series as a dedicated NetCDF variable.

    NetCDF structure
    ----------------
    dimensions:
        time
        variable
        point

    common variables:
        time(time)
        variable_index(variable)
        variable_name(variable)

        point_id(point)
        point_name(point)
        longitude(point)
        latitude(point)

    point variables:
        point_1(time, variable)
        point_2(time, variable)
        point_3(time, variable)

    Each point variable contains an array with shape:

        (n_time_steps, n_variables)

    Important
    ---------
    The columns in ``file_dframe`` are expected to already contain the final
    variable names.

    Example:

        file_fields = {
            "values_1": "rain",
            "values_2": "airt",
            "values_3": "sm"
        }

    The DataFrame must already contain:

        time, rain, airt, sm

    The file_fields mapping is used only to determine:

        variable_names = ["rain", "airt", "sm"]

    No DataFrame column renaming is performed.

    Example NetCDF data:

        variable_name[:] = ["rain", "airt", "sm"]

        point_1[:, 0] -> rain
        point_1[:, 1] -> airt
        point_1[:, 2] -> sm
    """

    # -------------------------------------------------------------------------------------
    # Validate output file
    if not isinstance(file_name, (str, os.PathLike)):
        raise TypeError(
            'The argument "file_name" must be a string or path-like object'
        )

    file_name = os.fspath(file_name)

    if not file_name:
        raise ValueError(
            'The argument "file_name" must be defined'
        )

    # -------------------------------------------------------------------------------------
    # Validate dataframe
    if not isinstance(file_dframe, pd.DataFrame):
        raise TypeError(
            'The object "file_dframe" must be a pandas DataFrame'
        )

    if file_dframe.empty:
        raise RuntimeError(
            f'No data available for point "{point_tag}"'
        )

    # -------------------------------------------------------------------------------------
    # Validate point tag
    if point_tag is None:
        raise RuntimeError(
            'The argument "point_tag" must be defined'
        )

    point_tag = str(point_tag).strip()

    if not point_tag:
        raise ValueError(
            'The argument "point_tag" must be a non-empty value'
        )

    if point_name is None:
        point_name = point_tag
    else:
        point_name = str(point_name)

    # -------------------------------------------------------------------------------------
    # Validate time options
    if not isinstance(time_index_label, str) or not time_index_label:
        raise ValueError(
            'The argument "time_index_label" must be a non-empty string'
        )

    if not isinstance(time_units, str) or not time_units:
        raise ValueError(
            'The argument "time_units" must be a non-empty string'
        )

    if not isinstance(time_calendar, str) or not time_calendar:
        raise ValueError(
            'The argument "time_calendar" must be a non-empty string'
        )

    # -------------------------------------------------------------------------------------
    # Validate compression
    if not isinstance(compression_level, int):
        raise TypeError(
            'The argument "compression_level" must be an integer'
        )

    if compression_level < 0 or compression_level > 9:
        raise ValueError(
            'The argument "compression_level" must be between 0 and 9'
        )

    # -------------------------------------------------------------------------------------
    # Validate mismatch threshold
    if mismatch_threshold_percentage is None:
        mismatch_threshold_percentage = 50.0

    mismatch_threshold_percentage = float(
        mismatch_threshold_percentage
    )

    if not 0.0 <= mismatch_threshold_percentage <= 100.0:
        raise ValueError(
            'The argument "mismatch_threshold_percentage" '
            'must be between 0 and 100'
        )

    # -------------------------------------------------------------------------------------
    # Validate no-data value
    try:
        file_no_data = float(file_no_data)

    except (TypeError, ValueError) as exc:

        raise ValueError(
            'The argument "file_no_data" must be numeric'
        ) from exc

    # -------------------------------------------------------------------------------------
    # Work on a copy
    point_dframe = file_dframe.copy()

    # -------------------------------------------------------------------------------------
    # Organize time index
    if time_index_label in point_dframe.columns:

        point_time_index = pd.DatetimeIndex(
            pd.to_datetime(
                point_dframe[time_index_label],
                errors="coerce"
            )
        )

        point_dframe = point_dframe.drop(
            columns=[time_index_label]
        )

        point_dframe.index = point_time_index

    elif isinstance(point_dframe.index, pd.DatetimeIndex):

        point_dframe.index = pd.DatetimeIndex(
            point_dframe.index
        )

    else:

        point_dframe.index = pd.DatetimeIndex(
            pd.to_datetime(
                point_dframe.index,
                errors="coerce"
            )
        )

    # -------------------------------------------------------------------------------------
    # Remove invalid timestamps
    invalid_time_mask = point_dframe.index.isna()

    if invalid_time_mask.any():

        invalid_time_steps = int(
            invalid_time_mask.sum()
        )

        log_stream.warning(
            ' ===> Removing %d invalid timestamp(s) for point "%s"',
            invalid_time_steps,
            point_tag
        )

        point_dframe = point_dframe.loc[
            ~invalid_time_mask
        ]

    if point_dframe.empty:
        raise RuntimeError(
            f'No valid timestamps available for point "{point_tag}"'
        )

    # -------------------------------------------------------------------------------------
    # Remove duplicate timestamps
    if point_dframe.index.has_duplicates:

        duplicated_time_mask = point_dframe.index.duplicated(
            keep="last"
        )

        duplicated_time_steps = int(
            duplicated_time_mask.sum()
        )

        log_stream.warning(
            ' ===> Removing %d duplicated timestamp(s) for point "%s"; '
            'the last occurrence is retained',
            duplicated_time_steps,
            point_tag
        )

        point_dframe = point_dframe.loc[
            ~duplicated_time_mask
        ]

    # -------------------------------------------------------------------------------------
    # Sort time index
    if sort_index:

        point_dframe = point_dframe.sort_index(
            ascending=ascending_index
        )

    point_time_index = pd.DatetimeIndex(
        point_dframe.index
    )

    # -------------------------------------------------------------------------------------
    # Organize final variable names
    if file_fields is None:

        variable_names = [
            str(column_name)
            for column_name in point_dframe.columns
        ]

    else:

        if not isinstance(file_fields, Mapping):
            raise TypeError(
                'The argument "file_fields" must be a mapping'
            )

        variable_names = [
            str(destination_name)
            for source_name, destination_name in file_fields.items()
            if source_name != time_index_label
        ]

    if not variable_names:
        raise RuntimeError(
            f'No data variables configured for point "{point_tag}"'
        )

    # -------------------------------------------------------------------------------------
    # Check duplicated variable names
    if len(variable_names) != len(set(variable_names)):
        raise ValueError(
            f'Duplicated final variable names detected: {variable_names}'
        )

    # -------------------------------------------------------------------------------------
    # Check final variables in DataFrame
    missing_variables = [
        variable_name
        for variable_name in variable_names
        if variable_name not in point_dframe.columns
    ]

    if missing_variables:
        raise KeyError(
            f'Missing final variables for point "{point_tag}": '
            f'{missing_variables}. Available fields are '
            f'{list(point_dframe.columns)}'
        )

    # Select final variables in configured order
    point_dframe = point_dframe[
        variable_names
    ].copy()

    # -------------------------------------------------------------------------------------
    # Convert all variables to numeric
    for variable_name in variable_names:

        point_dframe[variable_name] = pd.to_numeric(
            point_dframe[variable_name],
            errors="coerce"
        )

    # -------------------------------------------------------------------------------------
    # Validate coordinates
    if point_lon is None:

        point_lon = np.nan

    else:

        try:
            point_lon = float(point_lon)

        except (TypeError, ValueError):

            point_lon = np.nan

        if not np.isfinite(point_lon):
            point_lon = np.nan

    if point_lat is None:

        point_lat = np.nan

    else:

        try:
            point_lat = float(point_lat)

        except (TypeError, ValueError):

            point_lat = np.nan

        if not np.isfinite(point_lat):
            point_lat = np.nan

    # -------------------------------------------------------------------------------------
    # Create output folder
    file_folder = os.path.dirname(
        file_name
    )

    if file_folder:
        os.makedirs(
            file_folder,
            exist_ok=True
        )

    file_exists = os.path.exists(
        file_name
    )

    # Variable positions in the point matrix
    variable_indices = np.arange(
        len(variable_names),
        dtype=np.int32
    )

    # -------------------------------------------------------------------------------------
    # Create new NetCDF file
    if not file_exists:

        point_time_values = np.asarray(
            date2num(
                point_time_index.to_pydatetime(),
                units=time_units,
                calendar=time_calendar
            ),
            dtype=np.float64
        )

        with Dataset(
                file_name,
                mode="w",
                format="NETCDF4") as file_handle:

            # -------------------------------------------------------------------------
            # Dimensions
            file_handle.createDimension(
                "time",
                point_time_index.size
            )

            file_handle.createDimension(
                "variable",
                len(variable_names)
            )

            file_handle.createDimension(
                "point",
                None
            )

            # -------------------------------------------------------------------------
            # Time variable
            time_var = file_handle.createVariable(
                "time",
                "f8",
                ("time",)
            )

            time_var[:] = point_time_values

            time_var.units = time_units
            time_var.calendar = time_calendar
            time_var.standard_name = "time"
            time_var.long_name = "reference time"
            time_var.axis = "T"

            # Optional time attributes
            if time_fields is not None:

                if not isinstance(time_fields, Mapping):
                    raise TypeError(
                        'The argument "time_fields" must be a mapping'
                    )

                for attribute_name, attribute_value in time_fields.items():

                    if attribute_value is not None:

                        setattr(
                            time_var,
                            str(attribute_name),
                            attribute_value
                        )

            # -------------------------------------------------------------------------
            # Variable index
            variable_index_var = file_handle.createVariable(
                "variable_index",
                "i4",
                ("variable",)
            )

            variable_index_var[:] = variable_indices

            variable_index_var.long_name = (
                "zero-based variable position in each point data matrix"
            )

            # -------------------------------------------------------------------------
            # Variable name
            variable_name_var = file_handle.createVariable(
                "variable_name",
                str,
                ("variable",)
            )

            variable_name_var[:] = np.asarray(
                variable_names,
                dtype=object
            )

            variable_name_var.long_name = (
                "variable name associated with each matrix column"
            )

            # -------------------------------------------------------------------------
            # Point ID
            point_id_var = file_handle.createVariable(
                "point_id",
                str,
                ("point",)
            )

            point_id_var.long_name = "point identifier"
            point_id_var.cf_role = "timeseries_id"

            # -------------------------------------------------------------------------
            # Point name
            point_name_var = file_handle.createVariable(
                "point_name",
                str,
                ("point",)
            )

            point_name_var.long_name = "point descriptive name"

            # -------------------------------------------------------------------------
            # Longitude
            longitude_var = file_handle.createVariable(
                "longitude",
                "f8",
                ("point",),
                fill_value=np.nan
            )

            longitude_var.units = "degrees_east"
            longitude_var.standard_name = "longitude"
            longitude_var.long_name = "point longitude"
            longitude_var.axis = "X"

            # -------------------------------------------------------------------------
            # Latitude
            latitude_var = file_handle.createVariable(
                "latitude",
                "f8",
                ("point",),
                fill_value=np.nan
            )

            latitude_var.units = "degrees_north"
            latitude_var.standard_name = "latitude"
            latitude_var.long_name = "point latitude"
            latitude_var.axis = "Y"

            # -------------------------------------------------------------------------
            # Global attributes
            file_handle.title = "Point time-series collection"
            file_handle.featureType = "timeSeries"
            file_handle.Conventions = "CF-1.8"

            file_handle.time_coverage_start = (
                f"{point_time_index.min():%Y-%m-%d %H:%M:%S}"
            )

            file_handle.time_coverage_end = (
                f"{point_time_index.max():%Y-%m-%d %H:%M:%S}"
            )

            file_handle.number_of_variables = len(
                variable_names
            )

    # -------------------------------------------------------------------------------------
    # Append or update a point
    with Dataset(
            file_name,
            mode="a") as file_handle:

        # ---------------------------------------------------------------------------------
        # Validate required dimensions
        required_dimensions = [
            "time",
            "variable",
            "point"
        ]

        for dimension_name in required_dimensions:

            if dimension_name not in file_handle.dimensions:
                raise RuntimeError(
                    f'Dimension "{dimension_name}" is missing '
                    f'from NetCDF file "{file_name}"'
                )

        # ---------------------------------------------------------------------------------
        # Validate required variables
        required_variables = [
            "time",
            "variable_index",
            "variable_name",
            "point_id",
            "point_name",
            "longitude",
            "latitude"
        ]

        for required_variable in required_variables:

            if required_variable not in file_handle.variables:
                raise RuntimeError(
                    f'Variable "{required_variable}" is missing '
                    f'from NetCDF file "{file_name}"'
                )

        # ---------------------------------------------------------------------------------
        # Validate variable dimension
        file_variable_size = len(
            file_handle.dimensions["variable"]
        )

        if file_variable_size != len(variable_names):
            raise RuntimeError(
                f'Variable dimension mismatch for point "{point_tag}": '
                f'file has {file_variable_size} variables, '
                f'point has {len(variable_names)} variables'
            )

        # ---------------------------------------------------------------------------------
        # Validate stored variable indexes
        stored_variable_indices = np.asarray(
            file_handle.variables["variable_index"][:],
            dtype=np.int32
        )

        if not np.array_equal(
                stored_variable_indices,
                variable_indices):

            raise RuntimeError(
                f'Variable-index mismatch for point "{point_tag}": '
                f'stored={stored_variable_indices.tolist()}, '
                f'current={variable_indices.tolist()}'
            )

        # ---------------------------------------------------------------------------------
        # Validate stored variable names
        stored_variable_names = [
            _convert_nc_string(value)
            for value in file_handle.variables["variable_name"][:]
        ]

        if stored_variable_names != variable_names:
            raise RuntimeError(
                f'Variable-name mismatch for point "{point_tag}": '
                f'stored={stored_variable_names}, '
                f'current={variable_names}'
            )

        # ---------------------------------------------------------------------------------
        # Read saved time dimension
        time_var = file_handle.variables[
            "time"
        ]

        file_time_values = np.asarray(
            time_var[:],
            dtype=np.float64
        )

        file_time_units = getattr(
            time_var,
            "units",
            None
        )

        file_time_calendar = getattr(
            time_var,
            "calendar",
            "standard"
        )

        if file_time_units is None:
            raise RuntimeError(
                f'Time variable in "{file_name}" '
                f'does not define the "units" attribute'
            )

        file_time_values_datetime = num2date(
            file_time_values,
            units=file_time_units,
            calendar=file_time_calendar,
            only_use_cftime_datetimes=False,
            only_use_python_datetimes=True
        )

        file_time_index = pd.DatetimeIndex(
            pd.to_datetime(
                file_time_values_datetime
            )
        )

        if file_time_index.has_duplicates:
            raise RuntimeError(
                f'The saved NetCDF time dimension in "{file_name}" '
                f'contains duplicated timestamps'
            )

        file_n_steps = file_time_index.size
        point_n_steps = point_time_index.size

        # ---------------------------------------------------------------------------------
        # Compute time mismatch
        missing_time_index = file_time_index.difference(
            point_time_index
        )

        extra_time_index = point_time_index.difference(
            file_time_index
        )

        missing_n_steps = missing_time_index.size
        extra_n_steps = extra_time_index.size

        difference_n_steps = (
            missing_n_steps
            + extra_n_steps
        )

        if file_n_steps > 0:

            difference_percentage = (
                difference_n_steps
                / file_n_steps
                * 100.0
            )

        else:

            difference_percentage = 100.0

        # ---------------------------------------------------------------------------------
        # Warn about time mismatch
        if missing_n_steps > 0 or extra_n_steps > 0:

            log_stream.warning(
                ' ===> Time dimension mismatch for point "%s": '
                'file has %d step(s), point has %d step(s); '
                'missing=%d, extra=%d',
                point_tag,
                file_n_steps,
                point_n_steps,
                missing_n_steps,
                extra_n_steps
            )

            log_stream.warning(
                ' ===> Point "%s" will be aligned according to the '
                'saved NetCDF time dimension',
                point_tag
            )

            if difference_percentage > mismatch_threshold_percentage:

                log_stream.warning(
                    ' ===> Large time mismatch for point "%s": '
                    '%d different step(s), equal to %.2f%% of the '
                    'saved NetCDF time dimension. Threshold is %.2f%%',
                    point_tag,
                    difference_n_steps,
                    difference_percentage,
                    mismatch_threshold_percentage
                )

                log_stream.warning(
                    ' ===> Check input data for point "%s" and cut or fill '
                    'them according to the saved NetCDF dimensions',
                    point_tag
                )

            if extra_n_steps > 0:

                log_stream.warning(
                    ' ===> Cutting %d extra time step(s) for point "%s": '
                    '%s --> %s',
                    extra_n_steps,
                    point_tag,
                    extra_time_index.min().strftime(
                        "%Y-%m-%d %H:%M"
                    ),
                    extra_time_index.max().strftime(
                        "%Y-%m-%d %H:%M"
                    )
                )

            if missing_n_steps > 0:

                log_stream.warning(
                    ' ===> Filling %d missing time step(s) for point "%s" '
                    'using no-data value %s: %s --> %s',
                    missing_n_steps,
                    point_tag,
                    str(file_no_data),
                    missing_time_index.min().strftime(
                        "%Y-%m-%d %H:%M"
                    ),
                    missing_time_index.max().strftime(
                        "%Y-%m-%d %H:%M"
                    )
                )

        # ---------------------------------------------------------------------------------
        # Align point data to the saved time dimension
        point_dframe = point_dframe.reindex(
            file_time_index
        )

        # Convert DataFrame into matrix: (n_steps, n_variables)
        point_values = point_dframe[
            variable_names
        ].to_numpy(
            dtype=np.dtype(file_dtype),
            copy=True
        )

        point_values[
            ~np.isfinite(point_values)
        ] = file_no_data

        expected_shape = (
            file_n_steps,
            len(variable_names)
        )

        if point_values.shape != expected_shape:
            raise RuntimeError(
                f'Unexpected data shape for point "{point_tag}": '
                f'{point_values.shape}; expected {expected_shape}'
            )

        # ---------------------------------------------------------------------------------
        # Validate metadata dimensions
        metadata_dimensions = {
            "point_id": ("point",),
            "point_name": ("point",),
            "longitude": ("point",),
            "latitude": ("point",)
        }

        for variable_name, expected_dimensions in metadata_dimensions.items():

            current_dimensions = file_handle.variables[
                variable_name
            ].dimensions

            if current_dimensions != expected_dimensions:
                raise RuntimeError(
                    f'Invalid dimensions for variable "{variable_name}" '
                    f'in "{file_name}": {current_dimensions}; '
                    f'expected {expected_dimensions}'
                )

        # ---------------------------------------------------------------------------------
        # Find existing point
        existing_point_ids = [
            _convert_nc_string(value)
            for value in file_handle.variables["point_id"][:]
        ]

        if point_tag in existing_point_ids:

            point_index = existing_point_ids.index(
                point_tag
            )

            if not overwrite_point:
                raise RuntimeError(
                    f'Point "{point_tag}" already exists '
                    f'in NetCDF file "{file_name}"'
                )

            log_stream.warning(
                ' ===> Point "%s" already exists at index %d. '
                'Existing values will be overwritten',
                point_tag,
                point_index
            )

        else:

            point_index = len(
                existing_point_ids
            )

        # ---------------------------------------------------------------------------------
        # Write point metadata
        file_handle.variables[
            "point_id"
        ][point_index] = point_tag

        file_handle.variables[
            "point_name"
        ][point_index] = point_name

        file_handle.variables[
            "longitude"
        ][point_index] = point_lon

        file_handle.variables[
            "latitude"
        ][point_index] = point_lat

        # ---------------------------------------------------------------------------------
        # Create or validate point variable
        if point_tag not in file_handle.variables:

            point_var = file_handle.createVariable(
                point_tag,
                file_dtype,
                ("time", "variable"),
                fill_value=file_no_data,
                zlib=bool(compression),
                complevel=(
                    compression_level
                    if compression
                    else 0
                ),
                shuffle=bool(compression)
            )

            point_var.long_name = (
                f"time-series data for point {point_name}"
            )

            point_var.point_id = point_tag
            point_var.point_name = point_name
            point_var.point_index = point_index
            point_var.longitude = point_lon
            point_var.latitude = point_lat

            point_var.coordinates = (
                "time variable_index variable_name"
            )

        else:

            point_var = file_handle.variables[
                point_tag
            ]

            expected_dimensions = (
                "time",
                "variable"
            )

            if point_var.dimensions != expected_dimensions:
                raise RuntimeError(
                    f'Invalid dimensions for point variable "{point_tag}": '
                    f'{point_var.dimensions}; expected '
                    f'{expected_dimensions}'
                )

            if point_var.shape != expected_shape:
                raise RuntimeError(
                    f'Invalid shape for point variable "{point_tag}": '
                    f'{point_var.shape}; expected {expected_shape}'
                )

            point_var.point_id = point_tag
            point_var.point_name = point_name
            point_var.point_index = point_index
            point_var.longitude = point_lon
            point_var.latitude = point_lat

        # ---------------------------------------------------------------------------------
        # Write point matrix
        point_var[:, :] = point_values

        # ---------------------------------------------------------------------------------
        # Update global attributes
        file_handle.number_of_points = len(
            file_handle.dimensions["point"]
        )

        file_handle.number_of_variables = len(
            file_handle.dimensions["variable"]
        )

    log_stream.info(
        ' ---> Point "%s" saved in NetCDF file "%s" at index %d '
        'with array shape (%d, %d)',
        point_tag,
        file_name,
        point_index,
        point_values.shape[0],
        point_values.shape[1]
    )

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read datasets netcdf
def read_datasets_nc(
        file_name,
        file_fields=None,
        registry_fields=None,
        time_reference=None,
        time_start=None,
        time_end=None,
        time_rounding='H',
        time_frequency='H',
        ascending_index=False,
        sort_index=True,
        time_name='time',
        series_dimension_name='series',
        series_index_name='series_index',
        series_name_name='series_name',
        series_field_name='series_field',
        masked_to_nan=True,
        **kwargs) -> Dict[str, pd.DataFrame]:
    """
    Read all point time-series from a NetCDF collection.

    The NetCDF file is opened only once and all points are returned as a
    dictionary of pandas DataFrames.

    Expected NetCDF structure
    -------------------------
    Dimensions:
        time
        series

    Coordinate variable:
        time(time)

    Metadata variables:
        series_index(series)
        series_name(series)
        series_field(series), optional

    Point variables:
        point_1(time, series)
        point_2(time, series)
        ...

    Returned object
    ---------------
    {
        'point_1': DataFrame,
        'point_2': DataFrame,
        ...
    }

    Each DataFrame has the structure:

                         rain    airt      sm
        time
        2026-01-01 00:00      0.0    12.4    45.2
        2026-01-01 01:00      1.2    11.8    44.9

    Point selection is managed externally:

        datasets_obj = read_datasets_nc(...)
        point_dframe = datasets_obj['point_1']
    """

    # -------------------------------------------------------------------------
    # Initialize optional arguments
    if file_fields is None:
        file_fields = {}

    if registry_fields is None:
        registry_fields = {}

    # -------------------------------------------------------------------------
    # Check file path
    file_path_obj = Path(file_name)

    if not file_path_obj.exists():
        raise FileNotFoundError(
            f'NetCDF file not found: "{file_path_obj}"'
        )

    if not file_path_obj.is_file():
        raise FileNotFoundError(
            f'NetCDF path is not a file: "{file_path_obj}"'
        )

    # -------------------------------------------------------------------------
    # Initialize output dictionary
    points_obj: Dict[str, pd.DataFrame] = {}

    # -------------------------------------------------------------------------
    # Open NetCDF file only once
    try:

        with Dataset(str(file_path_obj), mode='r') as file_nc:

            # -----------------------------------------------------------------
            # Check required dimensions
            required_dimensions = [
                time_name,
                series_dimension_name,
            ]

            missing_dimensions = [
                dimension_name
                for dimension_name in required_dimensions
                if dimension_name not in file_nc.dimensions
            ]

            if missing_dimensions:
                raise RuntimeError(
                    f'Missing dimensions in "{file_path_obj}": '
                    f'{missing_dimensions}'
                )

            number_of_time_steps = len(
                file_nc.dimensions[time_name]
            )

            number_of_series = len(
                file_nc.dimensions[series_dimension_name]
            )

            if number_of_time_steps == 0:
                raise RuntimeError(
                    f'Dimension "{time_name}" is empty '
                    f'in "{file_path_obj}"'
                )

            if number_of_series == 0:
                raise RuntimeError(
                    f'Dimension "{series_dimension_name}" is empty '
                    f'in "{file_path_obj}"'
                )

            # -----------------------------------------------------------------
            # Check required variables
            required_variables = [
                time_name,
                series_index_name,
                series_name_name,
            ]

            missing_variables = [
                variable_name
                for variable_name in required_variables
                if variable_name not in file_nc.variables
            ]

            if missing_variables:
                raise RuntimeError(
                    f'Missing variables in "{file_path_obj}": '
                    f'{missing_variables}'
                )

            # -----------------------------------------------------------------
            # Read and decode time coordinate
            time_variable = file_nc.variables[time_name]

            if not hasattr(time_variable, 'units'):
                raise RuntimeError(
                    f'Time variable "{time_name}" does not have '
                    f'a "units" attribute in "{file_path_obj}"'
                )

            time_calendar = getattr(
                time_variable,
                'calendar',
                'standard'
            )

            time_values_raw = time_variable[:]

            try:
                time_values = num2date(
                    time_values_raw,
                    units=time_variable.units,
                    calendar=time_calendar,
                    only_use_cftime_datetimes=False,
                    only_use_python_datetimes=True,
                )

                time_values = pd.to_datetime(
                    np.asarray(time_values)
                )

            except Exception:

                # Fallback for calendars that cannot be directly converted
                time_values_cftime = num2date(
                    time_values_raw,
                    units=time_variable.units,
                    calendar=time_calendar,
                    only_use_cftime_datetimes=True,
                )

                time_values = pd.to_datetime([
                    value.strftime('%Y-%m-%d %H:%M:%S')
                    for value in time_values_cftime
                ])

            time_index = pd.DatetimeIndex(
                time_values,
                name='time'
            )

            if len(time_index) != number_of_time_steps:
                raise RuntimeError(
                    f'Time variable "{time_name}" contains '
                    f'{len(time_index)} values, but dimension '
                    f'"{time_name}" has length '
                    f'{number_of_time_steps}'
                )

            # -----------------------------------------------------------------
            # Read series indices
            series_indices = np.asarray(
                file_nc.variables[series_index_name][:],
                dtype=np.int32,
            ).reshape(-1)

            # -----------------------------------------------------------------
            # Read output series names
            series_names = _decode_netcdf_strings(
                file_nc.variables[series_name_name][:]
            )

            # -----------------------------------------------------------------
            # Read optional source fields
            if series_field_name in file_nc.variables:

                source_fields = _decode_netcdf_strings(
                    file_nc.variables[series_field_name][:]
                )

            else:

                source_fields = [None] * number_of_series

            # -----------------------------------------------------------------
            # Validate metadata lengths
            if len(series_indices) != number_of_series:
                raise RuntimeError(
                    f'Variable "{series_index_name}" contains '
                    f'{len(series_indices)} values, but dimension '
                    f'"{series_dimension_name}" has length '
                    f'{number_of_series}'
                )

            if len(series_names) != number_of_series:
                raise RuntimeError(
                    f'Variable "{series_name_name}" contains '
                    f'{len(series_names)} values, but dimension '
                    f'"{series_dimension_name}" has length '
                    f'{number_of_series}'
                )

            if len(source_fields) != number_of_series:
                raise RuntimeError(
                    f'Variable "{series_field_name}" contains '
                    f'{len(source_fields)} values, but dimension '
                    f'"{series_dimension_name}" has length '
                    f'{number_of_series}'
                )

            # -----------------------------------------------------------------
            # Validate series indices
            if np.any(series_indices < 0):
                raise RuntimeError(
                    f'Negative series indices found in '
                    f'"{series_index_name}": '
                    f'{series_indices.tolist()}'
                )

            if np.any(series_indices >= number_of_series):
                raise RuntimeError(
                    f'Series indices exceed dimension length: '
                    f'{series_indices.tolist()}. '
                    f'Maximum allowed index is '
                    f'{number_of_series - 1}'
                )

            if len(np.unique(series_indices)) != len(series_indices):
                raise RuntimeError(
                    f'Duplicated series indices found: '
                    f'{series_indices.tolist()}'
                )

            # -----------------------------------------------------------------
            # Validate series names
            series_names = [
                str(series_name).strip()
                for series_name in series_names
            ]

            empty_series_names = [
                series_position
                for series_position, series_name
                in enumerate(series_names)
                if not series_name
            ]

            if empty_series_names:
                raise RuntimeError(
                    f'Empty series names found at positions: '
                    f'{empty_series_names}'
                )

            if len(set(series_names)) != len(series_names):
                raise RuntimeError(
                    f'Duplicated series names found: '
                    f'{series_names}'
                )

            # -----------------------------------------------------------------
            # Optionally remap output series names using file_fields
            #
            # Example:
            #
            # file_fields = {
            #     'time': 'time',
            #     'values_1': 'rain',
            #     'values_2': 'airt',
            # }
            #
            # Stored NetCDF series name "rain" becomes "values_1".
            # -----------------------------------------------------------------
            series_output_names = []

            for series_name in series_names:

                output_name = series_name

                for field_name, file_field in file_fields.items():

                    if field_name == 'time':
                        continue

                    if file_field == series_name:
                        output_name = field_name
                        break

                series_output_names.append(output_name)

            if len(set(series_output_names)) != len(series_output_names):
                raise RuntimeError(
                    f'Duplicated output series names after applying '
                    f'"file_fields": {series_output_names}'
                )

            # -----------------------------------------------------------------
            # Define expected point variable structure
            excluded_variables = {
                time_name,
                series_index_name,
                series_name_name,
                series_field_name,
            }

            expected_dimensions = (
                time_name,
                series_dimension_name,
            )

            expected_shape = (
                number_of_time_steps,
                number_of_series,
            )

            # -----------------------------------------------------------------
            # Define optional requested time range
            if (time_start is not None) and (time_end is not None):

                time_start_expected = pd.Timestamp(
                    time_start
                ).floor(time_rounding.lower())

                time_end_expected = pd.Timestamp(
                    time_end
                ).floor(time_rounding.lower())

                if time_start_expected > time_end_expected:
                    raise RuntimeError(
                        f'Invalid time range: start '
                        f'"{time_start_expected}" is after end '
                        f'"{time_end_expected}"'
                    )

                time_range = pd.date_range(
                    start=time_start_expected,
                    end=time_end_expected,
                    freq=time_frequency.lower(),
                    name='time',
                )

            else:

                time_range = None

            # -----------------------------------------------------------------
            # Read all point variables
            for variable_name, variable_obj in file_nc.variables.items():

                # Skip coordinate and metadata variables
                if variable_name in excluded_variables:
                    continue

                # Point variables must have dimensions:
                # point_name(time, series)
                if variable_obj.dimensions != expected_dimensions:
                    continue

                # -------------------------------------------------------------
                # Read point values
                point_values = variable_obj[:]

                # Convert masked values to NaN
                if np.ma.isMaskedArray(point_values):

                    if masked_to_nan:
                        point_values = point_values.filled(np.nan)
                    else:
                        point_values = point_values.data

                # Make an independent and writable array
                point_values = np.array(
                    point_values,
                    dtype=np.float32,
                    copy=True,
                )

                # -------------------------------------------------------------
                # Validate point shape
                if point_values.shape != expected_shape:
                    raise RuntimeError(
                        f'Unexpected shape for point variable '
                        f'"{variable_name}": {point_values.shape}; '
                        f'expected {expected_shape}'
                    )

                # -------------------------------------------------------------
                # Convert explicit fill and missing values to NaN
                if masked_to_nan:

                    fill_values = []

                    if hasattr(variable_obj, '_FillValue'):
                        fill_values.append(
                            getattr(variable_obj, '_FillValue')
                        )

                    if hasattr(variable_obj, 'missing_value'):

                        missing_values = np.asarray(
                            getattr(variable_obj, 'missing_value')
                        ).reshape(-1)

                        fill_values.extend(
                            missing_values.tolist()
                        )

                    for fill_value in fill_values:

                        try:
                            fill_value_float = float(fill_value)
                        except (TypeError, ValueError):
                            continue

                        if np.isnan(fill_value_float):
                            continue

                        point_values[
                            np.isclose(
                                point_values,
                                fill_value_float,
                                equal_nan=False,
                            )
                        ] = np.nan

                # -------------------------------------------------------------
                # Organize the current point as a DataFrame
                point_data = {}

                for series_index, series_output_name in zip(
                        series_indices,
                        series_output_names):

                    point_data[series_output_name] = point_values[
                        :,
                        int(series_index),
                    ].copy()

                point_dframe = pd.DataFrame(
                    data=point_data,
                    index=time_index,
                )

                point_dframe.index.name = 'time'

                # -------------------------------------------------------------
                # Remove duplicated times, if present
                if point_dframe.index.has_duplicates:

                    log_stream.warning(
                        f' ===> Duplicated times found for point '
                        f'"{variable_name}" in "{file_path_obj}". '
                        f'Keeping the last occurrence'
                    )

                    point_dframe = point_dframe[
                        ~point_dframe.index.duplicated(
                            keep='last'
                        )
                    ]

                # -------------------------------------------------------------
                # Reindex over requested time interval
                if time_range is not None:
                    point_dframe = point_dframe.reindex(time_range)

                # -------------------------------------------------------------
                # Sort time index
                if sort_index:
                    point_dframe = point_dframe.sort_index(
                        ascending=ascending_index
                    )

                # -------------------------------------------------------------
                # Store metadata
                point_attrs = dict(registry_fields)

                point_attrs['time_reference'] = time_reference
                point_attrs['point_name'] = variable_name
                point_attrs['file_name'] = str(file_path_obj)
                point_attrs['series_names'] = series_output_names
                point_attrs['source_fields'] = source_fields

                # Add NetCDF variable attributes
                for attribute_name in variable_obj.ncattrs():

                    try:
                        point_attrs[attribute_name] = variable_obj.getncattr(
                            attribute_name
                        )
                    except Exception:
                        pass

                point_dframe.attrs = point_attrs

                # -------------------------------------------------------------
                # Store current point
                points_obj[variable_name] = point_dframe

    except OSError as exc:

        log_stream.error(
            f' ===> Error reading NetCDF file '
            f'"{file_path_obj}": {exc}'
        )

        raise RuntimeError(
            f'Unable to read NetCDF file "{file_path_obj}"'
        ) from exc

    # -------------------------------------------------------------------------
    # Check that at least one point was found
    if not points_obj:
        raise RuntimeError(
            f'No point variables with dimensions '
            f'{expected_dimensions} were found '
            f'in "{file_path_obj}"'
        )

    return points_obj


def _decode_netcdf_strings(values: Any) -> list[str]:
    """
    Decode NetCDF string variables into a list of Python strings.

    Supported formats:

    - variable-length strings
    - byte strings
    - fixed-width character arrays
    - NumPy object arrays
    """

    values_array = np.asarray(values)

    # -------------------------------------------------------------------------
    # Fixed-width character array
    #
    # Example:
    #
    # [
    #     [b'r', b'a', b'i', b'n'],
    #     [b'a', b'i', b'r', b't'],
    # ]
    # -------------------------------------------------------------------------
    if (
            values_array.ndim == 2
            and values_array.dtype.kind in {'S', 'U'}
    ):

        decoded_values = []

        for row in values_array:

            characters = []

            for value in row:

                if isinstance(value, bytes):

                    character = value.decode(
                        'utf-8',
                        errors='replace',
                    )

                elif isinstance(value, np.bytes_):

                    character = value.tobytes().decode(
                        'utf-8',
                        errors='replace',
                    )

                else:

                    character = str(value)

                characters.append(character)

            decoded_value = ''.join(characters)
            decoded_value = decoded_value.rstrip('\x00').strip()

            decoded_values.append(decoded_value)

        return decoded_values

    # -------------------------------------------------------------------------
    # One-dimensional string, byte or object array
    decoded_values = []

    for value in values_array.reshape(-1):

        if isinstance(value, bytes):

            decoded_value = value.decode(
                'utf-8',
                errors='replace',
            )

        elif isinstance(value, np.bytes_):

            decoded_value = value.tobytes().decode(
                'utf-8',
                errors='replace',
            )

        else:

            decoded_value = str(value)

        decoded_value = decoded_value.rstrip('\x00').strip()

        decoded_values.append(decoded_value)

    return decoded_values
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to check point in the file
def check_datasets_nc(file_name,
                      current_point,
                      first_point_expected=None,
                      last_point_expected=None):
    """Check if the required point variables exist in the NetCDF file."""

    if not os.path.exists(file_name):
        return (
            False,
            None if first_point_expected is None else False,
            None if last_point_expected is None else False,
        )

    with Dataset(file_name, mode="r") as file_handle:
        variables = file_handle.variables

        current_exists = current_point in variables

        first_exists = (
            None
            if first_point_expected is None
            else first_point_expected in variables
        )

        last_exists = (
            None
            if last_point_expected is None
            else last_point_expected in variables
        )

    return current_exists, first_exists, last_exists
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to write datasets in netcdf file (append mode for multiple points)
# ----------------------------------------------------------------------------------------------------------------------
# method to write datasets in NetCDF format, organized by point
def write_datasets_nc(
        file_path: str,
        point_name: str,
        dframe_combined: pd.DataFrame,
        file_fields: Mapping[str, str],
        point_id: str = None,
        point_long_name: str = None,
        point_lon: float = None,
        point_lat: float = None,
        time_name: str = "time",
        fill_value: float = -9999.0,
        dtype: str = "f4",
        compression_level: int = 4,
        overwrite_point: bool = False,
        mismatch_threshold_percentage: float = 50.0,
        **kwargs) -> None:
    """
    Save one point time series into a point-oriented NetCDF collection.

    NetCDF structure
    ----------------
    dimensions:
        time
        variable
        point

    shared variables:
        time(time)
        variable_index(variable)
        variable_name(variable)

        point_id(point)
        point_name(point)
        longitude(point)
        latitude(point)

    point variables:
        point_1(time, variable)
        point_2(time, variable)
        ...

    Field mapping
    -------------
    Example:

        file_fields = {
            "time": "time",
            "values_k1": "rain",
            "values_k2": "air_temperature",
            "values_k3": "soil_moisture"
        }

    The input DataFrame may contain either the source fields:

        values_k1
        values_k2
        values_k3

    or the final fields:

        rain
        air_temperature
        soil_moisture

    In both cases, the NetCDF variable names will be:

        variable_name = [
            "rain",
            "air_temperature",
            "soil_moisture"
        ]
    """

    # -------------------------------------------------------------------------------------
    # Validate dataframe
    if not isinstance(dframe_combined, pd.DataFrame):
        raise TypeError(
            '"dframe_combined" must be a pandas DataFrame'
        )

    if dframe_combined.empty:
        raise ValueError(
            f'DataFrame is empty for point "{point_name}"'
        )

    # -------------------------------------------------------------------------------------
    # Validate point variable name
    if not isinstance(point_name, str):
        raise TypeError(
            '"point_name" must be a string'
        )

    point_name = point_name.strip()

    if not point_name:
        raise ValueError(
            '"point_name" must be a non-empty string'
        )

    # -------------------------------------------------------------------------------------
    # Define point identifier
    if point_id is None:

        point_id = point_name

    else:

        point_id = str(
            point_id
        ).strip()

    if not point_id:
        raise ValueError(
            '"point_id" must be a non-empty value'
        )

    # -------------------------------------------------------------------------------------
    # Define descriptive point name
    if point_long_name is None:

        point_long_name = point_name

    else:

        point_long_name = str(
            point_long_name
        )

    # -------------------------------------------------------------------------------------
    # Validate time name
    if not isinstance(time_name, str):
        raise TypeError(
            '"time_name" must be a string'
        )

    time_name = time_name.strip()

    if not time_name:
        raise ValueError(
            '"time_name" must be a non-empty string'
        )

    # -------------------------------------------------------------------------------------
    # Validate field mapping
    if not isinstance(file_fields, Mapping):
        raise TypeError(
            '"file_fields" must be a mapping'
        )

    if not file_fields:
        raise ValueError(
            '"file_fields" must be a non-empty mapping'
        )

    # -------------------------------------------------------------------------------------
    # Validate compression
    if not isinstance(compression_level, int):
        raise TypeError(
            '"compression_level" must be an integer'
        )

    if compression_level < 0 or compression_level > 9:
        raise ValueError(
            '"compression_level" must be between 0 and 9'
        )

    # -------------------------------------------------------------------------------------
    # Validate mismatch threshold
    if mismatch_threshold_percentage is None:
        mismatch_threshold_percentage = 50.0

    try:

        mismatch_threshold_percentage = float(
            mismatch_threshold_percentage
        )

    except (TypeError, ValueError) as exc:

        raise ValueError(
            '"mismatch_threshold_percentage" must be numeric'
        ) from exc

    if not 0.0 <= mismatch_threshold_percentage <= 100.0:
        raise ValueError(
            '"mismatch_threshold_percentage" must be between 0 and 100'
        )

    # -------------------------------------------------------------------------------------
    # Validate fill value
    try:

        fill_value = float(
            fill_value
        )

    except (TypeError, ValueError) as exc:

        raise ValueError(
            '"fill_value" must be numeric'
        ) from exc

    # -------------------------------------------------------------------------------------
    # Validate output dtype
    try:

        output_dtype = np.dtype(
            dtype
        )

    except TypeError as exc:

        raise ValueError(
            f'Invalid NetCDF data type: "{dtype}"'
        ) from exc

    if not np.issubdtype(output_dtype, np.number):
        raise ValueError(
            f'NetCDF point data type must be numeric, not "{dtype}"'
        )

    # -------------------------------------------------------------------------------------
    # Define source and final variable names
    field_pairs = [
        (
            str(source_name),
            str(destination_name)
        )
        for source_name, destination_name in file_fields.items()
        if (
            str(source_name) != time_name
            and str(destination_name) != time_name
        )
    ]

    if not field_pairs:
        raise RuntimeError(
            f'No variables configured for point "{point_name}"'
        )

    source_variable_names = [
        source_name
        for source_name, destination_name in field_pairs
    ]

    final_variable_names = [
        destination_name
        for source_name, destination_name in field_pairs
    ]

    # -------------------------------------------------------------------------------------
    # Validate source and final names
    if len(source_variable_names) != len(set(source_variable_names)):
        raise ValueError(
            f'Duplicated source variable names: '
            f'{source_variable_names}'
        )

    if len(final_variable_names) != len(set(final_variable_names)):
        raise ValueError(
            f'Duplicated final variable names: '
            f'{final_variable_names}'
        )

    reserved_variable_names = {
        time_name,
        "variable_index",
        "variable_name",
        "point_id",
        "point_name",
        "longitude",
        "latitude"
    }

    invalid_final_names = [
        variable_name
        for variable_name in final_variable_names
        if variable_name in reserved_variable_names
    ]

    if invalid_final_names:
        raise ValueError(
            f'Final variable names use reserved NetCDF names: '
            f'{invalid_final_names}'
        )

    # -------------------------------------------------------------------------------------
    # Prepare point dataframe
    point_dframe = dframe_combined.copy()

    # -------------------------------------------------------------------------------------
    # Extract time from column or index
    if time_name in point_dframe.columns:

        point_time_index = pd.DatetimeIndex(
            pd.to_datetime(
                point_dframe[time_name],
                errors="coerce"
            )
        )

        point_dframe = point_dframe.drop(
            columns=[time_name]
        )

        point_dframe.index = point_time_index

    elif isinstance(point_dframe.index, pd.DatetimeIndex):

        point_dframe.index = pd.DatetimeIndex(
            point_dframe.index
        )

    else:

        converted_time_index = pd.to_datetime(
            point_dframe.index,
            errors="coerce"
        )

        if pd.isna(converted_time_index).all():
            raise TypeError(
                f'Time information not found for point "{point_name}". '
                f'Expected column "{time_name}" or a DatetimeIndex'
            )

        point_dframe.index = pd.DatetimeIndex(
            converted_time_index
        )

    # -------------------------------------------------------------------------------------
    # Remove invalid timestamps
    invalid_time_mask = point_dframe.index.isna()

    if invalid_time_mask.any():

        invalid_time_steps = int(
            invalid_time_mask.sum()
        )

        log_stream.warning(
            ' ===> Removing %d invalid timestamp(s) for point "%s"',
            invalid_time_steps,
            point_name
        )

        point_dframe = point_dframe.loc[
            ~invalid_time_mask
        ]

    if point_dframe.empty:
        raise RuntimeError(
            f'No valid timestamps available for point "{point_name}"'
        )

    # -------------------------------------------------------------------------------------
    # Remove duplicated timestamps
    if point_dframe.index.has_duplicates:

        duplicated_time_mask = point_dframe.index.duplicated(
            keep="last"
        )

        duplicated_time_steps = int(
            duplicated_time_mask.sum()
        )

        log_stream.warning(
            ' ===> Removing %d duplicated timestamp(s) for point "%s"; '
            'the last occurrence is retained',
            duplicated_time_steps,
            point_name
        )

        point_dframe = point_dframe.loc[
            ~duplicated_time_mask
        ]

    # -------------------------------------------------------------------------------------
    # Sort time
    point_dframe = point_dframe.sort_index(
        ascending=True
    )

    point_time_index = pd.DatetimeIndex(
        point_dframe.index
    )

    # -------------------------------------------------------------------------------------
    # Detect whether dataframe contains source or final variable names
    source_fields_available = all(
        variable_name in point_dframe.columns
        for variable_name in source_variable_names
    )

    final_fields_available = all(
        variable_name in point_dframe.columns
        for variable_name in final_variable_names
    )

    # Prefer source names when both layouts are available
    if source_fields_available:

        dataframe_variable_names = list(
            source_variable_names
        )

        log_stream.info(
            ' ===> Point "%s" uses source dataframe fields: %s',
            point_name,
            dataframe_variable_names
        )

    elif final_fields_available:

        dataframe_variable_names = list(
            final_variable_names
        )

        log_stream.info(
            ' ===> Point "%s" already uses final dataframe fields: %s',
            point_name,
            dataframe_variable_names
        )

    else:

        missing_source_variables = [
            variable_name
            for variable_name in source_variable_names
            if variable_name not in point_dframe.columns
        ]

        missing_final_variables = [
            variable_name
            for variable_name in final_variable_names
            if variable_name not in point_dframe.columns
        ]

        raise KeyError(
            f'Unable to identify DataFrame variable layout for point '
            f'"{point_name}". '
            f'Missing source fields: {missing_source_variables}; '
            f'missing final fields: {missing_final_variables}. '
            f'Available fields are {list(point_dframe.columns)}'
        )

    # -------------------------------------------------------------------------------------
    # Select input fields in configured order
    point_dframe = point_dframe[
        dataframe_variable_names
    ].copy()

    # -------------------------------------------------------------------------------------
    # Convert all fields to numeric
    for dataframe_variable_name in dataframe_variable_names:

        point_dframe[
            dataframe_variable_name
        ] = pd.to_numeric(
            point_dframe[dataframe_variable_name],
            errors="coerce"
        )

    # -------------------------------------------------------------------------------------
    # Rename source fields to final NetCDF variable names
    rename_fields = {
        dataframe_variable_name: final_variable_name
        for dataframe_variable_name, final_variable_name
        in zip(
            dataframe_variable_names,
            final_variable_names
        )
        if dataframe_variable_name != final_variable_name
    }

    if rename_fields:

        point_dframe = point_dframe.rename(
            columns=rename_fields
        )

    # Ensure final order
    point_dframe = point_dframe[
        final_variable_names
    ].copy()

    # Final names used in NetCDF
    variable_names = list(
        final_variable_names
    )

    # -------------------------------------------------------------------------------------
    # Organize longitude
    if point_lon is None:

        point_lon = np.nan

    else:

        try:

            point_lon = float(
                point_lon
            )

        except (TypeError, ValueError):

            point_lon = np.nan

    if not np.isfinite(point_lon):
        point_lon = np.nan

    # -------------------------------------------------------------------------------------
    # Organize latitude
    if point_lat is None:

        point_lat = np.nan

    else:

        try:

            point_lat = float(
                point_lat
            )

        except (TypeError, ValueError):

            point_lat = np.nan

    if not np.isfinite(point_lat):
        point_lat = np.nan

    # -------------------------------------------------------------------------------------
    # Prepare output path
    file_path_obj = Path(
        file_path
    )

    file_path_obj.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    file_exists = file_path_obj.exists()

    variable_indices = np.arange(
        len(variable_names),
        dtype=np.int32
    )

    time_units = "hours since 1970-01-01 00:00:00"
    time_calendar = "standard"

    # -------------------------------------------------------------------------------------
    # Create new NetCDF collection
    if not file_exists:

        point_time_values = np.asarray(
            date2num(
                point_time_index.to_pydatetime(),
                units=time_units,
                calendar=time_calendar
            ),
            dtype=np.float64
        )

        with Dataset(
                str(file_path_obj),
                mode="w",
                format="NETCDF4") as file_nc:

            # -------------------------------------------------------------------------
            # Dimensions
            file_nc.createDimension(
                time_name,
                len(point_time_index)
            )

            file_nc.createDimension(
                "variable",
                len(variable_names)
            )

            file_nc.createDimension(
                "point",
                None
            )

            # -------------------------------------------------------------------------
            # Time coordinate
            time_variable = file_nc.createVariable(
                time_name,
                "f8",
                (time_name,)
            )

            time_variable[:] = point_time_values

            time_variable.units = time_units
            time_variable.calendar = time_calendar
            time_variable.standard_name = "time"
            time_variable.long_name = "reference time"
            time_variable.axis = "T"

            # -------------------------------------------------------------------------
            # Variable index
            variable_index_obj = file_nc.createVariable(
                "variable_index",
                "i4",
                ("variable",)
            )

            variable_index_obj[:] = variable_indices

            variable_index_obj.long_name = (
                "zero-based variable position in each point matrix"
            )

            # -------------------------------------------------------------------------
            # Variable names
            variable_name_obj = file_nc.createVariable(
                "variable_name",
                str,
                ("variable",)
            )

            variable_name_obj[:] = np.asarray(
                variable_names,
                dtype=object
            )

            variable_name_obj.long_name = (
                "variable name associated with each point matrix column"
            )

            # -------------------------------------------------------------------------
            # Point identifiers
            point_id_obj = file_nc.createVariable(
                "point_id",
                str,
                ("point",)
            )

            point_id_obj.long_name = "point identifier"
            point_id_obj.cf_role = "timeseries_id"

            # -------------------------------------------------------------------------
            # Point descriptive names
            point_name_obj = file_nc.createVariable(
                "point_name",
                str,
                ("point",)
            )

            point_name_obj.long_name = "point descriptive name"

            # -------------------------------------------------------------------------
            # Longitude
            longitude_obj = file_nc.createVariable(
                "longitude",
                "f8",
                ("point",),
                fill_value=np.nan
            )

            longitude_obj.units = "degrees_east"
            longitude_obj.standard_name = "longitude"
            longitude_obj.long_name = "point longitude"
            longitude_obj.axis = "X"

            # -------------------------------------------------------------------------
            # Latitude
            latitude_obj = file_nc.createVariable(
                "latitude",
                "f8",
                ("point",),
                fill_value=np.nan
            )

            latitude_obj.units = "degrees_north"
            latitude_obj.standard_name = "latitude"
            latitude_obj.long_name = "point latitude"
            latitude_obj.axis = "Y"

            # -------------------------------------------------------------------------
            # Global attributes
            file_nc.title = "Point time-series dataset collection"
            file_nc.featureType = "timeSeries"
            file_nc.Conventions = "CF-1.8"

            file_nc.time_coverage_start = (
                f"{point_time_index.min():%Y-%m-%d %H:%M:%S}"
            )

            file_nc.time_coverage_end = (
                f"{point_time_index.max():%Y-%m-%d %H:%M:%S}"
            )

            file_nc.number_of_variables = len(
                variable_names
            )

            file_nc.number_of_points = 0

    # -------------------------------------------------------------------------------------
    # Append or update point
    with Dataset(
            str(file_path_obj),
            mode="a") as file_nc:

        # ---------------------------------------------------------------------------------
        # Validate required dimensions
        required_dimensions = [
            time_name,
            "variable",
            "point"
        ]

        missing_dimensions = [
            dimension_name
            for dimension_name in required_dimensions
            if dimension_name not in file_nc.dimensions
        ]

        if missing_dimensions:
            raise RuntimeError(
                f'Missing dimensions in "{file_path_obj}": '
                f'{missing_dimensions}. Available dimensions are '
                f'{list(file_nc.dimensions.keys())}'
            )

        # ---------------------------------------------------------------------------------
        # Validate required variables
        required_variables = [
            time_name,
            "variable_index",
            "variable_name",
            "point_id",
            "point_name",
            "longitude",
            "latitude"
        ]

        missing_variables = [
            variable_name
            for variable_name in required_variables
            if variable_name not in file_nc.variables
        ]

        if missing_variables:
            raise RuntimeError(
                f'Missing variables in "{file_path_obj}": '
                f'{missing_variables}. Available variables are '
                f'{list(file_nc.variables.keys())}'
            )

        # ---------------------------------------------------------------------------------
        # Validate dimension lengths
        stored_variable_count = len(
            file_nc.dimensions["variable"]
        )

        if stored_variable_count != len(variable_names):
            raise RuntimeError(
                f'Variable dimension mismatch for point "{point_name}": '
                f'file has {stored_variable_count} variable(s), '
                f'point has {len(variable_names)} variable(s)'
            )

        # ---------------------------------------------------------------------------------
        # Validate variable indices
        stored_variable_indices = np.asarray(
            file_nc.variables["variable_index"][:],
            dtype=np.int32
        ).reshape(-1)

        if not np.array_equal(
                stored_variable_indices,
                variable_indices):

            raise RuntimeError(
                f'Variable index mismatch for point "{point_name}": '
                f'stored={stored_variable_indices.tolist()}, '
                f'current={variable_indices.tolist()}'
            )

        # ---------------------------------------------------------------------------------
        # Validate variable names
        stored_variable_names = _decode_netcdf_strings(
            file_nc.variables["variable_name"][:]
        )

        stored_variable_names = [
            str(variable_name).strip()
            for variable_name in stored_variable_names
        ]

        if stored_variable_names != variable_names:
            raise RuntimeError(
                f'Variable name mismatch for point "{point_name}": '
                f'stored={stored_variable_names}, '
                f'current={variable_names}'
            )

        # ---------------------------------------------------------------------------------
        # Validate metadata dimensions
        metadata_dimensions = {
            "point_id": ("point",),
            "point_name": ("point",),
            "longitude": ("point",),
            "latitude": ("point",)
        }

        for metadata_name, expected_dimensions in metadata_dimensions.items():

            current_dimensions = file_nc.variables[
                metadata_name
            ].dimensions

            if current_dimensions != expected_dimensions:
                raise RuntimeError(
                    f'Invalid dimensions for metadata variable '
                    f'"{metadata_name}" in "{file_path_obj}": '
                    f'{current_dimensions}; expected '
                    f'{expected_dimensions}'
                )

        # ---------------------------------------------------------------------------------
        # Read saved time coordinate
        time_variable = file_nc.variables[
            time_name
        ]

        if not hasattr(time_variable, "units"):
            raise RuntimeError(
                f'Time variable "{time_name}" does not define '
                f'the "units" attribute in "{file_path_obj}"'
            )

        stored_time_values = num2date(
            time_variable[:],
            units=time_variable.units,
            calendar=getattr(
                time_variable,
                "calendar",
                "standard"
            ),
            only_use_cftime_datetimes=False,
            only_use_python_datetimes=True
        )

        stored_time_index = pd.DatetimeIndex(
            pd.to_datetime(
                np.asarray(stored_time_values)
            )
        )

        if stored_time_index.has_duplicates:
            raise RuntimeError(
                f'The saved time dimension in "{file_path_obj}" '
                f'contains duplicated timestamps'
            )

        stored_time_steps = len(
            stored_time_index
        )

        point_time_steps = len(
            point_time_index
        )

        # ---------------------------------------------------------------------------------
        # Compare point time with saved time
        missing_time_index = stored_time_index.difference(
            point_time_index
        )

        extra_time_index = point_time_index.difference(
            stored_time_index
        )

        missing_steps = len(
            missing_time_index
        )

        extra_steps = len(
            extra_time_index
        )

        different_steps = (
            missing_steps
            + extra_steps
        )

        if stored_time_steps > 0:

            difference_percentage = (
                different_steps
                / stored_time_steps
                * 100.0
            )

        else:

            difference_percentage = 100.0

        # ---------------------------------------------------------------------------------
        # Warn about time mismatch
        if missing_steps > 0 or extra_steps > 0:

            log_stream.warning(
                ' ===> Time dimension mismatch for point "%s": '
                'file has %d step(s), point has %d step(s); '
                'missing=%d, extra=%d',
                point_name,
                stored_time_steps,
                point_time_steps,
                missing_steps,
                extra_steps
            )

            log_stream.warning(
                ' ===> Point "%s" will be aligned according to the '
                'saved NetCDF time dimension',
                point_name
            )

            if difference_percentage > mismatch_threshold_percentage:

                log_stream.warning(
                    ' ===> Large time mismatch for point "%s": '
                    '%d different step(s), equal to %.2f%% of the '
                    'saved NetCDF time dimension. Threshold is %.2f%%',
                    point_name,
                    different_steps,
                    difference_percentage,
                    mismatch_threshold_percentage
                )

                log_stream.warning(
                    ' ===> Check input data for point "%s" and cut or fill '
                    'them according to the saved NetCDF dimensions',
                    point_name
                )

            if extra_steps > 0:

                log_stream.warning(
                    ' ===> Cutting %d extra time step(s) for point "%s": '
                    '%s --> %s',
                    extra_steps,
                    point_name,
                    extra_time_index.min().strftime(
                        "%Y-%m-%d %H:%M"
                    ),
                    extra_time_index.max().strftime(
                        "%Y-%m-%d %H:%M"
                    )
                )

            if missing_steps > 0:

                log_stream.warning(
                    ' ===> Filling %d missing time step(s) for point "%s" '
                    'using no-data value %s: %s --> %s',
                    missing_steps,
                    point_name,
                    str(fill_value),
                    missing_time_index.min().strftime(
                        "%Y-%m-%d %H:%M"
                    ),
                    missing_time_index.max().strftime(
                        "%Y-%m-%d %H:%M"
                    )
                )

        # ---------------------------------------------------------------------------------
        # Align point dataframe to saved time dimension
        point_dframe = point_dframe.reindex(
            stored_time_index
        )

        # ---------------------------------------------------------------------------------
        # Convert values to matrix
        point_data = point_dframe[
            variable_names
        ].to_numpy(
            dtype=output_dtype,
            copy=True
        )

        point_data[
            ~np.isfinite(point_data)
        ] = fill_value

        expected_shape = (
            stored_time_steps,
            len(variable_names)
        )

        if point_data.shape != expected_shape:
            raise RuntimeError(
                f'Unexpected shape for point "{point_name}": '
                f'{point_data.shape}; expected {expected_shape}'
            )

        # ---------------------------------------------------------------------------------
        # Read existing point identifiers
        existing_point_ids = _decode_netcdf_strings(
            file_nc.variables["point_id"][:]
        )

        existing_point_ids = [
            str(existing_point_id).strip()
            for existing_point_id in existing_point_ids
        ]

        if len(existing_point_ids) != len(set(existing_point_ids)):
            raise RuntimeError(
                f'Duplicated point IDs found in "{file_path_obj}": '
                f'{existing_point_ids}'
            )

        # ---------------------------------------------------------------------------------
        # Find point index
        if point_id in existing_point_ids:

            point_index = existing_point_ids.index(
                point_id
            )

            if not overwrite_point:
                raise RuntimeError(
                    f'Point "{point_id}" already exists in '
                    f'"{file_path_obj}"'
                )

            log_stream.warning(
                ' ===> Point "%s" already exists at index %d. '
                'Existing values will be overwritten',
                point_id,
                point_index
            )

        else:

            point_index = len(
                existing_point_ids
            )

        # ---------------------------------------------------------------------------------
        # Create or validate point variable
        expected_point_dimensions = (
            time_name,
            "variable"
        )

        if point_name in file_nc.variables:

            point_variable = file_nc.variables[
                point_name
            ]

            if point_variable.dimensions != expected_point_dimensions:
                raise RuntimeError(
                    f'Variable "{point_name}" already exists but has '
                    f'dimensions {point_variable.dimensions}; expected '
                    f'{expected_point_dimensions}'
                )

            if point_variable.shape != expected_shape:
                raise RuntimeError(
                    f'Variable "{point_name}" has shape '
                    f'{point_variable.shape}; expected {expected_shape}'
                )

            if not overwrite_point:
                raise RuntimeError(
                    f'Point variable "{point_name}" already exists in '
                    f'"{file_path_obj}"'
                )

        else:

            point_variable = file_nc.createVariable(
                point_name,
                dtype,
                expected_point_dimensions,
                fill_value=fill_value,
                zlib=compression_level > 0,
                complevel=compression_level,
                shuffle=compression_level > 0
            )

        # ---------------------------------------------------------------------------------
        # Write point metadata
        file_nc.variables[
            "point_id"
        ][point_index] = point_id

        file_nc.variables[
            "point_name"
        ][point_index] = point_long_name

        file_nc.variables[
            "longitude"
        ][point_index] = point_lon

        file_nc.variables[
            "latitude"
        ][point_index] = point_lat

        # ---------------------------------------------------------------------------------
        # Write point values
        point_variable[:, :] = point_data

        # ---------------------------------------------------------------------------------
        # Set point attributes
        point_variable.long_name = (
            f"Time series for {point_long_name}"
        )

        point_variable.point_id = point_id
        point_variable.point_name = point_long_name
        point_variable.point_index = point_index
        point_variable.longitude = point_lon
        point_variable.latitude = point_lat

        point_variable.coordinates = (
            f"{time_name} "
            "variable_index "
            "variable_name"
        )

        # ---------------------------------------------------------------------------------
        # Update global attributes
        file_nc.number_of_points = len(
            file_nc.dimensions["point"]
        )

        file_nc.number_of_variables = len(
            file_nc.dimensions["variable"]
        )

    log_stream.info(
        ' ---> Dataset point "%s" saved in "%s" at index %d '
        'with shape (%d, %d)',
        point_name,
        file_path_obj,
        point_index,
        point_data.shape[0],
        point_data.shape[1]
    )
# ----------------------------------------------------------------------------------------------------------------------
