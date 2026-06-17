"""
Library Features:

Name:          lib_data_io_csv
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import csv
import numpy as np
import pandas as pd

from copy import deepcopy

from lib_utils_generic import fill_tags2string, invert_dict
from lib_utils_obj import map_vars_dframe, sanitize_string, fill_tags_time
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
logging.getLogger('pandas').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read metrics
def read_metrics_csv(file_name,
                      file_fields,
                      registry_fields,
                      file_sep=',',
                      file_decimal='.',
                      **kwargs):

    try:
        metrics_data = pd.read_csv(
            file_name,
            sep=file_sep,
            decimal=file_decimal
        )
    except Exception as exc:
        log_stream.error(f' ===> Error reading file "{file_name}": {exc}')
        raise RuntimeError(f'Unable to read file "{file_name}"') from exc

    if file_fields is None:
        file_fields = {}

    return metrics_data
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read datasets in csv format
def read_datasets_csv(file_name,
                      file_fields,
                      registry_fields,
                      time_reference,
                      time_start=None,
                      time_end=None,
                      time_rounding='H',
                      time_frequency='H',
                      time_format='%Y%m%d%H%M',
                      file_sep=',',
                      file_decimal='.',
                      ascending_index=False,
                      sort_index=True,
                      point_tag=None,
                      point_name=None,
                      **kwargs):

    try:
        fields_data_raw = pd.read_csv(
            file_name,
            sep=file_sep,
            decimal=file_decimal
        )
    except Exception as exc:
        log_stream.error(f' ===> Error reading file "{file_name}": {exc}')
        raise RuntimeError(f'Unable to read file "{file_name}"') from exc

    if file_fields is None:
        file_fields = {}

    if registry_fields is None:
        registry_fields = {}

    time_col = file_fields.get('time', 'time')

    if time_col not in fields_data_raw.columns:
        raise RuntimeError(
            f'Time column "{time_col}" not found in file "{file_name}"'
        )

    # ---------------------------------------------------------------------
    # Case 1: read all points
    # Example input:
    # time | point_1 | point_2 | ... | point_1203
    # ---------------------------------------------------------------------
    if point_tag is None:

        fields_data_map = fields_data_raw.copy()

        if time_col != 'time':
            fields_data_map = fields_data_map.rename(columns={time_col: 'time'})

        columns_map = {
            file_col: field_name
            for field_name, file_col in file_fields.items()
            if field_name != 'time'
        }

        fields_data_map = fields_data_map.rename(columns=columns_map)

    # ---------------------------------------------------------------------
    # Case 2: read/extract one point
    # Example:
    # file_fields["values_k1"] = "point_{point_tag}"
    # point_tag = "1"
    # selected column = "point_1"
    #
    # Output:
    # time | values_k1
    # ---------------------------------------------------------------------
    else:

        columns_map = {}

        for field_name, file_col_tmpl in file_fields.items():

            if field_name == 'time':
                continue

            if isinstance(file_col_tmpl, str):
                file_col = file_col_tmpl.format(
                    point_tag=point_tag,
                    point_name=point_name
                )
            else:
                file_col = file_col_tmpl

            if file_col not in fields_data_raw.columns:
                raise RuntimeError(
                    f'Column "{file_col}" for field "{field_name}" '
                    f'not found in file "{file_name}"'
                )

            columns_map[file_col] = field_name

        columns_select = [time_col] + list(columns_map.keys())

        fields_data_map = fields_data_raw[columns_select].copy()

        fields_data_map = fields_data_map.rename(columns={
            time_col: 'time',
            **columns_map
        })

    fields_data_map['time'] = pd.to_datetime(
        fields_data_map['time'],
        format=time_format,
        errors='coerce'
    )

    fields_data_map = fields_data_map.dropna(subset=['time'])

    fields_data_map.index = pd.DatetimeIndex(fields_data_map['time'])
    fields_data_map.index.name = 'time'

    if (time_start is not None) and (time_end is not None):

        time_start = pd.Timestamp(time_start).floor(time_rounding.lower())
        time_end = pd.Timestamp(time_end).floor(time_rounding.lower())

        time_range = pd.date_range(
            time_start,
            time_end,
            freq=time_frequency.lower()
        )

        fields_data_select = pd.DataFrame(index=time_range)
        fields_data_select.index.name = 'time'
        fields_data_select = fields_data_select.join(fields_data_map)

    else:

        fields_data_select = deepcopy(fields_data_map)

    if sort_index:
        fields_data_select = fields_data_select.sort_index(
            ascending=ascending_index
        )

    fields_data_select.attrs = registry_fields
    fields_data_select.attrs['time_reference'] = time_reference

    return fields_data_select

# method to extract point from dframe
def extract_point_from_dframe(fields_data_all,
                              file_fields,
                              point_tag,
                              point_name=None,
                              file_name=None):

    if file_fields is None:
        file_fields = {}

    if point_tag is None:
        raise RuntimeError('Point tag must be defined to extract point data')

    value_fields = {
        field_name: file_col_tmpl
        for field_name, file_col_tmpl in file_fields.items()
        if field_name != 'time'
    }

    if not value_fields:
        raise RuntimeError('No value fields found in file_fields')

    columns_map = {}

    for field_name, file_col_tmpl in value_fields.items():

        file_col = file_col_tmpl.format(
            point_tag=point_tag,
            point_name=point_name
        )

        if file_col not in fields_data_all.columns:
            extra = f' in file "{file_name}"' if file_name else ''
            raise RuntimeError(
                f'Column "{file_col}" for point "{point_tag}" not found{extra}'
            )

        columns_map[file_col] = field_name

    fields_data_point = fields_data_all[list(columns_map.keys())].copy()

    fields_data_point = fields_data_point.rename(columns=columns_map)

    if 'time' in fields_data_all.columns:
        fields_data_point['time'] = fields_data_all['time']

    fields_data_point.index = fields_data_all.index
    fields_data_point.index.name = 'time'

    fields_data_point.attrs = deepcopy(fields_data_all.attrs)

    return fields_data_point
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read datasets in csv format
def read_datasets_csv_OLD(file_name,
                      file_fields, registry_fields,
                      time_reference, time_start=None, time_end=None,
                      time_rounding='H', time_frequency='Y', time_format='%Y%m%d%H%M',
                      file_sep=' ', file_decimal='.',
                      ascending_index=False, sort_index=True, **kwargs):
    # get file fields
    try:
        fields_data_raw = pd.read_csv(
            file_name, sep=file_sep, decimal=file_decimal, date_format=time_format)
    except Exception as exc:
        log_stream.warning(' ===> Library exception: ' + str(exc) + '. Try to use "date parser"')
        fields_data_raw = pd.read_csv(
            file_name, sep=file_sep, decimal=file_decimal, date_parser=time_format)

    if file_fields is None:
        file_fields = {}
    if registry_fields is None:
        registry_fields = {}

    # organize file fields
    tmp_fields = invert_dict(file_fields)
    fields_data_map = fields_data_raw.rename(columns=tmp_fields)
    fields_data_map.reset_index()
    fields_data_map.index = pd.DatetimeIndex(fields_data_map['time'])
    fields_data_map.index.name = 'time'

    # select file fields by time range
    if (time_start is not None) and (time_start is not None):

        time_start = pd.Timestamp(time_start).floor(time_rounding.lower())
        time_end = pd.Timestamp(time_end).floor(time_rounding.lower())
        time_range = pd.date_range(time_start, time_end, freq=time_frequency.lower())

        fields_data_select = pd.DataFrame(index=time_range)
        fields_data_select = fields_data_select.join(fields_data_map)

    else:
        fields_data_select = deepcopy(fields_data_map)

    # sort index
    if sort_index:
        if ascending_index:
            fields_data_select = fields_data_select.sort_index(ascending=True)
        else:
            fields_data_select = fields_data_select.sort_index(ascending=False)

    # add attributes
    if registry_fields is not None:
        fields_data_select.attrs = registry_fields
    fields_data_select.attrs['time_reference'] = time_reference

    return fields_data_select
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read parameters in csv format
def read_parameters_csv(file_name, file_fields, file_filters=None, file_sep=',', file_decimal='.'):

    # get file fields
    fields_data_raw = pd.read_table(file_name, sep=file_sep, decimal=file_decimal)
    fields_data_raw.columns = fields_data_raw.columns.str.strip()
    # map file fields
    fields_data_map = map_vars_dframe(fields_data_raw, file_fields)

    # get name fields
    var_data_name = fields_data_map['name'].values
    # parser tag
    if 'tag' not in list(fields_data_map.keys()):
        var_data_tag = []
        for string_name in var_data_name:
            string_tag = sanitize_string(string_name)
            var_data_tag.append(string_tag)

        fields_data_map['tag'] = var_data_tag

    # create fields dframe
    fields_dframe = pd.DataFrame(data=fields_data_map)

    if 'valid' in list(fields_dframe.columns):
        fields_dframe = fields_dframe[fields_dframe['valid'] == 1]

    if 'code' in list(fields_dframe.columns):
        fields_dframe['code'] = fields_dframe['code'].apply(str)
    else:
        log_stream.warning(' ===> Code field is not available in the registry file. Set default values.')
        rows_dframe = fields_dframe.shape[0]
        fields_dframe['code'] = [str(i) for i in range(1, rows_dframe + 1)]

    if 'tag' in list(fields_dframe.columns):
        fields_dframe['tag'] = fields_dframe['tag'].str.strip()

    if file_filters is not None:
        for filter_key, filter_value in file_filters.items():
            if filter_key in list(fields_dframe.columns):
                filter_data = fields_dframe[filter_key].values
                filter_id_list = []
                for filter_id_step, filter_row in enumerate(filter_data):

                    if not isinstance(filter_value, str):
                        filter_value = str(filter_value)
                    if not isinstance(filter_row, str):
                        filter_row = str(filter_row)

                    if filter_value in filter_row:
                        filter_id_list.append(filter_id_step)

                fields_dframe = fields_dframe.loc[filter_id_list]

    return fields_dframe
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read registry in csv format
def read_registry_csv(file_name, file_fields, file_filters=None, file_sep=',', file_decimal='.'):

    # get file fields
    fields_data_raw = pd.read_table(file_name, sep=file_sep, decimal=file_decimal)
    fields_data_raw.columns = fields_data_raw.columns.str.strip()
    # map file fields
    fields_data_map = map_vars_dframe(fields_data_raw, file_fields)

    # get name fields
    var_data_name = fields_data_map['name'].values
    # parser tag
    if 'tag' not in list(fields_data_map.keys()):
        var_data_tag = []
        for string_name in var_data_name:
            string_tag = sanitize_string(string_name)
            var_data_tag.append(string_tag)

        fields_data_map['tag'] = var_data_tag

    # create fields dframe
    fields_dframe = pd.DataFrame(data=fields_data_map)

    if 'valid' in list(fields_dframe.columns):
        fields_dframe = fields_dframe[fields_dframe['valid'] == 1]

    if 'code' in list(fields_dframe.columns):
        fields_dframe['code'] = fields_dframe['code'].apply(str)
    else:
        log_stream.warning(' ===> Code field is not available in the registry file. Set default values.')
        rows_dframe = fields_dframe.shape[0]
        fields_dframe['code'] = [str(i) for i in range(1, rows_dframe + 1)]

    if 'tag' in list(fields_dframe.columns):
        fields_dframe['tag'] = fields_dframe['tag'].str.strip()

    if file_filters is not None:
        for filter_key, filter_value in file_filters.items():
            if filter_value is not None:
                if filter_key in list(fields_dframe.columns):
                    filter_data = fields_dframe[filter_key].values
                    filter_id_list = []
                    for filter_id_step, filter_row in enumerate(filter_data):

                        if not isinstance(filter_value, str):
                            filter_value = str(filter_value)
                        if not isinstance(filter_row, str):
                            filter_row = str(filter_row)

                        if filter_value in filter_row:
                            filter_id_list.append(filter_id_step)

                    fields_dframe = fields_dframe.loc[filter_id_list]

    return fields_dframe
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write metrics csv
def write_metrics_csv(file_name, file_dframe, file_fields=None,
                      file_sep=',', file_decimal='.', file_float_format='%.2f',
                      file_index=False, file_header=True,
                      **kwargs):

    if file_fields is not None:
        file_fields = {}

    # dump file
    file_dframe.to_csv(
        file_name, mode='w',
        index=file_index, sep=file_sep, decimal=file_decimal,
        header=file_header, float_format=file_float_format,  quotechar='"')

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write datasets csv
def write_datasets_csv(file_name, file_dframe, file_fields=None, time_fields=None,
                       file_sep=',', file_decimal='.', file_float_format='%.2f',
                       file_index=True, file_header=True,
                       time_index_label='time', time_index_format='%Y-%m-%d %H:%M',
                       file_no_data=-9999,
                       ascending_index=False, sort_index=True, **kwargs):

    # organize file fields
    if file_fields is not None:
        file_dframe = file_dframe.rename(columns=file_fields)
    # set no data value
    if np.isfinite(file_no_data):
        file_dframe.fillna(file_no_data, inplace=True)
    # parse index
    if time_index_format is not None:
        file_dframe.index = file_dframe.index.strftime(time_index_format)

    # sort index
    if sort_index:
        if ascending_index:
            file_dframe = file_dframe.sort_index(ascending=True)
        else:
            file_dframe = file_dframe.sort_index(ascending=False)

    # remove time label if available in the columns
    if file_index:
        if time_index_label in list(file_dframe.columns):
            file_dframe = file_dframe.drop(columns=[time_index_label])

    # dump file
    file_dframe.to_csv(
        file_name, mode='w',
        index=file_index, sep=file_sep, decimal=file_decimal,
        index_label=time_index_label,
        header=file_header, float_format=file_float_format,  quotechar='"')

# ----------------------------------------------------------------------------------------------------------------------
