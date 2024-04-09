"""
Library Features:

Name:          lib_data_io_json
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
import numpy as np
import pandas as pd

from copy import deepcopy

from lib_utils_time import define_time_frequency
from lib_utils_obj import find_element_in_list
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert variable(s) to data
def convert_vars_to_data(obj_dframe_in, obj_fields=None, time_format='%Y-%m-%d %H:%M',
                         decimal_precision_data=2, decimal_precision_attrs=2,
                         no_data_generic=-9999, no_data_other=-9999, no_data_ref=-9998):

    # get time
    var_time = list(pd.DatetimeIndex(obj_dframe_in.index).strftime(time_format))

    # get data
    obj_data_tmp = {}
    for data_name, data_values in obj_dframe_in.items():

        if 'time' == data_name:
            var_data = list(pd.DatetimeIndex(data_values).strftime(time_format))
        elif 'ref' == data_name:
            data_values = np.around(data_values, decimal_precision_data)
            data_values[np.isnan(data_values)] = no_data_generic
            var_data = data_values.astype(str).tolist()
        elif 'other' in data_name:

            test_values = data_values.dropna()
            if not test_values.empty:
                data_values = np.around(data_values, decimal_precision_data)
                data_values[np.isnan(data_values)] = no_data_generic
                var_data = data_values.astype(str).tolist()
            else:
                data_null = np.zeros(shape=data_values.shape[0])
                data_null[:] = no_data_generic
                var_data = data_null.astype(str).tolist()
        else:
            log_stream.error(' ===> Field "' + data_name + '" was not expected in the object')
            raise NotImplemented('Case not implemented yet')

        obj_data_tmp[data_name] = var_data

    # get attributes
    attrs_data_in, attrs_list, attrs_null_no_data_other, attrs_null_no_data_ref = {}, None, None, None
    for data_name, data_values in obj_dframe_in.items():

        data_attrs = data_values.attrs
        if data_name in list(data_attrs.keys()):

            obj_attrs_tmp = data_attrs[data_name]

            if attrs_list is None:
                attrs_keys = list(obj_attrs_tmp.keys())
                attrs_values_no_data_other = [no_data_other] * attrs_keys.__len__()
                attrs_null_no_data_other = {k: str(v) for k, v in zip(attrs_keys, attrs_values_no_data_other)}
                attrs_values_no_data_ref = [no_data_ref] * attrs_keys.__len__()
                attrs_null_no_data_ref = {k: str(v) for k, v in zip(attrs_keys, attrs_values_no_data_ref)}

            attrs_workspace = {}
            for attrs_name, attrs_value in obj_attrs_tmp.items():

                if np.isnan(attrs_value):
                    attrs_value = no_data_generic

                if isinstance(attrs_value, float):
                    var_data = np.around(attrs_value, decimal_precision_attrs)
                    var_data = str(var_data)
                elif isinstance(attrs_value, int):
                    var_data = str(attrs_value)
                elif isinstance(attrs_value, str):
                    var_data = str(attrs_value)
                else:
                    log_stream.error(' ===> Field "' + data_name + '" was not expected in the object')
                    raise NotImplemented('Case not implemented yet')

                attrs_workspace[attrs_name] = var_data

            attrs_data_in[data_name] = attrs_workspace

        else:
            attrs_data_in[data_name] = None

    attrs_data_tmp = {}
    if attrs_null_no_data_other is not None:
        for attrs_name, attrs_values in attrs_data_in.items():
            if attrs_values is None:
                if attrs_name == 'ref':
                    attrs_data_tmp[attrs_name] = attrs_null_no_data_ref
                else:
                    attrs_data_tmp[attrs_name] = attrs_null_no_data_other
            else:
                attrs_data_tmp[attrs_name] = attrs_values

    # organize data and attributes
    obj_data_out, obj_attrs_out = {}, {}
    if obj_fields is not None:
        for data_name, data_value in obj_data_tmp.items():
            if data_name in list(obj_fields.keys()):
                field_name = obj_fields[data_name]
                obj_data_out[field_name] = data_value
                obj_attrs_out[field_name] = attrs_data_tmp[data_name]
            else:
                obj_data_out[data_name] = data_value
                obj_attrs_out[data_name] = attrs_data_tmp[data_name]
    else:
        for data_name, data_value in obj_dframe_in.items():
            obj_data_out[data_name] = data_value
            obj_attrs_out[data_name] = attrs_data_tmp[data_name]

    var_dframe_out = pd.DataFrame(data=obj_data_out, index=var_time)
    var_dframe_out.index.name = 'time'

    attrs_dframe_out = pd.DataFrame(data=obj_attrs_out)
    attrs_dframe_out.index.name = 'metrics'

    return var_dframe_out, attrs_dframe_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert data to variable(s)
def convert_data_to_vars(obj_data, obj_fields=None, delimiter_format=','):

    if obj_fields is not None:
        tmp_data = {}
        for data_name, data_value in obj_data.items():
            if data_name in list(obj_fields.values()):
                field_idx = list(obj_fields.values()).index(data_name)
                field_name = list(obj_fields.keys())[field_idx]
                tmp_data[field_name] = data_value

                if isinstance(data_value, pd.Series):
                    tmp_name = data_value.index.name
                    if tmp_name not in list(tmp_data.keys()):
                        tmp_data[tmp_name] = pd.DatetimeIndex(data_value.index.values)

            else:
                tmp_data[data_name] = data_value
    else:
        tmp_data = deepcopy(obj_data)

    dframe_data, dframe_attrs = {}, {}
    for field_name, field_data in tmp_data.items():

        if isinstance(field_data, pd.Series):
            field_parts = field_data.values
        elif isinstance(field_data, pd.DatetimeIndex):
            field_parts = deepcopy(field_data)
        elif isinstance(field_data, str):
            field_parts = field_data.split(delimiter_format)

            if not isinstance(field_parts, list):
                field_parts = [field_parts]
        else:
            log_stream.error(' ===> Field "' + field_name + '" format is not supported')
            raise NotImplemented('Case not implemented yet')

        if field_parts.__len__() > 1:
            if 'time' == field_name:
                if isinstance(field_parts, list):
                    var_data = pd.DatetimeIndex(field_parts)
                elif isinstance(field_parts, pd.DatetimeIndex):
                    var_data = deepcopy(field_parts)
                else:
                    log_stream.error(' ===> Values of "' + field_name + '" format are not supported')
                    raise NotImplemented('Case not implemented yet')
            elif 'ref' == field_name:
                if isinstance(field_parts, list):
                    var_data = np.array(field_parts).astype(float)
                elif isinstance(field_parts, np.ndarray):
                    var_data = field_parts.astype(float)
                else:
                    log_stream.error(' ===> Values of "' + field_name + '" format are not supported')
                    raise NotImplemented('Case not implemented yet')
            elif 'other' in field_name:
                if isinstance(field_parts, list):
                    var_data = np.array(field_parts).astype(float)
                elif isinstance(field_parts, np.ndarray):
                    var_data = field_parts.astype(float)
                else:
                    log_stream.error(' ===> Values of "' + field_name + '" format are not supported')
                    raise NotImplemented('Case not implemented yet')
            else:
                log_stream.error(' ===> Field "' + field_name + '" was not expected in the dataframe')
                raise NotImplemented('Case not implemented yet')

            dframe_data[field_name] = var_data

        elif field_parts.__len__() == 1:
            if re.search('\d', field_parts[0]):
                var_data = float(field_parts[0])
            else:
                var_data = field_parts[0]

            dframe_attrs[field_name] = var_data

        else:
            log_stream.warning(' ===> Field "' + field_name + '" was not converted by the procedure')

    dframe_time = dframe_data['time']
    dframe_data.pop('time')

    dframe_data = pd.DataFrame(data=dframe_data, index=dframe_time)
    dframe_data.index.name = 'time'
    dframe_data.attrs = dframe_attrs

    dframe_data = dframe_data.sort_index()

    return dframe_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert vars to dict
def convert_vars_to_dict(obj_data, time_format='%Y-%m-%d %H:%M', decimal_precision=2):

    obj_data = obj_data.round(decimal_precision)
    obj_data = obj_data.reset_index()

    obj_fields = list(obj_data.columns)

    dict_data = {}
    for field_name in obj_fields:

        field_data = obj_data[field_name]
        if 'time' == field_name:
            var_data = list(pd.DatetimeIndex(field_data.values).strftime(time_format))
        else:
            var_data = field_data.values
        dict_data[field_name] = var_data

    return dict_data

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to adjust data time
def adjust_data_time(dframe_data_generic, time_start=None, time_end=None):

    dframe_index = dframe_data_generic.index

    if time_start is None:
        time_start = dframe_index[0]
    if time_end is None:
        time_end = dframe_index[-1]

    time_frequency = define_time_frequency(dframe_index)

    return time_start, time_end, time_frequency

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to adjust data point
def adjust_data_point(dframe_data_generic, scale_factor=1, value_min=0, value_max=100, value_no_data=np.nan,
                      time_start=None, time_end=None, time_range=None, time_freq='D'):

    if time_range is None:
        time_range_expected = pd.date_range(start=time_start, end=time_end, freq=time_freq)
    else:
        time_range_expected = deepcopy(time_range)

    dframe_data_generic = dframe_data_generic.sort_index()
    time_range_dframe = dframe_data_generic.index.values

    if time_start in time_range_dframe:
        index_start = dframe_data_generic.index.get_loc(time_start)
    else:
        if time_start in dframe_data_generic.index:
            index_start = dframe_data_generic.index.get_loc(time_start, method='nearest')
        elif time_start < dframe_data_generic.index[0]:
            index_start = 0
        else:
            logging.error(' ===> Time start is not available in the dataframe')
            raise NotImplemented('Case not implemented yet')

    if time_end in time_range_dframe:
        index_end = dframe_data_generic.index.get_loc(time_end) + 1
    else:
        if time_end in dframe_data_generic.index:
            index_end = dframe_data_generic.index.get_loc(time_end, method='nearest') + 1
        elif time_end > dframe_data_generic.index[-1]:
            index_end = dframe_data_generic.index.size
        else:
            logging.error(' ===> Time end is not available in the dataframe')
            raise NotImplemented('Case not implemented yet')

    dframe_data_select = dframe_data_generic.iloc[index_start:index_end]
    # time_start_select, time_end_select = dframe_data_select.index[0], dframe_data_select.index[-1]

    dframe_data_expected = pd.DataFrame(index=time_range_expected)
    dframe_data_expected = dframe_data_expected.sort_index()

    dframe_data_expected = dframe_data_expected.join(dframe_data_select)

    if scale_factor is not None:
        dframe_data_expected = dframe_data_expected * scale_factor

    if value_min is not None:
        dframe_data_expected[dframe_data_expected < value_min] = value_no_data
    if value_max is not None:
        dframe_data_expected[dframe_data_expected > value_max] = value_no_data

    dframe_data_expected[dframe_data_expected == value_no_data] = np.nan

    return dframe_data_expected
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to fill data point
def fill_data_point(dframe_data_in,
                    fill_limit=2, fill_direction='forward',
                    fill_method='polynomial', fill_order=2):

    if fill_method == 'polynomial':
        dframe_data_out = dframe_data_in.interpolate(
            method='polynomial', order=fill_order,
            limit=fill_limit, limit_direction=fill_direction)
    else:
        dframe_data_out = deepcopy(dframe_data_in)

    return dframe_data_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data point
def resample_data_point(dframe_data_in,
                        resample_frequency=None, resample_method=None):

    if (resample_frequency is not None) and (resample_method is not None):

        var_time_list = dframe_data_in.index
        var_name_list = list(dframe_data_in.columns)

        if resample_method == 'average':

            dframe_data_out = pd.DataFrame(index=var_time_list)
            dframe_data_out = dframe_data_out.resample(resample_frequency).mean()
            dframe_data_out = dframe_data_out.sort_index()

            for var_name_step in var_name_list:

                series_data_in = dframe_data_in[var_name_step]
                series_data_out = series_data_in.resample(resample_frequency).mean()

                dframe_data_out = dframe_data_out.join(series_data_out)

        else:
            log_stream.error(' ===> Resample time operation "' + resample_method + '" is not supported')
            raise NotImplementedError('Case not implemented yet. Only "average" method is available')

    elif (resample_frequency is not None) and (resample_method is None):
        dframe_data_out = dframe_data_in.resample(resample_frequency)
    else:
        dframe_data_out = deepcopy(dframe_data_in)

    return dframe_data_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to join data point
def join_data_point(time_range, point_collections, point_registry,
                    tag_registry='tag', name_registry='name', reference_registry='tag'):

    # create a template dataframe
    point_dframe_tmpl = pd.DataFrame(index=time_range)
    point_dframe_tmpl = point_dframe_tmpl.sort_index()
    time_index_tmpl = point_dframe_tmpl.index

    # iterate over point(s)
    point_obj_collections = {}
    for point_id, point_row in point_registry.iterrows():

        # get point tag
        if tag_registry in list(point_row.keys()):
            point_tag = point_row[tag_registry]
        else:
            log_stream.error(' ===> Registry key "' + tag_registry + '" is not available in the dataframe')
            raise RuntimeError('Change reference key to correctly select the dataframe column')
        # get point name
        if name_registry in list(point_row.keys()):
            point_name = point_row[name_registry]
        else:
            log_stream.error(' ===> Registry key "' + name_registry + '" is not available in the dataframe')
            raise RuntimeError('Change reference key to correctly select the dataframe column')
        # get point reference
        if reference_registry in list(point_row.keys()):
            point_reference = point_row[reference_registry]
        else:
            log_stream.error(' ===> Registry key "' + reference_registry + '" is not available in the dataframe')
            raise RuntimeError('Change reference key to correctly select the dataframe column')

        # info point start
        log_stream.info(' ------> Compose time-series for point "' + point_reference + '" ... ')

        # iterate over datasets to join to common obj
        dset_obj_collections = {}
        for point_dset, point_dframe in point_collections.items():

            # info datasets start
            log_stream.info(' -------> Datasets "' + point_dset + '" ... ')

            # check datasets
            if point_dframe is not None:

                # get reference and other indexes
                point_index = point_dframe.index.values

                # find datasets key
                point_key = find_element_in_list([point_tag, point_name], list(point_dframe.columns))

                # check datasets key
                if point_key is not None:
                    if point_key in list(point_dframe.columns):

                        # get data
                        point_data = point_dframe[point_key].values
                        point_dframe_tmp = pd.DataFrame(data=point_data, index=point_index)
                        point_dframe_tmp = point_dframe_tmp.sort_index()
                        point_dframe_common = deepcopy(point_dframe_tmpl)

                        # join common and tmp dframe(s)
                        point_dframe_common = point_dframe_common.join(point_dframe_tmp)
                        # remove index duplicate(s)
                        point_dframe_common = point_dframe_common[~point_dframe_common.index.duplicated(keep='first')]

                        point_values_common = point_dframe_common.values[:, 0]

                    else:
                        log_stream.warning(' ===> Datasets key is not available in the source reference object')
                        point_values_common = np.zeros(shape=time_index_tmpl.shape[0])
                        point_values_common[:] = np.nan
                else:
                    log_stream.error(' ===> Datasets key is defined by NoneType')
                    raise RuntimeError('Key must be defined to correctly run the algorithm')
            else:
                log_stream.warning(' ===> Datasets are not available in the source reference object')
                point_values_common = np.zeros(shape=time_index_tmpl.shape[0])
                point_values_common[:] = np.nan

            # store datasets in a common object
            dset_obj_collections[point_dset] = point_values_common

            # info datasets end
            log_stream.info(' -------> Datasets "' + point_dset + '" ... DONE')

        # convert dict to dframe object
        dframe_obj = pd.DataFrame(data=dset_obj_collections, index=time_index_tmpl)
        dframe_obj.index.name = 'time'
        dframe_obj = dframe_obj.sort_index()

        # store dataframe to common obj
        point_obj_collections[point_reference] = dframe_obj

        log_stream.info(' ------> Compose time-series for point "' + point_reference + '" ... DONE')

    return point_obj_collections
# ----------------------------------------------------------------------------------------------------------------------
