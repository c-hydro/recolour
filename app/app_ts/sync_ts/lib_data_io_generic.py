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
# method to combine data over the expected time range
def combine_data_point_by_time(point_data_collections_src, point_geo,
                               time_start_expected, time_end_expected,
                               time_frequency_expected='D', time_reverse=True):
    # create time range expected
    time_range_expected = pd.date_range(start=time_start_expected, end=time_end_expected, freq=time_frequency_expected)

    # iterate over variable(s)
    point_data_collections_dst = {}
    if point_data_collections_src is not None:
        for point_var, point_data_raw in point_data_collections_src.items():

            # create null data
            null_data = np.zeros(shape=(time_range_expected.shape[0], point_geo.shape[0]))
            null_data[:, :] = np.nan
            null_dict = {}
            for point_id, point_label in enumerate(point_data_raw.columns):
                null_dict[point_label] = null_data[:, point_id]
            point_data_expected = pd.DataFrame(data=null_dict, index=time_range_expected)

            # update expected data with raw data
            point_data_expected.update(point_data_raw)
            # time reverse flag
            if time_reverse:
                point_data_expected = point_data_expected.sort_index(ascending=False)

            # store data in a common obj
            point_data_collections_dst[point_var] = point_data_expected

    else:
        log_stream.warning(' ===> No data available to combine over the expected time range')
        point_data_collections_dst = None

    return point_data_collections_dst

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to range data point
def range_data_point(point_collection,
                     time_start_reference, time_end_reference,
                     time_run_reference=None):

    if time_run_reference is not None:
        if isinstance(time_run_reference, str):
            time_run_reference = pd.Timestamp(time_run_reference)
        elif isinstance(time_run_reference, pd.Timestamp):
            pass
        else:
            log_stream.error(' ===> Time run reference is not correctly defined. Check your settings file.')
            raise RuntimeError('Time run reference is not correctly defined')

    if isinstance(time_start_reference, str):
        time_start_reference = pd.Timestamp(time_start_reference)
    elif isinstance(time_start_reference, pd.Timestamp):
        pass
    if isinstance(time_end_reference, str):
        time_end_reference = pd.Timestamp(time_end_reference)
    elif isinstance(time_end_reference, pd.Timestamp):
        pass

    time_start_collection, time_end_collection, frequency_collection, seconds_collections = [], [], [], []
    for point_name, point_data in point_collection.items():
        time_index_step = point_data.index
        time_start_step, time_end_step = time_index_step[0], time_index_step[-1]
        time_frequency_step = define_time_frequency(time_index_step)

        if time_frequency_step[0].isalpha():
            tmp_frequency_step = ''.join(['1', time_frequency_step])
        else:
            tmp_frequency_step = deepcopy(time_frequency_step)

        seconds_step = pd.to_timedelta(tmp_frequency_step).total_seconds()

        seconds_collections.append(seconds_step)
        frequency_collection.append(time_frequency_step)

        if time_start_reference is not None:
            if time_start_reference < time_start_step:
                time_start_step = deepcopy(time_start_reference)
        if time_end_reference is not None:
            if time_end_reference > time_end_step:
                time_end_step = deepcopy(time_end_reference)
            if time_run_reference > time_end_step:
                time_end_step = deepcopy(time_run_reference)

        time_start_collection.append(time_start_step)
        time_end_collection.append(time_end_step)

    idx_min = np.argmin(np.array(seconds_collections))

    time_frequency = frequency_collection[idx_min]
    time_start = pd.DatetimeIndex(time_start_collection).min()
    time_end = pd.DatetimeIndex(time_end_collection).max()

    if time_run_reference > time_end:
        time_end = deepcopy(time_run_reference)

    time_range = pd.date_range(start=time_start, end=time_end, freq=time_frequency)

    return time_frequency, time_start, time_end, time_range
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert registry to dict
def convert_registry_point_to_dict(reference_name, obj_registry, tag_registry='tag'):
    df_point = obj_registry.loc[obj_registry[tag_registry] == reference_name]
    dict_point = df_point.to_dict('records')[0]
    return dict_point
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert datasets to dict
def convert_datasets_point_to_dict(obj_data, time_format='%Y-%m-%d %H:%M', decimal_precision=2):

    obj_data = obj_data.round(decimal_precision)
    obj_data = obj_data.reset_index()

    dict_data = {}
    for column_name, column_data in obj_data.items():

        if 'time' == column_name:
            var_data = list(pd.DatetimeIndex(column_data.values).strftime(time_format))
        elif 'ref' == column_name:
            var_data = column_data.values.astype(str).tolist()
        elif 'other' in column_name:
            var_data = column_data.values.astype(str).tolist()
        else:
            log_stream.error(' ===> Column "' + column_name + '" was not expected in the dataframe')
            raise NotImplemented('Case not implemented yet')

        dict_data[column_name] = var_data

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
                      time_start=None, time_end=None, time_freq='H'):

    time_range_expected = pd.date_range(start=time_start, end=time_end, freq=time_freq)

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
    time_start_select, time_end_select = dframe_data_select.index[0], dframe_data_select.index[-1]

    dframe_data_expected = pd.DataFrame(index=time_range_expected)
    dframe_data_expected = dframe_data_expected.sort_index()

    dframe_data_expected = dframe_data_expected.join(dframe_data_select)

    if scale_factor is not None:
        dframe_data_expected = dframe_data_expected * scale_factor

    if value_min is not None:
        dframe_data_expected[dframe_data_expected < value_min] = value_no_data
    if value_max is not None:
        dframe_data_expected[dframe_data_expected > value_max] = value_no_data

    return dframe_data_expected
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to fill data point
def fill_data_point(dframe_data_in,
                    fill_limit=2, fill_direction='forward',
                    fill_method='polynomial', fill_order=2):

    if fill_method == 'polynomial':
        try:
            dframe_data_out = dframe_data_in.interpolate(
                method='polynomial', order=fill_order,
                limit=fill_limit, limit_direction=fill_direction)
        except BaseException as base_exp:
            log_stream.warning(' ===> Polynomial time-series interpolation failed: ' + str(base_exp))
            log_stream.warning(' ===> Try to use linear interpolation')

            dframe_data_out = dframe_data_in.interpolate(
                method='linear', limit=fill_limit, limit_direction=fill_direction)

    elif fill_method == 'linear':
        dframe_data_out = dframe_data_in.interpolate(
            method='linear', limit=fill_limit, limit_direction=fill_direction)
    else:
        dframe_data_out = deepcopy(dframe_data_in)

    return dframe_data_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to resample data point
def resample_data_point(dframe_data_in,
                        resample_frequency=None, resample_method=None, resample_direction='left'):

    if (resample_frequency is not None) and (resample_method is not None):

        var_time_list = dframe_data_in.index
        var_name_list = list(dframe_data_in.columns)

        if resample_method == 'average':

            data_frequency = dframe_data_in.index.freq.freqstr

            if not data_frequency[0].isnumeric():
                data_frequency = '1' + data_frequency
            if not resample_frequency[0].isnumeric():
                resample_frequency = '1' + resample_frequency

            seconds_data_frequency = int(pd.to_timedelta(data_frequency).total_seconds())
            seconds_resample_frequency = int(pd.to_timedelta(resample_frequency).total_seconds())
            str_data_frequency, str_resample_frequency = data_frequency[1:], resample_frequency[1:]

            # check resample and data frequency
            if seconds_data_frequency < seconds_resample_frequency:
                var_time_start, var_time_end = var_time_list[0], var_time_list[-1]

                flag_inverse_time = False
                if var_time_start > var_time_end:
                    var_time_list = sorted(var_time_list)
                    var_time_start, var_time_end = var_time_list[0], var_time_list[-1]
                    flag_inverse_time = True

                if resample_direction == 'right':
                    var_time_start = var_time_start - pd.Timedelta(1, str_data_frequency)
                elif resample_direction == 'left':
                    var_time_end = var_time_end + pd.Timedelta(1, str_data_frequency)
                else:
                    log_stream.error(' ===> Resample direction "' + resample_direction + '" is not allowed')
                    raise RuntimeError('Check the resample direction and choose "left" or "right"')

                var_time_list = pd.date_range(start=var_time_start, end=var_time_end, freq=str_resample_frequency)

                if flag_inverse_time:
                    var_time_list = var_time_list[::-1]

                if resample_direction == 'right':
                    dframe_data_out = pd.DataFrame(index=var_time_list[1:])
                    dframe_data_out = dframe_data_out.sort_index()
                    var_time_period = var_time_list
                elif resample_direction == 'left':
                    dframe_data_out = pd.DataFrame(index=var_time_list[:-1])
                    dframe_data_out = dframe_data_out.sort_index()
                    var_time_period = var_time_list
                else:
                    log_stream.error(' ===> Resample direction "' + resample_direction + '" is not allowed')
                    raise RuntimeError('Check the resample direction and choose "left" or "right"')

                for var_name_step in var_name_list:

                    series_data_in = dframe_data_in[var_name_step]

                    time_list, values_list = [], []
                    for id_period, time_period in enumerate(var_time_period):

                        time_check = False
                        if resample_direction == 'right':
                            time_check = id_period < len(var_time_period) - 1
                        elif resample_direction == 'left':
                            time_check = id_period < len(var_time_period) - 1
                        else:
                            log_stream.error(' ===> Resample direction "' + resample_direction + '" is not allowed')
                            raise RuntimeError('Check the resample direction and choose "left" or "right"')

                        # check time conditions (to avoid not defined index)
                        if time_check:

                            time_start, time_end = var_time_period[id_period], var_time_period[id_period + 1]
                            if resample_direction == 'right':
                                time_start = time_start + pd.Timedelta(1, str_data_frequency)
                            elif resample_direction == 'left':
                                time_end = time_end - pd.Timedelta(1, str_data_frequency)
                            else:
                                log_stream.error(' ===> Resample direction "' + resample_direction + '" is not allowed')
                                raise RuntimeError('Check the resample direction and choose "left" or "right"')

                            series_data_tmp = series_data_in[time_start:time_end]
                            value_tmp = series_data_tmp.mean()

                            values_list.append(value_tmp)
                            if resample_direction == 'left':
                                time_list.append(time_start)
                            elif resample_direction == 'right':
                                time_list.append(time_end)
                            else:
                                log_stream.error(' ===> Resample direction "' + resample_direction + '" is not allowed')
                                raise RuntimeError('Check the resample direction and choose "left" or "right"')

                    # create series out
                    series_data_out = pd.Series(data=values_list, index=time_list, name=var_name_step)

                    # join dframe out
                    dframe_data_out = dframe_data_out.join(series_data_out)

            elif seconds_data_frequency == seconds_resample_frequency:
                dframe_data_out = deepcopy(dframe_data_in)
            else:
                log_stream.error(' ===> Resample frequency operation is not supported')
                raise NotImplementedError('Case not implemented yet. Only "average" method is available')

        else:
            log_stream.error(' ===> Resample time operation "' + resample_method + '" is not supported')
            raise NotImplementedError('Case not implemented yet. Only "average" method is available')

    elif (resample_frequency is not None) and (resample_method is None):
        log_stream.error(' ===> Resample frequency and method are not defined')
        raise NotImplemented('Case not implemented yet.')
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
        dset_obj_collections, dset_obj_list = {}, []
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
                    log_stream.warning(' ===> Datasets key for "' + point_tag +
                                       '" is not defined or defined by NoneType')
                    point_values_common = np.zeros(shape=time_index_tmpl.shape[0])
                    point_values_common[:] = np.nan

            else:
                log_stream.warning(' ===> Datasets are not available in the source reference object')
                point_values_common = np.zeros(shape=time_index_tmpl.shape[0])
                point_values_common[:] = np.nan

            # store datasets in a common object
            if point_values_common is not None:
                # store data in the collections
                dset_obj_collections[point_dset] = point_values_common
                # info datasets end
                log_stream.info(' -------> Datasets "' + point_dset + '" ... DONE')
            else:
                # info datasets end
                log_stream.info(' -------> Datasets "' + point_dset +
                                '" ... SKIPPED. Common datasets is not defined; all datasets are not defined in the selected time-series')

        # convert dict to dframe object
        dframe_obj = pd.DataFrame(data=dset_obj_collections, index=time_index_tmpl)
        dframe_obj.index.name = 'time'
        dframe_obj = dframe_obj.sort_index()

        # iterate over columns to check if all elements are null or not
        for column_name in dframe_obj.columns:
            column_null = dframe_obj[column_name].isnull().all()
            if column_null:
                dframe_obj = None
                break

        # store dataframe to common obj
        point_obj_collections[point_reference] = dframe_obj

        log_stream.info(' ------> Compose time-series for point "' + point_reference + '" ... DONE')

    return point_obj_collections
# ----------------------------------------------------------------------------------------------------------------------
