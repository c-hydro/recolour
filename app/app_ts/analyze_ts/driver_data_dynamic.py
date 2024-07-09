"""
Class Features

Name:          driver_data_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd
from copy import deepcopy

from lib_data_io_pickle import read_obj, write_obj
from lib_data_io_csv import read_file_csv, write_file_csv
from lib_data_io_json import read_file_json, write_file_json
from lib_data_io_generic import (adjust_data_point, adjust_data_time,
                                 convert_data_to_vars, convert_vars_to_data, convert_vars_to_dict)

from lib_utils_analysis import (apply_time_series_scaling, apply_time_series_filter, apply_time_series_metrics,
                                join_time_series)
from lib_utils_obj import create_dict_from_list
from lib_utils_system import fill_tags2string, make_folder

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DriverData:

    # -------------------------------------------------------------------------------------
    # initialize class
    def __init__(self, time_run, time_range,
                 static_obj, source_dict, ancillary_dict, destination_dict=None,
                 flags_dict=None, template_dict=None, params_dict=None, methods_dict=None, tmp_dict=None):

        self.time_run = time_run
        self.time_range = time_range

        self.time_start = pd.DatetimeIndex([time_range[0], time_range[-1]]).min()
        self.time_end = pd.DatetimeIndex([time_range[0], time_range[-1]]).max()

        self.time_end, self.time_start = time_range[0], time_range[-1]
        self.static_obj = static_obj

        self.source_dict = source_dict
        self.ancillary_dict = ancillary_dict
        self.destination_dict = destination_dict

        self.flags_dict = flags_dict
        self.template_dict = template_dict
        self.params_dict = params_dict
        self.methods_dict = methods_dict
        self.tmp_dict = tmp_dict

        # data tag(s)
        self.file_name_tag, self.folder_name_tag, self.path_name_tag = 'file_name', 'folder_name', 'path_name'
        self.value_min_tag, self.value_max_tag, self.scale_factor_tag = 'value_min', 'value_max', 'scale_factor'
        self.value_nodata_tag = 'value_no_data'
        self.fields_tag, self.format_tag = 'fields', 'format'
        # time tag(s)
        self.time_start_tag, self.time_end_tag = 'time_start', 'time_end'
        self.time_reference_tag, self.time_period_tag = "time_reference", "time_period",
        self.time_freq_tag, self.time_rounding_tag = 'time_frequency', 'time_rounding'
        self.time_range_tag = 'time_range'

        # get reset flags
        self.reset_src = flags_dict['reset_dynamic_source']
        self.reset_dst = flags_dict['reset_dynamic_destination']

        # get static point
        self.point_obj = self.static_obj['point_obj']

        # get source data
        (folder_name_src, file_name_src, fields_src, format_src, scale_factor_src,
         value_min_src, value_max_src, value_no_data_src) = self.get_info_data(self.source_dict)
        file_path_src = os.path.join(folder_name_src, file_name_src)
        # zip source data
        self.dset_obj_src = self.zip_info_data(
            folder_name_src, file_name_src, file_path_src, fields_src, format_src, scale_factor_src,
            value_min_src, value_max_src, value_no_data_src)
        # get source time
        (time_ref_src, time_period_src, time_round_src, time_freq_src,
         time_start_src, time_end_src, time_range_src) = self.get_info_time(self.source_dict)
        # zip source time
        self.time_obj_src = self.zip_info_time(
            time_ref_src, time_period_src, time_round_src, time_freq_src,
            time_start_src, time_end_src, time_range_src)

        # get ancillary data
        (folder_name_anc, file_name_anc, _, _, _, _, _, _) = self.get_info_data(self.ancillary_dict)
        file_path_anc = os.path.join(folder_name_anc, file_name_anc)
        # zip ancillary data
        self.dset_obj_anc = self.zip_info_data(folder_name_anc, file_name_anc, file_path_anc)

        # get destination data
        (folder_name_dst_data, file_name_dst_data, fields_dst_data, format_dst_data, scale_factor_dst_data,
         value_min_dst_data, value_max_dst_data, value_no_data_dst_data) = self.get_info_data(
            self.destination_dict['data'])
        file_path_dst_data = os.path.join(folder_name_dst_data, file_name_dst_data)
        # zip destination data
        self.dset_obj_dst_data = self.zip_info_data(
            folder_name_dst_data, file_name_dst_data, file_path_dst_data,
            fields_dst_data, format_dst_data, scale_factor_dst_data,
            value_min_dst_data, value_max_dst_data, value_no_data_dst_data)

        # get destination metrics
        (folder_name_dst_met, file_name_dst_met, _, _, _, _, _, _) = self.get_info_data(
            self.destination_dict['metrics'])
        file_path_dst_met = os.path.join(folder_name_dst_met, file_name_dst_met)
        # zip destination metrics
        self.dset_obj_dst_met = self.zip_info_data(
            folder_name_dst_met, file_name_dst_met, file_path_dst_met)

        # get tmp data
        (folder_name_tmp, file_name_tmp, _, _, _, _, _, _) = self.get_info_data(self.tmp_dict)
        # zip tmp data
        self.dset_obj_tmp = self.zip_info_data(folder_name_tmp, file_name_tmp, None)

        # time-series name(s)
        self.ts_name_ref = 'ref'
        ts_name_tmp = list(fields_src.keys())

        if self.ts_name_ref not in ts_name_tmp:
            log_stream.error(' ===> Datasets reference tag "' + self.ts_name_ref +
                             '" must be in the source datasets')
            raise RuntimeError('Define the datasets reference tag in the source datasets')
        ts_name_other = []
        for ts_name_step in ts_name_tmp:
            if 'other' in ts_name_step:
                ts_name_other.append(ts_name_step)
        self.ts_name_other = sorted(ts_name_other)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to zip info time
    def zip_info_time(self, time_ref, time_period, time_round, time_freq, time_start, time_end, time_range):
        time_obj = {self.time_period_tag: time_period, self.time_reference_tag: time_ref,
                    self.time_freq_tag: time_freq, self.time_rounding_tag: time_round,
                    self.time_start_tag: time_start, self.time_end_tag: time_end,
                    self.time_range_tag: time_range}
        return time_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to zip info data
    def zip_info_data(self, folder_name, file_name, path_name, fields=None, file_format=None, scale_factor=None,
                      value_min=None, value_max=None, value_no_data=None):
        info_obj = {self.folder_name_tag: folder_name, self.file_name_tag: file_name, self.path_name_tag: path_name,
                    self.fields_tag: fields, self.format_tag: file_format, self.scale_factor_tag: scale_factor,
                    self.value_min_tag: value_min, self.value_max_tag: value_max, self.value_nodata_tag: value_no_data}
        return info_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get info data
    def get_info_data(self, obj_data):

        folder_name = check_key_of_obj(self.folder_name_tag, obj_data, value_data_default=None)
        file_name = check_key_of_obj(self.file_name_tag, obj_data, value_data_default=None)
        fields = check_key_of_obj(self.fields_tag, obj_data, value_data_default={})
        scale_factor = check_key_of_obj(self.scale_factor_tag, obj_data, value_data_default=1)
        file_format = check_key_of_obj(self.format_tag, obj_data, value_data_default=None)
        v_min = check_key_of_obj(self.value_min_tag, obj_data, value_data_default=None)
        v_max = check_key_of_obj(self.value_max_tag, obj_data, value_data_default=None)
        v_no_data = check_key_of_obj(self.value_nodata_tag, obj_data, value_data_default=np.nan)

        if v_no_data is None:
            v_no_data = np.nan

        return folder_name, file_name, fields, file_format, scale_factor, v_min, v_max, v_no_data
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get info time
    def get_info_time(self, obj_data):

        time_reference = check_key_of_obj(self.time_reference_tag, obj_data, value_data_default=self.time_run)
        time_period = check_key_of_obj(self.time_period_tag, obj_data, value_data_default=None)
        time_rounding = check_key_of_obj(self.time_rounding_tag, obj_data, value_data_default=None)
        time_freq = check_key_of_obj(self.time_freq_tag, obj_data, value_data_default=None)
        time_start = check_key_of_obj(self.time_start_tag, obj_data, value_data_default=None)
        time_end = check_key_of_obj(self.time_end_tag, obj_data, value_data_default=None)

        if time_start is None and time_end is None:
            if time_period is not None and time_freq is not None:
                if time_period == 'PERIOD_BY_TIME_SETTINGS':
                    time_range_by_set = self.time_range
                    time_start_by_set, time_end_by_set = time_range_by_set[0], time_range_by_set[-1]
                    if time_start_by_set > time_end_by_set:
                        time_end_by_set, time_start_by_set = time_range_by_set[0], time_range_by_set[-1]
                    time_range = pd.date_range(start=time_start_by_set, end=time_end_by_set, freq=time_freq)
                    time_start, time_end = time_range[0], time_range[-1]
                else:
                    time_range = pd.date_range(end=time_reference, periods=time_period, freq=time_freq)
                    time_start, time_end = time_range[0], time_range[-1]
            else:
                log_stream.error(' ===> The variables "time_period" and "time_frequency" are both undefined')
                raise RuntimeError('The variables "time_period" and "time_frequency" must be defined')
        elif time_start is not None and time_end is not None:
            time_start, time_end = pd.Timestamp(time_start), pd.Timestamp(time_end)
            time_range = pd.date_range(start=time_start, end=time_end, freq=time_freq)
        else:
            log_stream.error(' ===> The variables "time_start" and "time_end" must be both defined or undefined')
            raise NotImplemented('Case not implemented yet')

        return time_reference, time_period, time_rounding, time_freq, time_start, time_end, time_range

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to define file name
    def define_file_name(self, file_path_raw,
                         time_step=None, time_start=None, time_end=None,
                         point_name=None, point_code=None, point_tag=None):

        if time_step is None:
            time_step = deepcopy(time_end)

        template_time_generic = self.template_dict['time']
        template_datasets_generic = self.template_dict['datasets']

        generic_time_list = list(template_time_generic.keys())
        values_time_dict_generic = create_dict_from_list(generic_time_list, time_step)

        if time_end is not None:
            time_end_list = [key for key, value in values_time_dict_generic.items() if 'end' in key.lower()]
            values_dict_tmp = create_dict_from_list(time_end_list, time_end)
            values_time_dict_generic.update(**values_dict_tmp)

        if time_start is not None:
            time_start_list = [key for key, value in values_time_dict_generic.items() if 'start' in key.lower()]
            values_dict_tmp = create_dict_from_list(time_start_list, time_start)
            values_time_dict_generic.update(**values_dict_tmp)

        # generic_datasets_list = list(template_datasets_generic.keys())
        values_datasets_dict_generic = {}
        if point_name is not None:
            values_datasets_dict_generic['point_name'] = point_name
        if point_code is not None:
            values_datasets_dict_generic['point_code'] = str(point_code)
        if point_tag is not None:
            values_datasets_dict_generic['point_tag'] = point_tag

        file_path_def = fill_tags2string(file_path_raw, template_time_generic, values_time_dict_generic)[0]
        file_path_def = fill_tags2string(file_path_def, template_datasets_generic, values_datasets_dict_generic)[0]

        return file_path_def

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump data
    def dump_data(self, point_data_collections):

        # method start info
        log_stream.info(' ----> Dump time-series object(s) ... ')

        # get time(s)
        time_run = self.time_run
        # get static object
        point_registry = self.point_obj

        # get path(s)
        file_path_dst_data_raw = self.dset_obj_dst_data[self.path_name_tag]
        file_path_dst_met_raw = self.dset_obj_dst_met[self.path_name_tag]
        # get format(s)
        file_format_dst = self.dset_obj_dst_data[self.format_tag]

        # get flag(s)
        reset_dst = self.reset_dst

        # iterate over point(s)
        for point_fields in point_registry.to_dict(orient="records"):

            # get point information
            point_tag, point_name, point_code = point_fields['tag'], point_fields['name'], point_fields['code']
            # point_geo_x, point_geo_y = point_fields['longitude'], point_fields['latitude']

            # point info start
            log_stream.info(' -----> Point "' + point_tag + '" ... ')

            # define destination filename(s)
            file_path_dst_data_point = self.define_file_name(
                file_path_dst_data_raw,
                time_step=time_run, time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                point_name=point_tag, point_code=point_code, point_tag=point_tag)
            file_path_dst_met_point = self.define_file_name(
                file_path_dst_met_raw,
                time_step=time_run, time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            # check file destination reset
            if reset_dst:
                if os.path.exists(file_path_dst_data_point):
                    os.remove(file_path_dst_data_point)
                if os.path.exists(file_path_dst_met_point):
                    os.remove(file_path_dst_met_point)

            # check file destination availability
            if (not os.path.exists(file_path_dst_data_point)) or (not os.path.exists(file_path_dst_met_point)):

                # get destination data
                point_dframe_raw = point_data_collections[point_tag]

                # check destination data
                if point_dframe_raw is not None:

                    # method to convert data
                    point_dframe_converted, point_attrs_converted = convert_vars_to_data(
                        point_dframe_raw, obj_fields=self.dset_obj_src['fields'])

                    # check destination file format
                    if file_format_dst == 'csv':

                        # create datasets folder
                        folder_name_data_dst, file_name_data_dst = os.path.split(file_path_dst_data_point)
                        make_folder(folder_name_data_dst)
                        folder_name_met_dst, file_name_met_dst = os.path.split(file_path_dst_met_point)
                        make_folder(folder_name_met_dst)

                        # write data in csv format
                        write_file_csv(
                            file_path_dst_data_point, point_dframe_converted,
                            dframe_index_format='%Y-%m-%d %H:%M',
                            dframe_sep=',', dframe_decimal='.', dframe_float_format='%.2f',
                            dframe_index=True, dframe_header=True, dframe_index_label='time')
                        # write metrics in csv format
                        write_file_csv(
                            file_path_dst_met_point, point_attrs_converted,
                            dframe_index_format=None,
                            dframe_sep=',', dframe_decimal='.', dframe_float_format='%.2f',
                            dframe_index=True, dframe_header=True, dframe_index_label='metrics')

                    elif file_format_dst == 'json':

                        # convert dynamic information to dictionary format
                        point_datasets = convert_vars_to_dict(point_dframe_converted)

                        # merge point information
                        point_collections = {**point_fields, **point_datasets}

                        # create datasets folder
                        folder_name_data_dst, file_name_data_dst = os.path.split(file_path_dst_data_point)
                        make_folder(folder_name_data_dst)
                        folder_name_met_dst, file_name_met_dst = os.path.split(file_path_dst_met_point)
                        make_folder(folder_name_met_dst)

                        # write data in json format
                        write_file_json(file_path_dst_data_point, point_collections)
                        # write metrics in csv format
                        write_file_csv(
                            file_path_dst_met_point, point_attrs_converted,
                            dframe_index_format=None,
                            dframe_sep=',', dframe_decimal='.', dframe_float_format='%.2f',
                            dframe_index=True, dframe_header=True, dframe_index_label='metrics')

                    else:
                        log_stream.error(' ===> File data format is not supported')
                        raise NotImplemented('Case not implemented yet')

                    # point step end
                    log_stream.info(' -----> Point "' + point_tag + '"  ... DONE')

                else:
                    # point step end
                    log_stream.info(' -----> Point "' + point_tag + '"  ... SKIPPED. Datasets is defined by NoneType')

            else:
                # point step end
                log_stream.info(' -----> Point "' + point_tag + '"  ... SKIPPED. File(s) already exist(s)')

        # method end info
        log_stream.info(' ----> Dump time-series object(s) ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to analyze data
    def analyze_data(self, point_data_collections):

        # method start info
        log_stream.info(' ----> Analyze time-series object(s) ... ')

        # get method(s) object
        point_methods = self.methods_dict
        # get static object
        point_registry = self.point_obj

        # get time-series name(s)
        ts_name_ref, ts_name_other = self.ts_name_ref, self.ts_name_other

        # iterate over point(s)
        point_analysis_collections = {}
        for point_fields in point_registry.to_dict(orient="records"):

            # debug
            #point_fields = point_registry.to_dict(orient="records")[8]

            # get point information
            point_tag, point_name, point_code = point_fields['tag'], point_fields['name'], point_fields['code']
            # point_geo_x, point_geo_y = point_fields['longitude'], point_fields['latitude']

            # point info start
            log_stream.info(' -----> Point "' + point_tag + '" ... ')

            # get point data
            point_data_src = point_data_collections[point_tag]

            # get ts reference
            ts_ref = point_data_src[ts_name_ref]
            # check if ts reference is not None
            if ts_ref is not None:

                # remove nan(s)
                ts_ref = ts_ref.dropna()
                time_ref = ts_ref.index

                # check if ts reference is not empty
                if not ts_ref.empty:

                    # iterate over time-series other
                    ts_collections = {ts_name_ref: ts_ref}
                    for ts_name_step in ts_name_other:

                        # get point methods
                        ts_methods = point_methods[ts_name_step]
                        ts_other_src = point_data_src[ts_name_step]

                        # check if ts other is not None
                        if ts_other_src is not None:

                            # apply filter to time-series
                            if 'filter' in list(ts_methods.keys()):
                                # get parameters
                                filter_method_mode = ts_methods['filter']['mode']
                                filter_method_type = ts_methods['filter']['type']
                                filter_method_window = ts_methods['filter']['window']
                                # active method
                                if filter_method_mode:
                                    # method to filter time-series
                                    if filter_method_type is not None:
                                        ts_other_tmp = apply_time_series_filter(
                                            ts_other_src,
                                            ts_filter_type=filter_method_type, ts_filter_window=filter_method_window)
                                    else:
                                        ts_other_tmp = deepcopy(ts_other_src)
                                else:
                                    ts_other_tmp = deepcopy(ts_other_src)
                            else:
                                ts_other_tmp = deepcopy(ts_other_src)

                            # apply scaling to time-series
                            if 'scale' in list(ts_methods.keys()):
                                # get parameters
                                scale_method_mode = ts_methods['scale']['mode']
                                scale_method_type = ts_methods['scale']['type']
                                # active method
                                if scale_method_mode:

                                    # method to scale time-series
                                    if scale_method_type is not None:
                                        ts_other_scaled = apply_time_series_scaling(
                                            ts_ref, ts_other_tmp,
                                            ts_scale_method='cdf_beta_match')
                                    else:
                                        ts_other_scaled = None
                                else:
                                    ts_other_scaled = None
                            else:
                                ts_other_scaled = deepcopy(ts_other_src)

                            # compute time-series metrics
                            if ts_other_scaled is not None:
                                ts_metrics = apply_time_series_metrics(ts_ref, ts_other_scaled)
                                ts_other_scaled.attrs = ts_metrics
                            else:
                                log_stream.warning(' ===> Time series scaling for "' + ts_name_step +
                                                   '" is not applied due to NoneType object')

                            # collect time-series objects
                            if ts_other_scaled is not None:
                                ts_collections[ts_name_step] = ts_other_scaled
                            else:

                                # create null time-series
                                data_null = [np.nan] * time_ref.__len__()
                                obj_null = {'other': data_null}
                                ts_null = pd.Series(data=obj_null, index=time_ref)
                                ts_null.name = 'other'

                                ts_collections[ts_name_step] = ts_null
                                log_stream.warning(' ===> Time series selected collection for "' + ts_name_step +
                                                   '" is defined by NoneType object. Initialize with null series')
                        else:

                            # create null time-series
                            data_null = [np.nan] * time_ref.__len__()
                            obj_null = {'other': data_null}
                            ts_null = pd.Series(data=obj_null, index=time_ref)

                            # set time-series as None
                            ts_collections[ts_name_step] = ts_null

                            log_stream.warning(' ===> Time series raw collection for "' + ts_name_step +
                                               '" is defined by NoneType object. Initialize with null series')

                    # join time-series
                    ts_dframe = join_time_series(
                        ts_collections, ts_name_ref=ts_name_ref,
                        time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                        time_freq=self.time_obj_src['time_frequency'])
                    # store time-series
                    point_analysis_collections[point_tag] = ts_dframe

                else:
                    # store time-series as None (if ref datasets is empty)
                    log_stream.warning(' ===> Time series reference for point "' + point_tag +
                                       '" is empty. Time series collections are not analyzed')
                    point_analysis_collections[point_tag] = None

            # point info start
            log_stream.info(' -----> Point "' + point_tag + '" ... DONE')

        # method end info
        log_stream.info(' ----> Analyze time-series object(s) ... DONE.')

        return point_analysis_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize time-series object(s) ... ')

        # get time(s)
        time_run = self.time_run
        # get static object
        point_registry = self.point_obj

        # get path(s)
        file_path_src_raw = self.dset_obj_src[self.path_name_tag]
        file_path_anc_raw = self.dset_obj_anc[self.path_name_tag]
        file_path_dst_data_raw = self.dset_obj_dst_data[self.path_name_tag]
        file_path_dst_met_raw = self.dset_obj_dst_met[self.path_name_tag]
        # get format(s)
        file_format_src = self.dset_obj_src[self.format_tag]

        # get flag(s)
        reset_src, reset_dst = self.reset_src, self.reset_dst

        # iterate over point(s)
        point_data_collections = {}
        for point_fields in point_registry.to_dict(orient="records"):

            # debug
            # point_fields = point_registry.to_dict(orient="records")[8]

            # get point information
            point_tag, point_name, point_code = point_fields['tag'], point_fields['name'], point_fields['code']
            # point_geo_x, point_geo_y = point_fields['longitude'], point_fields['latitude']

            # point info start
            log_stream.info(' -----> Point "' + point_tag + '" ... ')

            # define source, ancillary and destination filename(s)
            file_path_src_point = self.define_file_name(
                file_path_src_raw,
                time_step=time_run, time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            file_path_anc_point = self.define_file_name(
                file_path_anc_raw,
                time_step=time_run, time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            file_path_dst_data_point = self.define_file_name(
                file_path_dst_data_raw,
                time_step=time_run, time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                point_name=point_tag, point_code=point_code, point_tag=point_tag)
            file_path_dst_met_point = self.define_file_name(
                file_path_dst_met_raw,
                time_step=time_run, time_start=self.time_obj_src['time_start'], time_end=self.time_obj_src['time_end'],
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            # reset source file
            if reset_src:
                if os.path.exists(file_path_anc_point):
                    os.remove(file_path_anc_point)
                if os.path.exists(file_path_dst_data_point):
                    os.remove(file_path_dst_data_point)
                if os.path.exists(file_path_dst_met_point):
                    os.remove(file_path_dst_met_point)

            # check ancillary file availability
            if not os.path.exists(file_path_anc_point):

                # get data start
                log_stream.info(' ------> Get datasets ... ')

                # check source reference file availability
                if os.path.exists(file_path_src_point):

                    # check source file format
                    if file_format_src == 'json':

                        # method to get data
                        point_data_raw = read_file_json(file_path_src_point)
                        # method to convert data
                        point_data_converted = convert_data_to_vars(
                            point_data_raw, obj_fields=self.dset_obj_src['fields'])
                        # method to select data
                        point_data_selected = adjust_data_point(
                            point_data_converted,
                            scale_factor=self.dset_obj_src[self.scale_factor_tag],
                            value_min=self.dset_obj_src[self.value_min_tag],
                            value_max=self.dset_obj_src[self.value_max_tag],
                            value_no_data=self.dset_obj_src[self.value_nodata_tag],
                            time_start=self.time_obj_src['time_start'],
                            time_end=self.time_obj_src['time_end'],
                            time_freq=self.time_obj_src['time_frequency'])

                    elif file_format_src == 'csv':

                        # method to get data
                        point_data_raw = read_file_csv(file_path_src_point, dframe_sep=';')
                        # method to convert data
                        point_data_converted = convert_data_to_vars(
                            point_data_raw, obj_fields=self.dset_obj_src['fields'])

                        # method to adjust data time
                        time_start_src, time_end_src, time_frequency_src = adjust_data_time(
                            point_data_converted,
                            time_start=self.time_obj_src['time_start'],
                            time_end=self.time_obj_src['time_end'])

                        # method to adjust data values
                        point_data_selected = adjust_data_point(
                            point_data_converted,
                            scale_factor=self.dset_obj_src[self.scale_factor_tag],
                            value_min=self.dset_obj_src[self.value_min_tag],
                            value_max=self.dset_obj_src[self.value_max_tag],
                            value_no_data=self.dset_obj_src[self.value_nodata_tag],
                            time_start=time_start_src,
                            time_end=time_end_src,
                            time_freq=time_frequency_src)

                    else:
                        log_stream.error(' ===> Source data type "' + file_format_src + '" is not supported.')
                        raise NotImplemented('Case not implemented yet')

                else:
                    log_stream.warning(' ===> Datasets file "' + file_path_src_point + '" was not available.')
                    point_data_selected = None

                # get data end
                log_stream.info(' ------> Get datasets ... DONE')

                # dump data start
                log_stream.info(' ------> Dump datasets ... ')
                # check data availability
                if point_data_selected is not None:

                    # method to dump data
                    folder_name_anc_point, file_name_anc_point = os.path.split(file_path_anc_point)
                    make_folder(folder_name_anc_point)
                    write_obj(file_path_anc_point, point_data_selected)

                    # dump data end
                    log_stream.info(' ------> Dump datasets ... DONE')

                else:
                    # dump data end
                    log_stream.info(' ------> Dump datasets ... SKIPPED. Datasets is not available')

                # point info end
                log_stream.info(' -----> Point "' + point_tag + '" ... DONE')

            else:

                # point info end
                point_data_selected = read_obj(file_path_anc_point)
                log_stream.info(' -----> Point "' + point_tag + '" ... SKIPPED. Data previously saved')

            point_data_collections[point_tag] = point_data_selected

        # method end info
        log_stream.info(' ----> Organize time-series object(s) ... DONE')

        return point_data_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to check key in a dictionary obj
def check_key_of_obj(key_data, obj_data, value_data_default=None):
    if key_data in list(obj_data.keys()):
        return obj_data[key_data]
    else:
        return value_data_default
# -------------------------------------------------------------------------------------
