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

from lib_data_io_csv import read_file_csv
from lib_data_io_pickle import read_obj, write_obj
from lib_data_io_json import write_file_json
from lib_data_io_generic import (adjust_data_point, adjust_data_time, range_data_point, combine_data_point_by_time,
                                 join_data_point, resample_data_point, fill_data_point,
                                 convert_registry_point_to_dict, convert_datasets_point_to_dict)
from lib_data_io_csv import write_file_csv

from lib_utils_obj import create_dict_from_list, remap_dict, remap_dframe
from lib_utils_system import fill_tags2string, make_folder

from lib_info_args import logger_name, time_format_algorithm

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
                 flags_dict=None, template_dict=None, params_dict=None, tmp_dict=None):

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
        self.tmp_dict = tmp_dict

        # define datasets tags and order them starting from reference tag
        self.dset_ref_tag = 'ref'
        dset_list_tmp = list(source_dict.keys())
        if self.dset_ref_tag not in dset_list_tmp:
            log_stream.error(' ===> Datasets reference tag "' + self.dset_ref_tag +
                             '" must be in the source datasets')
            raise RuntimeError('Define the datasets reference tag in the source datasets')
        dset_list_tmp.remove(self.dset_ref_tag)
        dset_list_tmp.sort()
        self.dset_list_tag = [self.dset_ref_tag] + dset_list_tmp

        # data tag(s)
        self.file_name_tag, self.folder_name_tag, self.path_name_tag = 'file_name', 'folder_name', 'path_name'
        self.value_min_tag, self.value_max_tag, self.scale_factor_tag = 'value_min', 'value_max', 'scale_factor'
        self.value_nodata_tag = 'value_no_data'
        self.fields_tag, self.format_tag, self.delimiter_tag = 'fields', 'format', 'delimiter'
        self.date_format_tag, self.decimal_precision_tag = 'date_format', 'decimal_precision'

        # time tag(s)
        self.time_period_tag, self.time_reference_tag = "time_period", 'time_reference'
        self.time_freq_tag, self.time_rounding_tag = 'time_frequency', 'time_rounding'
        self.time_start_tag, self.time_end_tag = 'time_start', 'time_end'

        # get reset flags
        self.reset_src = flags_dict['reset_dynamic_source']
        self.reset_dst = flags_dict['reset_dynamic_destination']

        # get parameter(s)
        self.par_resample_time_freq = self.params_dict['resample_time_frequency']
        self.par_resample_time_method = self.params_dict['resample_time_method']
        self.par_fill_time_method = self.params_dict['fill_time_method']
        self.par_fill_time_order = self.params_dict['fill_time_order']
        self.par_fill_time_limit = self.params_dict['fill_time_limit']
        self.par_fill_time_dir = self.params_dict['fill_time_direction']

        # get static point
        self.point_obj = self.static_obj['point_obj']

        # source object(s) iterating over datasets
        self.dset_obj_src, self.time_obj_src = {}, {}
        for dset_name in self.dset_list_tag:

            # get info data
            (folder_name_src, file_name_src, fields_src, format_src, scale_factor_src,
             value_min_src, value_max_src, value_no_data_src,
             delimiter_src, decimal_precision_src, date_format_src) = self.get_info_data(
                self.source_dict[dset_name])
            file_path_src = os.path.join(folder_name_src, file_name_src)
            # zip info data
            dset_obj_step = self.zip_info_data(folder_name_src, file_name_src, file_path_src,
                                               fields_src, format_src, scale_factor_src,
                                               value_min_src, value_max_src, value_no_data_src,
                                               delimiter_src, decimal_precision_src, date_format_src)
            # get info time
            (time_ref_src, time_period_src, time_round_src, time_freq_src,
             time_start_src, time_end_src) = self.get_info_time(self.source_dict[dset_name])
            # zip info time
            time_obj_step = self.zip_info_time(time_ref_src, time_period_src, time_round_src, time_freq_src,
                                               time_start_src, time_end_src)

            # store info data and time
            self.dset_obj_src[dset_name] = dset_obj_step
            self.time_obj_src[dset_name] = time_obj_step

        # ancillary object(s)
        (self.folder_name_anc, self.file_name_anc, _, _, _, _, _, _, _, _, _) = self.get_info_data(
            self.ancillary_dict)
        self.file_path_anc = os.path.join(self.folder_name_anc, self.file_name_anc)

        # destination object(s)
        (self.folder_name_dst, self.file_name_dst,
         self.fields_dst, self.format_dst, self.scale_factor_dst,
         self.value_min_dst, self.value_max_dst, self.value_no_data_dst,
         self.delimiter_dst, self.decimal_precision_dst, self.date_format_dst) = self.get_info_data(
            self.destination_dict)
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

        # tmp object(s)
        (self.folder_name_tmp, self.file_name_tmp, _, _, _, _, _, _, _, _, _) = self.get_info_data(self.tmp_dict)

        # check fields from source to destination
        for dset_name in self.dset_list_tag:
            if dset_name not in list(self.fields_dst.keys()):
                log_stream.warning(' ===> Datasets name "' + dset_name + '" is not defined in the destination object')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to zip info time
    def zip_info_time(self, time_ref, time_period, time_round, time_freq, time_start, time_end):
        time_obj = {self.time_period_tag: time_period, self.time_reference_tag: time_ref,
                    self.time_freq_tag: time_freq, self.time_rounding_tag: time_round,
                    self.time_start_tag: time_start, self.time_end_tag: time_end}
        return time_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to zip info data
    def zip_info_data(self, folder_name, file_name, path_name, fields, format, scale_factor,
                      value_min, value_max, value_no_data,
                      delimiter_data, decimal_precision, date_format):
        info_obj = {self.folder_name_tag: folder_name, self.file_name_tag: file_name, self.path_name_tag: path_name,
                    self.fields_tag: fields, self.format_tag: format, self.scale_factor_tag: scale_factor,
                    self.value_min_tag: value_min, self.value_max_tag: value_max,
                    self.value_nodata_tag: value_no_data, self.delimiter_tag: delimiter_data,
                    self.decimal_precision_tag: decimal_precision, self.date_format_tag: date_format}
        return info_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get info data
    def get_info_data(self, obj_data):

        folder_name = check_key_of_obj(self.folder_name_tag, obj_data, value_data_default=None)
        file_name = check_key_of_obj(self.file_name_tag, obj_data, value_data_default=None)
        fields = check_key_of_obj(self.fields_tag, obj_data, value_data_default={})
        scale_factor = check_key_of_obj(self.scale_factor_tag, obj_data, value_data_default=1)
        format = check_key_of_obj(self.format_tag, obj_data, value_data_default=None)
        vmin = check_key_of_obj(self.value_min_tag, obj_data, value_data_default=None)
        vmax = check_key_of_obj(self.value_max_tag, obj_data, value_data_default=None)
        vnodata = check_key_of_obj(self.value_nodata_tag, obj_data, value_data_default=np.nan)
        delimiter = check_key_of_obj(self.delimiter_tag, obj_data, value_data_default=',')
        decimal_precision = check_key_of_obj(self.decimal_precision_tag, obj_data, value_data_default=2)
        date_format = check_key_of_obj(self.date_format_tag, obj_data, value_data_default='%Y-%m-%d %H:%M')

        if vnodata is None:
            vnodata = np.nan

        return (folder_name, file_name, fields, format, scale_factor, vmin, vmax, vnodata,
                delimiter, decimal_precision, date_format)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get info time
    def get_info_time(self, obj_data):

        time_reference = check_key_of_obj(self.time_reference_tag, obj_data, value_data_default=None)
        time_period = check_key_of_obj(self.time_period_tag, obj_data, value_data_default=None)
        time_rounding = check_key_of_obj(self.time_rounding_tag, obj_data, value_data_default=None)
        time_freq = check_key_of_obj(self.time_freq_tag, obj_data, value_data_default=None)
        time_start = check_key_of_obj(self.time_start_tag, obj_data, value_data_default=None)
        time_end = check_key_of_obj(self.time_end_tag, obj_data, value_data_default=None)

        if time_start is not None:
            time_start = pd.Timestamp(time_start)
        if time_end is not None:
            time_end = pd.Timestamp(time_end)

        return time_reference, time_period, time_rounding, time_freq, time_start, time_end

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

        generic_datasets_list = list(template_datasets_generic.keys())
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
    def dump_data(self, dframe_point):

        # method start info
        log_stream.info(' ----> Organize destination object(s) ... ')

        # get time(s)
        time_run = self.time_run
        time_range_anls, time_start_anls, time_end_anls = self.time_range, self.time_start, self.time_end
        # get static object
        point_obj = self.point_obj

        # get file path(s)
        file_path_anc_raw, file_path_dst_raw = self.file_path_anc, self.file_path_dst
        # get file settings
        file_value_no_data_dst = self.value_no_data_dst
        file_delimiter_dst = self.delimiter_dst
        file_decimal_precision_dst = self.decimal_precision_dst
        file_date_format_dst = self.date_format_dst

        # iterate over dataframe
        if dframe_point is not None:

            # iterate over dataframe
            for point_key, point_dframe in dframe_point.items():

                # point step start
                log_stream.info(' -----> Point "' + point_key + '" ... ')

                # check point dataset availability (or exit if defined by NoneType)
                if point_dframe is not None:

                    # define destination file
                    file_path_dst_step = self.define_file_name(
                        file_path_dst_raw, time_step=time_run, time_start=time_start_anls, time_end=time_end_anls,
                        point_name=point_key, point_code=None, point_tag=point_key)

                    # define dataframe format for the output file
                    point_dframe = point_dframe.fillna(self.value_no_data_dst)
                    point_dframe = point_dframe.round(decimals=2)

                    # check destination file format
                    if self.format_dst == 'csv':

                        # remap fields name
                        point_dframe = remap_dframe(point_dframe, self.fields_dst)

                        # create datasets folder
                        folder_name_dst, file_name_dst = os.path.split(file_path_dst_step)
                        make_folder(folder_name_dst)

                        # write data in csv format
                        write_file_csv(
                            file_path_dst_step, point_dframe,
                            dframe_sep=file_delimiter_dst, dframe_decimal='.',
                            dframe_no_data=file_value_no_data_dst,
                            dframe_float_format='%.{:}f'.format(str(file_decimal_precision_dst)),
                            dframe_index=True, dframe_header=True,
                            dframe_index_label='time', dframe_index_format=file_date_format_dst)

                    elif self.format_dst == 'json':

                        # convert static information in dictionary format
                        point_registry = convert_registry_point_to_dict(point_key, point_obj, tag_registry='tag')
                        # convert dynamic information to dictionary format
                        point_datasets = convert_datasets_point_to_dict(point_dframe)
                        # remap fields name
                        point_datasets = remap_dict(point_datasets, self.fields_dst)

                        # merge point information
                        point_collections = {**point_registry, **point_datasets}

                        # create datasets folder
                        folder_name_dst, file_name_dst = os.path.split(file_path_dst_step)
                        make_folder(folder_name_dst)

                        # write data in json format
                        write_file_json(file_path_dst_step, point_collections)

                    else:
                        log_stream.error(' ===> File format is not supported')
                        raise NotImplemented('Case not implemented yet')

                    # point step end
                    log_stream.info(' -----> Point "' + point_key + '" ... DONE')

                else:
                    # point step end (no data available)
                    log_stream.info(' -----> Point "' + point_key + '" ... SKIPPED. Datasets are defined by NoneType')

            # method start info
            log_stream.info(' ----> Organize destination object(s) ... DONE')

        else:

            # method end info
            log_stream.info(' ----> Organize destination object(s) ... FAILED. All datasets are defined by NoneType')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to analyze data
    def analyze_data(self, dframe_point_in):

        # method start info
        log_stream.info(' ----> Analyze object(s) ... ')

        # get time(s)
        time_run = self.time_run
        time_range_anls = self.time_range
        time_start_anls, time_end_anls = self.time_start, self.time_end
        # get static object
        point_obj = self.point_obj

        # get path(s)
        file_path_anc_raw, file_path_dst_raw = self.file_path_anc, self.file_path_dst

        # define ancillary file
        file_path_anc_step = self.define_file_name(
            file_path_anc_raw, time_step=time_run, time_start=time_start_anls, time_end=time_end_anls)

        # iterate over dataframe
        if dframe_point_in is not None:

            # iterate over dataframe
            dframe_point_out = {}
            for point_key, point_dframe_tmp in dframe_point_in.items():

                # point step start
                log_stream.info(' -----> Point "' + point_key + '"  ... ')

                # check dataframe
                if point_dframe_tmp is not None:

                    # method to resample time-series
                    point_dframe_resampled = resample_data_point(
                        deepcopy(point_dframe_tmp),
                        resample_frequency=self.par_resample_time_freq, resample_method=self.par_resample_time_method)

                    # testing pre resampling
                    point_dframe_pre = point_dframe_resampled.dropna(how='all')

                    # method to fill time-series
                    point_dframe_filled = fill_data_point(
                        deepcopy(point_dframe_resampled),
                        fill_method=self.par_fill_time_method,
                        fill_order=self.par_fill_time_order,
                        fill_limit=self.par_fill_time_limit,
                        fill_direction=self.par_fill_time_dir)

                    # testing post resampling
                    point_dframe_post = point_dframe_filled.dropna(how='all')

                    # save in a common obj
                    dframe_point_out[point_key] = point_dframe_filled

                    # point step end
                    log_stream.info(' -----> Point "' + point_key + '"  ... DONE')

                else:
                    # point step end
                    log_stream.info(' -----> Point "' + point_key + '"  ... SKIPPED. Datasets are defined by NoneType')
                    # save in a common obj
                    dframe_point_out[point_key] = None

            # method end info
            log_stream.info(' ----> Analyze object(s) ... DONE.')

        else:
            # method end info
            dframe_point_out = None
            log_stream.info(' ----> Analyze object(s) ... FAILED. All datasets are not available')

        return dframe_point_out

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize source object(s) ... ')

        # get time(s)
        time_run = self.time_run
        time_range_anls = self.time_range
        time_start_anls, time_end_anls = self.time_start, self.time_end
        # get static object
        point_obj_registry = self.point_obj

        # get path(s)
        file_path_anc_raw = self.file_path_anc
        file_path_dst_raw = self.file_path_dst

        # get flag(s)
        reset_src, reset_dst = self.reset_src, self.reset_dst

        # define ancillary and destination files
        file_path_anc_step = self.define_file_name(
            file_path_anc_raw, time_step=time_run, time_start=time_start_anls, time_end=time_end_anls)
        file_path_dst_step = self.define_file_name(
            file_path_dst_raw, time_step=time_run, time_start=time_start_anls, time_end=time_end_anls)

        # reset source file
        if reset_src:
            if os.path.exists(file_path_anc_step):
                os.remove(file_path_anc_step)
            if os.path.exists(file_path_dst_step):
                os.remove(file_path_dst_step)

        # get data start
        log_stream.info(' -----> Get datasets ... ')

        # check destination file availability
        if not os.path.exists(file_path_anc_step):

            # iterate over datasets
            point_obj_collections_raw = {}
            for ((dset_name, dset_fields),
                 (time_name, time_fields)) in zip(self.dset_obj_src.items(), self.time_obj_src.items()):

                # define datasets tag
                assert dset_name == time_name, 'datasets and time names must be the same'
                dset_tag = list({dset_name, time_name})[0]

                # datasets start
                log_stream.info(' ------> Datasets "' + dset_tag + '" ... ')

                file_path_src_raw = dset_fields[self.path_name_tag]
                format_src, scale_factor_src = dset_fields[self.format_tag], dset_fields[self.scale_factor_tag]
                value_min_src, value_max_src = dset_fields[self.value_min_tag], dset_fields[self.value_max_tag]
                value_nodata_src, delimiter_src = dset_fields[self.value_nodata_tag], dset_fields[self.delimiter_tag]
                time_start_src, time_end_src = time_fields[self.time_start_tag], time_fields[self.time_end_tag]

                # define source file name
                file_path_src_step = self.define_file_name(
                    file_path_src_raw,
                    time_step=time_run, time_start=time_start_src, time_end=time_start_src)

                # check source reference file availability
                if os.path.exists(file_path_src_step):

                    # check source file format
                    if format_src == 'csv':

                        # method to get data
                        point_obj_data_raw = read_file_csv(file_path_src_step, dframe_sep=delimiter_src)

                        # method to adjust data time
                        time_start_src, time_end_src, time_frequency_src = adjust_data_time(
                            point_obj_data_raw, time_start=time_start_src, time_end=time_end_src)

                        # method to adjust data values
                        point_obj_data_sel = adjust_data_point(
                            point_obj_data_raw,
                            scale_factor=scale_factor_src,
                            value_min=value_min_src, value_max=value_max_src, value_no_data=value_nodata_src,
                            time_start=time_start_src, time_end=time_end_src, time_freq=time_frequency_src)
                    else:
                        log_stream.error(' ===> Source data type "' + format_src + '" is not supported.')
                        raise NotImplemented('Case not implemented yet')

                else:
                    log_stream.warning(' ===> Datasets file "' + file_path_src_step + '" was not available.')
                    point_obj_data_sel = None

                # collect data to a collections obj
                point_obj_collections_raw[dset_tag] = point_obj_data_sel

                # datasets start
                log_stream.info(' ------> Datasets "' + dset_tag + '" ... DONE')

            # get data end
            log_stream.info(' -----> Get datasets ... DONE')

            # check data start (ref and other_n)
            log_stream.info(' -----> Check datasets ... ')

            point_obj_check_ref, point_obj_check_other = False, False
            if point_obj_collections_raw is not None:

                if 'ref' in list(point_obj_collections_raw):
                    point_ref = point_obj_collections_raw['ref']
                    if point_ref is None:
                        log_stream.warning(' ===> Datasets reference is defined by NoneType. Procedure will be skipped')
                    else:
                        point_obj_check_ref = True
                else:
                    log_stream.warning(' ===> Datasets reference is not defined available. Procedure will be skipped')

                for key, value in list(point_obj_collections_raw.items()):
                    if key != 'ref':
                        point_data_other = point_obj_collections_raw[key]
                        if point_data_other is not None:
                            point_obj_check_other = True
                            break
                if not point_obj_check_other:
                    log_stream.warning(' ===> Datasets other_n are not defined available. Procedure will be skipped')
            else:
                log_stream.warning(' ===> All datasets are defined by NoneType. Procedure will be skipped')

            # check data start (ref and other_n)
            log_stream.info(' -----> Check datasets ... DONE')

            # apply check conditions of ref and other_n datasets
            if point_obj_check_ref and point_obj_check_other:

                # combine data start
                log_stream.info(' -----> Combine datasets by time ... ')
                # method to range data point
                time_frequency_expected, time_start_expected, time_end_expected, time_range_expected = range_data_point(
                    point_obj_collections_raw, time_run_reference=self.time_run,
                    time_start_reference=self.time_start, time_end_reference=self.time_end)
                # method to combine data point to the expected time range
                point_obj_collections_combined = combine_data_point_by_time(
                    point_obj_collections_raw, point_obj_registry,
                    time_start_expected=time_start_expected, time_end_expected=time_end_expected,
                    time_frequency_expected=time_frequency_expected, time_reverse=True)
                # combine data end
                log_stream.info(' -----> Combine datasets by time ... DONE')

                # join data start
                log_stream.info(' -----> Join datasets ... ')
                # method to join data point(s)
                point_obj_collections_joined = join_data_point(
                    time_range_expected, point_obj_collections_combined, point_obj_registry)
                # join data end
                log_stream.info(' -----> Join datasets ... DONE')

                # dump data start
                log_stream.info(' -----> Dump datasets ... ')
                # check data availability
                if point_obj_collections_joined is not None:

                    # method to dump data
                    folder_name_anc, file_name_anc = os.path.split(file_path_anc_step)
                    make_folder(folder_name_anc)
                    write_obj(file_path_anc_step, point_obj_collections_joined)

                    # dump data end
                    log_stream.info(' -----> Dump datasets ... DONE')

                else:
                    # dump data end (joined datasets are not available)
                    point_obj_collections_joined = None
                    log_stream.info(' -----> Dump datasets ... SKIPPED. Datasets is not available')

            else:
                # dump data end (reference or/and all other_kn datasets are defined by NoneType
                point_obj_collections_joined = None
                log_stream.info(' -----> Dump datasets ... SKIPPED. '
                                'Reference datasets is defined by NoneType or/and all '
                                'other_kn datasets are defined by NoneType')

        else:
            # read ancillary file
            point_obj_collections_joined = read_obj(file_path_anc_step)

            # get data end
            log_stream.info(' -----> Get datasets ... SKIPPED. Datasets were previously saved')

        # method end info
        log_stream.info(' ----> Organize source object(s) ... DONE')

        return point_obj_collections_joined

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
