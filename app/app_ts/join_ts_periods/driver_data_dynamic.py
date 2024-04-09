"""
Class Features

Name:          driver_data_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240301'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd
from copy import deepcopy

from lib_data_io_generic import join_data_point, range_data_point, combine_data_point_by_time
from lib_data_io_csv import read_file_csv, write_file_csv
from lib_data_io_pickle import read_obj, write_obj

from lib_utils_obj import create_dict_from_list, convert_dict_to_dframe
from lib_utils_system import fill_tags2string

from lib_info_args import logger_name, time_var_name

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt

time_format_algorithm = '%Y-%m'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class driver data
class DriverData:

    # ------------------------------------------------------------------------------------------------------------------
    # initialize class
    def __init__(self, time_reference, time_obj, static_obj,
                 source_dict, ancillary_dict, destination_dict,
                 flags_dict=None, template_dict=None,
                 tmp_dict=None):

        self.time_reference = pd.Timestamp(time_reference)
        self.time_obj = time_obj
        self.time_end, self.time_start = time_obj[0], time_obj[-1]

        self.static_obj = static_obj

        self.source_dict = source_dict
        self.ancillary_dict = ancillary_dict
        self.destination_dict = destination_dict

        self.flags_dict = flags_dict
        self.template_dict = template_dict
        self.tmp_dict = tmp_dict

        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.type_tag, self.no_data_tag, self.delimiter_tag = 'type', 'no_data', 'delimiter'
        self.date_format_tag, self.decimal_precision_tag = 'date_format', 'decimal_precision'

        self.reset_data_src = flags_dict['reset_dynamic_source']
        self.reset_data_dst = flags_dict['reset_dynamic_destination']

        self.point_static_collections = self.static_obj['point_obj']

        # source object(s)
        folder_name_src = self.source_dict[self.folder_name_tag]
        file_name_src = self.source_dict[self.file_name_tag]
        self.file_path_src = os.path.join(folder_name_src, file_name_src)
        self.type_src = self.source_dict[self.type_tag]
        # ancillary object(s)
        folder_name_anc = ancillary_dict[self.folder_name_tag]
        file_name_anc = ancillary_dict[self.file_name_tag]
        self.file_path_anc = os.path.join(folder_name_anc, file_name_anc)
        # destination object(s)
        folder_name_dst = destination_dict[self.folder_name_tag]
        file_name_dst = destination_dict[self.file_name_tag]
        self.file_path_dst = os.path.join(folder_name_dst, file_name_dst)
        self.type_dst = destination_dict[self.type_tag]
        self.date_format_dst = destination_dict[self.date_format_tag]
        self.delimiter_dst = destination_dict[self.delimiter_tag]
        self.decimal_precision_dst = destination_dict[self.decimal_precision_tag]
        self.no_data_dst = destination_dict[self.no_data_tag]

        # tmp object(s)
        self.folder_name_tmp = tmp_dict[self.folder_name_tag]
        self.file_name_tmp = tmp_dict[self.file_name_tag]

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to define file name
    def define_file_name(self, file_path_raw, time_step=None, time_start=None, time_end=None):

        if time_step is None:
            time_step = deepcopy(time_end)

        template_time_generic = self.template_dict

        generic_list = list(template_time_generic.keys())
        values_dict_generic = create_dict_from_list(generic_list, time_step)

        if time_end is not None:
            time_end_list = [key for key, value in values_dict_generic.items() if 'end' in key.lower()]
            values_dict_tmp = create_dict_from_list(time_end_list, time_end)
            values_dict_generic.update(**values_dict_tmp)

        if time_start is not None:
            time_start_list = [key for key, value in values_dict_generic.items() if 'start' in key.lower()]
            values_dict_tmp = create_dict_from_list(time_start_list, time_start)
            values_dict_generic.update(**values_dict_tmp)

        file_path_def = fill_tags2string(file_path_raw, template_time_generic, values_dict_generic)[0]
        return file_path_def

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to dump data
    def dump_data(self, point_collections):

        # method start info
        log_stream.info(' ----> Dump dynamic object(s) ... ')

        # get time(s)
        time_reference = self.time_reference

        # get path(s)
        file_path_dst_raw = self.file_path_dst

        # define source file name
        file_path_dst_def = self.define_file_name(file_path_dst_raw, time_step=time_reference)

        # get flag(s)
        reset_dst = self.reset_data_dst
        # reset destination file
        if reset_dst:
            if os.path.exists(file_path_dst_def):
                os.remove(file_path_dst_def)

        # check destination file availability
        if not os.path.exists(file_path_dst_def):

            # check destination file type
            if self.type_dst == 'csv_2d':

                # convert dictionary collections to dataframe
                point_dframe = convert_dict_to_dframe(point_collections)

                # create folder(s)
                folder_name, file_name = os.path.split(file_path_dst_def)
                os.makedirs(folder_name, exist_ok=True)

                # method to write data
                write_file_csv(file_path_dst_def, point_dframe,
                               dframe_sep=self.delimiter_dst, dframe_decimal='.',
                               dframe_float_format='%.{:}f'.format(self.decimal_precision_dst),
                               dframe_index=True, dframe_header=True,
                               dframe_index_label=time_var_name, dframe_index_format=self.date_format_dst,
                               dframe_no_data=self.no_data_dst)

            else:
                # file type not expected (error message)
                log_stream.error(' ===> File type "' + self.type_dst + "' not expected")
                raise NotImplementedError('File type not implemented yet')

        else:
            # file not found (warning message)
            log_stream.warning(' ===> File "' + file_path_dst_def + "' does not exists")

        # method end info
        log_stream.info(' ----> Dump dynamic object(s) ... DONE')

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize dynamic object(s) ... ')

        # get time(s)
        time_obj = self.time_obj
        # get static
        point_static_collections = self.point_static_collections

        # get path(s)
        file_path_src_raw, file_path_anc_raw = self.file_path_src, self.file_path_anc

        # get flag(s)
        reset_src = self.reset_data_src

        # define ancillary file
        file_path_anc_def = self.define_file_name(file_path_anc_raw, time_start=self.time_start, time_end=self.time_end)
        # reset ancillary file
        if reset_src:
            if os.path.exists(file_path_anc_def):
                os.remove(file_path_anc_def)

        # get data start
        log_stream.info(' -----> Get source datasets ... ')

        # check destination file availability
        if not os.path.exists(file_path_anc_def):

            # iterate over time step(s)
            point_dynamic_collections = None
            for time_step in time_obj:

                # time start info
                log_stream.info(' ------> Time "' + time_step.strftime(time_format_algorithm) + '" ... ')

                # define source file name
                file_path_src_step = self.define_file_name(file_path_src_raw, time_step=time_step)

                # check source file availability
                if os.path.exists(file_path_src_step):

                    # check source file type
                    if self.type_src == 'csv':

                        # method to get data
                        point_dynamic_data = read_file_csv(file_path_src_step)
                        # method to join data
                        point_dynamic_collections = join_data_point(
                            point_dynamic_data, point_dynamic_collections,
                            point_time=None, name_index='time')

                    else:
                        # file type not expected (error message)
                        log_stream.error(' ===> File type "' + self.type_src + "' not expected")
                        raise NotImplementedError('File type not implemented yet')

                else:
                    # file not found (warning message)
                    log_stream.warning(' ===> File "' + file_path_src_step + "' does not exists")

                # time end info
                log_stream.info(' ------> Time "' + time_step.strftime(time_format_algorithm) + '" ... DONE')

            # get data end
            log_stream.info(' -----> Get datasets ... DONE')

            # dump data start
            log_stream.info(' -----> Save datasets ... ')
            # check data availability
            if point_dynamic_collections is not None:

                # method to dump data
                folder_name_anc, file_name_anc = os.path.split(file_path_anc_def)
                os.makedirs(folder_name_anc, exist_ok=True)
                write_obj(file_path_anc_def, point_dynamic_collections)

                # dump data end
                log_stream.info(' -----> Save datasets ... DONE')

            else:
                # dump data end
                point_dynamic_collections = None
                log_stream.info(' -----> Save datasets ... SKIPPED. Datasets is not available')

        else:
            # read ancillary file
            point_dynamic_collections = read_obj(file_path_anc_def)

            # get data end
            log_stream.info(' -----> Get datasets ... SKIPPED. Datasets were previously saved')

        # method to range data point
        time_frequency, time_start, time_end = range_data_point(
            point_dynamic_collections, time_run_reference=self.time_reference,
            time_start_reference=self.time_start, time_end_reference=self.time_end)

        # method to combine data point to the expected time range
        point_dynamic_collections = combine_data_point_by_time(
            point_dynamic_collections, point_static_collections,
            time_start_expected=time_start, time_end_expected=time_end,
            time_frequency_expected=time_frequency, time_reverse=True)

        # method end info
        log_stream.info(' ----> Organize dynamic object(s) ... DONE')

        return point_dynamic_collections

    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
