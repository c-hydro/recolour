"""
Class Features

Name:          driver_data_dynamic_dst
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from copy import deepcopy

import pandas as pd

from lib_utils_obj import create_dict_from_list, map_vars_dframe
from lib_data_io_csv import write_file_csv
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
    def __init__(self, time_obj, static_obj, ancillary_dict, destination_dict,
                 flags_dict=None, template_dict=None, params_dict=None, tmp_dict=None):

        self.time_obj = time_obj

        self.time_start = pd.DatetimeIndex([time_obj[0], time_obj[-1]]).min()
        self.time_end = pd.DatetimeIndex([time_obj[0], time_obj[-1]]).max()

        self.time_end, self.time_start = time_obj[0], time_obj[-1]
        self.static_obj = static_obj

        self.ancillary_dict = ancillary_dict
        self.destination_dict = destination_dict

        self.flags_dict = flags_dict
        self.template_dict = template_dict
        self.params_dict = params_dict
        self.tmp_dict = tmp_dict

        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.fields_tag, self.no_data_tag = 'fields', 'no_data'

        self.format_dst = self.params_dict['format_destination']

        self.reset_dst = flags_dict['reset_dynamic_destination']

        self.grid_obj = self.static_obj['grid_obj']
        self.point_obj = self.static_obj['point_obj']

        # ancillary object(s)
        folder_name_anc = ancillary_dict[self.folder_name_tag]
        file_name_anc = ancillary_dict[self.file_name_tag]
        self.file_path_anc = os.path.join(folder_name_anc, file_name_anc)
        self.fields_anc = self.ancillary_dict[self.fields_tag]

        # destination object(s)
        folder_name_dst = self.destination_dict[self.folder_name_tag]
        file_name_dst = self.destination_dict[self.file_name_tag]
        self.file_path_dst = os.path.join(folder_name_dst, file_name_dst)
        self.fields_dst = self.destination_dict[self.fields_tag]
        self.no_data_dst = self.destination_dict[self.no_data_tag]

        # tmp object(s)
        self.folder_name_tmp = tmp_dict[self.folder_name_tag]
        self.file_name_tmp = tmp_dict[self.file_name_tag]

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to define file name
    def define_file_name(self, file_path_raw, time_step=None, time_start=None, time_end=None):

        if time_step is None:
            time_step = deepcopy(time_end)

        template_time_generic = self.template_dict['time']

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

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump datasets object
    def dump_obj_datasets(self, file_path, file_dframe_raw, file_fields=None):

        # info start method
        log_stream.info(' -----> Dump file datasets "' + file_path + '" ... ')

        if file_path.endswith('csv'):

            # map datasets fields
            file_dframe_map = map_vars_dframe(file_dframe_raw, file_fields)
            # fill nan(s) with no_data
            file_dframe_map = file_dframe_map.fillna(self.no_data_dst)

            # create datasets folder
            folder_name, file_name = os.path.split(file_path)
            make_folder(folder_name)

            # dump datasets dframe
            write_file_csv(
                file_path, file_dframe_map,
                dframe_sep=',', dframe_decimal='.', dframe_float_format='%.3f',
                dframe_index=True, dframe_header=True, dframe_index_label='time')

        else:
            log_stream.info(' -----> Dump file datasets "' + file_path + '" ... FAILED')
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Dump file datasets "' + file_path + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self, dframe_datasets):

        # method start info
        log_stream.info(' ----> Organize destination object(s) ... ')

        # get path(s)
        file_path_dst_raw, file_path_anc_raw = self.file_path_dst, self.file_path_anc

        # get flag(s)
        reset_dst = self.reset_dst

        # define ancillary file
        file_path_anc_def = self.define_file_name(file_path_anc_raw, time_start=self.time_start, time_end=self.time_end)
        # define destination file
        file_path_dst_def = self.define_file_name(file_path_dst_raw, time_start=self.time_start, time_end=self.time_end)

        # reset ancillary file
        if reset_dst:
            if os.path.exists(file_path_dst_def):
                os.remove(file_path_dst_def)

        # check destination file availability
        if not os.path.exists(file_path_dst_def):

            # check destination format
            if self.format_dst == 'csv_dr':

                # check availability of destination file(s)
                if dframe_datasets is not None:

                    # dump datasets dframe
                    self.dump_obj_datasets(file_path_dst_def, dframe_datasets, self.fields_dst)

                    # method end info
                    log_stream.info(' ----> Organize destination object(s) ... DONE')

                else:
                    # method end info
                    log_stream.info(' ----> Organize destination object(s) ... SKIPPED. '
                                    'Destination file were previously dumped.')
        else:
            # method end info
            log_stream.info(' ----> Organize destination object(s) ... FAILED.')
            log_stream.error(' ===> Destination format "' + self.format_dst + '" is not supported')
            raise NotImplemented('Case not implemented yet')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
