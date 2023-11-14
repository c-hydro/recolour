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

from lib_utils_fx import convert_obj_vars_to_point
from lib_utils_obj import create_dict_from_list, map_vars_dict
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
    def define_file_name(self, file_path_raw,
                         time_step=None, time_start=None, time_end=None,
                         point_name=None, point_code=None, point_tag=None, var_name=None):

        if time_step is None:
            time_step = deepcopy(time_end)

        template_time_generic = self.template_dict['time']
        template_datasets_generic = self.template_dict['datasets']

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

        generic_datasets_list = list(template_datasets_generic.keys())
        values_datasets_dict_generic = {}
        if point_name is not None:
            values_datasets_dict_generic['point_name'] = point_name
        if point_code is not None:
            values_datasets_dict_generic['point_code'] = str(point_code)
        if point_tag is not None:
            values_datasets_dict_generic['point_tag'] = point_tag
        if var_name is not None:
            values_datasets_dict_generic['var_name'] = var_name

        file_path_def = fill_tags2string(file_path_raw, template_time_generic, values_dict_generic)[0]
        file_path_def = fill_tags2string(file_path_def, template_datasets_generic, values_datasets_dict_generic)[0]

        return file_path_def

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump datasets object
    def dump_obj_datasets(self, file_path, file_obj_raw, file_fields=None):

        # info start method
        log_stream.info(' ------> Dump file datasets "' + file_path + '" ... ')

        if file_path.endswith('csv'):

            # map datasets fields
            file_obj_var = map_vars_dict(file_obj_raw, file_fields)
            # fill nan(s) with no_data
            file_obj_var = file_obj_var.fillna(self.no_data_dst)

            # create datasets folder
            folder_name, file_name = os.path.split(file_path)
            make_folder(folder_name)

            # dump datasets dframe
            write_file_csv(
                file_path, file_obj_var,
                dframe_sep=',', dframe_decimal='.', dframe_float_format='%.3f',
                dframe_index=True, dframe_header=True, dframe_index_label='time')

        else:
            log_stream.info(' ------> Dump file datasets "' + file_path + '" ... FAILED')
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Dump file datasets "' + file_path + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self, obj_datasets):

        # method start info
        log_stream.info(' ----> Organize destination object(s) ... ')

        # get path(s)
        file_path_dst_raw, file_path_anc_raw = self.file_path_dst, self.file_path_anc

        # get flag(s)
        reset_dst = self.reset_dst

        # define ancillary file
        file_path_anc_def = self.define_file_name(file_path_anc_raw, time_start=self.time_start, time_end=self.time_end)

        # check destination format
        if self.format_dst == 'csv_dr':

            # check availability of destination file(s)
            if obj_datasets is not None:

                # map datasets fields
                obj_datasets = map_vars_dict(obj_datasets, self.fields_dst)

                # iterate over variable(s)
                for var_name, var_dframe in obj_datasets.items():

                    # variable start info
                    log_stream.info(' -----> Variable "' + var_name + '" ... ')

                    # define destination file
                    file_path_dst_def = self.define_file_name(
                        file_path_dst_raw, time_start=self.time_start, time_end=self.time_end,
                        var_name=var_name)

                    # reset ancillary file
                    if reset_dst:
                        if os.path.exists(file_path_dst_def):
                            os.remove(file_path_dst_def)

                    # check destination file availability
                    if not os.path.exists(file_path_dst_def):

                        # dump datasets dframe
                        self.dump_obj_datasets(file_path_dst_def, var_dframe)

                        # variable end info
                        log_stream.info(' -----> Variable "' + var_name + '" ... DONE')

                    else:
                        # method end info
                        log_stream.info(' -----> Variable "' + var_name + '"  ... SKIPPED. '
                                        'Destination file were previously dumped.')

                # method end info
                log_stream.info(' ----> Organize destination object(s) ... DONE')

            else:
                # method end info
                log_stream.info(' ----> Organize destination object(s) ... SKIPPED. Datasets is defined by NoneType.')

        else:
            # method end info
            log_stream.info(' ----> Organize destination object(s) ... FAILED.')
            log_stream.error(' ===> Destination format "' + self.format_dst + '" is not supported')
            raise NotImplemented('Case not implemented yet')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
