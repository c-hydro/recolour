"""
Class Features

Name:          driver_data_dynamic_src
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd
from copy import deepcopy

from lib_data_io_nc import read_file_nc
from lib_data_io_tiff import read_file_tiff
from lib_data_io_pickle import read_obj, write_obj
from lib_data_io_generic import extract_data_grid2point, join_data_point, combine_data_point_by_time

from lib_utils_fx import convert_dataset_to_rzsm
from lib_utils_obj import create_dict_from_list, create_dataset, filter_dataset, convert_dataset_to_data_array
from lib_utils_system import fill_tags2string, make_folder
from lib_utils_time import define_time_frequency

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
    def __init__(self, time_obj, static_obj, source_dict, ancillary_dict,
                 flags_dict=None, template_dict=None, params_dict=None, tmp_dict=None):

        self.time_obj = time_obj

        self.time_start = pd.DatetimeIndex([time_obj[0], time_obj[-1]]).min()
        self.time_end = pd.DatetimeIndex([time_obj[0], time_obj[-1]]).max()
        self.time_frequency = define_time_frequency(time_obj)

        self.time_end, self.time_start = time_obj[0], time_obj[-1]
        self.static_obj = static_obj

        self.source_dict = source_dict
        self.ancillary_dict = ancillary_dict

        self.flags_dict = flags_dict
        self.template_dict = template_dict
        self.params_dict = params_dict
        self.tmp_dict = tmp_dict

        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.fields_tag = 'fields'

        self.reset_src = flags_dict['reset_dynamic_source']

        self.format_src = self.params_dict['format_source']
        self.geo_method_search = self.params_dict['geo_method_search']
        self.geo_radius_influence = self.params_dict['geo_radius_influence']
        self.geo_neighbours = self.params_dict['geo_neighbours']
        self.geo_spatial_operation = self.params_dict['geo_spatial_operation']
        self.geo_spatial_window = self.params_dict['geo_spatial_window']
        self.geo_spatial_mask = self.params_dict['geo_spatial_mask']

        self.grid_obj = self.static_obj['grid_obj']
        self.point_obj = self.static_obj['point_obj']

        # source object(s)
        folder_name_src = self.source_dict[self.folder_name_tag]
        file_name_src = self.source_dict[self.file_name_tag]
        self.file_path_src = os.path.join(folder_name_src, file_name_src)
        self.fields_src = self.source_dict[self.fields_tag]

        # ancillary object(s)
        folder_name_anc = ancillary_dict[self.folder_name_tag]
        file_name_anc = ancillary_dict[self.file_name_tag]
        self.file_path_anc = os.path.join(folder_name_anc, file_name_anc)
        self.fields_anc = self.ancillary_dict[self.fields_tag]

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
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize source object(s) ... ')

        # get time(s)
        time_obj = self.time_obj
        # get static
        grid_obj_geo, point_obj_geo = self.grid_obj, self.point_obj

        # get path(s)
        file_path_src_raw, file_path_anc_raw = self.file_path_src, self.file_path_anc

        # get flag(s)
        reset_src = self.reset_src

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
            point_collections = None
            for time_step in time_obj:

                # time start info
                log_stream.info(' ------> Time "' + time_step.strftime(time_format_algorithm) + '" ... ')

                # define source file name
                file_path_src_step = self.define_file_name(file_path_src_raw, time_step=time_step)

                # check source file availability
                if os.path.exists(file_path_src_step):

                    # check source file format
                    if self.format_src == 'ecmwf_grid_nc':

                        # method to get data
                        grid_data, grid_attrs, grid_geo_x, grid_geo_y, geo_attrs = read_file_nc(
                            file_path_src_step, file_variables_selected=list(self.fields_src.values()))

                        # check source data availability
                        if grid_data is not None:

                            # method to create dataset
                            grid_obj_dset_raw = create_dataset(
                                grid_data, grid_geo_x, grid_geo_y,
                                data_time=time_step, data_attrs=grid_attrs, common_attrs=geo_attrs)
                            # method to select dataset
                            grid_obj_dset_def = filter_dataset(grid_obj_dset_raw, dset_vars_filter=self.fields_src)
                            # method to convert dset to darray collection
                            grid_obj_da_vars = convert_dataset_to_data_array(grid_obj_dset_def)

                            # method to compute destination datasets
                            grid_obj_da_rzsm = convert_dataset_to_rzsm(grid_obj_da_vars)
                        else:
                            # destination datasets not available
                            grid_obj_da_rzsm = None

                    else:
                        log_stream.error(' ===> Source data type "' + self.format_src + '" is not supported.')
                        raise NotImplemented('Case not implemented yet')

                    # check destination data availability
                    if grid_obj_da_rzsm is not None:
                        # method to get data point(s)
                        point_data = extract_data_grid2point(
                            grid_obj_da_rzsm, grid_obj_geo, point_obj_geo,
                            method_spatial_operation=self.geo_spatial_operation,
                            method_spatial_mask=self.geo_spatial_mask)
                    else:
                        # data point(s) not available
                        point_data = None

                    # method to join data point(s)
                    point_collections = join_data_point(time_step, point_data, point_collections)

                else:
                    # file not found (warning message)
                    log_stream.warning(' ===> File "' + file_path_src_step + "' does not exists")

                # time end info
                log_stream.info(' ------> Time "' + time_step.strftime(time_format_algorithm) + '" ... DONE')

            # get data end
            log_stream.info(' -----> Get datasets ... DONE')

            # dump data start
            log_stream.info(' -----> Dump datasets ... ')
            # check data availability
            if point_collections is not None:

                # method to dump data
                folder_name_anc, file_name_anc = os.path.split(file_path_anc_def)
                make_folder(folder_name_anc)
                write_obj(file_path_anc_def, point_collections)

                # dump data end
                log_stream.info(' -----> Dump datasets ... DONE')

            else:
                # dump data end
                point_collections = None
                log_stream.info(' -----> Dump datasets ... SKIPPED. Datasets is not available')

        else:
            # read ancillary file
            point_collections = read_obj(file_path_anc_def)

            # get data end
            log_stream.info(' -----> Get datasets ... SKIPPED. Datasets were previously saved')

        # method to combine data point to the expected time range
        point_collections = combine_data_point_by_time(
            point_collections, point_obj_geo,
            time_start_expected=self.time_start, time_end_expected=self.time_end,
            time_frequency_expected=self.time_frequency, time_reverse=True)

        # method end info
        log_stream.info(' ----> Organize source object(s) ... DONE')

        return point_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
