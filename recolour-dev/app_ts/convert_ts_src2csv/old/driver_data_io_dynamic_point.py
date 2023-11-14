"""
Class Features

Name:          driver_data_io_dynamic_point
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '202200502'
Version:       '1.5.0'
"""

######################################################################################
# Library
import logging
import os
import re
import glob
import pandas as pd

from copy import deepcopy

from lib_data_io_pickle import write_obj

from lib_utils_data_point import get_data_mat, concatenate_data_series
from lib_utils_system import unzip_filename
from lib_utils_time import define_time_range
from lib_utils_system import fill_tags2string, make_folder, change_extension

from lib_info_args import logger_name, zip_extension

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverData for point soil moisture datasets
class DriverData:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, time_step, time_reference, src_dict, anc_dict, alg_dict=None,
                 geo_dict=None, time_dict=None, tmp_dict=None,
                 template_tags_dict=None,
                 flag_data_src='point',
                 flag_data_updating=True):

        self.time_step = pd.Timestamp(time_step)
        self.time_reference = pd.Timestamp(time_reference)

        self.src_dict = src_dict
        self.anc_dict = anc_dict
        self.alg_dict = alg_dict
        self.geo_dict = geo_dict
        self.tmp_dict = tmp_dict

        self.flag_data_src = flag_data_src

        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'
        self.file_type_tag = 'obj_type'
        self.file_compression_tag = 'obj_compression'

        self.grid_terrain_tag = 'terrain'
        self.grid_cn_tag = 'cn'
        self.grid_cnet_tag = 'channels_network'
        self.points_registry_tag = 'stations_registry'

        self.geo_x_tag = 'Longitude'
        self.geo_y_tag = 'Latitude'

        self.template_tags_dict = template_tags_dict

        # algorithm object(s)
        self.alg_catchment_name = alg_dict['catchment_name']
        self.alg_point_geo_method_search = alg_dict['geo_method_search']
        self.alg_point_geo_radius_influence = alg_dict['geo_radius_influence']
        self.alg_point_geo_neighbours = alg_dict['geo_neighbours']

        # time object(s)
        self.time_dict = time_dict[self.flag_data_src]
        self.time_range, self.time_start, self.time_end = define_time_range(self.time_dict)

        # source object(s)
        self.file_type_src = src_dict[self.flag_data_src][self.file_type_tag]
        self.folder_name_src_raw = self.src_dict[self.flag_data_src][self.folder_name_tag]
        self.file_name_src_raw = self.src_dict[self.flag_data_src][self.file_name_tag]

        self.file_path_src = self.collect_file_list(
            self.folder_name_src_raw, self.file_name_src_raw, file_time_range=self.time_range)

        # ancillary object(s)
        self.file_name_anc_raw = anc_dict[self.flag_data_src][self.file_name_tag]
        self.folder_name_anc_raw = anc_dict[self.flag_data_src][self.folder_name_tag]

        self.file_path_anc = self.collect_file_list(
            self.folder_name_anc_raw, self.file_name_anc_raw, file_time_range=pd.DatetimeIndex([self.time_reference]))

        # tmp object(s)
        self.folder_name_tmp_raw = tmp_dict[self.folder_name_tag]
        self.file_name_tmp_raw = tmp_dict[self.file_name_tag]

        self.file_extension_zip = zip_extension
        self.file_extension_unzip = 'bin'

        self.flag_data_updating = flag_data_updating

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to collect ancillary file
    def collect_file_list(self, folder_name_raw, file_name_raw, file_time_range=None):

        catchment_name = self.alg_catchment_name

        if '*' in folder_name_raw:
            log_stream.error(' ===> Special character "*" is not supported in "folder_name" definition')
            raise NotImplementedError('Case not implemented yet')

        if '*' in file_name_raw:

            template_values_step = {'catchment_name': catchment_name,
                                    'point_code': None, 'point_name': None,
                                    'time_start': None, 'time_end': None}

            folder_obj_def = fill_tags2string(
                folder_name_raw, self.template_tags_dict, template_values_step)
            folder_name_def = folder_obj_def[0]
            file_obj_def = fill_tags2string(
                file_name_raw, self.template_tags_dict, template_values_step)
            file_name_def = file_obj_def[0]

            file_name_obj = os.path.join(folder_name_def, file_name_def)

        else:

            file_name_list = []
            for time_step in file_time_range:

                template_values_step = {'source_sub_path_time_grid': time_step,
                                        'source_datetime_grid': time_step,
                                        'source_sub_path_time_point': time_step,
                                        'source_datetime_point': time_step,
                                        'ancillary_sub_path_time_grid': time_step,
                                        'ancillary_datetime_grid': time_step,
                                        'ancillary_sub_path_time_point': time_step,
                                        'ancillary_datetime_point': time_step,
                                        'destination_sub_path_time': time_step,
                                        'destination_datetime': time_step,
                                        'catchment_name': catchment_name,
                                        'time_start': None,
                                        'time_end': None,
                                        'point_code': None,
                                        'point_name': None}

                folder_tags_def = fill_tags2string(
                    folder_name_raw, self.template_tags_dict, template_values_step)
                folder_name_def = folder_tags_def[0]

                if file_name_raw is not None:

                    file_tags_def = fill_tags2string(
                        file_name_raw, self.template_tags_dict, template_values_step)
                    file_name_def = file_tags_def[0]

                    file_path_def = os.path.join(folder_name_def, file_name_def)
                else:
                    file_path_def = folder_name_def

                file_name_list.append(file_path_def)

            if file_name_list.__len__() == 1:
                file_name_obj = file_name_list[0]
            elif file_name_list.__len__() > 1:
                file_name_obj = deepcopy(file_name_list)
            else:
                log_stream.error(' ===> File list is empty')
                raise RuntimeError('File list must be defined by one or more elements')

        return file_name_obj

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to search
    @staticmethod
    def search_file_list(file_path_undefined):

        if '*' in file_path_undefined:
            file_path_list = glob.glob(file_path_undefined)
            file_path_list = sorted(file_path_list)
        else:
            file_path_list = deepcopy(file_path_undefined)

        if isinstance(file_path_list, str):
            file_path_list = [file_path_list]

        return file_path_list

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to search file time period
    @staticmethod
    def search_file_time(file_path_list):

        file_time_start, file_time_end = [], []
        for file_id, file_path_step in enumerate(file_path_list):

            file_time_tmp = re.search(r'\d{12}_\d{12}', file_path_step)
            file_time_str_start, file_time_str_end = file_time_tmp.group().split('_')

            file_time_str_start = str(pd.Timestamp(file_time_str_start))
            file_time_str_end = str(pd.Timestamp(file_time_str_end))

            file_time_start.append(file_time_str_start)
            file_time_end.append(file_time_str_end)

        return file_time_start, file_time_end
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize data
    def organize_data(self, var_name='soil_moisture', var_min=0.0, var_max=100.0, var_scale_factor=1.0):

        log_stream.info(' ----> Organize soil moisture point ... ')

        time_range = self.time_range
        file_type_src = self.file_type_src

        file_path_src_raw = self.file_path_src
        file_path_anc_raw = self.file_path_anc

        geo_da_terrain = self.geo_dict[self.grid_terrain_tag]
        geo_da_cn = self.geo_dict[self.grid_cn_tag]
        geo_da_cnet = self.geo_dict[self.grid_cnet_tag]
        point_dframe_registry = self.geo_dict[self.points_registry_tag]

        flag_data_updating = self.flag_data_updating

        for point_id, point_row in point_dframe_registry.iterrows():

            point_name, point_code = point_row['point_name'], point_row['point_code']

            template_values_dict = {'point_code': point_code, 'point_name': point_name}
            file_path_src_def = fill_tags2string(file_path_src_raw, self.template_tags_dict, template_values_dict)[0]
            file_path_anc_def = fill_tags2string(file_path_anc_raw, self.template_tags_dict, template_values_dict)[0]

            if flag_data_updating:
                if os.path.exists(file_path_anc_def):
                    os.remove(file_path_anc_def)

            log_stream.info(' -----> Point "' + point_name + '" ... ')

            if not os.path.exists(file_path_anc_def):

                file_path_list_src = self.search_file_list(file_path_src_def)
                file_time_start_src, file_time_end_src = self.search_file_time(file_path_list_src)

                log_stream.info(' ------> Get datasets ... ')

                sm_point_series = None
                for time_start_step, time_end_step, file_src_step in zip(
                        file_time_start_src, file_time_end_src, file_path_list_src):

                    log_stream.info(' --------> Time "' + str(time_start_step) + '" - "'
                                    + str(time_end_step) + '"... ')

                    if os.path.exists(file_src_step):

                        if file_src_step.endswith(self.file_extension_zip):
                            file_tmp_step = change_extension(file_src_step, self.file_extension_unzip)
                            unzip_filename(file_src_step, file_tmp_step)
                        else:
                            file_tmp_step = file_src_step

                        if file_type_src == 'point_mat':

                            sm_point_base = get_data_mat(
                                file_tmp_step, var_name=point_name,
                                var_value_min=var_min, var_value_max=var_max, var_scale_factor=var_scale_factor)

                        else:
                            log_stream.error(' ===> Source data type "' + file_type_src + '" is not supported.')
                            raise NotImplementedError('Only "point_mat" type is available.')

                        # Concatenate series of different period
                        sm_point_series = concatenate_data_series(sm_point_base, sm_point_series)

                        log_stream.info(' --------> Time "' + str(time_start_step) + '" - "'
                                        + str(time_end_step) + '"... DONE')

                    else:
                        log_stream.info(' --------> Time "' + str(time_start_step) + '" - "'
                                        + str(time_end_step) + '"... FAILED')
                        log_stream.warning(' ===> File: "' + file_src_step + '" does not exist')

                log_stream.info(' ------> Get datasets ... DONE')

                # Save point time series
                log_stream.info(' ------> Save datasets ... ')
                if sm_point_series is not None:
                    folder_name_anc, file_name_anc = os.path.split(file_path_anc_def)
                    make_folder(folder_name_anc)

                    write_obj(file_path_anc_def, sm_point_series)

                    log_stream.info(' ------> Save datasets ... DONE')

                else:
                    log_stream.info(' ------> Save datasets ... FAILED')
                    log_stream.error(' ===> Datasets are defined by NoneType')
                    raise IOError('Datasets must be defined to correctly run the algorithm')

                log_stream.info(' -----> Point "' + point_name + '" ... DONE')

            else:

                log_stream.info(' -----> Point "' + point_name + '" ... SKIPPED. Dataframe previously computed')

        log_stream.info(' ----> Organize soil moisture point ... DONE')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
