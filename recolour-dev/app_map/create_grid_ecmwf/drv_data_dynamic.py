"""
Class Features

Name:          drv_data_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from copy import deepcopy

from lib_info_args import time_format_algorithm
from lib_info_args import (geo_dim_name_x, geo_dim_name_y, geo_coord_name_x, geo_coord_name_y,
                           geo_var_name_x, geo_var_name_y)

from lib_utils_time import set_time_file
from lib_utils_generic import make_folder, reset_folder
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_nc import read_file_nc, organize_file_nc, write_file_nc
from lib_data_io_tiff import organize_file_tiff, write_file_tiff

from lib_utils_io import fill_string_with_time

from lib_fx_utils import organize_data, convert_data, check_data

import lib_fx_methods_datasets as fx_methods_datasets
import lib_fx_methods_common as fx_methods_common

from drv_data_fx import DrvFx, get_fx_method, get_fx_settings
# debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_time_reference, alg_time_datasets,
                 alg_static, alg_settings,
                 tag_section_flags='flags', tag_section_template='template',
                 tag_section_methods_datasets='methods_datasets', tag_section_methods_common='methods_common',
                 tag_section_datasets='datasets',
                 tag_section_log='log', tag_section_tmp='tmp'):

        self.alg_time_reference = alg_time_reference
        self.alg_time_datasets = alg_time_datasets

        self.alg_static = alg_static

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_template = alg_settings[tag_section_template]
        self.alg_methods_dset = alg_settings[tag_section_methods_datasets]
        self.alg_methods_common = alg_settings[tag_section_methods_common]

        self.alg_datasets_src = alg_settings[tag_section_datasets]['dynamic']['source']
        self.alg_datasets_anc_raw = alg_settings[tag_section_datasets]['dynamic']['ancillary']['raw']
        self.alg_datasets_anc_def = alg_settings[tag_section_datasets]['dynamic']['ancillary']['def']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_settings[tag_section_log]
        self.alg_tmp = alg_settings[tag_section_tmp]

        self.tag_folder_name, self.tag_file_name = 'folder_name', 'file_name'
        self.tag_time_frequency, self.tag_time_rounding = 'time_frequency', 'time_rounding'
        self.tag_time_period = 'time_period'
        self.tag_format = 'format'

        self.reset_datasets_anc_raw = self.alg_flags['reset_datasets_ancillary_raw']
        self.reset_datasets_anc_def = self.alg_flags['reset_datasets_ancillary_def']
        self.reset_datasets_dst = self.alg_flags['reset_datasets_destination']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_src = self.alg_datasets_src[self.tag_folder_name]
        self.file_name_src = self.alg_datasets_src[self.tag_file_name]
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)
        self.time_frequency_src = self.alg_datasets_src[self.tag_time_frequency]
        self.time_rounding_src = self.alg_datasets_src[self.tag_time_rounding]
        self.time_period_src = self.alg_datasets_src[self.tag_time_period]

        self.folder_name_anc_raw = self.alg_datasets_anc_raw[self.tag_folder_name]
        self.file_name_anc_raw = self.alg_datasets_anc_raw[self.tag_file_name]
        self.file_path_anc_raw = os.path.join(self.folder_name_anc_raw, self.file_name_anc_raw)

        self.folder_name_anc_def = self.alg_datasets_anc_def[self.tag_folder_name]
        self.file_name_anc_def = self.alg_datasets_anc_def[self.tag_file_name]
        self.file_path_anc_def = os.path.join(self.folder_name_anc_def, self.file_name_anc_def)

        self.folder_name_dst = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_datasets_dst[self.tag_file_name]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)
        self.time_frequency_dst = self.alg_datasets_dst[self.tag_time_frequency]
        self.time_rounding_dst = self.alg_datasets_dst[self.tag_time_rounding]
        self.time_period_dst = self.alg_datasets_dst[self.tag_time_period]
        self.format_dst = self.alg_datasets_dst[self.tag_format]

        self.grid_data, self.grid_attrs = self.alg_static['grid_data'], self.alg_static['grid_attrs']

        self.fx_settings_dset_compute = self.alg_methods_dset['fx_compute_data']
        self.fx_settings_common_crop = self.alg_methods_common['fx_crop_data']
        self.fx_settings_common_resample = self.alg_methods_common['fx_resample_data']
        self.fx_settings_common_mask = self.alg_methods_common['fx_mask_data']

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize dynamic datasets ... ')

        # get file path
        file_path_src_tmpl = self.file_path_src
        file_path_anc_raw_tmpl, file_path_anc_def_tmpl = self.file_path_anc_raw, self.file_path_anc_def
        file_path_dst_tmpl = self.file_path_dst

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            logging.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # set file time(s)
            alg_time_subset = set_time_file(
                time_start=alg_time_group, time_end=None,
                time_frequency=self.time_frequency_src, time_rounding=self.time_rounding_src,
                time_period=self.time_period_src, time_reverse=True)

            # iterate over time file(s)
            for alg_time_step in alg_time_subset:

                # info start time step
                logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) + '" ... ')

                # method to fill the filename(s)
                file_path_src_step = fill_string_with_time(file_path_src_tmpl, alg_time_step, self.alg_template)
                file_path_anc_raw_step = fill_string_with_time(file_path_anc_raw_tmpl, alg_time_step, self.alg_template)
                file_path_anc_def_step = fill_string_with_time(file_path_anc_def_tmpl, alg_time_step, self.alg_template)
                file_path_dst_step = fill_string_with_time(file_path_dst_tmpl, alg_time_step, self.alg_template)

                # clean ancillary datasets (if ancillary flag(s) is activated)
                if self.reset_datasets_anc_raw or self.reset_datasets_anc_def:
                    if os.path.exists(file_path_anc_raw_step):
                        os.remove(file_path_anc_raw_step)
                    if os.path.exists(file_path_anc_def_step):
                        os.remove(file_path_anc_def_step)
                    if os.path.exists(file_path_dst_step):
                        os.remove(file_path_dst_step)
                # clean destination datasets (if ancillary flag(s) is activated)
                if self.reset_datasets_dst:
                    if os.path.exists(file_path_dst_step):
                        os.remove(file_path_dst_step)
                # clean ancillary and destination datasets if are not available together
                if (not os.path.exists(file_path_anc_raw_step)) and (not os.path.exists(file_path_dst_step)):
                    if os.path.exists(file_path_anc_raw_step):
                        os.remove(file_path_anc_raw_step)
                    if os.path.exists(file_path_anc_def_step):
                        os.remove(file_path_anc_def_step)
                    if os.path.exists(file_path_dst_step):
                        os.remove(file_path_dst_step)

                # check file ancillary availability
                if not os.path.exists(file_path_anc_raw_step):
                    # check file source availability
                    if os.path.exists(file_path_src_step):

                        # get dataset source
                        dset_src_step = read_file_nc(file_path_src_step, file_variables_selected=None)

                        # organize variable(s) source
                        vars_src_step = organize_data(dset_src_step,
                                                      var_name_upd={'lon': geo_var_name_x, 'lat': geo_var_name_y})

                        # save variable(s) obj to ancillary file
                        folder_name_anc_raw_step, file_name_anc_raw_step = os.path.split(file_path_anc_raw_step)
                        make_folder(folder_name_anc_raw_step)
                        write_file_obj(file_path_anc_raw_step, vars_src_step)

                        # info end time step
                        logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                     '" ... DONE')
                    else:
                        # info end time step
                        logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                     '" ... FAILED. File "' + file_path_src_step + '" not available.')
                else:
                    # info end time step
                    logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                 '" ... SKIPPED. Datasets previously saved.')

            # info end time
            logging.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... DONE')

        # info end method
        logging.info(' ---> Organize dynamic datasets ... DONE')

    # method to analyze data
    def analyze_data(self):

        # info start method
        logging.info(' ---> Analyze dynamic data ... ')

        # get file path
        file_path_anc_raw_tmpl, file_path_anc_def_tmpl = self.file_path_anc_raw, self.file_path_anc_def
        # get grid info
        grid_data = self.grid_data

        # fx settings
        fx_settings_dset_cmp = self.fx_settings_dset_compute
        fx_settings_common_crop = self.fx_settings_common_crop
        fx_settings_common_res = self.fx_settings_common_resample
        fx_settings_common_mask = self.fx_settings_common_mask

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            logging.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # set file time(s)
            alg_time_subset = set_time_file(
                time_start=alg_time_group, time_end=None,
                time_frequency=self.time_frequency_src, time_rounding=self.time_rounding_src,
                time_period=self.time_period_src, time_reverse=True)

            # iterate over time file(s)
            for alg_time_step in alg_time_subset:

                # info start time step
                logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) + '" ... ')

                # method to fill the filename(s)
                file_path_anc_raw_step = fill_string_with_time(file_path_anc_raw_tmpl, alg_time_step, self.alg_template)
                file_path_anc_def_step = fill_string_with_time(file_path_anc_def_tmpl, alg_time_step, self.alg_template)

                # check file ancillary availability
                if not os.path.exists(file_path_anc_def_step):
                    # check file source availability
                    if os.path.exists(file_path_anc_raw_step):

                        # get variable(s) obj
                        fx_obj_data_src = read_file_obj(file_path_anc_raw_step)

                        # info start compute datasets
                        logging.info(' ------> (1) Compute datasets ... ')

                        # get fx dset settings
                        (fx_active_dset_cmp, fx_name_dset_cmp,
                            fx_vars_dset_cmp, fx_params_dset_cmp) = get_fx_settings(fx_settings_dset_cmp)
                        # get fx dset handle
                        fx_handle_dset_cmp = get_fx_method(fx_name_dset_cmp, fx_methods_datasets)

                        # flag to activate compute part
                        if fx_active_dset_cmp:
                            # initialize fx driver
                            driver_fx_dset = DrvFx(fx_obj_data_src,
                                                   fx_name=fx_name_dset_cmp, fx_handle=fx_handle_dset_cmp,
                                                   fx_params=fx_params_dset_cmp,
                                                   fx_vars_in=fx_vars_dset_cmp['in'], fx_vars_out=fx_vars_dset_cmp['out'])

                            # organize fx data source
                            fx_kwargs_data, fx_kwargs_geo = driver_fx_dset.organize_fx_src()
                            # execute fx method
                            fx_result_data = driver_fx_dset.execute_fx(fx_kwargs_data)
                            # organize fx data destination
                            fx_destination_data = driver_fx_dset.organize_fx_dst(fx_result_data, fx_kwargs_geo)

                            # method to convert data
                            fx_obj_data_cmp = convert_data(fx_destination_data)

                            # info end compute datasets
                            logging.info(' ------> (1) Compute datasets ... DONE')

                        else:
                            logging.info(' ------> (1) Compute datasets ... FAILED')
                            logging.error(' ===> Compute datasets part is needed by the algorithm to correctly run')
                            raise RuntimeError('Active this part in the configuration file')

                        # info start crop datasets
                        logging.info(' ------> (2) Crop datasets ... ')
                        # get fx crop settings
                        (fx_active_common_crop, fx_name_common_crop,
                            fx_vars_common_crop, fx_params_common_crop) = get_fx_settings(fx_settings_common_crop)
                        # get fx crop handle
                        fx_handle_common_crop = get_fx_method(fx_name_common_crop, fx_methods_common)

                        # flag to activate crop part
                        if fx_active_common_crop:
                            # call crop method
                            fx_obj_data_crop = fx_handle_common_crop(
                                obj_in=fx_obj_data_cmp,
                                geo_y_lower=grid_data.bb_bottom, geo_y_upper=grid_data.bb_top,
                                geo_x_left=grid_data.bb_left, geo_x_right=grid_data.bb_right)
                            # info end crop datasets
                            logging.info(' ------> (2) Crop datasets ... DONE')
                        else:
                            # info end crop datasets
                            fx_obj_data_crop = deepcopy(fx_obj_data_cmp)
                            logging.info(' ------> (2) Crop datasets ... SKIPPED. Flag is not activated')

                        # info start resample datasets
                        logging.info(' ------> (3) Resample datasets ... ')
                        # get fx resample settings
                        (fx_active_common_res, fx_name_common_res,
                            fx_vars_common_res, fx_params_common_res) = get_fx_settings(fx_settings_common_res)
                        # get fx resample handle
                        fx_handle_common_res = get_fx_method(fx_name_common_res, fx_methods_common)

                        # flag to activate resample part
                        if fx_active_common_res:
                            # call resample method
                            fx_obj_data_res = fx_handle_common_res(
                                obj_in=fx_obj_data_crop,
                                geo_x_out=grid_data[geo_var_name_x].values,
                                geo_y_out=grid_data[geo_var_name_y].values,
                                var_name_geo_x=geo_var_name_x, var_name_geo_y=geo_var_name_y,
                                coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                                dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y,
                                resampling_method=fx_params_common_res['method'],
                                resampling_max_distance=fx_params_common_res['max_distance'],
                                resampling_min_neighbours=fx_params_common_res['min_neighbours'],
                                resampling_neighbours=fx_params_common_res['neighbours'],
                                resampling_extend_data=fx_params_common_res['extend_data']
                            )

                            # info end resample datasets
                            logging.info(' ------> (3) Resample datasets ... DONE')
                        else:
                            # info end resample datasets
                            fx_obj_data_res = deepcopy(fx_obj_data_crop)
                            logging.info(' ------> (3) Resample datasets ... SKIPPED. Flag is not activated')

                        # info start mask datasets
                        logging.info(' ------> (4) Mask datasets ... ')
                        # get fx mask settings
                        (fx_active_common_mask, fx_name_common_mask,
                            fx_vars_common_mask, fx_params_common_mask) = get_fx_settings(fx_settings_common_mask)
                        # get fx mask handle
                        fx_handle_common_mask = get_fx_method(fx_name_common_mask, fx_methods_common)

                        # flag to activate resample part
                        if fx_active_common_res and fx_active_common_mask:
                            # call mask method
                            fx_obj_data_mask = fx_handle_common_mask(
                                fx_obj_data_res, grid_data,
                                mask_value_min=fx_params_common_mask['min_value'],
                                mask_value_max=fx_params_common_mask['max_value'])
                            # info end mask datasets
                            logging.info(' ------> (4) Mask datasets ... DONE')
                        else:
                            # info end mask datasets
                            fx_obj_data_mask = deepcopy(fx_obj_data_res)
                            logging.info(' ------> (4) Mask datasets ... SKIPPED. '
                                         'Flag resample or/and flag mask is/are not activated')

                        # check data
                        if check_data(fx_obj_data_mask, thr_data='all'):

                            # save weights collection in pickle format
                            folder_name_anc_def_step, file_name_anc_def_step = os.path.split(file_path_anc_def_step)
                            make_folder(folder_name_anc_def_step)
                            write_file_obj(file_path_anc_def_step, fx_obj_data_mask)

                            # info end time group
                            logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                         '" ... DONE')
                        else:
                            logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                         '" ... FAILED. Some/All datasets are defined by NoneType')
                    else:
                        # info end time group
                        logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                     '" ... FAILED. Source datasets not available')

                else:
                    # info end time group
                    logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                 '" ... SKIPPED. Destination datasets previously computed')

            # info end time
            logging.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... DONE')

        # info end method
        logging.info(' ---> Analyze dynamic data ... DONE')

    # method to save data
    def dump_data(self):

        # info start method
        logging.info(' ---> Dump dynamic datasets ... ')

        # get file path
        file_path_anc_def_tmpl, file_path_dst_tmpl = self.file_path_anc_def, self.file_path_dst
        # get grid info
        grid_data, grid_attrs = self.grid_data, self.grid_attrs

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            logging.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # set file time(s)
            alg_time_subset = set_time_file(
                time_start=alg_time_group, time_end=None,
                time_frequency=self.time_frequency_src, time_rounding=self.time_rounding_src,
                time_period=self.time_period_src, time_reverse=True)

            # iterate over time file(s)
            for alg_time_step in alg_time_subset:

                # info start time step
                logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) + '" ... ')

                # method to fill the filename(s)
                file_path_anc_def_step = fill_string_with_time(file_path_anc_def_tmpl, alg_time_step, self.alg_template)
                file_path_dst_step = fill_string_with_time(file_path_dst_tmpl, alg_time_step, self.alg_template)

                # check file destination availability
                if not os.path.exists(file_path_dst_step):
                    # check file source availability
                    if os.path.exists(file_path_anc_def_step):

                        # method to get data obj
                        variable_collection = read_file_obj(file_path_anc_def_step)

                        # check destination format
                        if self.format_dst == 'netcdf':
                            # method to organize netcdf dataset
                            variable_dset = organize_file_nc(obj_variable=variable_collection, obj_time=alg_time_step)
                            # method to write netcdf dataset
                            folder_name_dst_step, file_name_dst_step = os.path.split(file_path_dst_step)
                            make_folder(folder_name_dst_step)
                            write_file_nc(file_path_dst_step, variable_dset)

                        elif (self.format_dst == 'tiff') or (self.format_dst == 'tif'):
                            # method to organize tiff dataset
                            variable_data, variable_attrs = organize_file_tiff(
                                variable_collection, alg_time_step,
                                obj_transform=grid_attrs['transform'], obj_proj=grid_attrs['crs'])
                            # method to write tiff dataset
                            folder_name_dst_step, file_name_dst_step = os.path.split(file_path_dst_step)
                            make_folder(folder_name_dst_step)
                            write_file_tiff(
                                file_name=file_path_dst_step, file_data=variable_data,
                                **variable_attrs)

                        else:
                            logging.error(' ===> Destination format "' + self.format_dst + '" is not supported')
                            raise NotImplemented('Case not implemented yet')

                        # info end time group
                        logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                     '" ... DONE')
                    else:
                        # info end time group
                        logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                     '" ... FAILED. Source datasets not available')

                else:
                    # info end time group
                    logging.info(' -----> Time Step "' + alg_time_step.strftime(time_format_algorithm) +
                                 '" ... SKIPPED. Destination datasets previously computed')

            # info end time
            logging.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... DONE')

        # info end method
        logging.info(' ---> Dump dynamic datasets ... DONE')

# -------------------------------------------------------------------------------------
