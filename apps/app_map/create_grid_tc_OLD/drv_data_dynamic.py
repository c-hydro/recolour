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

import pandas as pd

from copy import deepcopy

from lib_info_args import time_format_algorithm
from lib_info_args import (geo_dim_name_x, geo_dim_name_y, geo_coord_name_x, geo_coord_name_y,
                           geo_var_name_x, geo_var_name_y)

from lib_utils_generic import make_folder, reset_folder
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_nc import read_file_nc, organize_file_nc, write_file_nc
from lib_data_io_tiff import organize_file_tiff, write_file_tiff

from lib_info_args import logger_name
from lib_utils_io import fill_path_with_tags

from lib_fx_utils import get_fx_method, get_fx_settings, remove_nans, check_data

import lib_fx_methods_datasets as fx_methods_datasets
import lib_fx_methods_common as fx_methods_common

from cpl_data_dynamic import select_datasets_fields, join_dataset_obj
from cpl_data_dynamic import CouplerDataset, CouplerAncillary

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
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
        self.alg_template_dset = alg_settings[tag_section_template]['datasets']
        self.alg_template_time = alg_settings[tag_section_template]['time']
        self.alg_methods_dset = alg_settings[tag_section_methods_datasets]
        self.alg_methods_common = alg_settings[tag_section_methods_common]

        self.alg_datasets_src = alg_settings[tag_section_datasets]['dynamic']['source']
        self.alg_datasets_anc_pnt_raw = alg_settings[tag_section_datasets]['dynamic']['ancillary']['points']['raw']
        self.alg_datasets_anc_pnt_def = alg_settings[tag_section_datasets]['dynamic']['ancillary']['points']['def']
        self.alg_datasets_anc_map = alg_settings[tag_section_datasets]['dynamic']['ancillary']['maps']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_settings[tag_section_log]
        self.alg_tmp = alg_settings[tag_section_tmp]

        self.tag_dset_data_name, self.tag_dset_metrics_name = 'soil_moisture', 'metrics'
        self.tag_dset_ref, self.tag_dset_k1, self.tag_dset_k2 = 'ref', 'k1', 'k2'

        self.tag_folder_name, self.tag_file_name = 'folder_name', 'file_name'
        self.tag_variable = 'variable'
        self.tag_time_frequency, self.tag_time_rounding = 'time_frequency', 'time_rounding'
        self.tag_time_period = 'time_period'
        self.tag_format = 'format'

        self.reset_datasets_anc_pnt_raw = self.alg_flags['reset_datasets_ancillary_points_raw']
        self.reset_datasets_anc_pnt_def = self.alg_flags['reset_datasets_ancillary_points_def']
        self.reset_datasets_anc_map = self.alg_flags['reset_datasets_ancillary_maps']
        self.reset_datasets_dst = self.alg_flags['reset_datasets_destination']
        self.reset_logs = self.alg_flags['reset_logs']

        self.alg_datasets_src_data = self.alg_datasets_src[self.tag_dset_data_name]
        self.alg_datasets_src_metrics = self.alg_datasets_src[self.tag_dset_metrics_name]

        self.alg_datasets_src_ref = self.alg_datasets_src_data[self.tag_dset_ref]
        self.folder_name_src_ref_tmpl = self.alg_datasets_src_ref[self.tag_folder_name]
        self.file_name_src_ref_tmpl = self.alg_datasets_src_ref[self.tag_file_name]
        self.file_path_src_ref_tmpl = os.path.join(self.folder_name_src_ref_tmpl, self.file_name_src_ref_tmpl)
        self.variable_src_ref = self.alg_datasets_src_ref[self.tag_variable]
        self.time_frequency_src_ref = self.alg_datasets_src_ref[self.tag_time_frequency]
        self.time_rounding_src_ref = self.alg_datasets_src_ref[self.tag_time_rounding]
        self.time_period_src_ref = self.alg_datasets_src_ref[self.tag_time_period]

        self.alg_datasets_src_k1 = self.alg_datasets_src_data[self.tag_dset_k1]
        self.folder_name_src_k1_tmpl = self.alg_datasets_src_k1[self.tag_folder_name]
        self.file_name_src_k1_tmpl = self.alg_datasets_src_k1[self.tag_file_name]
        self.file_path_src_k1_tmpl = os.path.join(self.folder_name_src_k1_tmpl, self.file_name_src_k1_tmpl)
        self.variable_src_k1 = self.alg_datasets_src_k1[self.tag_variable]
        self.time_frequency_src_k1 = self.alg_datasets_src_k1[self.tag_time_frequency]
        self.time_rounding_src_k1 = self.alg_datasets_src_k1[self.tag_time_rounding]
        self.time_period_src_k1 = self.alg_datasets_src_k1[self.tag_time_period]

        self.alg_datasets_src_k2 = self.alg_datasets_src_data[self.tag_dset_k2]
        self.folder_name_src_k2_tmpl = self.alg_datasets_src_k2[self.tag_folder_name]
        self.file_name_src_k2_tmpl = self.alg_datasets_src_k2[self.tag_file_name]
        self.file_path_src_k2_tmpl = os.path.join(self.folder_name_src_k2_tmpl, self.file_name_src_k2_tmpl)
        self.variable_src_k2 = self.alg_datasets_src_k2[self.tag_variable]
        self.time_frequency_src_k2 = self.alg_datasets_src_k2[self.tag_time_frequency]
        self.time_rounding_src_k2 = self.alg_datasets_src_k2[self.tag_time_rounding]
        self.time_period_src_k2 = self.alg_datasets_src_k2[self.tag_time_period]

        self.folder_name_src_metrics_tmpl = self.alg_datasets_src_metrics[self.tag_folder_name]
        self.file_name_src_metrics_tmpl = self.alg_datasets_src_metrics[self.tag_file_name]
        self.file_path_src_metrics_tmpl = os.path.join(self.folder_name_src_metrics_tmpl, self.file_name_src_metrics_tmpl)
        self.variable_src_metrics = self.alg_datasets_src_metrics[self.tag_variable]

        self.folder_name_anc_pnt_raw_tmpl = self.alg_datasets_anc_pnt_raw[self.tag_folder_name]
        self.file_name_anc_pnt_raw_tmpl = self.alg_datasets_anc_pnt_raw[self.tag_file_name]
        self.file_path_anc_pnt_raw_tmpl = os.path.join(self.folder_name_anc_pnt_raw_tmpl, self.file_name_anc_pnt_raw_tmpl)

        self.folder_name_anc_pnt_def_tmpl = self.alg_datasets_anc_pnt_def[self.tag_folder_name]
        self.file_name_anc_pnt_def_tmpl = self.alg_datasets_anc_pnt_def[self.tag_file_name]
        self.file_path_anc_pnt_def_tmpl = os.path.join(self.folder_name_anc_pnt_def_tmpl, self.file_name_anc_pnt_def_tmpl)

        self.folder_name_anc_map_tmpl = self.alg_datasets_anc_map[self.tag_folder_name]
        self.file_name_anc_map_tmpl = self.alg_datasets_anc_map[self.tag_file_name]
        self.file_path_anc_map_tmpl = os.path.join(self.folder_name_anc_map_tmpl, self.file_name_anc_map_tmpl)

        self.folder_name_dst_tmpl = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst_tmpl = self.alg_datasets_dst[self.tag_file_name]
        self.file_path_dst_tmpl = os.path.join(self.folder_name_dst_tmpl, self.file_name_dst_tmpl)
        self.variable_dst = self.alg_datasets_dst[self.tag_variable]
        self.time_frequency_dst = self.alg_datasets_dst[self.tag_time_frequency]
        self.time_rounding_dst = self.alg_datasets_dst[self.tag_time_rounding]
        self.time_period_dst = self.alg_datasets_dst[self.tag_time_period]
        self.format_dst = self.alg_datasets_dst[self.tag_format]

        self.grid_geo_data, self.grid_geo_attrs = self.alg_static['geo']['data'], self.alg_static['geo']['attrs']
        self.grid_ref, self.grid_k1, self.grid_k2 = self.alg_static['ref'], self.alg_static['k1'], self.alg_static['k2']
        self.grid_cells = self.alg_static['cells']

        self.alg_dset_list = [self.tag_dset_ref, self.tag_dset_k1, self.tag_dset_k2]

        self.fx_settings_dset_scale = self.alg_methods_dset['fx_scale_data']
        self.fx_settings_dset_weigh = self.alg_methods_dset['fx_weigh_data']
        self.fx_settings_common_resample = self.alg_methods_common['fx_resample_data']
        self.fx_settings_common_filter = self.alg_methods_common['fx_filter_data']
        self.fx_settings_common_mask = self.alg_methods_common['fx_mask_data']

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ---> Organize dynamic datasets ... ')

        # get datasets src data and metrics
        alg_datasets_src_data = self.alg_datasets_src_data
        alg_datasets_src_metrics = self.alg_datasets_src_metrics
        # get cells
        alg_grid_cells = self.grid_cells
        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # get file path
        file_path_anc_pnt_raw_tmpl = self.file_path_anc_pnt_raw_tmpl
        file_path_anc_pnt_def_tmpl = self.file_path_anc_pnt_def_tmpl
        file_path_anc_map_tmpl = self.file_path_anc_map_tmpl
        file_path_dst_tmpl = self.file_path_dst_tmpl

        # set ancillary flag
        reset_datasets_anc_pnt_raw = self.reset_datasets_anc_pnt_raw
        reset_datasets_anc_pnt_def = self.reset_datasets_anc_pnt_def

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary map file name
            file_path_anc_map_step = fill_path_with_tags(
                file_path_anc_map_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)
            # define destination file name
            file_path_dst_step = fill_path_with_tags(
                file_path_dst_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)

            # iterate over cell(s)
            for alg_cell_group in alg_grid_cells:

                # info start cell group
                alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) + '" ... ')

                # define ancillary file name
                file_path_anc_pnt_raw_step = fill_path_with_tags(
                    file_path_anc_pnt_raw_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)
                file_path_anc_pnt_def_step = fill_path_with_tags(
                    file_path_anc_pnt_def_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                # reset ancillary and destination file(s) (if needed)
                if reset_datasets_anc_pnt_raw or reset_datasets_anc_pnt_def:
                    if os.path.exists(file_path_anc_pnt_def_step):
                        os.remove(file_path_anc_pnt_def_step)
                    if os.path.exists(file_path_anc_pnt_raw_step):
                        os.remove(file_path_anc_pnt_raw_step)
                    if os.path.exists(file_path_anc_map_step):
                        os.remove(file_path_anc_map_step)
                    if os.path.exists(file_path_dst_step):
                        os.remove(file_path_dst_step)

                # check availability of ancillary file
                if not os.path.exists(file_path_anc_pnt_raw_step):

                    # info start metrics group
                    alg_logger.info(' ------> (1) Metrics Group ... ')

                    # method to include metrics
                    alg_metrics_args = select_datasets_fields(alg_datasets_src_metrics)
                    # initialize ancillary class
                    coupler_ancillary = CouplerAncillary(
                        cell_data=alg_cell_group, dataset_template=alg_template_dset,
                        **alg_metrics_args)
                    # method to get metrics
                    alg_metrics_obj = coupler_ancillary.get_ancillary()
                    # method to organize ancillary
                    alg_metrics_collections = coupler_ancillary.organize_ancillary(alg_metrics_obj)

                    # info end metrics group
                    alg_logger.info(' ------> (1) Metrics Group ... DONE')

                    # info start dataset group
                    alg_logger.info(' ------> (2) Dataset Group ... ')

                    # iterate over data name(s)
                    alg_dset_collections, alg_time_collections = {}, {}
                    for alg_dset_name in self.alg_dset_list:

                        # info start dataset group
                        alg_logger.info(' -------> Data "' + alg_dset_name + '" ... ')

                        # get datasets field(s)
                        alg_dset_fields = alg_datasets_src_data[alg_dset_name]

                        # prepare datasets args
                        alg_dset_args = select_datasets_fields(alg_dset_fields)

                        # initialize dataset class
                        coupler_dataset = CouplerDataset(
                            time_data=alg_time_group, cell_data=alg_cell_group,
                            dataset_template=alg_template_dset, time_template=alg_template_time,
                            **alg_dset_args)

                        # method to get data
                        alg_dset_obj, alg_dset_time = coupler_dataset.get_data()
                        # method to organize data
                        alg_dset_obj = coupler_dataset.organize_data(alg_dset_obj, alg_dset_time)

                        # store data in the collections object
                        alg_dset_collections[alg_dset_name] = alg_dset_obj
                        alg_time_collections[alg_dset_name] = alg_dset_time

                        # info start dataset group
                        alg_logger.info(' -------> Data "' + alg_dset_name + '" ... DONW')

                    # info start dataset group
                    alg_logger.info(' ------> (2) Dataset Group ... DONE')

                    # info start collection group
                    alg_logger.info(' ------> (3) Collections group ... ')

                    # method to join dataset(s)
                    alg_dframe_collections = join_dataset_obj(
                        alg_dset_collections, alg_time_collections, alg_metrics_collections)

                    # save ancillary datasets in pickle format
                    if alg_dframe_collections is not None:
                        folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_raw_step)
                        make_folder(folder_name_anc)
                        write_file_obj(file_path_anc_pnt_raw_step, alg_dframe_collections)

                        # info end collection group
                        alg_logger.info(' ------> (3) Collections group ... DONE')

                    else:

                        # info end collection group
                        alg_logger.info(' ------> (3) Collections group ... SKIPPED. Data are not available')

                    # info end cell group
                    alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) +
                                    '" ... DONE')

                else:

                    # info end cell group
                    alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) +
                                    '" ... SKIPPED. Datasets previously saved')

            # info end time
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... DONE')

        # info end method
        alg_logger.info(' ---> Organize dynamic datasets ... DONE')

    # method to analyze data points
    def analyze_data_points(self):

        # info start method
        alg_logger.info(' ---> Analyze dynamic data points ... ')

        # get file path
        file_path_anc_pnt_raw_tmpl = self.file_path_anc_pnt_raw_tmpl
        file_path_anc_pnt_def_tmpl = self.file_path_anc_pnt_def_tmpl
        # get cells
        alg_grid_cells = self.grid_cells
        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # set ancillary flag
        reset_datasets_anc_pnt_def = self.reset_datasets_anc_pnt_def

        # fx settings
        fx_settings_dset_scale = self.fx_settings_dset_scale
        fx_settings_dset_weigh = self.fx_settings_dset_weigh

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # iterate over cell(s)
            for alg_cell_group in alg_grid_cells:

                # info start cell group
                alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) + '" ... ')

                # define ancillary file name raw and dwf
                file_path_anc_pnt_raw_step = fill_path_with_tags(
                    file_path_anc_pnt_raw_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)
                file_path_anc_pnt_def_step = fill_path_with_tags(
                    file_path_anc_pnt_def_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                # reset ancillary file(s) (if needed)
                if reset_datasets_anc_pnt_def:
                    if os.path.exists(file_path_anc_pnt_def_step):
                        os.remove(file_path_anc_pnt_def_step)

                # check file ancillary availability
                if not os.path.exists(file_path_anc_pnt_def_step):
                    # check file source availability
                    if os.path.exists(file_path_anc_pnt_raw_step):

                        # get variable(s) obj
                        obj_data_cell_src = read_file_obj(file_path_anc_pnt_raw_step)

                        # info start scale datasets
                        alg_logger.info(' ------> (1) Scale datasets ... ')
                        # get fx scale settings
                        (fx_active_dset_scale, fx_name_dset_scale,
                            fx_vars_dset_scale, fx_params_dset_scale) = get_fx_settings(fx_settings_dset_scale)
                        # get fx scale method
                        fx_handle_dset_scale = get_fx_method(fx_name_dset_scale, fx_methods_datasets)

                        # flag to activate scale part
                        if fx_active_dset_scale:
                            # call scale method
                            obj_data_cell_scaled = fx_handle_dset_scale(
                                obj_data=obj_data_cell_src,
                                variables_data=fx_vars_dset_scale, parameters_data=fx_params_dset_scale)
                            # info end scale datasets
                            alg_logger.info(' ------> (1) Scale datasets ... DONE')
                        else:
                            # info end scale datasets
                            alg_logger.info(' ------> (1) Scale datasets ... FAILED')
                            alg_logger.error(' ===> Scale datasets part is needed by the algorithm to correctly run')
                            raise RuntimeError('Active this part in the configuration file')

                        # info start weigh datasets
                        alg_logger.info(' ------> (2) Weigh datasets ... ')
                        # get fx weigh settings
                        (fx_active_dset_weigh, fx_name_dset_weigh,
                            fx_vars_dset_weigh, fx_params_dset_weigh) = get_fx_settings(fx_settings_dset_weigh)
                        # get fx weigh method
                        fx_handle_dset_weigh = get_fx_method(fx_name_dset_weigh, fx_methods_datasets)

                        # flag to activate weigh part
                        if fx_active_dset_weigh:
                            # call weigh method
                            obj_data_cell_weighted = fx_handle_dset_weigh(
                                obj_data=obj_data_cell_scaled,
                                variables_data=fx_vars_dset_weigh,
                                parameters_data=fx_params_dset_weigh['lut'],
                                active_ref_k1=fx_params_dset_weigh['flags']['active_ref_k1'],
                                active_ref_k2=fx_params_dset_weigh['flags']['active_ref_k2'],
                                active_ref=fx_params_dset_weigh['flags']['active_ref'],
                                active_k1=fx_params_dset_weigh['flags']['active_k1'],
                                active_k2=fx_params_dset_weigh['flags']['active_k2'])
                            # info end weigh datasets
                            alg_logger.info(' ------> (2) Weigh datasets ... DONE')
                        else:
                            # info end weigh datasets
                            alg_logger.info(' ------> (2) Weigh datasets ... FAILED')
                            alg_logger.error(' ===> Weigh datasets part is needed by the algorithm to correctly run')
                            raise RuntimeError('Active this part in the configuration file')

                        # save ancillary datasets in pickle format
                        folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_def_step)
                        make_folder(folder_name_anc)
                        write_file_obj(file_path_anc_pnt_def_step, obj_data_cell_weighted)

                        # info end cell group
                        alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) + '" ... DONE')

                    else:
                        # info end cell group
                        alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) +
                                        '" ... FAILED. Datasets are not available')

                else:
                    # info end cell group
                    alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) +
                                    '" ... SKIPPED. Datasets previously saved')

            # info end time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... DONE')

        # info end method
        alg_logger.info(' ---> Analyze dynamic data points ... DONE')

    # method to analyze data maps
    def analyze_data_maps(self):

        # info start method
        alg_logger.info(' ---> Analyze dynamic data maps ... ')

        # get file path
        file_path_anc_pnt_def_tmpl = self.file_path_anc_pnt_def_tmpl
        file_path_anc_map_tmpl = self.file_path_anc_map_tmpl

        # get grid data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs

        # get cells
        alg_grid_cells = self.grid_cells
        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # fx settings
        fx_settings_common_resample = self.fx_settings_common_resample
        fx_settings_common_filter = self.fx_settings_common_filter
        fx_settings_common_mask = self.fx_settings_common_mask

        # flags
        reset_datasets_anc_map = self.reset_datasets_anc_map

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary file name
            file_path_anc_map_step = fill_path_with_tags(
                file_path_anc_map_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)

            # flag to reset ancillary datasets
            if reset_datasets_anc_map:
                if os.path.exists(file_path_anc_map_step):
                    os.remove(file_path_anc_map_step)

            # check file ancillary availability
            if not os.path.exists(file_path_anc_map_step):

                # info start merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... ')

                # iterate over cell(s)
                obj_data_cell_merge, obj_data_cell_list = None, []
                for alg_cell_group in alg_grid_cells:

                    # info start cell group
                    alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... ')

                    # define ancillary file name
                    file_path_anc_pnt_def_step = fill_path_with_tags(
                        file_path_anc_pnt_def_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                        tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                    # check file source availability
                    if os.path.exists(file_path_anc_pnt_def_step):

                        # get variable(s) obj
                        obj_data_cell_src = read_file_obj(file_path_anc_pnt_def_step)
                        # remove nan(s) obj
                        obj_data_cell_src = remove_nans(obj_data_cell_src, keep_finite=True)

                        # merge variable(s) obj
                        if obj_data_cell_merge is None:
                            obj_data_cell_merge = deepcopy(obj_data_cell_src)
                        else:
                            obj_data_cell_merge = pd.concat([obj_data_cell_merge, obj_data_cell_src],
                                                            axis=0, ignore_index=True)

                        # append cells idx
                        obj_data_cell_list.append(alg_cell_group)

                        # info end cell group
                        alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... DONE')

                    else:
                        # info end cell group
                        alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) +
                                        '" ... FAILED. Datasets are not available')

                # info end merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... DONE')

                # info start resample datasets
                alg_logger.info(' -----> (2) Resample datasets ... ')
                # get fx resample settings
                (fx_active_common_resample, fx_name_common_resample,
                 fx_vars_common_resample, fx_params_common_resample) = get_fx_settings(fx_settings_common_resample)
                # get fx resample handle
                fx_handle_common_resample = get_fx_method(fx_name_common_resample, fx_methods_common)

                # flag to activate resample part
                if fx_active_common_resample:
                    # call resample method
                    obj_data_cell_resample = fx_handle_common_resample(
                        obj_in=obj_data_cell_merge,
                        geo_x_out=grid_geo_data[geo_var_name_x].values,
                        geo_y_out=grid_geo_data[geo_var_name_y].values,
                        var_mapping_in=fx_vars_common_resample['in'],
                        var_mapping_out=fx_vars_common_resample['out'],
                        var_name_geo_x=geo_var_name_x, var_name_geo_y=geo_var_name_y,
                        coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                        dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y,
                        resampling_max_distance=fx_params_common_resample['max_distance'],
                        resampling_min_neighbours=fx_params_common_resample['min_neighbours'],
                        resampling_neighbours=fx_params_common_resample['neighbours']
                    )

                    # info end resample datasets
                    alg_logger.info(' -----> (2) Resample datasets ... DONE')
                else:
                    # info end resample datasets
                    alg_logger.info(' -----> (2) Resample datasets ... FAILED')
                    alg_logger.error(' ===> Resample datasets part is needed by the algorithm to correctly run')
                    raise RuntimeError('Active this part in the configuration file')

                # info start mask datasets
                alg_logger.info(' -----> (3) Filter datasets ... ')
                # get fx mask settings
                (fx_active_common_filter, fx_name_common_filter,
                 fx_vars_common_filter, fx_params_common_filter) = get_fx_settings(fx_settings_common_filter)
                # get fx mask handle
                fx_handle_common_filter = get_fx_method(fx_name_common_filter, fx_methods_common)

                # flag to activate filter part
                if fx_active_common_filter and fx_active_common_filter:
                    # call mask method
                    obj_data_cell_filter = fx_handle_common_filter(
                        obj_data_cell_resample,
                        var_mapping_in=fx_vars_common_filter['in'],
                        var_mapping_out=fx_vars_common_filter['out'],
                        filter_type=fx_params_common_filter['kernel_type'],
                        filter_pixels_std=fx_params_common_filter['kernel_pixels_std'],
                        filter_mode=fx_params_common_filter['kernel_mode'])
                    # info end filter datasets
                    alg_logger.info(' -----> (3) Filter datasets ... DONE')
                else:
                    # info end filter datasets
                    obj_data_cell_filter = deepcopy(obj_data_cell_resample)
                    alg_logger.info(' -----> (3) Filter datasets ... SKIPPED. Flag filter is not activated')

                # info start mask datasets
                alg_logger.info(' -----> (4) Mask datasets ... ')
                # get fx mask settings
                (fx_active_common_mask, fx_name_common_mask,
                 fx_vars_common_mask, fx_params_common_mask) = get_fx_settings(fx_settings_common_mask)
                # get fx mask handle
                fx_handle_common_mask = get_fx_method(fx_name_common_mask, fx_methods_common)

                # flag to activate mask part
                if fx_active_common_mask and fx_active_common_mask:
                    # call mask method
                    obj_data_cell_mask = fx_handle_common_mask(
                        obj_data_cell_filter, grid_geo_data,
                        var_mapping_in=fx_vars_common_mask['in'],
                        var_mapping_out=fx_vars_common_mask['out'],
                        mask_value_min=fx_params_common_mask['min_value'],
                        mask_value_max=fx_params_common_mask['max_value'])
                    # info end mask datasets
                    alg_logger.info(' -----> (4) Mask datasets ... DONE')
                else:
                    # info end mask datasets
                    obj_data_cell_mask = deepcopy(obj_data_cell_resample)
                    alg_logger.info(' -----> (4) Mask datasets ... SKIPPED. Flag mask is not activated')

                # check data
                if check_data(obj_data_cell_mask, thr_data='any'):

                    # save data collection in pickle format
                    folder_name_anc, file_name_anc = os.path.split(file_path_anc_map_step)
                    make_folder(folder_name_anc)
                    write_file_obj(file_path_anc_map_step, obj_data_cell_mask)

                    # info end time
                    alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                    '" ... DONE')
                else:
                    # info end time
                    alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                    '" ... FAILED. Datasets are defined by NoneType')

        # info end method
        alg_logger.info(' ---> Analyze dynamic data maps ... DONE')

    # method to dump data
    def dump_data(self):

        # info start method
        alg_logger.info(' ---> Dump dynamic datasets ... ')

        # get file path
        file_path_anc_map_tmpl, file_path_dst_tmpl = self.file_path_anc_map_tmpl, self.file_path_dst_tmpl
        # get grid info
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        # get variable(s)
        variable_dst = self.variable_dst

        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # flags
        reset_datasets_dst = self.reset_datasets_dst

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary file name
            file_path_anc_map_step = fill_path_with_tags(
                file_path_anc_map_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)
            # define destination file name
            file_path_dst_step = fill_path_with_tags(
                file_path_dst_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)

            # reset destination file (if needed)
            if reset_datasets_dst:
                if os.path.exists(file_path_dst_step):
                    os.remove(file_path_dst_step)

            # check file ancillary availability
            if not os.path.exists(file_path_dst_step):

                # check file source availability
                if os.path.exists(file_path_anc_map_step):

                    # method to get data obj
                    obj_data_map = read_file_obj(file_path_anc_map_step)

                    # check destination format
                    if self.format_dst == 'netcdf':

                        # method to organize netcdf dataset
                        variable_dset = organize_file_nc(obj_variable=obj_data_map, obj_time=alg_time_group)
                        # method to write netcdf dataset
                        folder_name_dst_step, file_name_dst_step = os.path.split(file_path_dst_step)
                        make_folder(folder_name_dst_step)
                        write_file_nc(file_path_dst_step, variable_dset)

                    elif (self.format_dst == 'tiff') or (self.format_dst == 'tif'):

                        # method to organize tiff dataset
                        variable_data, variable_attrs = organize_file_tiff(
                            obj_data_map, alg_time_group, obj_variable=variable_dst,
                            obj_transform=grid_geo_attrs['transform'], obj_proj=grid_geo_attrs['crs'])
                        # method to write tiff dataset
                        folder_name_dst_step, file_name_dst_step = os.path.split(file_path_dst_step)
                        make_folder(folder_name_dst_step)
                        write_file_tiff(
                            file_name=file_path_dst_step, file_data=variable_data,
                            **variable_attrs)

                    else:
                        alg_logger.error(' ===> Destination format "' + self.format_dst + '" is not supported')
                        raise NotImplemented('Case not implemented yet')

                    # info end time
                    alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                    '" ... DONE')

                else:
                    # info end time
                    alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                    '" ... FAILED. Datasets are not available')

            else:
                # info end time
                alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                '" ... SKIPPED. Datasets previously computed')

        # info end method
        alg_logger.info(' ---> Dump dynamic datasets ... DONE')

# -------------------------------------------------------------------------------------
