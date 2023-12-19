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

import numpy as np
import pandas as pd

from copy import deepcopy

from lib_info_args import time_format_algorithm
from lib_info_args import (geo_dim_name_x, geo_dim_name_y, geo_coord_name_x, geo_coord_name_y,
                           geo_var_name_x, geo_var_name_y)

from lib_utils_generic import make_folder
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_cell import get_grid_data
from lib_data_io_nc import organize_file_nc, write_file_nc
from lib_data_io_tiff import organize_file_tiff, write_file_tiff

from lib_info_args import logger_name
from lib_utils_io import fill_path_with_tags
from lib_utils_time import split_time_parts
from lib_utils_geo import resample_grid_to_points

from lib_fx_utils import get_fx_method, get_fx_settings, remove_nans, check_data

import lib_fx_methods_datasets as fx_methods_datasets
import lib_fx_methods_common as fx_methods_common

from cpl_data_dynamic import select_datasets_fields, join_dataset_obj
from cpl_data_dynamic import CouplerDataset, CouplerAncillary

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
from lib_utils_plot import plot_data_2d
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
                 tag_section_datasets='datasets', tag_section_cells='cells',
                 tag_section_log='log', tag_section_tmp='tmp'):

        self.alg_time_reference = alg_time_reference
        self.alg_time_datasets = alg_time_datasets

        self.alg_static = alg_static

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_template_dset = alg_settings[tag_section_template]['datasets']
        self.alg_template_time = alg_settings[tag_section_template]['time']
        self.alg_method_cells = alg_settings[tag_section_cells]
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

        self.max_distance_ref_ancillary = self.alg_method_cells['max_distance']['ancillary']
        self.max_distance_ref = self.alg_method_cells['max_distance']['ref']
        self.max_distance_ref_k1 = self.alg_method_cells['max_distance']['k1']
        self.max_distance_ref_k2 = self.alg_method_cells['max_distance']['k2']

        self.dset_max_distance = self.alg_method_cells['max_distance']
        self.dset_max_timedelta = self.alg_method_cells['max_timedelta']
        self.dset_remap_flag = self.alg_method_cells['remap_datasets_flag']
        self.dset_remap_name = self.alg_method_cells['remap_datasets_name']

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ---> Organize dynamic datasets ... ')

        # get grid data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        geo_x_1d = grid_geo_data['longitude'].values
        geo_y_1d = grid_geo_data['latitude'].values
        geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
        geo_mask = grid_geo_data.values

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
                    alg_dset_collections, alg_time_collections, alg_map_collection = {}, {}, {}
                    var_data_map = None
                    for alg_dset_name in self.alg_dset_list:

                        # info start dataset group
                        alg_logger.info(' -------> Data "' + alg_dset_name + '" ... ')

                        if alg_dset_name == 'ref':
                            dset_time_reference = alg_time_group
                        else:
                            dset_time_reference = alg_time_collections['ref']

                        # get datasets field(s)
                        alg_dset_fields = alg_datasets_src_data[alg_dset_name]

                        # prepare datasets args
                        alg_dset_args = select_datasets_fields(alg_dset_fields)

                        # initialize dataset class
                        coupler_dataset = CouplerDataset(
                            time_data=alg_time_group, cell_data=alg_cell_group, grid_ref=self.grid_ref,
                            dataset_name=alg_dset_name,
                            dataset_template=alg_template_dset, time_template=alg_template_time,
                            **alg_dset_args)

                        # method to get data
                        alg_dset_src, alg_time_src = coupler_dataset.get_data()

                        # get time tolerance(s)
                        dset_time_period, dset_time_frequency = split_time_parts(self.dset_max_timedelta[alg_dset_name])

                        # method to organize time
                        alg_dset_select, alg_time_select = coupler_dataset.organize_time(
                            alg_dset_src,
                            time_reference_dset=dset_time_reference, time_reference_group=alg_time_group,
                            time_tolerance_period=dset_time_period, time_tolerance_frequency=dset_time_frequency)

                        # method to organize data
                        alg_dset_dst = coupler_dataset.organize_data(alg_dset_select, alg_time_select)

                        # store data in the collections object
                        alg_dset_collections[alg_dset_name] = alg_dset_dst
                        alg_time_collections[alg_dset_name] = alg_time_select

                        var_data_map = np.zeros(shape=(geo_x_2d.shape[0], geo_y_2d.shape[1]))
                        var_data_map[:] = np.nan
                        if alg_dset_name == 'k1':
                            var_name = 'soil_moisture_k1'
                        elif alg_dset_name == 'k2':
                            var_name = 'soil_moisture_k2'
                        elif alg_dset_name == 'ref':
                            var_name = 'soil_moisture_ref'
                        else:
                            alg_logger.error(' ===> Dataset "' + alg_dset_name + '" is not available')
                            raise NotImplemented('Case not implemented yet')

                        var_data_grid, idx_data_grid = get_grid_data(alg_dset_dst,
                                                                     geo_mask, geo_x_2d, geo_y_2d,
                                                                     var_name=var_name)

                        if var_data_grid is not None:
                            var_data_map[idx_data_grid[:, 0], idx_data_grid[:, 1]] = var_data_grid[idx_data_grid[:, 0], idx_data_grid[: ,1]]
                            # plot_data_2d(var_data_map, geo_x_2d, geo_y_2d)
                            alg_map_collection[alg_dset_name] = deepcopy(var_data_map)

                        # info start dataset group
                        alg_logger.info(' -------> Data "' + alg_dset_name + '" ... DONE')

                    # info start dataset group
                    alg_logger.info(' ------> (2) Dataset Group ... DONE')

                    # info start collection group
                    alg_logger.info(' ------> (3) Collections group ... ')

                    # method to join dataset(s)
                    alg_dframe_collections = join_dataset_obj(
                        alg_dset_collections, alg_time_collections, alg_metrics_collections,
                        max_dist_k1=self.dset_max_distance['k1'], max_dist_k2=self.dset_max_distance['k2'],
                        max_dist_anc=self.dset_max_distance['ancillary'])

                    ''' debug
                    import numpy as np
                    from repurpose.resample import resample_to_grid

                    var_geo_x = alg_dframe_collections['longitude'].values
                    var_geo_y = alg_dframe_collections['latitude'].values

                    var_data_1 = alg_dframe_collections['soil_moisture_ref'].values
                    var_data_2 = alg_dframe_collections['soil_moisture_k1'].values
                    var_data_3 = alg_dframe_collections['soil_moisture_k2'].values

                    geo_x_1d = grid_geo_data['longitude'].values
                    geo_y_1d = grid_geo_data['latitude'].values
                    geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
                    geo_mask = grid_geo_data.values

                    values_obj = resample_to_grid(
                        {'data': var_data_2},
                        var_geo_x, var_geo_y, geo_x_2d, geo_y_2d,
                        search_rad=self.max_distance_ref_k1, fill_values=np.nan,
                        min_neighbours=1, neighbours=8)
                    var_data = values_obj['data']
                    var_data[geo_mask == 0] = np.nan

                    plot_data_2d(var_data, geo_x_2d, geo_y_2d)
                    '''

                    # save ancillary datasets in pickle format
                    if alg_dframe_collections is not None:
                        folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_raw_step)
                        make_folder(folder_name_anc)

                        alg_obj_collection = {'data': alg_dframe_collections,
                                              'time': alg_time_collections,
                                              'map': alg_map_collection}
                        write_file_obj(file_path_anc_pnt_raw_step, alg_obj_collection)

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

        # get grid data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        geo_x_1d = grid_geo_data['longitude'].values
        geo_y_1d = grid_geo_data['latitude'].values
        geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
        geo_mask = grid_geo_data.values

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

            # prepare datasets
            file_path_anc_pnt_raw_list, file_path_anc_pnt_def_list = [], []
            file_ancillary_update = False
            for alg_cell_group in alg_grid_cells:

                # define ancillary file name raw and dwf
                file_path_anc_pnt_raw_step = fill_path_with_tags(
                    file_path_anc_pnt_raw_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                file_path_anc_pnt_def_step = fill_path_with_tags(
                    file_path_anc_pnt_def_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                file_path_anc_pnt_raw_list.append(file_path_anc_pnt_raw_step)
                file_path_anc_pnt_def_list.append(file_path_anc_pnt_def_step)

            for file_path_anc_pnt_def_tmp in file_path_anc_pnt_def_list:
                if reset_datasets_anc_pnt_def:
                    if os.path.exists(file_path_anc_pnt_def_tmp):
                        os.remove(file_path_anc_pnt_def_tmp)
                    file_ancillary_update = True
                else:
                    if not os.path.exists(file_path_anc_pnt_def_tmp):
                        file_ancillary_update = True
                        break

            # flag to update datasets
            if file_ancillary_update:

                # info start merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... ')

                # iterate over cell(s)
                obj_data_cell_merged, obj_time_cell_merged, obj_data_cell_list = None, None, []
                obj_data_map_merge, obj_data_time_merged = {}, {}
                for alg_cell_group, file_path_anc_pnt_raw_step in zip(alg_grid_cells, file_path_anc_pnt_raw_list):

                    # info start cell group
                    alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... ')

                    # check file source availability
                    if os.path.exists(file_path_anc_pnt_raw_step):

                        # get variable(s) obj
                        obj_cell_src = read_file_obj(file_path_anc_pnt_raw_step)
                        obj_data_cell_src, obj_time_cell_src = obj_cell_src['data'], obj_cell_src['time']
                        obj_data_map_src = obj_cell_src['map']

                        # info start points
                        alg_logger.info(' -------> Points ... ')

                        # merge variable(s) obj
                        if obj_data_cell_merged is None:
                            obj_data_cell_merged = deepcopy(obj_data_cell_src)
                        else:
                            obj_data_cell_merged = pd.concat([obj_data_cell_merged, obj_data_cell_src],
                                                             axis=0, ignore_index=True)
                        # merge time obj
                        if obj_time_cell_merged is None:
                            obj_time_cell_merged = deepcopy(obj_time_cell_src)
                        else:
                            pass

                        # append cells idx
                        obj_data_cell_list.append(alg_cell_group)
                        # info start points
                        alg_logger.info(' -------> Points ... DONE')

                        # info start time
                        alg_logger.info(' -------> Time ... ')
                        obj_data_time_merged[alg_cell_group] = obj_time_cell_src
                        # info end time
                        alg_logger.info(' -------> Time ... DONE')

                        # info start maps
                        alg_logger.info(' -------> Maps ... ')
                        # iterate over dataset(s)
                        for alg_dset_name in self.alg_dset_list:

                            # info start dataset group
                            alg_logger.info(' --------> Data "' + alg_dset_name + '" ... ')

                            if alg_dset_name in list(obj_data_map_src.keys()):

                                obj_data_map_step = obj_data_map_src[alg_dset_name]

                                if obj_data_map_step is not None:
                                    if alg_dset_name not in list(obj_data_map_merge.keys()):
                                        var_data_map = np.zeros(shape=(geo_x_2d.shape[0], geo_y_2d.shape[1]))
                                        var_data_map[:] = np.nan
                                    else:
                                        var_data_map = deepcopy(obj_data_map_merge[alg_dset_name])

                                    var_data_map[np.isfinite(obj_data_map_step)] = obj_data_map_step[np.isfinite(obj_data_map_step)]
                                    obj_data_map_merge[alg_dset_name] = var_data_map

                                    #plot_data_2d(var_data_map, geo_x_2d, geo_y_2d)

                                    # info end dataset group
                                    alg_logger.info(' --------> Data "' + alg_dset_name +
                                                    '" ... DONE')

                                else:
                                    # info end dataset group
                                    alg_logger.info(' --------> Data "' + alg_dset_name +
                                                    '" ... SKIPPED. Map is defined by NoneType')

                            else:
                                # info end dataset group
                                alg_logger.info(' --------> Data "' + alg_dset_name +
                                                '" ... SKIPPED. All maps are defined by NoneType')

                        # info start maps
                        alg_logger.info(' -------> Maps ... DONE')

                        # info end cell group
                        alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... DONE')

                    else:
                        # info end cell group
                        alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) +
                                        '" ... FAILED. Datasets are not available')

                # info end merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... DONE')

                ''' debug
                plot_data_2d(obj_data_map_merge['ref'], geo_x_2d, geo_y_2d)
                plot_data_2d(obj_data_map_merge['k1'], geo_x_2d, geo_y_2d)
                '''

                # info start remap datasets
                alg_logger.info(' -----> (2) Remap datasets ... ')

                # iterate over dataset(s)
                for alg_dset_name in self.alg_dset_list:

                    # info start dataset remap
                    alg_logger.info(' ------> Data "' + alg_dset_name + '" ... ')

                    alg_remap_flag = False
                    if alg_dset_name in list(self.dset_remap_flag.keys()):
                        alg_remap_flag = self.dset_remap_flag[alg_dset_name]
                    else:
                        alg_logger.warning(' ===> Remap flag for dataset "' + alg_dset_name + '" is not available')
                    alg_remap_name = False
                    if alg_dset_name in list(self.dset_remap_name.keys()):
                        alg_remap_name = self.dset_remap_name[alg_dset_name]
                    else:
                        alg_logger.warning(' ===> Remap name for dataset "' + alg_dset_name + '" is not available')

                    alg_max_distance = 25000
                    if alg_dset_name in list(self.dset_max_distance.keys()):
                        alg_max_distance = self.dset_max_distance[alg_dset_name]
                    else:
                        alg_logger.warning(' ===> Max distance for dataset "' + alg_dset_name + '" is not available')

                    # remap cells using map datasets and merged points
                    if alg_remap_flag:
                        obj_data_cell_merged = resample_grid_to_points(
                            obj_data_map_merge, obj_data_cell_merged,
                            geo_mask, geo_x_2d, geo_y_2d,
                            var_name_dset=alg_dset_name, var_name_data=alg_remap_name,
                            search_rad=alg_max_distance, debug=False)

                        # info end dataset remap
                        alg_logger.info(' ------> Data "' + alg_dset_name + '" ... DONE')
                    else:
                        # info end dataset remap
                        alg_logger.info(' ------> Data "' + alg_dset_name + '" ... SKIPPED. Remap is not active')

                # info end remap datasets
                alg_logger.info(' -----> (2) Remap datasets ... DONE')

                # info start scale datasets
                alg_logger.info(' ------> (3) Scale datasets ... ')
                # get fx scale settings
                (fx_active_dset_scale, fx_name_dset_scale,
                 fx_vars_dset_scale, fx_params_dset_scale) = get_fx_settings(fx_settings_dset_scale)
                # get fx scale method
                fx_handle_dset_scale = get_fx_method(fx_name_dset_scale, fx_methods_datasets)

                # flag to activate scale part
                if fx_active_dset_scale:
                    # call scale method
                    obj_data_cell_scaled = fx_handle_dset_scale(
                        obj_data=obj_data_cell_merged,
                        variables_data=fx_vars_dset_scale, parameters_data=fx_params_dset_scale)
                    # info end scale datasets
                    alg_logger.info(' ------> (3) Scale datasets ... DONE')
                else:
                    # info end scale datasets
                    alg_logger.info(' ------> (3) Scale datasets ... FAILED')
                    alg_logger.error(' ===> Scale datasets part is needed by the algorithm to correctly run')
                    raise RuntimeError('Active this part in the configuration file')

                # info start weigh datasets
                alg_logger.info(' ------> (4) Weigh datasets ... ')
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
                    alg_logger.info(' ------> (4) Weigh datasets ... DONE')
                else:
                    # info end weigh datasets
                    alg_logger.info(' ------> (4) Weigh datasets ... FAILED')
                    alg_logger.error(' ===> Weigh datasets part is needed by the algorithm to correctly run')
                    raise RuntimeError('Active this part in the configuration file')

                # info start save datasets
                alg_logger.info(' ------> (5) Save datasets ... ')
                # iterate over cell(s)
                for alg_cell_group, file_path_anc_pnt_def_step in zip(alg_grid_cells, file_path_anc_pnt_def_list):

                    # info start cell group
                    alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) + '" ... ')

                    # check file def availability
                    if not os.path.exists(file_path_anc_pnt_def_step):

                        if alg_cell_group in list(obj_data_time_merged.keys()):

                            # get data and time information
                            obj_data_cell_select = obj_data_cell_weighted.loc[
                                obj_data_cell_weighted['cell'] == alg_cell_group]
                            obj_data_time_select = obj_data_time_merged[alg_cell_group]

                            # save ancillary datasets in pickle format
                            folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_def_step)
                            make_folder(folder_name_anc)

                            obj_cell_collection = {'data': obj_data_cell_select, 'time': obj_data_time_select}
                            write_file_obj(file_path_anc_pnt_def_step, obj_cell_collection)

                            # info end cell group
                            alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) + '" ... DONE')

                        else:
                            # info end cell group
                            alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) +
                                            '" ... SKIPPED. Datasets is not available')

                    else:
                        # info end cell group
                        alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) +
                                        '" ... FAILED. Datasets are not available')

                    # info end save datasets
                    alg_logger.info(' ------> (5) Save datasets ... DONE')

                # info end time group
                alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                '" ... DONE')

            else:
                # info start time group
                alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                '" ... SKIPPED. All datasets are already available')

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
                obj_data_cell_merge, obj_time_cell_merge, obj_data_cell_list = None, None, []
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
                        obj_cell_src = read_file_obj(file_path_anc_pnt_def_step)
                        obj_data_cell_src, obj_time_cell_src = obj_cell_src['data'], obj_cell_src['time']

                        # remove nan(s) obj
                        obj_data_cell_src = remove_nans(obj_data_cell_src, keep_finite=True)

                        # merge variable(s) obj
                        if obj_data_cell_merge is None:
                            obj_data_cell_merge = deepcopy(obj_data_cell_src)
                        else:
                            obj_data_cell_merge = pd.concat([obj_data_cell_merge, obj_data_cell_src],
                                                            axis=0, ignore_index=True)
                        # merge time obj
                        if obj_time_cell_merge is None:
                            obj_time_cell_merge = deepcopy(obj_time_cell_src)
                        else:
                            pass

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
                    obj_data_map_resample = fx_handle_common_resample(
                        var_data_obj=obj_data_cell_merge,
                        geo_mask_out=grid_geo_data.values,
                        geo_x_out=grid_geo_data[geo_var_name_x].values,
                        geo_y_out=grid_geo_data[geo_var_name_y].values,
                        var_mapping_in=fx_vars_common_resample['in'],
                        var_mapping_out=fx_vars_common_resample['out'],
                        var_name_geo_x=geo_var_name_x, var_name_geo_y=geo_var_name_y,
                        coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                        dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y,
                        resampling_max_distance=fx_params_common_resample['max_distance'],
                        resampling_min_neighbours=fx_params_common_resample['min_neighbours'],
                        resampling_neighbours=fx_params_common_resample['neighbours'],
                        resampling_remove_artifacts=fx_params_common_resample['remove_artifacts'],
                        resampling_datasets_artifacts=fx_params_common_resample['datasets_artifacts']
                    )

                    # info end resample datasets
                    alg_logger.info(' -----> (2) Resample datasets ... DONE')
                else:
                    # info end resample datasets
                    alg_logger.info(' -----> (2) Resample datasets ... FAILED')
                    alg_logger.error(' ===> Resample datasets part is needed by the algorithm to correctly run')
                    raise RuntimeError('Active this part in the configuration file')

                ''' debug
                var_data = obj_data_map_resample['soil_moisture_resampled']
                plot_data_2d(var_data)
                '''

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
                    obj_data_map_filter = fx_handle_common_filter(
                        obj_data_map_resample,
                        var_mapping_in=fx_vars_common_filter['in'],
                        var_mapping_out=fx_vars_common_filter['out'],
                        filter_type=fx_params_common_filter['kernel_type'],
                        filter_pixels_std=fx_params_common_filter['kernel_pixels_std'],
                        filter_mode=fx_params_common_filter['kernel_mode'])
                    # info end filter datasets
                    alg_logger.info(' -----> (3) Filter datasets ... DONE')
                else:
                    # info end filter datasets
                    obj_data_map_filter = deepcopy(obj_data_map_resample)
                    alg_logger.info(' -----> (3) Filter datasets ... SKIPPED. Flag filter is not activated')

                ''' debug
                var_data = obj_data_map_filter['soil_moisture_filtered']
                plot_data_2d(var_data)
                '''

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
                    obj_data_map_mask = fx_handle_common_mask(
                        obj_data_map_filter, grid_geo_data,
                        var_mapping_in=fx_vars_common_mask['in'],
                        var_mapping_out=fx_vars_common_mask['out'],
                        mask_value_min=fx_params_common_mask['min_value'],
                        mask_value_max=fx_params_common_mask['max_value'])
                    # info end mask datasets
                    alg_logger.info(' -----> (4) Mask datasets ... DONE')
                else:
                    # info end mask datasets
                    obj_data_map_mask = deepcopy(obj_data_map_resample)
                    alg_logger.info(' -----> (4) Mask datasets ... SKIPPED. Flag mask is not activated')

                ''' debug
                var_data_mf = obj_data_map_mask['soil_moisture_masked_filtered']
                var_data_mr = obj_data_map_mask['soil_moisture_masked_resampled']
                var_flags_m = obj_data_map_mask['flags_masked']
                var_flags_r = obj_data_map_mask['flags_resampled']
                plot_data_2d(var_data_mf)
                plot_data_2d(var_data_mr)
                plot_data_2d(var_flags_m)
                plot_data_2d(var_flags_r)
                '''

                # check data
                if check_data(obj_data_map_mask, thr_data='any'):

                    # save data collection in pickle format
                    folder_name_anc, file_name_anc = os.path.split(file_path_anc_map_step)
                    make_folder(folder_name_anc)

                    obj_map_collection = {'data': obj_data_map_mask, 'time': obj_time_cell_merge}
                    write_file_obj(file_path_anc_map_step, obj_map_collection)

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
                    obj_map = read_file_obj(file_path_anc_map_step)
                    obj_data_map, obj_time_map = obj_map['data'], obj_map['time']

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
