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
import xarray as xr

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
from lib_utils_geo import (resample_grid_to_points, resample_points_to_grid,
                           filter_points_nans, convert_grid_to_swath)

from lib_fx_utils import get_fx_method, get_fx_settings, remove_nans, check_data

import lib_fx_methods_datasets as fx_methods_datasets
import lib_fx_methods_common as fx_methods_common

from repurpose.resample import resample_to_grid
from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel

from cpl_data_dynamic import select_datasets_fields, join_dataset_obj
from cpl_data_dynamic import CouplerDataset, CouplerAncillary

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
from lib_utils_plot import plot_data_2d
from lib_utils_debug import convert_cell_to_grid

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
        self.alg_datasets_anc_pnt_src = alg_settings[tag_section_datasets]['dynamic']['ancillary']['points_source']
        self.alg_datasets_anc_pnt_dst_data = alg_settings[tag_section_datasets]['dynamic']['ancillary']['points_destination']['data']
        self.alg_datasets_anc_pnt_dst_mets = alg_settings[tag_section_datasets]['dynamic']['ancillary']['points_destination']['metrics']
        self.alg_datasets_anc_map = alg_settings[tag_section_datasets]['dynamic']['ancillary']['maps']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_settings[tag_section_log]
        self.alg_tmp = alg_settings[tag_section_tmp]

        self.tag_dset_data_name, self.tag_dset_metrics_name = 'soil_moisture', 'metrics'
        self.tag_dset_ref, self.tag_dset_k1, self.tag_dset_k2 = 'ref', 'k1', 'k2'

        self.tag_folder_name, self.tag_file_name = 'folder_name', 'file_name'
        self.tag_variable = 'variable'
        self.tag_time_frequency, self.tag_time_rounding = 'time_frequency', 'time_rounding'
        self.tag_time_period, self.tag_time_stamp = 'time_period', 'time_stamp'
        self.tag_format = 'format'

        self.reset_datasets_anc_pnt_src = self.alg_flags['reset_datasets_ancillary_points_src']
        self.reset_datasets_anc_pnt_dst_data = self.alg_flags['reset_datasets_ancillary_points_dst_data']
        self.reset_datasets_anc_pnt_dst_mets = self.alg_flags['reset_datasets_ancillary_points_dst_metrics']
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
        self.file_path_src_metrics_tmpl = os.path.join(
            self.folder_name_src_metrics_tmpl, self.file_name_src_metrics_tmpl)
        self.variable_src_metrics = self.alg_datasets_src_metrics[self.tag_variable]

        self.folder_name_anc_pnt_src_tmpl = self.alg_datasets_anc_pnt_src[self.tag_folder_name]
        self.file_name_anc_pnt_src_tmpl = self.alg_datasets_anc_pnt_src[self.tag_file_name]
        self.file_path_anc_pnt_src_tmpl = os.path.join(
            self.folder_name_anc_pnt_src_tmpl, self.file_name_anc_pnt_src_tmpl)

        self.folder_name_anc_pnt_dst_data_tmpl = self.alg_datasets_anc_pnt_dst_data[self.tag_folder_name]
        self.file_name_anc_pnt_dst_data_tmpl = self.alg_datasets_anc_pnt_dst_data[self.tag_file_name]
        self.file_path_anc_pnt_dst_data_tmpl = os.path.join(
            self.folder_name_anc_pnt_dst_data_tmpl, self.file_name_anc_pnt_dst_data_tmpl)
        self.folder_name_anc_pnt_dst_mets_tmpl = self.alg_datasets_anc_pnt_dst_mets[self.tag_folder_name]
        self.file_name_anc_pnt_dst_mets_tmpl = self.alg_datasets_anc_pnt_dst_mets[self.tag_file_name]
        self.file_path_anc_pnt_dst_mets_tmpl = os.path.join(
            self.folder_name_anc_pnt_dst_mets_tmpl, self.file_name_anc_pnt_dst_mets_tmpl)

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
        if self.tag_time_stamp in list(self.alg_datasets_dst.keys()):
            self.time_stamp_dst = self.alg_datasets_dst[self.tag_time_stamp]
        else:
            self.time_stamp_dst = 'cday'

        self.format_dst = self.alg_datasets_dst[self.tag_format]

        self.grid_geo_data, self.grid_geo_attrs = self.alg_static['geo']['data'], self.alg_static['geo']['attrs']
        self.grid_mask_data, self.grid_mask_attrs = self.alg_static['mask']['data'], self.alg_static['mask']['attrs']
        self.grid_rm_data, self.grid_rm_attrs = self.alg_static['removed']['data'], self.alg_static['removed']['attrs']
        self.grid_ref, self.grid_k1, self.grid_k2 = self.alg_static['ref'], self.alg_static['k1'], self.alg_static['k2']
        self.grid_cells = self.alg_static['cells']

        self.alg_dset_list = [self.tag_dset_ref, self.tag_dset_k1, self.tag_dset_k2]

        self.fx_settings_dset_scale = self.alg_methods_dset['fx_scale_data']
        self.fx_settings_dset_weigh = self.alg_methods_dset['fx_weigh_data']
        self.fx_settings_common_resample = self.alg_methods_common['fx_resample_data']
        self.fx_settings_common_filter = self.alg_methods_common['fx_filter_data']
        self.fx_settings_common_update = self.alg_methods_common['fx_update_data']
        self.fx_settings_common_mask = self.alg_methods_common['fx_mask_data']

        self.max_distance_ref_ancillary = self.alg_method_cells['max_distance']['ancillary']
        self.max_distance_ref = self.alg_method_cells['max_distance']['ref']
        self.max_distance_ref_k1 = self.alg_method_cells['max_distance']['k1']
        self.max_distance_ref_k2 = self.alg_method_cells['max_distance']['k2']

        self.dset_max_distance = self.alg_method_cells['max_distance']
        self.dset_max_timedelta = self.alg_method_cells['max_timedelta']
        self.dset_remap_flag = self.alg_method_cells['remap_datasets_flag']
        self.dset_remap_name = self.alg_method_cells['remap_datasets_name']

    # method to organize data points
    def organize_data_points(self):

        # info start method
        alg_logger.info(' ---> Organize dynamic points ... ')

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
        file_path_anc_pnt_src_tmpl = self.file_path_anc_pnt_src_tmpl
        file_path_anc_pnt_dst_data_tmpl = self.file_path_anc_pnt_dst_data_tmpl
        file_path_anc_pnt_dst_mets_tmpl = self.file_path_anc_pnt_dst_mets_tmpl
        file_path_anc_map_tmpl = self.file_path_anc_map_tmpl
        file_path_dst_tmpl = self.file_path_dst_tmpl

        # set ancillary flag
        reset_datasets_anc_pnt_src = self.reset_datasets_anc_pnt_src
        reset_datasets_anc_pnt_dst_data = self.reset_datasets_anc_pnt_dst_data
        reset_datasets_anc_pnt_dst_mets = self.reset_datasets_anc_pnt_dst_mets

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary points file name
            file_path_anc_pnt_dst_data_step = fill_path_with_tags(
                file_path_anc_pnt_dst_data_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                tmpl_tags_dset=alg_template_dset)
            file_path_anc_pnt_dst_mets_step = fill_path_with_tags(
                file_path_anc_pnt_dst_mets_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                tmpl_tags_dset=alg_template_dset)

            # define ancillary map file name
            file_path_anc_map_step = fill_path_with_tags(
                file_path_anc_map_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)
            # define destination file name
            file_path_dst_step = fill_path_with_tags(
                file_path_dst_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)

            # reset ancillary points data and destination file(s) (if needed)
            if reset_datasets_anc_pnt_dst_data:
                if os.path.exists(file_path_anc_pnt_dst_data_step):
                    os.remove(file_path_anc_pnt_dst_data_step)
                if os.path.exists(file_path_anc_map_step):
                    os.remove(file_path_anc_map_step)
                if os.path.exists(file_path_dst_step):
                    os.remove(file_path_dst_step)
            # reset ancillary points metrics file(s) (if needed)
            if reset_datasets_anc_pnt_dst_mets:
                if os.path.exists(file_path_anc_pnt_dst_mets_step):
                    os.remove(file_path_anc_pnt_dst_mets_step)

            # iterate over cell(s)
            for alg_cell_group in alg_grid_cells:

                # info start cell group
                alg_logger.info(' -----> Cell Group "' + str(alg_cell_group) + '" ... ')

                # define ancillary file name
                file_path_anc_pnt_src_step = fill_path_with_tags(
                    file_path_anc_pnt_src_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                # reset ancillary points raw file(s)
                if reset_datasets_anc_pnt_src:
                    if os.path.exists(file_path_anc_pnt_src_step):
                        os.remove(file_path_anc_pnt_src_step)

                # check availability of ancillary file
                if not os.path.exists(file_path_anc_pnt_src_step):

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

                        # get reference time (using the reference datasets selection)
                        if alg_dset_name == 'ref':
                            dset_time_reference = alg_time_group
                        else:
                            dset_time_reference = alg_time_collections['ref']

                        # check reference time format
                        if dset_time_reference is not None:

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

                            # check the selected datasets and time
                            if alg_dset_dst is not None and alg_time_select is not None:

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
                                    var_data_map[idx_data_grid[:, 0], idx_data_grid[:, 1]] = var_data_grid[
                                        idx_data_grid[:, 0], idx_data_grid[:, 1]]
                                    # plot_data_2d(var_data_map, geo_x_2d, geo_y_2d)
                                    alg_map_collection[alg_dset_name] = deepcopy(var_data_map)

                                # info end dataset group
                                alg_logger.info(
                                    ' -------> Data "' + alg_dset_name + '" ... DONE')
                            else:
                                # info end dataset group
                                alg_logger.info(
                                    ' -------> Data "' + alg_dset_name +
                                    '" ... SKIPPED. Datasets is not available due to the selected dataset time')
                                # store data in the collections object
                                alg_dset_collections[alg_dset_name] = None
                                alg_time_collections[alg_dset_name] = None

                        else:

                            # info end dataset group
                            alg_logger.info(
                                ' -------> Data "' + alg_dset_name +
                                '" ... SKIPPED. Datasets is not available due to the reference dataset time')
                            # store data in the collections object
                            alg_dset_collections[alg_dset_name] = None
                            alg_time_collections[alg_dset_name] = None

                    # info start dataset group
                    alg_logger.info(' ------> (2) Dataset Group ... DONE')

                    # info start collection group
                    alg_logger.info(' ------> (3) Collections group ... ')

                    # method to join dataset(s)
                    alg_dframe_collections = join_dataset_obj(
                        alg_dset_collections, alg_time_collections, alg_metrics_collections,
                        max_dist_k1=self.dset_max_distance['k1'], max_dist_k2=self.dset_max_distance['k2'],
                        max_dist_anc=self.dset_max_distance['ancillary'],
                        remap_k1_2_ref=self.dset_remap_flag['k1'], remap_k2_2_ref=self.dset_remap_flag['k2'],
                        tag_name_ref='ref', tag_name_k1='k1', tag_name_k2='k2',
                        )

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
                        folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_src_step)
                        make_folder(folder_name_anc)

                        alg_obj_collection = {'data': alg_dframe_collections,
                                              'time': alg_time_collections,
                                              'map': alg_map_collection}
                        write_file_obj(file_path_anc_pnt_src_step, alg_obj_collection)

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
        alg_logger.info(' ---> Organize dynamic points ... DONE')

    # method to analyze ancillary data points metrics
    def analyze_data_points_metrics(self):

        # info start method
        alg_logger.info(' ---> Analyze dynamic points metrics ... ')

        # get grid data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        geo_x_1d = grid_geo_data['longitude'].values
        geo_y_1d = grid_geo_data['latitude'].values
        geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
        geo_mask = grid_geo_data.values

        # get file path
        file_path_anc_pnt_src_tmpl = self.file_path_anc_pnt_src_tmpl
        file_path_anc_pnt_dst_mets_tmpl = self.file_path_anc_pnt_dst_mets_tmpl
        # get cells
        alg_grid_cells = self.grid_cells
        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # set ancillary flag
        reset_datasets_anc_pnt_src = self.reset_datasets_anc_pnt_src
        reset_datasets_anc_pnt_dst_mets = self.reset_datasets_anc_pnt_dst_mets

        # metrics variable name(s)
        variable_metrics_list = list(self.variable_src_metrics.values())

        # filter settings
        filter_type = 'box' # or 'gauss'
        box_width = 5
        gauss_stddev, gauss_mode = 2, 'center'

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary file name points destination metrics
            file_path_anc_pnt_dst_mets_step = fill_path_with_tags(
                file_path_anc_pnt_dst_mets_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                tmpl_tags_dset=alg_template_dset)

            # prepare datasets
            file_path_anc_pnt_src_list, alg_cell_list = [], []
            for alg_cell_group in alg_grid_cells:

                # define ancillary file name points source
                file_path_anc_pnt_src_step = fill_path_with_tags(
                    file_path_anc_pnt_src_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                # check availability of source file(s)
                if os.path.exists(file_path_anc_pnt_src_step):
                    file_path_anc_pnt_src_list.append(file_path_anc_pnt_src_step)
                    alg_cell_list.append(alg_cell_group)
                else:
                    # message of missing datasets
                    alg_logger.warning(' ===> Cell Group "' + str(alg_cell_group) +
                                       '" ... SKIPPED. Datasets are not available.')

            # check updating of ancillary file(s)
            if reset_datasets_anc_pnt_dst_mets or reset_datasets_anc_pnt_src:
                if os.path.exists(file_path_anc_pnt_dst_mets_step):
                    os.remove(file_path_anc_pnt_dst_mets_step)

            # check file availability
            if not os.path.exists(file_path_anc_pnt_dst_mets_step):

                # info start merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... ')

                # iterate over cell(s)
                obj_data_cell_merged, obj_time_cell_merged, obj_data_cell_list = None, None, []
                for alg_cell_group, file_path_anc_pnt_src_step in zip(alg_cell_list, file_path_anc_pnt_src_list):

                    # info start cell group
                    alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... ')

                    # check file source availability
                    if os.path.exists(file_path_anc_pnt_src_step):

                        # get variable(s) obj
                        obj_cell_src = read_file_obj(file_path_anc_pnt_src_step)
                        obj_data_cell_src, obj_time_cell_src = obj_cell_src['data'], obj_cell_src['time']

                        # info start points
                        alg_logger.info(' -------> Points ... ')

                        # add information to the obj
                        obj_cell_src['cell'] = alg_cell_group

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

                        # info end cell group
                        alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... DONE')

                    else:
                        # info end cell group
                        alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) +
                                        '" ... FAILED. Datasets are not available')

                # info end merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... DONE')

                # info start filter datasets
                alg_logger.info(' -----> (2) Filter datasets ... ')

                # create mask reference
                mask_2d, mask_geox, mask_geoy = resample_points_to_grid(
                    geo_mask.ravel(),
                    var_geox_1d_in=geo_x_2d.ravel(),
                    var_geoy_1d_in=geo_y_2d.ravel(),
                    var_geox_1d_out=obj_data_cell_merged['longitude'].values,
                    var_geoy_1d_out=obj_data_cell_merged['latitude'].values)
                mask_2d = np.flipud(mask_2d)

                # convert grid to swath
                mask_1d = convert_grid_to_swath(
                    mask_2d, mask_geox, mask_geoy,
                    obj_data_cell_merged['longitude'].values, obj_data_cell_merged['latitude'].values,
                    search_rad=50000, neighbours=1)

                # create swath dataframe
                obj_data_cell_metrics = pd.DataFrame(
                    data={'mask': mask_1d,
                          'longitude': obj_data_cell_merged['longitude'].values,
                          'latitude': obj_data_cell_merged['latitude'].values})

                ''' debug
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_metrics, search_rad=50000, var_name_data='mask')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                '''

                # iterate over metrics variable(s)
                variable_metrics_list = variable_metrics_list + ['cell']
                for variable_metrics_name in variable_metrics_list:

                    # info start variable
                    alg_logger.info(' ------> Variable "' + variable_metrics_name + '" ... ')

                    if variable_metrics_name not in ['longitude', 'latitude', 'location_id']:

                        # get variable data
                        cell_data_1d_raw = obj_data_cell_merged[variable_metrics_name].values
                        cell_geo_x_1d_raw = obj_data_cell_merged['longitude'].values
                        cell_geo_y_1d_raw = obj_data_cell_merged['latitude'].values

                        cell_data_1d_finite, cell_geo_x_1d_finite, cell_geo_y_1d_finite = filter_points_nans(
                            cell_data_1d_raw, cell_geo_x_1d_raw, cell_geo_y_1d_raw)

                        if variable_metrics_name == 'cell':
                            cell_data_1d_raw = cell_data_1d_raw.astype(int)

                        # resample points to grid
                        cell_data_2d_resampled, _, _ = resample_points_to_grid(
                            cell_data_1d_raw,
                            var_geox_1d_in=cell_geo_x_1d_raw,
                            var_geoy_1d_in=cell_geo_y_1d_raw,
                            var_geox_1d_out=obj_data_cell_merged['longitude'].values,
                            var_geoy_1d_out=obj_data_cell_merged['latitude'].values)
                        cell_data_2d_resampled = np.flipud(cell_data_2d_resampled)

                        if variable_metrics_name != 'cell':
                            # box values and create a filtered grid
                            if filter_type == 'box':
                                obj_kernel = Box2DKernel(box_width)
                            elif filter_type == 'gauss':
                                obj_kernel = Gaussian2DKernel(x_stddev=gauss_stddev, mode=gauss_mode)
                            else:
                                alg_logger.error(' ===> Filter type "' + filter_type + '" is not available')
                                raise NotImplemented('Case not implemented yet')

                            # filter grid using convolution
                            cell_data_2d_filtered = convolve(cell_data_2d_resampled, obj_kernel)
                            cell_data_2d_filtered[mask_2d == 0] = np.nan

                        else:
                            # skip filter for cell variable
                            cell_data_2d_filtered = deepcopy(cell_data_2d_resampled)
                            cell_data_2d_filtered[mask_2d == 0] = np.nan

                        # convert grid to swath
                        cell_data_1d_filtered = convert_grid_to_swath(
                            cell_data_2d_filtered, mask_geox, mask_geoy,
                            obj_data_cell_merged['longitude'].values, obj_data_cell_merged['latitude'].values,
                            search_rad=50000, neighbours=1)
                        if variable_metrics_name == 'cell':
                            cell_data_1d_filtered = cell_data_1d_filtered.astype(int)
                            cell_data_1d_filtered[cell_data_1d_filtered < 0] = -9999

                        # update swath dataframe
                        obj_data_cell_metrics[variable_metrics_name] = cell_data_1d_filtered

                        ''' debug
                        var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                            obj_data_cell_metrics, search_rad=50000, var_name_data=variable_metrics_name)
                        var_data_2d[mask_2d == 0] = np.nan
                        plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                        '''

                        # info end variable (done
                        alg_logger.info(' ------> Variable "' + variable_metrics_name + '" ... DONE')

                    else:
                        # info end variable (skipped)
                        alg_logger.info(' ------> Variable "' + variable_metrics_name + '" ... SKIPPED')

                # info start filter datasets
                alg_logger.info(' -----> (2) Filter datasets ... DONE')

                # info start save datasets
                alg_logger.info(' -----> (3) Save datasets ... ')

                # save ancillary datasets in pickle format
                folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_dst_mets_step)
                make_folder(folder_name_anc)

                # check data availability
                if not obj_data_cell_metrics.empty:

                    # organize data
                    obj_cell_collection = {'metrics': obj_data_cell_metrics}
                    write_file_obj(file_path_anc_pnt_dst_mets_step, obj_cell_collection)

                    # info end save datasets
                    alg_logger.info(' -----> (3) Save datasets ... DONE')

                else:

                    # info end save datasets
                    alg_logger.info(' -----> (3) Save datasets ... SKIPPED. Data are not available')

            else:
                # info start time group
                alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                '" ... SKIPPED. All datasets are already available')

        # info end method
        alg_logger.info(' ---> Analyze dynamic points metrics ... DONE')

    # method to analyze data points variables
    def analyze_data_points_variables(self):

        # info start method
        alg_logger.info(' ---> Analyze dynamic points variables ... ')

        # get grid data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        geo_x_1d = grid_geo_data['longitude'].values
        geo_y_1d = grid_geo_data['latitude'].values
        geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
        geo_mask = grid_geo_data.values

        # get file path
        file_path_anc_pnt_src_tmpl = self.file_path_anc_pnt_src_tmpl
        file_path_anc_pnt_dst_data_tmpl = self.file_path_anc_pnt_dst_data_tmpl
        file_path_anc_pnt_dst_mets_tmpl = self.file_path_anc_pnt_dst_mets_tmpl
        # get cells
        alg_grid_cells = self.grid_cells
        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # set ancillary flag
        reset_datasets_anc_pnt_src = self.reset_datasets_anc_pnt_src
        reset_datasets_anc_pnt_dst_data = self.reset_datasets_anc_pnt_dst_data

        # fx settings
        fx_settings_dset_scale = self.fx_settings_dset_scale
        fx_settings_dset_weigh = self.fx_settings_dset_weigh

        # metrics variable name(s)
        variable_metrics_list = list(self.variable_src_metrics.values())

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary file name points destination data and metrics
            file_path_anc_pnt_dst_data_step = fill_path_with_tags(
                file_path_anc_pnt_dst_data_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                tmpl_tags_dset=alg_template_dset)

            file_path_anc_pnt_dst_mets_step = fill_path_with_tags(
                file_path_anc_pnt_dst_mets_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                tmpl_tags_dset=alg_template_dset)

            # prepare datasets
            file_path_anc_pnt_src_list, alg_cell_list = [], []
            for alg_cell_group in alg_grid_cells:

                # define ancillary file name raw and dwf
                file_path_anc_pnt_src_step = fill_path_with_tags(
                    file_path_anc_pnt_src_tmpl, alg_time_group, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(alg_cell_group)}, tmpl_tags_dset=alg_template_dset)

                # check availability of source file(s)
                if os.path.exists(file_path_anc_pnt_src_step):
                    file_path_anc_pnt_src_list.append(file_path_anc_pnt_src_step)
                    alg_cell_list.append(alg_cell_group)
                else:
                    # message of missing datasets
                    alg_logger.warning(' ===> Cell Group "' + str(alg_cell_group) +
                                       '" ... SKIPPED. Datasets are not available.')

            # check updating of ancillary file(s)
            if reset_datasets_anc_pnt_dst_data or reset_datasets_anc_pnt_src:
                if os.path.exists(file_path_anc_pnt_dst_data_step):
                    os.remove(file_path_anc_pnt_dst_data_step)

            # check file availability
            if not os.path.exists(file_path_anc_pnt_dst_data_step):

                # info start merge datasets
                alg_logger.info(' -----> (1) Merge datasets ... ')

                # iterate over cell(s)
                obj_data_cell_merged, obj_time_cell_merged, obj_data_cell_list = None, None, []
                obj_data_map_merge, obj_data_time_merged = {}, {}
                for alg_cell_group, file_path_anc_pnt_src_step in zip(alg_cell_list, file_path_anc_pnt_src_list):

                    # info start cell group
                    alg_logger.info(' ------> Cell Group "' + str(alg_cell_group) + '" ... ')

                    # check file source availability
                    if os.path.exists(file_path_anc_pnt_src_step):

                        # get variable(s) obj
                        obj_cell_src = read_file_obj(file_path_anc_pnt_src_step)
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

                                    # plot_data_2d(var_data_map, geo_x_2d, geo_y_2d)

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
                if 'k1' in list(obj_data_map_merge.keys()):
                    plot_data_2d(obj_data_map_merge['k1'], geo_x_2d, geo_y_2d)
                if 'k2' in list(obj_data_map_merge.keys()):
                    plot_data_2d(obj_data_map_merge['k2'], geo_x_2d, geo_y_2d)
                '''

                # info start collect metrics
                alg_logger.info(' -----> (2) Collect metrics ... ')

                # check file source availability
                if os.path.exists(file_path_anc_pnt_dst_mets_step):
                    # get metrics obj
                    obj_cell_src = read_file_obj(file_path_anc_pnt_dst_mets_step)
                    obj_mets_cell_merged = obj_cell_src['metrics']

                    # update data cell merged
                    for variable_name in obj_mets_cell_merged.columns:
                        obj_data_cell_merged[variable_name] = obj_mets_cell_merged[variable_name].values

                else:
                    # info end cell group
                    alg_logger.error(' ===> File metrics "' + file_path_anc_pnt_dst_mets_step + '" is not available')
                    raise FileExistsError('File must be available to correctly run the algorithm')

                # create mask reference
                mask_2d, mask_geox, mask_geoy = resample_points_to_grid(
                    geo_mask.ravel(),
                    var_geox_1d_in=geo_x_2d.ravel(),
                    var_geoy_1d_in=geo_y_2d.ravel(),
                    var_geox_1d_out=obj_data_cell_merged['longitude'].values,
                    var_geoy_1d_out=obj_data_cell_merged['latitude'].values)
                mask_2d = np.flipud(mask_2d)

                # convert grid to swath
                mask_1d = convert_grid_to_swath(
                    mask_2d, mask_geox, mask_geoy,
                    obj_data_cell_merged['longitude'].values, obj_data_cell_merged['latitude'].values,
                    search_rad=50000, neighbours=1)

                # info end collect metrics
                alg_logger.info(' -----> (2) Collect metrics ... DONE')

                # info start remap datasets
                alg_logger.info(' -----> (3) Remap datasets ... ')

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

                ''' debug
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_mets_cell_merged, search_rad=50000, var_name_data='mean_ref')
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_mets_cell_merged, search_rad=50000, var_name_data='mean_k1')
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_mets_cell_merged, search_rad=50000, var_name_data='mean_k2')
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_ref')
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_k1')
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_k2')
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                '''

                # info end remap datasets
                alg_logger.info(' -----> (3) Remap datasets ... DONE')

                # info start scale datasets
                alg_logger.info(' -----> (4) Scale datasets ... ')
                # get fx scale settings
                (fx_active_dset_scale, fx_name_dset_scale,
                 fx_vars_dset_scale, fx_params_dset_scale) = get_fx_settings(fx_settings_dset_scale)
                # get fx scale method
                fx_handle_dset_scale = get_fx_method(fx_name_dset_scale, fx_methods_datasets)

                # flag to activate scale part
                if fx_active_dset_scale:
                    # call scale method
                    obj_data_cell_scaled = fx_handle_dset_scale(
                        obj_data=obj_data_cell_merged, obj_metrics=obj_mets_cell_merged,
                        variables_data=fx_vars_dset_scale, parameters_data=fx_params_dset_scale)
                    # info end scale datasets
                    alg_logger.info(' -----> (4) Scale datasets ... DONE')
                else:
                    # info end scale datasets
                    alg_logger.info(' -----> (4) Scale datasets ... FAILED')
                    alg_logger.error(' ===> Scale datasets part is needed by the algorithm to correctly run')
                    raise RuntimeError('Active this part in the configuration file')

                ''' debug
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_k1_scaled')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_k2_scaled')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                '''

                # info start weigh datasets
                alg_logger.info(' -----> (5) Weigh datasets ... ')
                # get fx weigh settings
                (fx_active_dset_weigh, fx_name_dset_weigh,
                 fx_vars_dset_weigh, fx_params_dset_weigh) = get_fx_settings(fx_settings_dset_weigh)
                # get fx weigh method
                fx_handle_dset_weigh = get_fx_method(fx_name_dset_weigh, fx_methods_datasets)

                # flag to activate weigh part
                if fx_active_dset_weigh:
                    # call weigh method
                    obj_data_cell_weighted = fx_handle_dset_weigh(
                        obj_data=obj_data_cell_scaled, obj_metrics=obj_mets_cell_merged,
                        variables_data=fx_vars_dset_weigh,
                        parameters_data=fx_params_dset_weigh['lut'],
                        active_ref_k1=fx_params_dset_weigh['flags']['active_ref_k1'],
                        active_ref_k2=fx_params_dset_weigh['flags']['active_ref_k2'],
                        active_ref=fx_params_dset_weigh['flags']['active_ref'],
                        active_k1=fx_params_dset_weigh['flags']['active_k1'],
                        active_k2=fx_params_dset_weigh['flags']['active_k2'])
                    # info end weigh datasets
                    alg_logger.info(' -----> (5) Weigh datasets ... DONE')
                else:
                    # info end weigh datasets
                    alg_logger.info(' -----> (5) Weigh datasets ... FAILED')
                    alg_logger.error(' ===> Weigh datasets part is needed by the algorithm to correctly run')
                    raise RuntimeError('Active this part in the configuration file')

                ''' debug
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_weighted')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_ref_weighted')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_k1_weighted')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                var_data_2d, geo_x_2d, geo_y_2d = convert_cell_to_grid(
                    obj_data_cell_merged, search_rad=50000, var_name_data='soil_moisture_k2_weighted')
                var_data_2d[mask_2d == 0] = np.nan
                plot_data_2d(var_data_2d, geo_x_2d, geo_y_2d)
                '''

                # info start save datasets
                alg_logger.info(' -----> (6) Save datasets ... ')

                # save ancillary datasets in pickle format
                folder_name_anc, file_name_anc = os.path.split(file_path_anc_pnt_dst_data_step)
                make_folder(folder_name_anc)

                # check data availability
                if not obj_data_cell_weighted.empty:

                    # organize data
                    obj_cell_collection = {'data': obj_data_cell_weighted, 'time': obj_data_time_merged}
                    write_file_obj(file_path_anc_pnt_dst_data_step, obj_cell_collection)

                    # info end save datasets
                    alg_logger.info(' -----> (6) Save datasets ... DONE')

                else:

                    # info end save datasets
                    alg_logger.info(' -----> (6) Save datasets ... SKIPPED. Data are not available')

                # info end time group
                alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                '" ... DONE')

            else:
                # info start time group
                alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                '" ... SKIPPED. All datasets are already available')

        # info end method
        alg_logger.info(' ---> Analyze dynamic points variables ... DONE')

    # method to analyze data maps
    def analyze_data_maps(self):

        # info start method
        alg_logger.info(' ---> Analyze dynamic maps ... ')

        # get grid data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        geo_x_1d = grid_geo_data['longitude'].values
        geo_y_1d = grid_geo_data['latitude'].values
        geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
        geo_mask = grid_geo_data.values

        # get file path
        file_path_anc_pnt_dst_data_tmpl = self.file_path_anc_pnt_dst_data_tmpl
        file_path_anc_map_tmpl = self.file_path_anc_map_tmpl

        # get geo data and attributes
        grid_geo_data, grid_geo_attrs = self.grid_geo_data, self.grid_geo_attrs
        # get mask data and attributes
        grid_mask_data, grid_mask_attrs = self.grid_mask_data, self.grid_mask_attrs
        # get rm data and attributes
        grid_rm_data, grid_rm_attrs = self.grid_rm_data, self.grid_rm_attrs

        # get cells
        alg_grid_cells = self.grid_cells
        # get dataset and time templates
        alg_template_dset, alg_template_time = self.alg_template_dset, self.alg_template_time

        # fx settings
        fx_settings_common_resample = self.fx_settings_common_resample
        fx_settings_common_filter = self.fx_settings_common_filter
        fx_settings_common_update = self.fx_settings_common_update
        fx_settings_common_mask = self.fx_settings_common_mask

        # flags
        reset_datasets_anc_map = self.reset_datasets_anc_map

        # iterate over time step(s)
        for alg_time_group in self.alg_time_datasets:

            # info start time group
            alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) + '" ... ')

            # define ancillary points dataset file name
            file_path_anc_pnt_dst_data_step = fill_path_with_tags(
                file_path_anc_pnt_dst_data_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)

            # define ancillary map file name
            file_path_anc_map_step = fill_path_with_tags(
                file_path_anc_map_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)

            # flag to reset ancillary datasets
            if reset_datasets_anc_map:
                if os.path.exists(file_path_anc_map_step):
                    os.remove(file_path_anc_map_step)

            # check file ancillary availability
            if not os.path.exists(file_path_anc_map_step):

                # info start collect datasets
                alg_logger.info(' -----> (1) Collect datasets ... ')

                # check file source availability
                if os.path.exists(file_path_anc_pnt_dst_data_step):
                    # get datasets obj
                    obj_cell_src = read_file_obj(file_path_anc_pnt_dst_data_step)
                    obj_data_cell_merged, obj_time_cell_merged = obj_cell_src['data'], obj_cell_src['time']
                    # remove nan(s) obj
                    obj_data_cell_merged = remove_nans(obj_data_cell_merged, keep_finite=True)
                else:
                    # info end cell group
                    alg_logger.error(' ===> File datasets "' + file_path_anc_pnt_dst_data_step + '" is not available')
                    raise FileExistsError('File must be available to correctly run the algorithm')

                # create mask reference
                mask_2d, mask_geox, mask_geoy = resample_points_to_grid(
                    geo_mask.ravel(),
                    var_geox_1d_in=geo_x_2d.ravel(),
                    var_geoy_1d_in=geo_y_2d.ravel(),
                    var_geox_1d_out=obj_data_cell_merged['longitude'].values,
                    var_geoy_1d_out=obj_data_cell_merged['latitude'].values)
                mask_2d = np.flipud(mask_2d)

                # convert grid to swath
                mask_1d = convert_grid_to_swath(
                    mask_2d, mask_geox, mask_geoy,
                    obj_data_cell_merged['longitude'].values, obj_data_cell_merged['latitude'].values,
                    search_rad=50000, neighbours=1)

                # info end collect metrics
                alg_logger.info(' -----> (1) Collect datasets ... DONE')


                ''' debug
                import numpy as np
                from repurpose.resample import resample_to_grid

                var_geo_x = obj_data_cell_merge['longitude'].values
                var_geo_y = obj_data_cell_merge['latitude'].values

                var_data_1 = obj_data_cell_merge['soil_moisture_ref'].values
                var_data_2 = obj_data_cell_merge['soil_moisture_k1'].values
                var_data_3 = obj_data_cell_merge['soil_moisture_k2'].values

                geo_x_1d = grid_geo_data['longitude'].values
                geo_y_1d = grid_geo_data['latitude'].values
                geo_x_2d, geo_y_2d = np.meshgrid(geo_x_1d, geo_y_1d)
                geo_mask = grid_geo_data.values

                max_distance = self.max_distance_ref_k1
                max_distance = 50000
                values_obj = resample_to_grid(
                    {'data': var_data_1},
                    var_geo_x, var_geo_y, geo_x_2d, geo_y_2d,
                    search_rad=max_distance, fill_values=np.nan,
                    min_neighbours=1, neighbours=8)
                var_data = values_obj['data']
                var_data[geo_mask == 0] = np.nan

                plot_data_2d(var_data, geo_x_2d, geo_y_2d)

                var_data_filled = deepcopy(var_data)
                '''

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
                    obj_data_map_resampled = fx_handle_common_resample(
                        var_data_obj=obj_data_cell_merged,
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
                var_data = obj_data_map_resampled['soil_moisture_resampled'].values
                var_data[geo_mask == 0] = np.nan
                plot_data_2d(var_data)
                '''

                # info start filter datasets
                alg_logger.info(' -----> (3) Filter datasets ... ')
                # get fx filter settings
                (fx_active_common_filter, fx_name_common_filter,
                 fx_vars_common_filter, fx_params_common_filter) = get_fx_settings(fx_settings_common_filter)
                # get fx filter handle
                fx_handle_common_filter = get_fx_method(fx_name_common_filter, fx_methods_common)

                # flag to activate filter part
                if fx_active_common_filter and fx_active_common_filter:
                    # call filter method
                    obj_data_map_filtered = fx_handle_common_filter(
                        obj_data_map_resampled,
                        var_mapping_in=fx_vars_common_filter['in'],
                        var_mapping_out=fx_vars_common_filter['out'],
                        filter_type=fx_params_common_filter['kernel_type'],
                        filter_pixels_std=fx_params_common_filter['kernel_pixels_std'],
                        filter_mode=fx_params_common_filter['kernel_mode'])
                    # info end filter datasets
                    alg_logger.info(' -----> (3) Filter datasets ... DONE')

                    # stored benchmark datasets
                    obj_data_map_filtered['soil_moisture_filtered_benchmark'] = deepcopy(
                        obj_data_map_filtered['soil_moisture_filtered'])
                    obj_data_map_filtered['soil_moisture_resampled_benchmark'] = deepcopy(
                        obj_data_map_filtered['soil_moisture_resampled'])

                else:
                    # save obj and stored benchmark datasets (only for compatibility with the next steps)
                    obj_data_map_filtered = deepcopy(obj_data_map_resampled)
                    obj_data_map_filtered['soil_moisture_filtered_benchmark'] = deepcopy(
                        obj_data_map_filtered['soil_moisture_filtered'])
                    obj_data_map_filtered['soil_moisture_resampled_benchmark'] = deepcopy(
                        obj_data_map_filtered['soil_moisture_resampled'])
                    # info end filter datasets
                    alg_logger.info(' -----> (3) Filter datasets ... SKIPPED. Flag filter is not activated')

                ''' debug
                var_data_f = obj_data_map_filtered['soil_moisture_filtered'].values
                var_data_f[geo_mask == 0] = np.nan
                plot_data_2d(var_data_f)
                var_data_r = obj_data_map_filtered['soil_moisture_resampled'].values
                var_data_r[geo_mask == 0] = np.nan
                plot_data_2d(var_data_r)
                '''

                # info start update datasets
                alg_logger.info(' -----> (4) Update datasets ... ')
                # get fx update settings
                (fx_active_common_update, fx_name_common_update,
                 fx_vars_common_update, fx_params_common_update) = get_fx_settings(fx_settings_common_update)
                # get fx update handle
                fx_handle_common_update = get_fx_method(fx_name_common_update, fx_methods_common)

                # flag to activate update part
                if fx_active_common_update and fx_active_common_update:
                    # call update method
                    obj_data_map_updated = fx_handle_common_update(
                        obj_data_map_filtered,
                        cell_obj=obj_data_cell_merged,
                        geo_obj=grid_geo_data,
                        mask_obj=grid_mask_data,
                        removed_obj=grid_rm_data,
                        var_mapping_in=fx_vars_common_update['in'],
                        var_mapping_out=fx_vars_common_update['out'])
                    # info end update datasets
                    alg_logger.info(' -----> (4) Update datasets ... DONE')
                else:
                    # info end update datasets
                    obj_data_map_updated = deepcopy(obj_data_map_filtered)
                    alg_logger.info(' -----> (4) Update datasets ... SKIPPED. Flag update is not activated')

                ''' debug
                var_data_f = obj_data_map_updated['soil_moisture_filtered'].values
                var_data_f[geo_mask == 0] = np.nan
                plot_data_2d(var_data_f)
                var_data_r = obj_data_map_updated['soil_moisture_resampled'].values
                var_data_r[geo_mask == 0] = np.nan
                plot_data_2d(var_data_r)
                var_data_rb = obj_data_map_updated['soil_moisture_resampled_benchmark'].values
                var_data_rb[geo_mask == 0] = np.nan
                plot_data_2d(var_data_rb)
                '''

                # info start mask datasets
                alg_logger.info(' -----> (5) Mask datasets ... ')
                # get fx mask settings
                (fx_active_common_mask, fx_name_common_mask,
                 fx_vars_common_mask, fx_params_common_mask) = get_fx_settings(fx_settings_common_mask)
                # get fx mask handle
                fx_handle_common_mask = get_fx_method(fx_name_common_mask, fx_methods_common)

                # flag to activate mask part
                if fx_active_common_mask and fx_active_common_mask:
                    # call mask method
                    obj_data_map_masked = fx_handle_common_mask(
                        obj_data_map_updated, grid_geo_data,
                        var_mapping_in=fx_vars_common_mask['in'],
                        var_mapping_out=fx_vars_common_mask['out'],
                        mask_value_min=fx_params_common_mask['min_value'],
                        mask_value_max=fx_params_common_mask['max_value'])
                    # info end mask datasets
                    alg_logger.info(' -----> (5) Mask datasets ... DONE')
                else:
                    # info end mask datasets
                    obj_data_map_masked = deepcopy(obj_data_map_updated)
                    alg_logger.info(' -----> (5) Mask datasets ... SKIPPED. Flag mask is not activated')

                ''' debug
                var_data_mf = obj_data_map_masked['soil_moisture_masked_filtered']
                plot_data_2d(var_data_mf)
                var_data_mr = obj_data_map_masked['soil_moisture_masked_resampled']
                plot_data_2d(var_data_mr)
                var_flags_m = obj_data_map_masked['flags_masked']
                plot_data_2d(var_flags_m)
                var_flags_r = obj_data_map_masked['flags_resampled']
                plot_data_2d(var_flags_r)
                #'''

                # check data
                if check_data(obj_data_map_masked, thr_data='any'):

                    # save data collection in pickle format
                    folder_name_anc, file_name_anc = os.path.split(file_path_anc_map_step)
                    make_folder(folder_name_anc)

                    obj_map_collection = {'data': obj_data_map_masked, 'time': obj_time_cell_merged}
                    write_file_obj(file_path_anc_map_step, obj_map_collection)

                    # info end time
                    alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                    '" ... DONE')
                else:
                    # info end time
                    alg_logger.info(' ----> Time Group "' + alg_time_group.strftime(time_format_algorithm) +
                                    '" ... FAILED. Datasets are defined by NoneType')

        # info end method
        alg_logger.info(' ---> Analyze dynamic maps ... DONE')

    # method to dump data maps
    def dump_data_maps(self):

        # info start method
        alg_logger.info(' ---> Dump dynamic maps ... ')

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

            # reformat time stamp of the destination file
            if self.time_stamp_dst == 'cday':
                alg_time_dump = deepcopy(alg_time_group)
                alg_logger.info(' -----> Time Dump is defined by "cday" tag "' +
                                alg_time_dump.strftime(time_format_algorithm) + '" ... DONE')
            elif self.time_stamp_dst == 'pday':
                alg_time_dump = alg_time_group - pd.Timedelta(days=1)
                alg_logger.info(' -----> Time Dump is defined by "pday" tag "' +
                                alg_time_dump.strftime(time_format_algorithm) + '" ... DONE')
            elif self.time_stamp_dst == 'fday':
                alg_time_dump = alg_time_group + pd.Timedelta(days=1)
                alg_logger.info(' -----> Time Dump is defined by "fday" tag "' +
                                alg_time_dump.strftime(time_format_algorithm) + '" ... DONE')
            else:
                alg_logger.error(' ===> Time stamp modifier "' + self.time_stamp_dst +
                                 '" is not expected by the algorithm (pday, cday, fday')
                raise NotImplemented('Case not implemented yet')

            # define ancillary file name
            file_path_anc_map_step = fill_path_with_tags(
                file_path_anc_map_tmpl, alg_time_group, tmpl_tags_time=alg_template_time)
            # define destination file name
            file_path_dst_step = fill_path_with_tags(
                file_path_dst_tmpl, alg_time_dump, tmpl_tags_time=alg_template_time)

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
        alg_logger.info(' ---> Dump dynamic maps ... DONE')

# -------------------------------------------------------------------------------------
