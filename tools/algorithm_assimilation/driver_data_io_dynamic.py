"""
Class Features

Name:          driver_data_io_dynamic_datasets_obs
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210225'
Version:       '1.0.0'
"""

######################################################################################
# Library
import logging
import os
import re

import pandas as pd
import numpy as np

from copy import deepcopy
from datetime import datetime
from osgeo import gdalconst
import scipy.ndimage as ndimage

from lib_data_io_nc import read_file_model
from lib_data_io_nc import read_file_datasets, split_file_datasets, wrap_fx_datasets, merge_fx_datasets

from lib_utils_system import fill_tags2string, make_folder
from lib_data_io_tiff import write_file_tiff, read_file_tiff

from lib_utils_time import set_period
from lib_data_analysis import wrap_fx_analysis, merge_fx_analysis
from lib_data_io_generic import read_obj, write_obj

from lib_utils_io import create_filename_tmp
from lib_utils_geo import create_map_idx

# Debug
import matplotlib.pylab as plt
from matplotlib import colors
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverDynamic
class DriverDynamic:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, time_run=None, time_range=None, time_chunks=None,
                 src_dict=None, anc_dict=None, anl_dict=None, dest_dict=None, tmp_dict=None,
                 parameters_dict=None,
                 tag_datasets_terrain='terrain', tag_datasets_map='map', tag_datasets_channel_network='channel_network',
                 tag_datasets_generic='datasets', tag_model_generic='model',
                 tag_datasets_sm='soil_moisture', tag_model_sm='soil_moisture',
                 tag_analysis_sm='soil_moisture',
                 static_data_collection=None,
                 flag_active_variables=None,
                 flag_cleaning_ancillary_model_collections=True,
                 flag_cleaning_ancillary_dset_obj=True,
                 flag_cleaning_ancillary_dset_ref=True,
                 flag_cleaning_ancillary_dset_grp=True,
                 flag_cleaning_ancillary_dset_collections=True,
                 flag_cleaning_analysis=True,
                 flag_cleaning_destination=True,
                 alg_template_tags=None):

        self.time_run = time_run
        self.time_range = time_range
        self.time_chunks = time_chunks

        self.src_dict = src_dict
        self.anc_dict = anc_dict
        self.anl_dict = anl_dict
        self.dest_dict = dest_dict
        self.tmp_dict = tmp_dict
        self.parameters_dict = parameters_dict

        self.alg_template_tags = alg_template_tags

        self.tag_datasets_terrain = tag_datasets_terrain
        self.tag_datasets_map = tag_datasets_map
        self.tag_datasets_channel_network = tag_datasets_channel_network
        self.tag_datasets_generic = tag_datasets_generic
        self.tag_datasets_sm = tag_datasets_sm
        self.tag_model_generic = tag_model_generic
        self.tag_model_sm = tag_model_sm
        self.tag_analysis_sm = tag_analysis_sm

        self.tag_file_name = 'file_name'
        self.tag_folder_name = 'folder_name'
        self.tag_file_compression = 'file_compression'
        self.tag_var_name = 'var_name'
        self.tag_var_scale_factor = 'var_scale_factor'
        self.tag_var_mask_cnet = 'var_mask_channel_network'
        self.tag_time_start = 'time_start'
        self.tag_time_end = 'time_end'
        self.tag_time_rounding = 'time_rounding'
        self.tag_time_frequency = 'time_frequency'

        self.tag_timeseries_obj = 'time_series_obj'
        self.tag_timeseries_ref = 'time_series_reference'
        self.tag_timeseries_grp = 'time_series_group'
        self.tag_timeseries_clt = 'time_series_collections'

        self.tag_anl_timeseries_grp = 'time_series_group'
        self.tag_anl_statistics_grp = 'statistics_group'
        self.tag_anl_timeseries_clt = 'time_series_collections'
        self.tag_anl_statistics_clt = 'statistics_collections'
        self.tag_anl_statistics_map = 'statistics_map'

        self.terrain_data = static_data_collection[self.tag_datasets_terrain]
        self.map_data = static_data_collection[self.tag_datasets_map]
        self.channel_network_data = static_data_collection[self.tag_datasets_channel_network]

        self.var_filter_list = [6, 12, 32]
        self.var_filter_format = 'swi_t{0:02d}'

        self.var_filter_name = [self.var_filter_format.format(var_name) for var_name in self.var_filter_list]

        folder_name_src_sm = self.src_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_folder_name]
        file_name_src_sm = self.src_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_file_name]
        var_name_src_sm = self.src_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_var_name]
        file_path_src_sm = os.path.join(folder_name_src_sm, file_name_src_sm)

        folder_name_src_model = self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_folder_name]
        file_name_src_model = self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_file_name]
        var_name_src_model = self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_var_name]
        file_path_src_model = os.path.join(folder_name_src_model, file_name_src_model)
        file_compression_flag_src_model = self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_file_compression]

        file_obj_src_model = self.define_file_string(
            file_path_src_model,
            time_start=self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_time_start],
            time_end=self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_time_end],
            time_frequency=self.src_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_time_frequency]
        )

        folder_name_anc_sm_obj = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_obj][self.tag_folder_name]
        file_name_anc_sm_obj = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_obj][self.tag_file_name]
        var_name_anc_sm_obj = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_obj][self.tag_var_name]
        file_path_anc_sm_obj = os.path.join(folder_name_anc_sm_obj, file_name_anc_sm_obj)

        folder_name_anc_sm_ref = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_ref][self.tag_folder_name]
        file_name_anc_sm_ref = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_ref][self.tag_file_name]
        var_name_anc_sm_ref = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_ref][self.tag_var_name]
        file_path_anc_sm_ref = os.path.join(folder_name_anc_sm_ref, file_name_anc_sm_ref)

        folder_name_anc_sm_grp = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_grp][self.tag_folder_name]
        file_name_anc_sm_grp = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_grp][self.tag_file_name]
        var_name_anc_sm_grp = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_grp][self.tag_var_name]
        file_path_anc_sm_grp = os.path.join(folder_name_anc_sm_grp, file_name_anc_sm_grp)

        folder_name_anc_sm_clt = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_clt][self.tag_folder_name]
        file_name_anc_sm_clt = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_clt][self.tag_file_name]
        var_name_anc_sm_clt = self.anc_dict[self.tag_datasets_generic][self.tag_datasets_sm][self.tag_timeseries_clt][self.tag_var_name]
        file_path_anc_sm_clt = os.path.join(folder_name_anc_sm_clt, file_name_anc_sm_clt)

        self.var_datasets_name = [var_name_anc_sm_clt] + self.var_filter_name

        file_obj_anc_sm_clt = self.define_file_string(
            file_path_anc_sm_clt,
            time_start=None,
            time_end=None,
            extra_args={'variable_name_list': self.var_datasets_name}
        )

        folder_name_anc_model_clt = self.anc_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_timeseries_clt][self.tag_folder_name]
        file_name_anc_model_clt = self.anc_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_timeseries_clt][self.tag_file_name]
        var_name_anc_model_clt = self.anc_dict[self.tag_model_generic][self.tag_datasets_sm][self.tag_timeseries_clt][self.tag_var_name]
        file_path_anc_model_clt = os.path.join(folder_name_anc_model_clt, file_name_anc_model_clt)

        folder_name_anl_data_grp = self.anl_dict[self.tag_analysis_sm][self.tag_anl_timeseries_grp][self.tag_folder_name]
        file_name_anl_data_grp = self.anl_dict[self.tag_analysis_sm][self.tag_anl_timeseries_grp][self.tag_file_name]
        file_path_anl_data_grp = os.path.join(folder_name_anl_data_grp, file_name_anl_data_grp)

        file_obj_anl_data_grp = self.define_file_string(
            file_path_anl_data_grp,
            time_start=None,
            time_end=None,
            extra_args={'variable_name_list': self.var_datasets_name}
        )

        folder_name_anl_stats_grp = self.anl_dict[self.tag_analysis_sm][self.tag_anl_statistics_grp][self.tag_folder_name]
        file_name_anl_stats_grp = self.anl_dict[self.tag_analysis_sm][self.tag_anl_statistics_grp][self.tag_file_name]
        file_path_anl_stats_grp = os.path.join(folder_name_anl_stats_grp, file_name_anl_stats_grp)

        file_obj_anl_stats_grp = self.define_file_string(
            file_path_anl_stats_grp,
            time_start=None,
            time_end=None,
            extra_args={'variable_name_list': self.var_datasets_name}
        )

        folder_name_anl_data_clt = self.anl_dict[self.tag_analysis_sm][self.tag_anl_timeseries_clt][self.tag_folder_name]
        file_name_anl_data_clt = self.anl_dict[self.tag_analysis_sm][self.tag_anl_timeseries_clt][self.tag_file_name]
        file_path_anl_data_clt = os.path.join(folder_name_anl_data_clt, file_name_anl_data_clt)

        file_obj_anl_data_clt = self.define_file_string(
            file_path_anl_data_clt,
            time_start=None,
            time_end=None,
            extra_args={'variable_name_list': self.var_datasets_name}
        )

        folder_name_anl_stats_clt = self.anl_dict[self.tag_analysis_sm][self.tag_anl_statistics_clt][self.tag_folder_name]
        file_name_anl_stats_clt = self.anl_dict[self.tag_analysis_sm][self.tag_anl_statistics_clt][self.tag_file_name]
        file_path_anl_stats_clt = os.path.join(folder_name_anl_stats_clt, file_name_anl_stats_clt)

        file_obj_anl_stats_clt = self.define_file_string(
            file_path_anl_stats_clt,
            time_start=None,
            time_end=None,
            extra_args={'variable_name_list': self.var_datasets_name}
        )

        folder_name_anl_stats_map = self.anl_dict[self.tag_analysis_sm][self.tag_anl_statistics_map][self.tag_folder_name]
        file_name_anl_stats_map = self.anl_dict[self.tag_analysis_sm][self.tag_anl_statistics_map][self.tag_file_name]
        file_path_anl_stats_map = os.path.join(folder_name_anl_stats_map, file_name_anl_stats_map)

        file_obj_anl_stats_map = self.define_file_string(
            file_path_anl_stats_map,
            time_start=None,
            time_end=None,
            extra_args={'variable_name_list': self.var_datasets_name}
        )

        folder_name_dest_map = self.dest_dict[self.tag_datasets_sm][self.tag_folder_name]
        file_name_dest_map = self.dest_dict[self.tag_datasets_sm][self.tag_file_name]
        var_name_dest_map = self.dest_dict[self.tag_datasets_sm][self.tag_var_name]
        var_scale_factor_dest_map = self.dest_dict[self.tag_datasets_sm][self.tag_var_scale_factor]
        var_mask_cnet_dest_map = self.dest_dict[self.tag_datasets_sm][self.tag_var_mask_cnet]
        file_path_dest_map = os.path.join(folder_name_dest_map, file_name_dest_map)
        file_compression_flag_dest_map = self.dest_dict[self.tag_datasets_sm][self.tag_file_compression]

        file_obj_dest_map = self.define_file_string(
            file_path_dest_map,
            time_start=self.dest_dict[self.tag_datasets_sm][self.tag_time_start],
            time_end=self.dest_dict[self.tag_datasets_sm][self.tag_time_end],
            time_frequency=self.dest_dict[self.tag_datasets_sm][self.tag_time_frequency]
        )

        self.time_period_dest = set_period(
            time_start=self.dest_dict[self.tag_datasets_sm][self.tag_time_start],
            time_end=self.dest_dict[self.tag_datasets_sm][self.tag_time_end],
            time_freq=self.dest_dict[self.tag_datasets_sm][self.tag_time_frequency])

        self.file_path_src_dst = {self.tag_datasets_sm: file_path_src_sm}
        self.var_name_src_dst = {self.tag_datasets_sm: var_name_src_sm}

        self.file_path_anc_dst_obj = {self.tag_datasets_sm: file_path_anc_sm_obj}
        self.var_name_anc_dst_obj = {self.tag_datasets_sm: var_name_anc_sm_obj}
        self.file_path_anc_dst_ref = {self.tag_datasets_sm: file_path_anc_sm_ref}
        self.var_name_anc_dst_ref = {self.tag_datasets_sm: var_name_anc_sm_ref}
        self.file_path_anc_dst_grp = {self.tag_datasets_sm: file_path_anc_sm_grp}
        self.var_name_anc_dst_grp = {self.tag_datasets_sm: var_name_anc_sm_grp}
        self.file_path_anc_dst_clt = {self.tag_datasets_sm: file_obj_anc_sm_clt}
        self.var_name_anc_dst_clt = {self.tag_datasets_sm: self.var_datasets_name}

        self.file_path_anl_data_grp = {self.tag_analysis_sm: file_obj_anl_data_grp}
        self.file_path_anl_stats_grp = {self.tag_analysis_sm: file_obj_anl_stats_grp}
        self.file_path_anl_data_clt = {self.tag_analysis_sm: file_obj_anl_data_clt}
        self.file_path_anl_stats_clt = {self.tag_analysis_sm: file_obj_anl_stats_clt}
        self.file_path_anl_stats_map = {self.tag_analysis_sm: file_obj_anl_stats_map}

        self.file_path_dest_map = {self.tag_datasets_sm: file_obj_dest_map}
        self.var_name_dset_map = {self.tag_datasets_sm: var_name_dest_map}
        self.var_scale_factor_dest_map = {self.tag_datasets_sm: var_scale_factor_dest_map}
        self.var_mask_cnet_dest_map = {self.tag_datasets_sm: var_mask_cnet_dest_map}

        self.file_path_src_model = {self.tag_model_sm: file_obj_src_model}
        self.var_name_src_model = {self.tag_model_sm: var_name_src_model}
        self.file_compression_flag_src_model = file_compression_flag_src_model
        self.file_compression_ext_src_model = '.gz'

        self.file_path_anc_model_clt = {self.tag_model_sm: file_path_anc_model_clt}
        self.var_name_anc_model_clt = {self.tag_model_sm: var_name_anc_model_clt}

        self.file_compression_flag_dest_map = file_compression_flag_dest_map
        self.file_compression_ext_dest_map = '.gz'

        self.var_geo_x = 'longitude'
        self.var_geo_y = 'latitude'

        self.flag_active_variables = flag_active_variables

        self.flag_cleaning_anc_model_collections = flag_cleaning_ancillary_model_collections
        self.flag_cleaning_anc_dset_obj = flag_cleaning_ancillary_dset_obj
        self.flag_cleaning_anc_dset_ref = flag_cleaning_ancillary_dset_ref
        self.flag_cleaning_anc_dset_grp = flag_cleaning_ancillary_dset_grp
        self.flag_cleaning_anc_dset_collections = flag_cleaning_ancillary_dset_collections
        self.flag_cleaning_analysis = flag_cleaning_analysis
        self.flag_cleaning_dest = flag_cleaning_destination

        self.terrain_attributes = {'ncols': self.terrain_data['values'].shape[1],
                                   'nrows': self.terrain_data['values'].shape[0],
                                   'nodata_value': -9999.0,
                                   'xllcorner': self.terrain_data['bb_left'],
                                   'yllcorner': self.terrain_data['bb_bottom'],
                                   'cellsize': np.mean([self.terrain_data['res_lon'], self.terrain_data['res_lat']])}

        self.terrain_idxs = create_map_idx(self.terrain_data['values'])

        self.var_chunks = self.parameters_dict['chunks']
        self.var_cpu_datasets = self.parameters_dict['cpu_datasets']
        self.var_cpu_analysis = self.parameters_dict['cpu_analysis']

        self.format_chunks_element = '{0:06d}'
        self.format_chunks_group = '{0:06d}'

        self.folder_tmp = self.tmp_dict[self.tag_folder_name]

        self.pearson_thr = 0.5
        self.pearson_nodata = np.nan
        self.isolated_pixels_thr = 50

        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define  filename
    def define_file_string(self, file_path_raw, time_start=None, time_end=None, time_frequency='H',
                           extra_args=None):

        if (time_start is not None) and (time_end is not None):
            time_period = pd.date_range(start=time_start, end=time_end, freq=time_frequency)
        else:
            time_period = ['NA']

        alg_template_tags = self.alg_template_tags

        var_name_list = None
        if extra_args is not None:
            if 'variable_name_list' in extra_args:
                var_name_list = extra_args['variable_name_list']

        file_path_list = []
        for time_step in time_period:

            alg_template_values_raw = {
                'source_dset_datetime': time_step,
                'source_dset_sub_path_time': time_step,
                'source_model_datetime': time_step,
                'source_model_sub_path_time': time_step,
                'ancillary_datetime': time_step,
                'ancillary_sub_path_time': time_step,
                'destination_datetime': time_step,
                'destination_sub_path_time': time_step,
            }

            if not var_name_list:
                file_path_obj = fill_tags2string(file_path_raw, alg_template_tags, alg_template_values_raw)
                file_path_list.append(file_path_obj)
            else:
                for var_name_step in var_name_list:
                    alg_template_values_extra = {'dataset_var_name':var_name_step,
                                                 'model_var_name': var_name_step}
                    alg_template_values_step = {**alg_template_values_raw, **alg_template_values_extra}

                    file_path_obj = fill_tags2string(file_path_raw, alg_template_tags, alg_template_values_step)
                    file_path_list.append(file_path_obj)

        return file_path_list

    # -------------------------------------------------------------------------------------
    '''
    # -------------------------------------------------------------------------------------
    # Method to check time step
    @staticmethod
    def validate_time(var_time_str, file_name):

        file_time = re.search(r'\d{4}\d{2}\d{2}\d{2}\d{2}', file_name)
        file_time_stamp = pd.Timestamp(datetime.strptime(file_time.group(), '%Y%m%d%H%M'))
        var_time_stamp = pd.Timestamp(var_time_str)

        if var_time_stamp == file_time_stamp:
            time_check = True
        else:
            time_check = False

        return time_check

    # -------------------------------------------------------------------------------------
    '''
    # -------------------------------------------------------------------------------------
    # Method to dump dynamic data
    def dump_dynamic_data(self, var_name='soil_moisture', var_scale='datasets_scaled_lr',
                          var_limit_min=0, var_limit_max=1, var_nodata=np.nan):

        logging.info(' ---> Dump dynamic data ... ')

        time_period = self.time_period_dest

        terrain_values = self.terrain_data['values']
        terrain_transform = self.terrain_data['transform']
        cnet_values = self.channel_network_data['values']

        terrain_attrs = self.terrain_attributes

        file_obj_data_clt = self.file_path_anl_data_clt
        file_obj_stats_map = self.file_path_anl_stats_map
        file_obj_data_map = self.file_path_dest_map

        var_dump_idx = self.var_datasets_name.index(var_name)

        flag_cleaning_destination = self.flag_cleaning_dest

        for var_name_active in self.flag_active_variables:

            logging.info(' ----> Variable: ' + var_name_active + ' ... ')

            file_path_data_clt = file_obj_data_clt[var_name_active]
            file_path_stats_map = file_obj_stats_map[var_name_active]
            file_path_data_map = file_obj_data_map[var_name_active]

            file_path_data = file_path_data_clt[var_dump_idx]
            file_path_stats = file_path_stats_map[var_dump_idx]

            if os.path.exists(file_path_stats):
                var_stats_2d = read_file_tiff(file_path_stats)[0]
                var_stats_2d = var_stats_2d[:, :, 0]
            else:
                var_stats_2d = None

            if os.path.exists(file_path_data):

                file_obj_data = read_obj(file_path_data)
                file_obj_attrs = file_obj_data.attrs

                var_model_id = np.asarray(file_obj_attrs['model_id'], dtype=int)

                if var_name_active in list(self.var_scale_factor_dest_map.keys()):
                    var_scale_factor = self.var_scale_factor_dest_map[var_name_active]
                else:
                    var_scale_factor = 1

                if var_name_active in list(self.var_mask_cnet_dest_map.keys()):
                    var_mask_cnet = self.var_mask_cnet_dest_map[var_name_active]
                else:
                    var_mask_cnet = False

                for time_step, file_path_step in zip(time_period, file_path_data_map):

                    logging.info(' -----> Time: ' + time_step.strftime('%Y-%m-%d %H:%M') + ' ... ')

                    if flag_cleaning_destination:
                        if os.path.exists(file_path_step):
                            os.remove(file_path_step)

                    if not os.path.exists(file_path_step):
                        var_data_2d_init = np.zeros([terrain_values.shape[0], terrain_values.shape[1]])
                        var_data_2d_init[:, :] = np.nan

                        logging.info(' ------> Select data ... ')
                        var_data_1d_in = file_obj_data.sel(time=time_step)[var_scale].values
                        var_data_1d_in = np.asarray(var_data_1d_in, dtype=float).flatten()

                        # var_data_1d_tmp1 = deepcopy(var_data_2d_init.flatten())
                        # var_data_1d_tmp1[var_model_id] = var_data_1d_in
                        # var_data_2d_tmp1 = np.reshape(var_data_1d_tmp1, [terrain_values.shape[0], terrain_values.shape[1]])

                        if var_limit_min is not None:
                            var_data_1d_in[var_data_1d_in < var_limit_min] = var_nodata
                        if var_limit_max is not None:
                            var_data_1d_in[var_data_1d_in > var_limit_max] = var_nodata

                        var_data_1d_out = deepcopy(var_data_2d_init.flatten())
                        var_data_1d_out[var_model_id] = var_data_1d_in
                        var_data_2d_out = np.reshape(var_data_1d_out, [terrain_values.shape[0], terrain_values.shape[1]])

                        # apply variable scale factor
                        var_data_2d_out = var_data_2d_out / var_scale_factor
                        # filter over domain
                        var_data_2d_out[terrain_values < 0] = var_nodata

                        if var_stats_2d is not None:
                            var_data_2d_out[np.isnan(var_stats_2d)] = var_nodata
                        logging.info(' ------> Select data ... DONE')

                        logging.info(' ------> Mask data ... ')
                        if var_mask_cnet:
                            var_data_2d_out[cnet_values == 1] = var_nodata
                            logging.info(' ------> Mask data ... DONE')
                        else:
                            logging.info(' ------> Mask data ... SKIPPED. Channel network mask not activated')

                        logging.info(' ------> Save data ... ')
                        if not np.all(np.isnan(var_data_2d_out)):

                            folder_name, file_name = os.path.split(file_path_step)
                            make_folder(folder_name)

                            write_file_tiff(file_path_step, var_data_2d_out, None,
                                            terrain_values.shape[1], terrain_values.shape[0], terrain_transform,
                                            file_epsg=4326,
                                            file_metadata={
                                                'field_description': var_name,
                                                'field_analysis': 'pearson_filtered with thr ' + str(self.pearson_thr)},
                                            file_scale_factor=1,
                                            file_format=gdalconst.GDT_Float32)

                            logging.info(' ------> Save data ... DONE')
                        else:
                            logging.info(' ------> Save data ... SKIPPED. All values are NANs')

                        logging.info(' -----> Time: ' + time_step.strftime('%Y-%m-%d %H:%M') +
                                     ' ... DONE')

                    else:

                        logging.info(' -----> Time: ' + time_step.strftime('%Y-%m-%d %H:%M') +
                                     ' ... SKIPPED. Datasets already done')

                logging.info(' ----> Variable: ' + var_name_active + ' ... DONE')

            else:

                logging.info(' ----> Variable: ' + var_name_active + ' ... FAILED. Datasets is not available')

        logging.info(' ---> Dump dynamic data ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to view dynamic statistics
    @staticmethod
    def view_dynamic_stats(map, file_name):

        plt.figure()
        cmap = colors.ListedColormap(['gray', 'red', 'orange', 'yellow', 'blue', 'green'])
        boundaries = [-1, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

        plt.imshow(map, cmap=cmap, norm=norm)
        plt.clim(-9999, 1)
        plt.colorbar()

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title('Pearson - Variable "' + file_analysis_variable + '" scaled Vs HMC-SM')

        # plt.show()
        fig_dpi = 120
        plt.savefig('pearson_' + file_analysis_variable + '.png', dpi=fig_dpi)

        print('ciao')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to map dynamic data
    def filter_dynamic_stats(self):

        logging.info(' ---> Filter dynamic statistics ... ')

        terrain_transform = self.terrain_data['transform']
        terrain_idxs = self.terrain_idxs

        flag_cleaning_analysis = self.flag_cleaning_analysis

        for var_name in self.flag_active_variables:

            logging.info(' ----> Variable: ' + var_name + ' ... ')

            file_analysis_stats_clt = self.file_path_anl_stats_clt[var_name]
            file_analysis_stats_map = self.file_path_anl_stats_map[var_name]

            for file_analysis_variable, file_stats_clt_step, file_stats_map_step in zip(
                    self.var_datasets_name, file_analysis_stats_clt, file_analysis_stats_map):

                logging.info(' -----> Filter: ' + file_analysis_variable + ' ... ')

                if os.path.exists(file_stats_clt_step):

                    if flag_cleaning_analysis:
                        if os.path.exists(file_stats_map_step):
                            os.remove(file_stats_map_step)

                    if not os.path.exists(file_stats_map_step):

                        logging.info(' ------> Read statistics data ... ')
                        obj_analysis_stats = read_obj(file_stats_clt_step)
                        logging.info(' ------> Read statistics data ... DONE')

                        logging.info(' ------> Filter statistics data ... ')
                        file_data_score_2d = np.zeros([terrain_idxs.shape[0], terrain_idxs.shape[1]])
                        file_data_score_1d = file_data_score_2d.flatten()

                        pearson_values_1d = np.asarray(obj_analysis_stats['datasets_pearson_lr'], dtype=float)
                        model_id = np.asarray(obj_analysis_stats['model_id'], dtype=int)

                        file_data_score_1d[model_id] = pearson_values_1d
                        file_data_score_2d = np.reshape(file_data_score_1d, [terrain_idxs.shape[0], terrain_idxs.shape[1]])

                        file_mask_idxs_1d = np.argwhere(file_data_score_2d.flatten() >= self.pearson_thr)

                        file_mask_score_1d = np.zeros([file_data_score_2d.shape[0], file_data_score_2d.shape[1]]).flatten()
                        file_mask_score_1d[file_mask_idxs_1d] = 1
                        file_mask_score_1d = file_mask_score_1d.astype(int)
                        file_mask_score_2d = np.reshape(file_mask_score_1d,
                                                        [file_data_score_2d.shape[0], file_data_score_2d.shape[1]])

                        struct_label = ndimage.generate_binary_structure(2, 2)
                        file_mask_filtered_2d = np.copy(file_mask_score_2d)
                        id_regions, num_ids = ndimage.label(file_mask_filtered_2d, structure=struct_label)
                        size_ids = np.array(ndimage.sum(file_mask_score_2d, id_regions, range(num_ids + 1)))

                        file_mask_idx_1d = np.zeros([file_data_score_2d.shape[0], file_data_score_2d.shape[1]]).flatten()
                        file_mask_filtered_1d = np.copy(file_mask_score_2d).flatten()
                        for num_id, size_id in zip(range(num_ids + 1), size_ids):
                            if size_id <= self.isolated_pixels_thr:
                                idx_pixel = np.argwhere(id_regions.flatten() == num_id)
                                file_mask_filtered_1d[idx_pixel] = 0
                                file_mask_idx_1d[idx_pixel] = 1

                        file_mask_filtered_2d = np.reshape(
                            file_mask_filtered_1d, [file_data_score_2d.shape[0], file_data_score_2d.shape[1]])
                        file_mask_idx_2d = np.reshape(
                            file_mask_idx_1d, [file_data_score_2d.shape[0], file_data_score_2d.shape[1]])

                        file_data_score_2d[file_mask_filtered_2d == 0] = self.pearson_nodata
                        file_data_score_2d[file_data_score_2d < self.pearson_thr] = self.pearson_nodata
                        logging.info(' ------> Filter statistics data ... DONE')

                        logging.info(' ------> Save statistics data ... ')
                        folder_name, file_name = os.path.split(file_stats_map_step)
                        make_folder(folder_name)

                        write_file_tiff(file_stats_map_step, file_data_score_2d, None,
                                        file_data_score_2d.shape[1], file_data_score_2d.shape[0], terrain_transform,
                                        file_epsg=4326,
                                        file_metadata={
                                            'field_description': file_analysis_variable,
                                            'field_analysis': 'pearson_filtered with thr ' + str(self.pearson_thr)},
                                        file_scale_factor=1,
                                        file_format=gdalconst.GDT_Float32)

                        logging.info(' ------> Save statistics data ... DONE')

                        logging.info(' -----> Filter: ' + file_analysis_variable + ' ... DONE')

                    else:
                        logging.info(' -----> Filter: ' + file_analysis_variable +
                                     ' ... SKIPPED. Filter already applied')
                else:
                    logging.info(' -----> Filter: ' + file_analysis_variable +
                                 ' ... FAILED. Statistics datasets was not available')

            logging.info(' ---> Filter dynamic statistics ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to analyze dynamic data
    def analyze_dynamic_data(self):

        logging.info(' ---> Analyze dynamic data ... ')

        var_chunks = self.var_chunks
        var_cpu = self.var_cpu_analysis

        terrain_data = self.terrain_data
        terrain_idxs = self.terrain_idxs
        map_data = self.map_data

        flag_cleaning_analysis = self.flag_cleaning_analysis

        for var_name in self.flag_active_variables:

            logging.info(' ----> Variable: ' + var_name + ' ... ')

            file_ancillary_dst_clt = self.file_path_anc_dst_clt[var_name]
            file_ancillary_model_clt = self.file_path_anc_model_clt[var_name]

            file_analysis_data_grp = self.file_path_anl_data_grp[var_name]
            file_analysis_stats_grp = self.file_path_anl_stats_grp[var_name]
            file_analysis_data_clt = self.file_path_anl_data_clt[var_name]
            file_analysis_stats_clt = self.file_path_anl_stats_clt[var_name]

            var_name_dst = self.var_name_anc_dst_clt[var_name]
            var_name_model = self.var_name_anc_model_clt[var_name]

            if not all([os.path.exists(f) for f in file_analysis_data_clt]):

                # Method to wrap fx analysis method
                wrap_fx_analysis(file_ancillary_model_clt, file_ancillary_dst_clt,
                                 file_analysis_data_grp, file_analysis_stats_grp,
                                 terrain_data=terrain_data, map_data=map_data,
                                 var_model_name=var_name_model, var_dst_name=var_name_dst,
                                 var_chunks=var_chunks,
                                 format_chunks_element=self.format_chunks_element,
                                 format_chunk_group=self.format_chunks_group,
                                 cpu_n=var_cpu,
                                 flag_cleaning_file=flag_cleaning_analysis)

                # Method to merge fx analysis outcome
                merge_fx_analysis(file_analysis_data_grp, file_analysis_stats_grp,
                                  file_analysis_data_clt, file_analysis_stats_clt,
                                  var_datasets_name='datasets_scaled_lr',
                                  terrain_data=terrain_data, map_data=map_data,
                                  var_model_name=var_name_model, var_dst_name=var_name_dst,
                                  var_chunks=var_chunks,
                                  format_chunks_element=self.format_chunks_element,
                                  format_chunk_group=self.format_chunks_group)

            logging.info(' ----> Variable: ' + var_name + ' ... DONE')

        logging.info(' ---> Analyze dynamic data ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize dynamic data
    def organize_dynamic_data_model(self):

        logging.info(' ---> Organize dynamic data MODEL ... ')

        folder_tmp = self.folder_tmp
        terrain_idxs = self.terrain_idxs

        flag_cleaning_anc_clt = self.flag_cleaning_anc_model_collections

        for var_name in self.flag_active_variables:

            logging.info(' ----> Variable: ' + var_name + ' ... ')

            file_source_list = self.file_path_src_model[var_name]
            file_ancillary_path_clt = self.file_path_anc_model_clt[var_name]

            file_var_name_in = self.var_name_src_model[var_name]
            file_var_name_out = var_name

            if not isinstance(file_source_list, list):
                file_source_list = [file_source_list]

            if os.path.exists(file_source_list[0]):
                read_file_model(
                    file_source_list, folder_tmp=folder_tmp,
                    file_idxs=terrain_idxs,
                    var_file_ancillary=file_ancillary_path_clt, var_upd_ancillary=flag_cleaning_anc_clt,
                    var_name_in=file_var_name_in, var_name_out=file_var_name_out)

                logging.info(' ----> Variable: ' + var_name + ' ... DONE')

            else:
                logging.warning(' ----> Variable: ' + var_name + ' ... FAILED.')
                raise IOError('File ' + file_source_list[0] + ' not found!')

        logging.info(' ---> Organize dynamic data MODEL ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize dynamic data
    def organize_dynamic_data_obs(self):

        logging.info(' ---> Organize dynamic data OBS ... ')

        var_chunks = self.var_chunks
        var_cpu = self.var_cpu_datasets

        flag_cleaning_anc_obj = self.flag_cleaning_anc_dset_obj
        flag_cleaning_anc_ref = self.flag_cleaning_anc_dset_ref
        flag_cleaning_anc_grp = self.flag_cleaning_anc_dset_grp
        flag_cleaning_anc_collections = self.flag_cleaning_anc_dset_collections

        for var_name in self.flag_active_variables:

            logging.info(' ----> Variable: ' + var_name + ' ... ')

            file_source_path = self.file_path_src_dst[var_name]
            file_ancillary_path_obj = self.file_path_anc_dst_obj[var_name]
            file_ancillary_path_ref = self.file_path_anc_dst_ref[var_name]
            file_ancillary_path_grp = self.file_path_anc_dst_grp[var_name]
            file_ancillary_obj_clt = self.file_path_anc_dst_clt[var_name]

            file_var_name_in = self.var_name_src_dst[var_name]
            file_var_name_out = var_name

            var_name_clt_list = self.var_name_anc_dst_clt[var_name]

            if flag_cleaning_anc_collections:
                for file_ancillary_path_clt in file_ancillary_obj_clt:
                    if os.path.exists(file_ancillary_path_clt):
                        os.remove(file_ancillary_path_clt)

            logging.info(' -----> Get and organize datasets ... ')

            if not all([os.path.exists(f) for f in file_ancillary_obj_clt]):

                if os.path.exists(file_source_path):
                    file_obj = read_file_datasets(
                        file_source_path,
                        file_obj=file_ancillary_path_obj, file_upd=flag_cleaning_anc_obj,
                        var_name_in=file_var_name_in, var_name_out=file_var_name_out)

                    file_chunks = split_file_datasets(
                        file_obj, file_ref=file_ancillary_path_ref, file_upd=flag_cleaning_anc_ref,
                        var_chunks=var_chunks, format_chunks_element=self.format_chunks_element)

                    wrap_fx_datasets(file_obj[file_var_name_out],  file_obj['time'],
                                     var_file_ancillary=file_ancillary_path_grp, var_upd_ancillary=flag_cleaning_anc_grp,
                                     var_filter_list=self.var_filter_list, format_filter_list=self.var_filter_format,
                                     file_chunks=file_chunks,
                                     format_chunks_element=self.format_chunks_element,
                                     format_chunk_group=self.format_chunks_group,
                                     cpu_n=var_cpu)

                    merge_fx_datasets(
                        var_name_clt_list,
                        file_ancillary_path_grp, file_ancillary_obj_clt,
                        file_chunks=file_chunks,
                        format_chunks_element=self.format_chunks_element,
                        format_chunk_group=self.format_chunks_group)

                else:
                    logging.info(' -----> Get and organize datasets ... FAILED. Source data is not available')
                    raise IOError('File ' + file_source_path + ' not found!')

                logging.info(' -----> Get and organize datasets ... DONE')

            else:
                logging.info(' -----> Get and organize datasets ... SKIPPED. Datasets are previously organized')

            logging.info(' ----> Variable: ' + var_name + ' ... DONE')

        logging.info(' ---> Organize dynamic data OBS ... DONE')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
