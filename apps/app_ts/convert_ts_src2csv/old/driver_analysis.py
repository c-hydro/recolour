"""
Class Features

Name:          driver_analysis
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
import pytesmo.scaling as scaling

from copy import deepcopy

from lib_data_io_pickle import read_obj, write_obj

from lib_utils_system import fill_tags2string, make_folder

from lib_info_args import logger_name, zip_extension

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
import matplotlib.pylab as plt
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverAnalysis
class DriverAnalysis:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, time_step, time_reference,
                 anc_dict, dst_dict, alg_dict=None,
                 geo_dict=None, time_dict=None, tmp_dict=None,
                 template_tags_dict=None,
                 flag_data_anc_point='point', flag_data_anc_grid='grid',
                 flag_data_dst='analysis',
                 flag_data_updating=True, flag_data_season=True):

        self.time_step = pd.Timestamp(time_step)
        self.time_reference = pd.Timestamp(time_reference)

        self.anc_dict = anc_dict
        self.dst_dict = dst_dict
        self.alg_dict = alg_dict
        self.geo_dict = geo_dict
        self.tmp_dict = tmp_dict

        self.flag_data_anc_point = flag_data_anc_point
        self.flag_data_anc_grid = flag_data_anc_grid
        self.flag_data_dst = flag_data_dst

        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'

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
        self.alg_point_geo_spatial_operation = alg_dict['geo_spatial_operation']
        self.alg_point_geo_spatial_mask = alg_dict['geo_spatial_mask']

        self.alg_datasets_reference = alg_dict['datasets_reference']
        self.alg_datasets_other = alg_dict['datasets_other']

        # ancillary object(s)
        self.file_name_anc_point_raw = anc_dict[self.flag_data_anc_point][self.file_name_tag]
        self.folder_name_anc_point_raw = anc_dict[self.flag_data_anc_point][self.folder_name_tag]

        self.file_path_anc_point = self.collect_file_list(
            self.folder_name_anc_point_raw, self.file_name_anc_point_raw,
            file_time_range=pd.DatetimeIndex([self.time_reference]))

        self.file_name_anc_grid_raw = anc_dict[self.flag_data_anc_grid][self.file_name_tag]
        self.folder_name_anc_grid_raw = anc_dict[self.flag_data_anc_grid][self.folder_name_tag]

        self.file_path_anc_grid = self.collect_file_list(
            self.folder_name_anc_grid_raw, self.file_name_anc_grid_raw,
            file_time_range=pd.DatetimeIndex([self.time_reference]))

        # destination object(s)
        self.folder_name_dst_raw = self.dst_dict[self.flag_data_dst][self.folder_name_tag]
        self.file_name_dst_raw = self.dst_dict[self.flag_data_dst][self.file_name_tag]

        self.file_path_dst = self.collect_file_list(
            self.folder_name_dst_raw, self.file_name_dst_raw, file_time_range=pd.DatetimeIndex([self.time_reference]))

        # tmp object(s)
        self.folder_name_tmp_raw = tmp_dict[self.folder_name_tag]
        self.file_name_tmp_raw = tmp_dict[self.file_name_tag]

        self.file_extension_zip = zip_extension
        self.file_extension_unzip = 'bin'

        self.flag_data_updating = flag_data_updating
        self.flag_data_season = flag_data_season

        # season filter
        if self.flag_data_season:
            self.lut_data_season = {
                1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA',
                7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON', 12: 'DJF'}
        else:
            self.lut_data_season = {
                1: 'ALL', 2: 'ALL', 3: 'ALL', 4: 'ALL', 5: 'ALL', 6: 'ALL',
                7: 'ALL', 8: 'ALL', 9: 'ALL', 10: 'ALL', 11: 'ALL', 12: 'ALL'}

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
    # Method to fill the time series
    @staticmethod
    def fill_ts(ts_reference, ts_other, fill_limit=2, fill_direction='forward',
                fill_method='polynomial', fill_order=2):

        if fill_method == 'polynomial':
            ts_reference_filled = ts_reference.interpolate(
                method='polynomial', order=fill_order, limit=fill_limit, limit_direction=fill_direction)

            ts_other_filled = ts_other.interpolate(
                method='polynomial', order=fill_order, limit=fill_limit, limit_direction=fill_direction)

        else:
            log_stream.error(' ===> Fill method "' + fill_method + '" is not supported')
            raise NotImplementedError('Case not implemented yet. Only "polynomial" method is available')

        return ts_reference_filled, ts_other_filled
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to synchronize the time series
    @staticmethod
    def sync_ts(ts_reference, ts_other_in, temporal_window='30min', temporal_operation='average'):

        name_reference = ts_reference.name
        time_reference = ts_reference.index
        time_other = ts_other_in.index

        ts_other_out = None
        for time_stamp_center in time_reference:

            time_stamp_left = time_stamp_center - pd.Timedelta(temporal_window)
            time_stamp_right = time_stamp_center + pd.Timedelta(temporal_window)

            if time_stamp_left < time_other[0] - pd.Timedelta(temporal_window):
                time_stamp_left = time_other[0]
            if time_stamp_right > time_other[-1] + pd.Timedelta(temporal_window):
                time_stamp_right = time_other[-1]

            time_stamp_period = pd.date_range(start=time_stamp_left, end=time_stamp_right, freq=temporal_window)

            ts_other_subset = ts_other_in.loc[time_stamp_period]

            if temporal_operation == 'average':
                value_other_cmp = ts_other_subset.mean()
            else:
                log_stream.error(' ===> Temporal operation "' + temporal_operation + '" is not supported')
                raise NotImplementedError('Case not implemented yet. Only "average" method is available')

            if ts_other_out is None:
                ts_other_out = pd.Series(data=value_other_cmp, index=[time_stamp_center], name=name_reference)
            else:
                # Create tmp series
                ts_other_tmp = pd.Series(data=value_other_cmp, index=[time_stamp_center], name=name_reference)
                # Append new line to series
                ts_other_out = pd.concat([ts_other_out, ts_other_tmp])

        return ts_other_out

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to join the time series
    @staticmethod
    def join_ts(ts_reference, ts_other, ts_name='ts', suffix_reference='_reference', suffix_other='_other'):

        ts_reference.name = ts_name + suffix_reference
        ts_other.name = ts_name + suffix_other

        dframe_joined = pd.concat((ts_reference, ts_other), axis=1).dropna()

        return dframe_joined
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to add seasonal description
    @staticmethod
    def add_dframe_seasons(dframe_data, column_season='season', lut_season=None):

        dframe_time = dframe_data.index

        grp_season = [lut_season.get(pd.Timestamp(t_stamp).month) for t_stamp in dframe_time]
        dframe_data[column_season] = grp_season

        return dframe_data
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize data
    def organize_data(self):

        log_stream.info(' ----> Organize soil moisture analysis ... ')

        file_path_anc_point_raw = self.file_path_anc_point
        file_path_anc_grid_raw = self.file_path_anc_grid
        file_path_dst_raw = self.file_path_dst

        geo_da_terrain = self.geo_dict[self.grid_terrain_tag]
        geo_da_cn = self.geo_dict[self.grid_cn_tag]
        geo_da_cnet = self.geo_dict[self.grid_cnet_tag]
        point_dframe_registry = self.geo_dict[self.points_registry_tag]

        flag_data_updating = self.flag_data_updating

        sm_analysis_collection = {}
        for point_id, point_row in point_dframe_registry.iterrows():

            point_name, point_code = point_row['point_name'], point_row['point_code']

            log_stream.info(' -----> Point "' + point_name + '" ... ')

            template_values_dict = {'point_code': point_code, 'point_name': point_name}
            file_path_anc_point_def = fill_tags2string(
                file_path_anc_point_raw,  self.template_tags_dict, template_values_dict)[0]
            file_path_anc_grid_def = fill_tags2string(
                file_path_anc_grid_raw, self.template_tags_dict, template_values_dict)[0]
            file_path_dst_def = fill_tags2string(
                file_path_dst_raw, self.template_tags_dict, template_values_dict)[0]

            if flag_data_updating:
                if os.path.exists(file_path_dst_def):
                    os.remove(file_path_dst_def)

            if not os.path.exists(file_path_dst_def):

                log_stream.info(' ------> Get datasets time-series ... ')
                if os.path.exists(file_path_anc_point_def):
                    point_ts = read_obj(file_path_anc_point_def)
                else:
                    log_stream.error(' ===> Point datasets is not available for point "' + point_name + '"')
                    raise FileNotFoundError('File "' + file_path_anc_point_def + '" not found')

                if os.path.exists(file_path_anc_grid_def):
                    grid_ts = read_obj(file_path_anc_grid_def)
                else:
                    log_stream.error(' ===> Grid datasets is not available for point "' + point_name + '"')
                    raise FileNotFoundError('File "' + file_path_anc_grid_def + '" not found')

                log_stream.info(' ------> Get datasets time-series ... DONE')

                log_stream.info(' ------> Order datasets time-series ... ')
                if self.alg_datasets_reference == 'grid' and self.alg_datasets_other == 'point':
                    ref_ts = deepcopy(grid_ts)
                    other_ts = deepcopy(point_ts)
                elif self.alg_datasets_reference == 'point' and self.alg_datasets_other == 'grid':
                    other_ts = deepcopy(grid_ts)
                    ref_ts = deepcopy(point_ts)
                else:
                    log_stream.error(' ===> Datasets order composed by reference: "' + self.alg_datasets_reference +
                                     '" and other "' + self.alg_datasets_reference + '" is not supported')
                    raise NotImplementedError('Check the datasets order; supported flag are "grid" and "point"')
                log_stream.info(' ------> Order datasets time-series ... DONE')

                log_stream.info(' ------> Adjust datasets time-series ... ')

                other_ts_synchronized = self.sync_ts(ref_ts, other_ts)
                ref_ts_filled, other_ts_filled = self.fill_ts(ref_ts, other_ts_synchronized)
                dframe_joined = self.join_ts(ref_ts_filled, other_ts_filled)

                dframe_seasonal = self.add_dframe_seasons(
                    dframe_joined, column_season='season', lut_season=self.lut_data_season)

                log_stream.info(' ------> Adjust datasets time-series ... DONE')

                '''
                # DEBUG START
                dframe_scaled = scaling.scale(dframe_joined, method='cdf_beta_match', reference_index=1)

                # now the scaled ascat data and insitu_sm are in the same space
                fig2, ax2 = plt.subplots()
                dframe_scaled.plot(figsize=(15, 5), title='scaled data', ax=ax2)
                plt.show()
                # DEBUG END
                '''

                # Save point collections
                log_stream.info(' ------> Save datasets time-series ... ')
                if dframe_seasonal is not None:
                    folder_name_dst, file_name_dst = os.path.split(file_path_dst_def)
                    make_folder(folder_name_dst)

                    write_obj(file_path_dst_def, dframe_seasonal)

                    log_stream.info(' ------> Save datasets time-series ... DONE')
                else:
                    log_stream.info(' ------> Save datasets time-series ... FAILED')
                    log_stream.error(' ===> Datasets are defined by NoneType')
                    raise IOError('Datasets must be defined to correctly run the algorithm')

                log_stream.info(' -----> Point "' + point_name + '" ... DONE')

            else:
                log_stream.info(' -----> Point "' + point_name + '" ... SKIPPED. Data previously computed')

        log_stream.info(' ----> Organize soil moisture analysis ... DONE')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
