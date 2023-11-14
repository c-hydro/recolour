"""
Class Features

Name:          driver_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import os
from copy import deepcopy

import pandas as pd

from lib_data_io_geo import read_grid_data, read_point_data, check_obj_data
from lib_data_io_pickle import read_obj, write_obj

from lib_utils_geo_OLD import get_grid_value_from_xy, get_grid_idx_from_xy, get_idx_by_win
from lib_utils_system import fill_tags2string, make_folder

from lib_info_args import logger_name
from lib_info_args import (geo_coord_name_x, geo_coord_name_y,
                           geo_var_name_x, geo_var_name_y,
                           geo_dim_name_x, geo_dim_name_y)

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class DriverData
class DriverData:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, src_dict, dst_dict, tmp_dict=None,
                 template_dict=None, flags_dict=None, params_dict=None):

        self.src_dict = src_dict
        self.dst_dict = dst_dict
        self.tmp_dict = tmp_dict

        self.template_dict = template_dict
        self.flags_dict = flags_dict
        self.params_dict = params_dict

        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.fields_tag = 'fields'

        self.grid_reference_tag, self.point_reference_tag = 'grid_reference', 'point_reference'

        self.reset_static = flags_dict['reset_static']

        self.geo_method_search = self.params_dict['geo_method_search']
        self.geo_radius_influence = self.params_dict['geo_radius_influence']
        self.geo_neighbours = self.params_dict['geo_neighbours']
        self.geo_spatial_operation = self.params_dict['geo_spatial_operation']
        self.geo_spatial_window = self.params_dict['geo_spatial_window']
        self.geo_spatial_mask = self.params_dict['geo_spatial_mask']

        self.field_value_point_tag = 'point_value'
        self.field_idx_1d_point_tag = 'point_idx_1d'
        self.field_idx_2d_x_point_tag = 'point_idx_2d_x'
        self.field_idx_2d_y_point_tag = 'point_idx_2d_y'
        self.field_longitude_grid_tag = 'grid_longitude'
        self.field_latitude_grid_tag = 'grid_latitude'
        self.field_distance_grid_tag = 'grid_distance'

        # source object(s)
        folder_name_src_grid = self.src_dict[self.grid_reference_tag][self.folder_name_tag]
        file_name_src_grid = self.src_dict[self.grid_reference_tag][self.file_name_tag]
        self.file_path_src_grid = os.path.join(folder_name_src_grid, file_name_src_grid)
        self.fields_src_grid = self.src_dict[self.grid_reference_tag][self.fields_tag]

        folder_name_src_point = self.src_dict[self.point_reference_tag][self.folder_name_tag]
        file_name_src_point = self.src_dict[self.point_reference_tag][self.file_name_tag]
        self.file_path_src_point = os.path.join(folder_name_src_point, file_name_src_point)
        self.fields_src_point = self.src_dict[self.point_reference_tag][self.fields_tag]

        # destination object(s)
        folder_name_dst = self.dst_dict[self.folder_name_tag]
        file_name_dst = self.dst_dict[self.file_name_tag]
        self.file_path_dst = os.path.join(folder_name_dst, file_name_dst)

        # tmp object(s)
        self.folder_name_tmp_raw = tmp_dict[self.folder_name_tag]
        self.file_name_tmp_raw = tmp_dict[self.file_name_tag]

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to join geo obj
    def join_geo_obj(self, grid_obj, point_obj):

        log_stream.info(' -----> Join the point and the grid information ... ')

        grid_x_1d = grid_obj[geo_coord_name_x].values
        grid_y_1d = grid_obj[geo_coord_name_y].values
        grid_x_2d, grid_y_2d = np.meshgrid(grid_x_1d, grid_y_1d)

        point_dframe_out = pd.DataFrame()
        for point_id, point_row_in in point_obj.iterrows():

            # copy the point obj
            point_row_out = deepcopy(point_row_in)

            # get point name
            point_name, point_tag = point_row_in['name'], point_row_in['tag']
            # get point location
            point_x_center = point_row_in['longitude']
            point_y_center = point_row_in['latitude']

            # info start point
            log_stream.info(' ------> Point "' + point_tag + '" ... ')

            grid_value_center, grid_x_center, grid_y_center = get_grid_value_from_xy(
                grid_obj, point_x_center, point_y_center, select_method=self.geo_method_search)

            idx_cont_center, idx_x_center, idx_y_center, distance_center = get_grid_idx_from_xy(
                grid_x_center, grid_y_center, grid_x_2d, grid_y_2d,
                geo_radius_influence=self.geo_radius_influence, geo_neighbours=self.geo_neighbours)

            point_idx_window, point_xy_window = get_idx_by_win(
                idx_x_center, idx_y_center, grid_x_2d, grid_y_2d,
                geo_win_x=self.geo_spatial_window, geo_win_y=self.geo_spatial_window)

            # Iterate over a window around the center point
            grid_list_value = []
            grid_list_x, grid_list_y = [], []
            idx_list_cont, idx_list_x, idx_list_y, distance_list = [], [], [], []
            for point_idx_step, point_xy_step in zip(point_idx_window, point_xy_window):

                point_x_step, point_y_step = point_xy_step[0], point_xy_step[1]

                grid_value_window, grid_x_window, grid_y_window = get_grid_value_from_xy(
                    grid_obj, point_x_step, point_y_step, select_method=self.geo_method_search)

                idx_cont_window, idx_x_window, idx_y_window, distance_window = get_grid_idx_from_xy(
                    grid_x_window, grid_y_window, grid_x_2d, grid_y_2d,
                    geo_radius_influence=self.geo_radius_influence,
                    geo_neighbours=self.geo_neighbours)

                grid_list_value.append(float(grid_value_window))
                grid_list_x.append(float(grid_x_window))
                grid_list_y.append(float(grid_y_window))

                idx_list_cont.append(int(idx_cont_window))
                idx_list_x.append(int(idx_x_window))
                idx_list_y.append(int(idx_y_window))
                distance_list.append(float(distance_window))

            assert point_row_in['longitude'] == point_x_center, \
                'coordinate x of geo point "' + point_name + '" is not the same found in the source file'
            assert point_row_in['latitude'] == point_y_center, \
                'coordinate y of geo point "' + point_name + '" is not the same found in the source file'

            point_row_out[self.field_value_point_tag] = grid_list_value
            point_row_out[self.field_idx_1d_point_tag] = idx_list_cont
            point_row_out[self.field_idx_2d_x_point_tag] = idx_list_x
            point_row_out[self.field_idx_2d_y_point_tag] = idx_list_y
            point_row_out[self.field_distance_grid_tag] = distance_list
            point_row_out[self.field_longitude_grid_tag] = grid_list_x
            point_row_out[self.field_latitude_grid_tag] = grid_list_y

            point_dframe_out = point_dframe_out.append(point_row_out, ignore_index=True)

            log_stream.info(' ------> Point "' + point_tag + '" ... DONE')

        # geo_point_dframe_out = geo_point_dframe_out.set_index('point_name')

        log_stream.info(' -----> Join the point and the grid information ... DONE')

        return point_dframe_out

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get geo grid file
    @staticmethod
    def get_geo_grid(file_name):

        log_stream.info(' -----> Read file grid "' + file_name + '" ... ')
        if ((file_name.endswith('txt') or file_name.endswith('asc')) or
                (file_name.endswith('tiff') or file_name.endswith('tif'))):
            grid_obj = read_grid_data(file_name, output_format='data_array')
            check_obj_data(file_name, grid_obj, data_mandatory=True)
            log_stream.info(' -----> Read file grid "' + file_name + '" ... DONE')
        else:
            log_stream.info(' -----> Read file grid "' + file_name + '" ... FAILED')
            log_stream.error(' ===> File is mandatory to run the application')
            raise FileNotFoundError('File not found')

        return grid_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get geo point file
    def get_geo_point(self, file_name):

        log_stream.info(' -----> Read file point "' + file_name + '" ... ')

        if file_name.endswith('csv') or file_name.endswith('txt'):
            point_obj = read_point_data(file_name,
                                        file_columns_remap=self.fields_src_point)
            log_stream.info(' -----> Read file point "' + file_name + '" ... DONE')
        else:
            log_stream.info(' -----> Read file point "' + file_name + '" ... FAILED')
            log_stream.error(' ===> File is mandatory to run the application')
            raise FileNotFoundError('File not found')

        return point_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize static information ... ')

        # get file path(s)
        file_path_src_grid, file_path_src_point = self.file_path_src_grid, self.file_path_src_point
        file_path_dst = self.file_path_dst
        # get flag(s)
        reset_static = self.reset_static

        # reset destination file (if needed)
        if reset_static:
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)

        # check destination file availability
        if not os.path.exists(file_path_dst):

            # get grid object
            obj_grid_base, obj_grid_attrs = self.get_geo_grid(file_path_src_grid)
            # get point obj
            obj_point_base = self.get_geo_point(file_path_src_point)

            # join grid and point object(s)
            obj_point_extended = self.join_geo_obj(obj_grid_base, obj_point_base)

            # merge objects to collections
            obj_collections = {
                'grid_obj': obj_grid_base, 'grid_attrs': obj_grid_attrs, 'point_obj': obj_point_extended}

            # dump collections to destination file
            folder_name_anc, file_name_anc = os.path.split(file_path_dst)
            make_folder(folder_name_anc)

            write_obj(file_path_dst, obj_collections)

        else:

            # get collections from destination file
            obj_collections = read_obj(file_path_dst)

        # method end info
        log_stream.info(' ----> Organize static information ... DONE')

        return obj_collections
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
