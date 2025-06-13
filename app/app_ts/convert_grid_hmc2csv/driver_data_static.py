"""
Class Features

Name:          driver_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240702'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import os
from copy import deepcopy

import pandas as pd

from lib_data_io_geo import read_grid_data, read_point_data, check_obj_data
from lib_data_io_pickle import read_obj, write_obj

from lib_utils_geo import get_grid_value_from_xy, get_grid_idx_from_xy, get_idx_by_win
from lib_utils_geo import compute_volume_max, mask_volume_max
from lib_utils_system import fill_tags2string, make_folder

from lib_utils_obj import create_darray_2d

from lib_info_args import logger_name
from lib_info_args import (geo_coord_name_x, geo_coord_name_y,
                           geo_var_name_x, geo_var_name_y,
                           geo_dim_name_x, geo_dim_name_y)

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class DriverData
class DriverData:

    # ------------------------------------------------------------------------------------------------------------------
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
        self.format_tag, self.fields_tag = 'format', 'fields'
        self.delimiter_tag = 'delimiter'

        self.point_reference_tag = 'point_reference'
        self.grid_reference_tag, self.soil_type_reference_tag = 'grid_reference', 'soil_type_reference'
        self.river_network_reference_tag = 'river_network_reference'

        self.reset_static = flags_dict['reset_static']

        if 'format_static_point' in list(params_dict.keys()):
            self.format_src_point = params_dict['format_static_point']
        else:
            self.format_src_point = 'csv'
        if 'format_static_grid' in list(params_dict.keys()):
            self.format_src_grid = params_dict['format_static_grid']
        else:
            self.format_src_grid = 'tiff'
        self.format_src_pnt = self.params_dict['format_static_point']
        self.format_src_grid = self.params_dict['format_static_grid']
        if 'format_static_soil_type' in list(params_dict.keys()):
            self.format_src_soil_type = self.params_dict['format_static_soil_type']
        else:
            self.format_src_soil_type = None
        if 'format_static_river_network' in list(params_dict.keys()):
            self.format_src_river_net = self.params_dict['format_static_river_network']
        else:
            self.format_src_river_net = None
        self.geo_method_search = self.params_dict['geo_method_search']
        self.geo_radius_influence = self.params_dict['geo_radius_influence']
        self.geo_neighbours = self.params_dict['geo_neighbours']
        self.geo_spatial_operation = self.params_dict['geo_spatial_operation']
        self.geo_spatial_window = self.params_dict['geo_spatial_window']
        self.geo_spatial_mask = self.params_dict['geo_spatial_mask']

        self.field_value_point_tag = 'point_value'
        self.field_soil_type_point_tag = 'point_soil_type'
        self.field_soil_vmax_point_tag = 'point_soil_vmax'
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

        if self.soil_type_reference_tag in list(self.src_dict.keys()):
            folder_name_src_soil_type = self.src_dict[self.soil_type_reference_tag][self.folder_name_tag]
            file_name_src_soil_type = self.src_dict[self.soil_type_reference_tag][self.file_name_tag]
            self.file_path_src_soil_type = os.path.join(folder_name_src_soil_type, file_name_src_soil_type)
            self.fields_src_soil_type = self.src_dict[self.soil_type_reference_tag][self.fields_tag]
        else:
            self.file_path_src_soil_type, self.fields_src_soil_type = None, None

        if self.river_network_reference_tag in list(self.src_dict.keys()):
            folder_name_src_river_net = self.src_dict[self.river_network_reference_tag][self.folder_name_tag]
            file_name_src_river_net = self.src_dict[self.river_network_reference_tag][self.file_name_tag]
            self.file_path_src_river_net = os.path.join(folder_name_src_river_net, file_name_src_river_net)
            self.fields_src_river_net = self.src_dict[self.river_network_reference_tag][self.fields_tag]
        else:
            self.file_path_src_river_net, self.fields_src_river_net = None, None

        folder_name_src_point = self.src_dict[self.point_reference_tag][self.folder_name_tag]
        file_name_src_point = self.src_dict[self.point_reference_tag][self.file_name_tag]
        self.file_path_src_point = os.path.join(folder_name_src_point, file_name_src_point)
        self.fields_src_point = self.src_dict[self.point_reference_tag][self.fields_tag]

        if self.delimiter_tag in list(self.src_dict[self.point_reference_tag].keys()):
            self.delimiter_src_point = self.src_dict[self.point_reference_tag][self.delimiter_tag]
        else:
            self.delimiter_src_point = ','

        # destination object(s)
        folder_name_dst = self.dst_dict[self.folder_name_tag]
        file_name_dst = self.dst_dict[self.file_name_tag]
        self.file_path_dst = os.path.join(folder_name_dst, file_name_dst)

        # tmp object(s)
        self.folder_name_tmp_raw = tmp_dict[self.folder_name_tag]
        self.file_name_tmp_raw = tmp_dict[self.file_name_tag]

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to join geo grid obj
    def join_geo_grid_obj(self, soil_type_obj, soil_type_attrs, river_net_obj, river_net_attrs):

        log_stream.info(' -----> Join the grid and the soil information ... ')

        if (soil_type_obj is not None) and (river_net_obj is not None):

            grid_soil_type_values = soil_type_obj.values
            grid_x_1d = soil_type_obj[geo_coord_name_x].values
            grid_y_1d = soil_type_obj[geo_coord_name_y].values
            grid_x_2d, grid_y_2d = np.meshgrid(grid_x_1d, grid_y_1d)

            grid_river_net_values = river_net_obj.values

            soil_vmax_attrs = deepcopy(soil_type_attrs)

            grid_soil_vmax_values = compute_volume_max(grid_soil_type_values,
                                                       no_data_value=-9999, fill_value=np.nan)
            grid_soil_vmax_values = mask_volume_max(grid_soil_vmax_values, grid_river_net_values,
                                                    river_value=1, fill_value=np.nan)

            ''' debug
            import matplotlib.pylab as plt
            plt.figure()
            plt.imshow(grid_soil_type_values)
            plt.colorbar()
            plt.figure()
            plt.imshow(grid_soil_vmax_values)
            plt.colorbar()
            plt.show()
            '''

            soil_vmax_obj = create_darray_2d(
                grid_soil_vmax_values, grid_x_2d[0, :], grid_y_2d[:, 0],
                coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y)

            soil_vmax_obj.attrs = soil_vmax_attrs

            log_stream.info(' -----> Join the grid and the soil information ... DONE')

        else:
            log_stream.info(' -----> Join the grid and the soil information ... SKIPPED. '
                            'Soil type information not available')
            soil_vmax_obj, soil_vmax_attrs = None, None

        return soil_vmax_obj, soil_vmax_attrs

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to join geo point obj
    def join_geo_point_obj(self, grid_obj, soil_type_obj, soil_vmax_obj, point_obj):

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
            grid_list_soil_type, grid_list_soil_vmax = [], []
            grid_list_x, grid_list_y = [], []
            idx_list_cont, idx_list_x, idx_list_y, distance_list = [], [], [], []
            for point_idx_step, point_xy_step in zip(point_idx_window, point_xy_window):

                point_x_step, point_y_step = point_xy_step[0], point_xy_step[1]

                grid_value_window, grid_x_window, grid_y_window = get_grid_value_from_xy(
                    grid_obj, point_x_step, point_y_step, select_method=self.geo_method_search)

                if (soil_type_obj is not None) and (soil_vmax_obj is not None):
                    soil_type_window, _, _ = get_grid_value_from_xy(
                        soil_type_obj, point_x_step, point_y_step, select_method=self.geo_method_search)
                    soil_vmax_window, _, _ = get_grid_value_from_xy(
                        soil_vmax_obj, point_x_step, point_y_step, select_method=self.geo_method_search)
                else:
                    soil_type_window, soil_vmax_window = np.nan, np.nan

                idx_cont_window, idx_x_window, idx_y_window, distance_window = get_grid_idx_from_xy(
                    grid_x_window, grid_y_window, grid_x_2d, grid_y_2d,
                    geo_radius_influence=self.geo_radius_influence,
                    geo_neighbours=self.geo_neighbours)

                grid_list_value.append(float(grid_value_window))
                grid_list_soil_type.append(float(soil_type_window))
                grid_list_soil_vmax.append(float(soil_vmax_window))
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
            point_row_out[self.field_soil_type_point_tag] = grid_list_soil_type
            point_row_out[self.field_soil_vmax_point_tag] = grid_list_soil_vmax
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

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to get geo grid file
    @staticmethod
    def get_geo_grid(file_name, data_mandatory=True):

        log_stream.info(' -----> Read file grid "' + file_name + '" ... ')
        if ((file_name.endswith('txt') or file_name.endswith('asc')) or
                (file_name.endswith('tiff') or file_name.endswith('tif'))):
            grid_obj = read_grid_data(file_name, output_format='data_array')
            check_obj_data(file_name, grid_obj, data_mandatory=data_mandatory)
            log_stream.info(' -----> Read file grid "' + file_name + '" ... DONE')
        else:
            log_stream.info(' -----> Read file grid "' + file_name + '" ... FAILED')
            log_stream.error(' ===> File is mandatory to run the application')
            raise FileNotFoundError('File not found')

        return grid_obj
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to get geo point file
    def get_geo_point(self, file_name):

        log_stream.info(' -----> Read file point "' + file_name + '" ... ')

        if file_name.endswith('csv') or file_name.endswith('txt'):
            point_obj = read_point_data(
                file_name, file_delimiter=self.delimiter_src_point,
                file_columns_remap=self.fields_src_point)
            log_stream.info(' -----> Read file point "' + file_name + '" ... DONE')
        else:
            log_stream.info(' -----> Read file point "' + file_name + '" ... FAILED')
            log_stream.error(' ===> File is mandatory to run the application')
            raise FileNotFoundError('File not found')

        return point_obj
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize static information ... ')

        # get file path(s)
        file_path_src_grid, file_path_src_soil_type = self.file_path_src_grid, self.file_path_src_soil_type
        file_path_src_river_net = self.file_path_src_river_net
        file_path_src_point = self.file_path_src_point
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
            if self.format_src_grid == 'tiff':
                obj_grid_base, obj_grid_attrs = self.get_geo_grid(file_path_src_grid, data_mandatory=True)
            elif self.format_src_grid == 'ascii':
                obj_grid_base, obj_grid_attrs = self.get_geo_grid(file_path_src_grid, data_mandatory=True)
            else:
                log_stream.error(' ===> File grid format "' + self.format_src_grid + '" not expected')
                raise NotImplementedError('Case not implemented yet')

            if obj_grid_base is None:
                log_stream.error(' ===> File grid "' + str(file_path_src_grid) + '" is not correctly read')
                raise RuntimeError('Base grid is defined by NoneType. Grid must be defined by object')

            # get soil type object
            if self.format_src_soil_type is None:
                log_stream.warning(' ===> File "soil_type" format is not defined. Datasets will not available')
                obj_soil_type_base, obj_soil_type_attrs = None, None
            elif self.format_src_soil_type == 'ascii':
                obj_soil_type_base, obj_soil_type_attrs = self.get_geo_grid(
                    file_path_src_soil_type, data_mandatory=False)
            else:
                log_stream.error(' ===> File grid format "' + self.format_src_soil_type + '" not expected')
                raise NotImplementedError('Case not implemented yet')

            # get soil object
            if self.format_src_river_net is None:
                log_stream.warning(' ===> File "river_network" is not defined. Datasets will not available')
                obj_river_net_base, obj_river_net_attrs = None, None
            elif self.format_src_soil_type == 'ascii':
                obj_river_net_base, obj_river_net_attrs = self.get_geo_grid(
                    file_path_src_river_net, data_mandatory=False)
            else:
                log_stream.error(' ===> File grid format "' + self.format_src_river_net + '" not expected')
                raise NotImplementedError('Case not implemented yet')

            # get point obj
            if self.format_src_point == 'csv':
                obj_point_base = self.get_geo_point(file_path_src_point)
            else:
                log_stream.error(' ===> File point format "' + self.format_src_point + '" not expected')
                raise NotImplementedError('Case not implemented yet')

            # join grid and soil object(s)
            obj_soil_vmax_base, obj_soil_vmax_attrs = self.join_geo_grid_obj(
                obj_soil_type_base, obj_soil_type_attrs, obj_river_net_base, obj_river_net_attrs)

            # join grid and point object(s)
            obj_point_extended = self.join_geo_point_obj(
                obj_grid_base, obj_soil_type_base, obj_soil_vmax_base, obj_point_base)

            # merge objects to collections
            obj_collections = {
                'grid_obj': obj_grid_base, 'grid_attrs': obj_grid_attrs,
                'soil_type_obj': obj_soil_type_base, 'soil_type_attrs': obj_soil_type_attrs,
                'soil_vmax_obj': obj_soil_vmax_base, 'soil_vmax_attrs': obj_soil_vmax_attrs,
                'river_network_obj': obj_river_net_base, 'river_network_attrs': obj_river_net_attrs,
                'point_obj': obj_point_extended}

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
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
