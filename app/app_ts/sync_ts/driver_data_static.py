"""
Class Features

Name:          driver_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20200515'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_io_geo import read_point_data
from lib_data_io_pickle import read_obj, write_obj

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
    def __init__(self, src_dict, dst_dict=None, tmp_dict=None,
                 template_dict=None, flags_dict=None, params_dict=None):

        self.src_dict = src_dict
        self.dst_dict = dst_dict
        self.tmp_dict = tmp_dict

        self.template_dict = template_dict
        self.flags_dict = flags_dict
        self.params_dict = params_dict

        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.fields_tag, self.delimiter_tag = 'fields', 'delimiter'

        self.reset_static = flags_dict['reset_static']

        # source object(s)
        folder_name_src = self.src_dict[self.folder_name_tag]
        file_name_src = self.src_dict[self.file_name_tag]
        self.file_path_src = os.path.join(folder_name_src, file_name_src)
        self.fields_src = self.src_dict[self.fields_tag]

        if self.delimiter_tag in list(self.src_dict.keys()):
            self.delimiter_src = self.src_dict[self.delimiter_tag]
        else:
            self.delimiter_src = ','

        # destination object(s)
        folder_name_dst = self.dst_dict[self.folder_name_tag]
        file_name_dst = self.dst_dict[self.file_name_tag]
        self.file_path_dst = os.path.join(folder_name_dst, file_name_dst)

        # tmp object(s)
        self.folder_name_tmp_raw = tmp_dict[self.folder_name_tag]
        self.file_name_tmp_raw = tmp_dict[self.file_name_tag]

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get geo point file
    def get_geo_point(self, file_name):

        log_stream.info(' -----> Read file point "' + file_name + '" ... ')

        if file_name.endswith('csv') or file_name.endswith('txt'):
            point_obj = read_point_data(
                file_name, file_delimiter=self.delimiter_src,
                file_columns_remap=self.fields_src)
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
        file_path_src, file_path_dst = self.file_path_src, self.file_path_dst
        # get flag(s)
        reset_static = self.reset_static

        # reset destination file (if needed)
        if reset_static:
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)

        # check destination file availability
        if not os.path.exists(file_path_dst):

            # get point obj
            obj_point = self.get_geo_point(file_path_src)

            # merge objects to collections
            obj_collections = {'point_obj': obj_point}

            # dump object to destination file
            folder_name_anc, file_name_anc = os.path.split(file_path_dst)
            make_folder(folder_name_anc)

            write_obj(file_path_dst, obj_collections)

        else:

            # get object from destination file
            obj_collections = read_obj(file_path_dst)

        # method end info
        log_stream.info(' ----> Organize static information ... DONE')

        return obj_collections
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
