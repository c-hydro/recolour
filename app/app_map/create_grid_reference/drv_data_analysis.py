"""
Class Features

Name:          drv_data_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_analysis import create_file_name, dump_file_tiff, warp_grid_data
from lib_data_io_geo import read_grid_data
from lib_data_io_pickle import write_file_obj, read_file_obj
from lib_utils_generic import make_folder
from lib_info_args import geo_coord_name_x, geo_coord_name_y, geo_dim_name_x, geo_dim_name_y
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_settings,
                 tag_section_flags='flags',
                 tag_section_grid='grid',
                 tag_section_parameters='parameters',
                 tag_section_log='log', tag_section_tmp='tmp'):

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_params = alg_settings[tag_section_parameters]

        self.alg_grid_src = alg_settings[tag_section_grid]['source']
        self.alg_grid_anc_raw = alg_settings[tag_section_grid]['ancillary']['raw']
        self.alg_grid_anc_def = alg_settings[tag_section_grid]['ancillary']['def']
        self.alg_grid_dst = alg_settings[tag_section_grid]['destination']

        self.alg_log = alg_settings[tag_section_log]
        self.alg_tmp = alg_settings[tag_section_tmp]

        self.alg_params_mask = self.alg_params['mask']
        self.alg_params_interp = self.alg_params['interpolate']

        self.reset_datasets_anc_raw = self.alg_flags['reset_datasets_ancillary_raw']
        self.reset_datasets_anc_def = self.alg_flags['reset_datasets_ancillary_def']
        self.reset_datasets_dst = self.alg_flags['reset_datasets_destination']
        self.reset_tmp = self.alg_flags['reset_tmp']
        self.reset_logs = self.alg_flags['reset_logs']

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'

        self.folder_name_anc_def = self.alg_grid_anc_def[self.tag_folder_name]
        self.file_name_anc_def = self.alg_grid_anc_def[self.tag_file_name]
        self.file_path_anc_def = os.path.join(self.folder_name_anc_def, self.file_name_anc_def)

        self.folder_name_tmp = self.alg_tmp[self.tag_folder_name]
        self.file_name_tmp = self.alg_tmp[self.tag_file_name]

    # method to remove tmp file(s)
    @staticmethod
    def __remove_tmp(file_list):
        for file_step in file_list:
            if os.path.exists(file_step):
                os.remove(file_step)

    # method to process data
    def process_data(self, grid_data_src, grid_attrs_src):

        # info start method
        logging.info(' ---> Process datasets ... ')
        # get parameters
        alg_params_mask, alg_params_interp = self.alg_params_mask, self.alg_params_interp

        # reset ancillary file (if flag is true)
        if self.reset_datasets_anc_def:
            if os.path.exists(self.file_path_anc_def):
                os.remove(self.file_path_anc_def)

        # check source data
        if (grid_data_src is not None) and (grid_attrs_src is not None):

            # check ancillary file availability
            if not os.path.exists(self.file_path_anc_def):

                # info start get datasets
                logging.info(' ----> (1) Get grid data ... ')
                # define tiff source file name
                file_path_src_tmp = create_file_name(folder_name=self.folder_name_tmp)
                # dump tiff source file name
                dump_file_tiff(file_path_src_tmp, grid_data_src, grid_attrs_src,
                               grid_no_data=alg_params_mask['grid_mask_no_data'])
                # info end get datasets
                logging.info(' ----> (1) Get grid data ... DONE')

                # info start warp datasets
                logging.info(' ----> (2) Warp grid data ... ')
                # define tiff warping file name
                file_path_dst_tmp = create_file_name(folder_name=self.folder_name_tmp)
                # warp tiff
                warp_grid_data(file_path_src_tmp, file_path_dst_tmp,
                               grid_warp_method=alg_params_interp['grid_interpolation_method'],
                               grid_warp_resolution=alg_params_interp['grid_interpolation_resolution'],
                               grid_warp_active=alg_params_interp['grid_interpolation_active'],
                               grid_warp_copy=True)

                # read grid da and attributes
                grid_da_dst, grid_attrs_dst = read_grid_data(
                    file_path_dst_tmp,
                    coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                    dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y)

                # info end warp datasets
                logging.info(' ----> (2) Warp grid data ... DONE')

                # info start warp datasets
                logging.info(' ----> (3) Organize grid data ... ')
                # organize grid obj
                grid_obj = {'grid_data': grid_da_dst, 'grid_attrs': grid_attrs_dst}

                # save grid obj
                folder_name_anc_def, file_name_anc_def = os.path.split(self.file_path_anc_def)
                make_folder(folder_name_anc_def)
                write_file_obj(self.file_path_anc_def, grid_obj)

                # remove tmp file(s)
                if self.reset_tmp:
                    self.__remove_tmp([file_path_src_tmp, file_path_dst_tmp])

                # info start warp datasets
                logging.info(' ----> (3) Organize grid data ... DONE')

                # info end method
                logging.info(' ---> Process datasets ... DONE')

            else:
                # read grid obj
                grid_obj = read_file_obj(self.file_path_anc_def)

                # info end method
                logging.info(' ---> Process datasets ... DONE. Datasets previously saved')

        else:

            # info end method
            grid_obj = {'grid_data': None, 'grid_attrs': None}
            logging.info(' ---> Process datasets ... SKIPPED. Destination datasets previously saved')

        return grid_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
