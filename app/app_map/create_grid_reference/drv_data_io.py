"""
Class Features

Name:          drv_data_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_analysis import dump_file_tiff

from lib_data_io_geo import read_grid_data, check_grid_data
from lib_data_io_pickle import write_file_obj, read_file_obj
from lib_utils_generic import make_folder
from lib_info_args import geo_coord_name_x, geo_coord_name_y, geo_dim_name_x, geo_dim_name_y

from lib_data_analysis import compute_grid_mask

# debugging
import matplotlib.pylab as plt
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

        self.folder_name_src = self.alg_grid_src[self.tag_folder_name]
        self.file_name_src = self.alg_grid_src[self.tag_file_name]
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)

        self.folder_name_anc_raw = self.alg_grid_anc_raw[self.tag_folder_name]
        self.file_name_anc_raw = self.alg_grid_anc_raw[self.tag_file_name]
        self.file_path_anc_raw = os.path.join(self.folder_name_anc_raw, self.file_name_anc_raw)

        self.folder_name_anc_def = self.alg_grid_anc_def[self.tag_folder_name]
        self.file_name_anc_def = self.alg_grid_anc_def[self.tag_file_name]
        self.file_path_anc_def = os.path.join(self.folder_name_anc_def, self.file_name_anc_def)

        self.folder_name_dst = self.alg_grid_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_grid_dst[self.tag_file_name]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

        self.folder_name_tmp = self.alg_tmp[self.tag_folder_name]
        self.file_name_tmp = self.alg_tmp[self.tag_file_name]
        self.file_path_tmp = os.path.join(self.folder_name_tmp, self.file_name_tmp)

    @staticmethod
    def __remove_tmp(file_list):
        for file_step in file_list:
            if os.path.exists(file_step):
                os.remove(file_step)

    # method to dump data
    def dump_data(self, grid_data_dst, grid_attrs_dst):

        # info start method
        logging.info(' ---> Dump destination datasets ... ')
        # get parameters
        alg_params_mask = self.alg_params_mask

        # check destination data
        if (grid_data_dst is not None) and (grid_attrs_dst is not None):
            # check destination file availability
            if not os.path.exists(self.file_path_dst):

                # dump tiff destination file name
                dump_file_tiff(self.file_path_dst, grid_data_dst, grid_attrs_dst,
                               grid_no_data=alg_params_mask['grid_mask_no_data'])

                # remove tmp file(s)
                if self.reset_tmp:
                    if self.reset_datasets_anc_raw:
                        self.__remove_tmp([self.file_path_anc_raw])
                    if self.reset_datasets_anc_def:
                        self.__remove_tmp([self.file_path_anc_def])

                # info end method
                logging.info(' ---> Dump destination datasets ... DONE')
            else:
                # info end method
                logging.info(' ---> Dump destination datasets ... DONE. Datasets previously saved')
        else:
            # info end method
            logging.info(' ---> Dump destination datasets ... DONE. Datasets previously saved')

    # method to get data
    def get_data(self):

        # info start method
        logging.info(' ---> Get source datasets ... ')
        # get parameters
        alg_params_mask = self.alg_params_mask

        # reset ancillary file (if flag is true)
        if self.reset_datasets_anc_raw or self.reset_datasets_anc_def:
            if os.path.exists(self.file_path_anc_raw):
                os.remove(self.file_path_anc_raw)
            if os.path.exists(self.file_path_anc_def):
                os.remove(self.file_path_anc_def)
            if os.path.exists(self.file_path_dst):
                os.remove(self.file_path_dst)
        # reset destination file (if flag is true)
        if self.reset_datasets_dst:
            if os.path.exists(self.file_path_dst):
                os.remove(self.file_path_dst)

        # check destination file availability
        if not os.path.exists(self.file_path_dst):
            # check ancillary file availability
            if not os.path.exists(self.file_path_anc_raw):

                # read grid da and attributes
                grid_da, grid_attrs = read_grid_data(
                    self.file_path_src,
                    coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                    dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y)
                # check grid da and attributes
                check_grid_data(self.file_path_src, grid_da, grid_attrs, data_mandatory=True)
                # create grid mask
                grid_da = compute_grid_mask(
                    grid_da,
                    grid_mask_min=alg_params_mask['grid_mask_min'],
                    grid_mask_max=alg_params_mask['grid_mask_max'],
                    grid_mask_data=alg_params_mask['grid_mask_data'],
                    grid_mask_no_data=alg_params_mask['grid_mask_no_data'],
                    grid_mask_format=alg_params_mask['grid_mask_format'],
                    grid_mask_active=alg_params_mask['grid_mask_active'])

                '''
                # debug
                plt.figure()
                plt.imshow(grid_da.values)
                plt.colorbar()
                plt.show()
                '''

                # organize grid obj
                grid_obj = {'grid_data': grid_da, 'grid_attrs': grid_attrs}

                # save grid obj
                folder_name_anc_raw, file_name_anc_raw = os.path.split(self.file_path_anc_raw)
                make_folder(folder_name_anc_raw)
                write_file_obj(self.file_path_anc_raw, grid_obj)

            else:

                # info end method
                logging.info(' ---> Get source datasets ... DONE. Datasets previously saved')
                # read grid obj
                grid_obj = read_file_obj(self.file_path_anc_raw)

        else:

            # info end method
            grid_obj = {'grid_data': None, 'grid_attrs': None}
            logging.info(' ---> Get source datasets ... SKIPPED. Destination datasets previously saved')

        return grid_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
