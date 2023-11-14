"""
Class Features

Name:          drv_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_io_geo import read_grid_data, check_grid_data, write_grid_data
from lib_data_io_pickle import write_file_obj, read_file_obj
from lib_info_args import geo_coord_name_x, geo_coord_name_y, geo_dim_name_x, geo_dim_name_y

from lib_utils_generic import make_folder
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_settings,
                 tag_section_flags='flags',
                 tag_section_grid='grid_reference', tag_section_cells='cells_reference',
                 tag_section_datasets='datasets',
                 tag_section_log='log'):

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_datasets_src = alg_settings[tag_section_datasets]['static']['source']
        self.alg_datasets_anc = alg_settings[tag_section_datasets]['static']['ancillary']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['static']['destination']

        self.alg_grid_ref_src = self.alg_datasets_src[tag_section_grid]
        self.alg_grid_ref_anc = self.alg_datasets_anc[tag_section_grid]
        self.alg_grid_ref_dst = self.alg_datasets_dst[tag_section_grid]
        self.alg_cells_ref_dst = self.alg_datasets_dst[tag_section_cells]

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name, self.tag_file_name = 'folder_name', 'file_name'

        self.reset_geo_anc = self.alg_flags['reset_geo_ancillary']

        self.folder_name_grid_ref_src = self.alg_grid_ref_src[self.tag_folder_name]
        self.file_name_grid_ref_src = self.alg_grid_ref_src[self.tag_file_name]
        self.file_path_grid_ref_src = os.path.join(self.folder_name_grid_ref_src, self.file_name_grid_ref_src)

        self.folder_name_grid_ref_anc = self.alg_grid_ref_anc[self.tag_folder_name]
        self.file_name_grid_ref_anc = self.alg_grid_ref_anc[self.tag_file_name]
        self.file_path_grid_ref_anc = os.path.join(self.folder_name_grid_ref_anc, self.file_name_grid_ref_anc)

        self.folder_name_grid_ref_dst = self.alg_grid_ref_dst[self.tag_folder_name]
        self.file_name_grid_ref_dst = self.alg_grid_ref_dst[self.tag_file_name]
        self.file_path_grid_ref_dst = os.path.join(self.folder_name_grid_ref_dst, self.file_name_grid_ref_dst)

        self.folder_name_cells_ref_dst = self.alg_cells_ref_dst[self.tag_folder_name]
        self.file_name_cells_ref_dst = self.alg_cells_ref_dst[self.tag_file_name]
        self.file_path_cells_ref_dst = os.path.join(self.folder_name_cells_ref_dst, self.file_name_cells_ref_dst)

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize static datasets ... ')

        # reset ancillary and destination file
        if self.reset_geo_anc:
            if os.path.exists(self.file_path_grid_ref_anc):
                os.remove(self.file_path_grid_ref_anc)
            if os.path.exists(self.file_path_grid_ref_dst):
                os.remove(self.file_path_grid_ref_dst)

        # check file ancillary availability
        if not os.path.exists(self.file_path_grid_ref_anc):

            # read grid data
            grid_data_ref_src, grid_attrs_ref_src = read_grid_data(
                self.file_path_grid_ref_src,
                coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y)
            # check grid data
            check_grid_data(self.file_path_grid_ref_src, grid_data_ref_src, grid_attrs_ref_src, data_mandatory=True)
            # write grid data
            write_grid_data(self.file_path_grid_ref_dst, grid_data_ref_src)

            # organize grid obj
            grid_obj = {'grid_data': grid_data_ref_src, 'grid_attrs': grid_attrs_ref_src}

            # save grid obj
            folder_name_grid_ref_anc, file_name_grid_ref_anc = os.path.split(self.file_path_grid_ref_anc)
            make_folder(folder_name_grid_ref_anc)
            write_file_obj(self.file_path_grid_ref_anc, grid_obj)

        else:
            # read grid obj
            grid_obj = read_file_obj(self.file_path_grid_ref_anc)

        # info end method
        logging.info(' ---> Organize static datasets ... DONE')

        return grid_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
