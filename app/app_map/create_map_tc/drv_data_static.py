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

from lib_data_io_cell import read_grid_cell, find_grid_fields, find_grid_cells
from lib_data_io_geo import read_grid_data, check_grid_data, remap_grid_data
from lib_data_io_pickle import read_file_obj, write_file_obj

from lib_info_args import logger_name
from lib_utils_generic import make_folder

# set logger obj
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_settings,
                 tag_section_flags='flags', tag_section_cells='cells',
                 tag_section_datasets='datasets',
                 tag_section_cell_grid='cell_grid',
                 tag_section_geo_grid='geo_grid',
                 tag_section_mask_grid='mask_grid', tag_section_removed_grid='removed_area_grid',
                 tag_section_log='log'):

        self.alg_flags = alg_settings[tag_section_flags]

        self.alg_cells = alg_settings[tag_section_cells]
        self.alg_datasets = alg_settings[tag_section_datasets]['static']

        self.alg_src_cell_grid = self.alg_datasets['source'][tag_section_cell_grid]
        self.alg_src_geo_grid = self.alg_datasets['source'][tag_section_geo_grid]
        self.alg_src_mask_grid = self.alg_datasets['source'][tag_section_mask_grid]
        self.alg_src_rm_grid = self.alg_datasets['source'][tag_section_removed_grid]
        self.alg_dst_workspace = self.alg_datasets['destination']

        self.alg_log = alg_settings[tag_section_log]

        self.tag_name_dset_ref, self.tag_name_dset_k1, self.tag_name_dset_k2 = 'ref', 'k1', 'k2'

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'

        self.reset_datasets_geo = self.alg_flags['reset_datasets_geo']

        self.folder_name_src_cell_grid_ref = self.alg_src_cell_grid[self.tag_name_dset_ref][self.tag_folder_name]
        self.file_name_src_cell_grid_ref = self.alg_src_cell_grid[self.tag_name_dset_ref][self.tag_file_name]
        self.file_path_src_cell_grid_ref = os.path.join(
            self.folder_name_src_cell_grid_ref, self.file_name_src_cell_grid_ref)

        self.folder_name_src_cell_grid_k1 = self.alg_src_cell_grid[self.tag_name_dset_k1][self.tag_folder_name]
        self.file_name_src_cell_grid_k1 = self.alg_src_cell_grid[self.tag_name_dset_k1][self.tag_file_name]
        self.file_path_src_cell_grid_k1 = os.path.join(
            self.folder_name_src_cell_grid_k1, self.file_name_src_cell_grid_k1)

        self.folder_name_src_cell_grid_k2 = self.alg_src_cell_grid[self.tag_name_dset_k2][self.tag_folder_name]
        self.file_name_src_cell_grid_k2 = self.alg_src_cell_grid[self.tag_name_dset_k2][self.tag_file_name]
        self.file_path_src_cell_grid_k2 = os.path.join(
            self.folder_name_src_cell_grid_k2, self.file_name_src_cell_grid_k1)

        self.folder_name_src_geo_grid = self.alg_src_geo_grid[self.tag_folder_name]
        self.file_name_src_geo_grid = self.alg_src_geo_grid[self.tag_file_name]
        self.file_path_src_geo_grid = os.path.join(self.folder_name_src_geo_grid, self.file_name_src_geo_grid)

        self.folder_name_src_mask_grid = self.alg_src_mask_grid[self.tag_folder_name]
        self.file_name_src_mask_grid = self.alg_src_mask_grid[self.tag_file_name]
        self.file_path_src_mask_grid = os.path.join(self.folder_name_src_mask_grid, self.file_name_src_mask_grid)

        self.folder_name_src_rm_grid = self.alg_src_rm_grid[self.tag_folder_name]
        self.file_name_src_rm_grid = self.alg_src_rm_grid[self.tag_file_name]
        self.file_path_src_rm_grid = os.path.join(self.folder_name_src_rm_grid, self.file_name_src_rm_grid)

        self.folder_name_dst_workspace = self.alg_dst_workspace[self.tag_folder_name]
        self.file_name_dst_workspace = self.alg_dst_workspace[self.tag_file_name]
        self.file_path_dst_workspace = os.path.join(self.folder_name_dst_workspace, self.file_name_dst_workspace)

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ---> Organize static datasets ... ')

        # check reset flag
        if self.reset_datasets_geo:
            if os.path.exists(self.file_path_dst_workspace):
                os.remove(self.file_path_dst_workspace)

        # check obj workspace
        if not os.path.exists(self.file_path_dst_workspace):

            # read geo grid
            obj_geo_data, obj_geo_attrs, obj_geo_lons, obj_geo_lats = read_grid_data(
                self.file_path_src_geo_grid, var_idx=0)
            check_grid_data(self.file_path_src_geo_grid, obj_geo_data, obj_geo_attrs, data_mandatory=True)

            # read mask grid (areas to be update)
            tmp_mask_data, tmp_mask_attrs, tmp_mask_lons, tmp_mask_lats = read_grid_data(
                self.file_path_src_mask_grid, var_idx=2)
            check_grid_data(self.file_path_src_mask_grid, tmp_mask_data, tmp_mask_attrs, data_mandatory=False)
            # remap mask grid (if dimensions are not the same)
            if tmp_mask_data is not None:
                obj_mask_data = remap_grid_data(
                    tmp_mask_data, tmp_mask_lons, tmp_mask_lats,
                    obj_geo_lons, obj_geo_lats, obj_geo_attrs,
                    interp_method='nearest', interp_no_data=0)
                obj_mask_attrs = obj_mask_data.attrs
            else:
                alg_logger.warning(' ===> Mask grid is not available. Check the settings')
                obj_mask_data, obj_mask_attrs = None, None

            # read remove grid (areas to removed)
            tmp_rm_data, tmp_rm_attrs, tmp_rm_lons, tmp_rm_lats = read_grid_data(
                self.file_path_src_rm_grid, var_idx=2)
            check_grid_data(self.file_path_src_rm_grid, tmp_rm_data, tmp_rm_attrs, data_mandatory=False)
            # remap removed grid (if dimensions are not the same)
            if tmp_rm_data is not None:
                obj_rm_data = remap_grid_data(
                    tmp_rm_data, tmp_rm_lons, tmp_rm_lats,
                    obj_geo_lons, obj_geo_lats, obj_geo_attrs,
                    interp_method='nearest', interp_no_data=0)
                obj_rm_attrs = obj_rm_data.attrs
            else:
                alg_logger.warning(' ===> Mask grid is not available. Check the settings')
                obj_rm_data, obj_rm_attrs = None, None

            # find grid cells
            obj_grid_cells = find_grid_cells(
                cells_list=self.alg_cells['cell_list'],
                cell_start=self.alg_cells['cell_start'], cell_end=self.alg_cells['cell_end'])

            # read grid obj
            obj_grid_ref = read_grid_cell(self.file_path_src_cell_grid_ref)
            obj_grid_k1 = read_grid_cell(self.file_path_src_cell_grid_k1)
            obj_grid_k2 = read_grid_cell(self.file_path_src_cell_grid_k2)
            # find grid fields
            obj_fields_ref, obj_fields_k1, obj_fields_k2 = find_grid_fields(
                obj_grid_ref, obj_grid_k1, obj_grid_k2,
                max_distance_k1=self.alg_cells['max_distance']['k1'],
                max_distance_k2=self.alg_cells['max_distance']['k2'])

            # organize object workspace
            obj_workspace = {
                'geo': {'data': obj_geo_data, 'attrs': obj_geo_attrs},
                'mask': {'data': obj_mask_data, 'attrs': obj_mask_attrs},
                'removed': {'data': obj_rm_data, 'attrs': obj_rm_attrs},
                "ref": obj_fields_ref, 'k1': obj_fields_k1, 'k2': obj_fields_k2,
                'cells': obj_grid_cells
            }

            # save object workspace
            folder_name_dst, file_name_dst = os.path.split(self.file_path_dst_workspace)
            make_folder(folder_name_dst)
            write_file_obj(self.file_path_dst_workspace, obj_workspace)

        else:
            # read object workspace
            obj_workspace = read_file_obj(self.file_path_dst_workspace)

        # info end method
        alg_logger.info(' ---> Organize static datasets ... DONE')

        return obj_workspace

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
