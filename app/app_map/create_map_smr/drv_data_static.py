"""
Class Features

Name:          drv_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230822'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_info_args import logger_name
from lib_data_io_geo import read_file_raster
from lib_utils_grid import get_grid_cells, read_grid_file, get_grid_obj

# set logger
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_obj_settings):

        self.alg_flags = alg_obj_settings['flags']
        self.alg_info = alg_obj_settings['info']
        self.alg_cells = alg_obj_settings['cells']
        self.alg_parameters = alg_obj_settings['parameters']
        self.alg_log = alg_obj_settings['log']

        self.alg_obj_src_ref_nrt = alg_obj_settings['datasets']['static']['source'][
            'reference']['cell_datasets_data_nrt']
        self.alg_obj_src_ref_dr = alg_obj_settings['datasets']['static']['source'][
            'reference']['cell_datasets_data_dr']
        self.alg_obj_src_other_dr = alg_obj_settings['datasets']['static']['source'][
            'other']['cell_datasets_data_dr']
        self.alg_obj_dst = alg_obj_settings['datasets']['static']['destination']

        self.tag_folder_name, self.tag_file_name = 'folder_name', 'file_name'

        self.folder_name_cell_src_ref_nrt = self.alg_obj_src_ref_nrt[self.tag_folder_name]
        self.file_name_cell_src_ref_nrt = self.alg_obj_src_ref_nrt[self.tag_file_name]
        self.file_path_cell_src_ref_nrt = os.path.join(
            self.folder_name_cell_src_ref_nrt, self.file_name_cell_src_ref_nrt)

        self.folder_name_cell_src_ref_dr = self.alg_obj_src_ref_dr[self.tag_folder_name]
        self.file_name_cell_src_ref_dr = self.alg_obj_src_ref_dr[self.tag_file_name]
        self.file_path_cell_src_ref_dr = os.path.join(
            self.folder_name_cell_src_ref_dr, self.file_name_cell_src_ref_dr)

        self.folder_name_cell_src_other_dr = self.alg_obj_src_other_dr[self.tag_folder_name]
        self.file_name_cell_src_other_dr = self.alg_obj_src_other_dr[self.tag_file_name]
        self.file_path_cell_src_other_dr = os.path.join(
            self.folder_name_cell_src_other_dr, self.file_name_cell_src_other_dr)

        self.folder_name_cell_dst = self.alg_obj_dst[self.tag_folder_name]
        self.file_name_cell_dst = self.alg_obj_dst[self.tag_file_name]
        self.file_path_cell_dst = os.path.join(self.folder_name_cell_dst, self.file_name_cell_dst)

        self.cell_start, self.cell_end = self.alg_cells['cell_start'], self.alg_cells['cell_end']
        self.cell_list = self.alg_cells['cell_list']

    # method to wrap file grid
    def wrap_file_grid(self, file_path):

        # check file availability
        if os.path.exists(file_path):
            # read grid file
            obj_data = read_grid_file(file_path)
            # get grid object
            grid_data = get_grid_obj(obj_data)
            # read reference product grid and create grid cells and gpis
            cell_data, gpis_data = get_grid_cells(
                cell_start=self.cell_start, cell_end=self.cell_end, cells_list=self.cell_list,
                cell_grid_reference=grid_data)
        else:
            # message error
            alg_logger.error(' ===> File grid ' + file_path + ' not found')
            raise IOError('File grid not found. Check your path.')

        return obj_data, grid_data, cell_data, gpis_data

    # method to wrap file raster
    def wrap_file_raster(self, file_path):

        # check file availability
        if os.path.exists(file_path):
            # read grid file
            obj_data = read_file_raster(file_path)
        else:
            # message error
            alg_logger.error(' ===> File raster ' + file_path + ' not found')
            raise IOError('File raster not found. Check your path.')

        return obj_data

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ----> Organize static datasets ... ')

        # get grid cell source reference datasets nrt
        (obj_cell_src_ref_nrt, grid_cell_src_ref_nrt,
         list_cell_src_ref_nrt, list_gpis_src_ref_nrt) = self.wrap_file_grid(self.file_path_cell_src_ref_nrt)
        # get grid cell source reference datasets dr
        (obj_cell_src_ref_dr, grid_cell_src_ref_dr,
         list_cell_src_ref_dr, list_gpis_src_ref_dr) = self.wrap_file_grid(self.file_path_cell_src_ref_dr)
        # get grid cell source other datasets dr
        (obj_cell_src_other_dr, grid_cell_src_other_dr,
         list_cell_src_other_dr, list_gpis_src_other_dr) = self.wrap_file_grid(self.file_path_cell_src_other_dr)

        # get grid destination datasets
        obj_grid_dst = self.wrap_file_raster(self.file_path_cell_dst)

        # organize obj collections
        obj_collections = {
            'grid_src_ref_nrt': grid_cell_src_ref_nrt, 'grid_src_ref_dr': grid_cell_src_ref_dr,
            'grid_src_other_dr': grid_cell_src_other_dr,
            'cell_src_ref_nrt': list_cell_src_ref_nrt, 'gpis_src_ref_nrt': list_gpis_src_ref_nrt,
            'cell_src_ref_dr': list_cell_src_ref_dr, 'gpis_src_ref_dr': list_gpis_src_ref_dr,
            'cell_src_other_dr': list_cell_src_other_dr, 'gpis_src_other_dr': list_gpis_src_other_dr,
            'grid_dst': obj_grid_dst}

        # info end method
        alg_logger.info(' ----> Organize static datasets ... DONE')

        return obj_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
