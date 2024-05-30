"""
Class Features

Name:          drv_cell_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230822'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_info_args import logger_name
from lib_utils_grid import get_grid_cells, read_grid_file, get_grid_obj

# set logger
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver cell
class DrvCell:

    # method to initialize class
    def __init__(self, alg_obj_settings):

        self.alg_flags = alg_obj_settings['flags']
        self.alg_info = alg_obj_settings['info']
        self.alg_cells = alg_obj_settings['cells']
        self.alg_parameters = alg_obj_settings['parameters']
        self.alg_log = alg_obj_settings['log']

        self.alg_cells_data = alg_obj_settings['datasets']['static']['cell_datasets_data']
        self.alg_cells_ref = alg_obj_settings['datasets']['static']['cell_datasets_ref']

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'

        self.folder_name_cell_data = self.alg_cells_data[self.tag_folder_name]
        self.file_name_cell_data = self.alg_cells_data[self.tag_file_name]
        self.file_path_cell_data = os.path.join(self.folder_name_cell_data, self.file_name_cell_data)

        self.folder_name_cell_ref = self.alg_cells_ref[self.tag_folder_name]
        self.file_name_cell_ref = self.alg_cells_ref[self.tag_file_name]
        self.file_path_cell_ref = os.path.join(self.folder_name_cell_ref, self.file_name_cell_ref)

        self.cell_start, self.cell_end = self.alg_cells['cell_start'], self.alg_cells['cell_end']
        self.cell_list = self.alg_cells['cell_list']

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ----> Organize static datasets ... ')

        # get grid cell data datasets
        obj_cell_data = read_grid_file(self.file_path_cell_data)
        grid_cell_data = get_grid_obj(obj_cell_data)
        # get grid cell reference datasets
        obj_cell_ref = read_grid_file(self.file_path_cell_ref)
        grid_cell_ref = get_grid_obj(obj_cell_ref)

        # read data product grid and create grid cells and gpis
        list_cell_data, list_gpis_data = get_grid_cells(
            cell_start=self.cell_start, cell_end=self.cell_end, cells_list=self.cell_list,
            cell_grid_reference=grid_cell_data)

        # read reference product grid and create grid cells and gpis
        list_cell_ref, list_gpis_ref = get_grid_cells(
            cell_start=self.cell_start, cell_end=self.cell_end, cells_list=self.cell_list,
            cell_grid_reference=grid_cell_ref)

        # organize obj collections
        obj_collections = {
            'grid_cell_data': grid_cell_data, 'grid_cell_ref': grid_cell_ref,
            'list_cell_data': list_cell_data, 'list_gpis_data': list_gpis_data,
            'list_cell_ref': list_cell_ref, 'list_gpis_ref': list_gpis_ref}

        # info end method
        alg_logger.info(' ----> Organize static datasets ... DONE')

        return obj_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
