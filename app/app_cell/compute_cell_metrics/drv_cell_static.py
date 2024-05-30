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
    def __init__(self, alg_settings,
                 tag_section_flags='flags', tag_section_info='info', tag_section_cells='cells',
                 tag_section_grid_product='grid_product',
                 tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_log='log'):

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_info = alg_settings[tag_section_info]
        self.alg_cells = alg_settings[tag_section_cells]
        self.alg_parameters = alg_settings[tag_section_params]
        self.alg_datasets = alg_settings[tag_section_datasets]['static']

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'

        self.folder_name_datasets = self.alg_datasets[self.tag_folder_name]
        self.file_name_datasets = self.alg_datasets[self.tag_file_name]
        self.file_path_datasets = os.path.join(self.folder_name_datasets, self.file_name_datasets)

        self.cell_start, self.cell_end = self.alg_cells['cell_start'], self.alg_cells['cell_end']
        self.cell_list = self.alg_cells['cell_list']

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ----> Organize static datasets ... ')

        # grid reference product
        grid_reference = read_grid_file(self.file_path_datasets)

        # read product grid and create grid cells and gpis
        grid_cells, grid_gpis = get_grid_cells(
            cell_start=self.cell_start, cell_end=self.cell_end, cells_list=self.cell_list,
            path_grid=self.folder_name_datasets, file_grid=self.file_name_datasets)

        # method to get grid object
        grid_obj = get_grid_obj(grid_reference)

        # organize grid collections
        grid_collections = {
            'grid_cells': grid_cells, 'grid_gpis': grid_gpis,
            'grid_reference': grid_reference, 'grid_obj': grid_obj}

        # info end method
        alg_logger.info(' ----> Organize static datasets ... DONE')

        return grid_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
