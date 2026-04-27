"""
Class Features

Name:          drv_map_static
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
from lib_utils_grid import get_grid_cells, read_grid_file

# set logger
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver map
class DrvMap:

    # method to initialize class
    def __init__(self, alg_settings,
                 tag_section_flags='flags', tag_section_info='info', tag_section_cells='cells',
                 tag_section_grid_product='grid_product',
                 tag_section_grid_reference='grid_reference',
                 tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_log='log'):

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_info = alg_settings[tag_section_info]
        self.alg_cells = alg_settings[tag_section_cells]
        self.alg_parameters = alg_settings[tag_section_params]
        self.alg_datasets = alg_settings[tag_section_datasets]['static']

        self.alg_grid_product = self.alg_datasets[tag_section_grid_product]
        self.alg_grid_reference = self.alg_datasets[tag_section_grid_reference]

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'

        self.folder_name_grid_product = self.alg_grid_product[self.tag_folder_name]
        self.file_name_grid_product = self.alg_grid_product[self.tag_file_name]
        self.file_path_grid_product = os.path.join(self.folder_name_grid_product, self.file_name_grid_product)

        self.folder_name_grid_reference = self.alg_grid_reference[self.tag_folder_name]
        self.file_name_grid_reference = self.alg_grid_reference[self.tag_file_name]
        self.file_path_grid_reference = os.path.join(self.folder_name_grid_reference, self.file_name_grid_reference)

        self.cell_start, self.cell_end = self.alg_cells['cell_start'], self.alg_cells['cell_end']
        self.cell_list = self.alg_cells['cell_list']

    # method to organize data
    def organize_data(self):

        # info start method
        alg_logger.info(' ----> Organize static datasets ... ')

        # grid reference product
        grid_reference_product = read_grid_file(self.file_path_grid_product)
        # read reference grid
        grid_reference_domain = read_file_raster(self.file_path_grid_reference)

        # read product grid and create grid cells and gpis
        grid_product_cells, grid_product_gpis = get_grid_cells(
            cell_start=self.cell_start, cell_end=self.cell_end, cells_list=self.cell_list,
            path_grid=self.folder_name_grid_product, file_grid=self.file_name_grid_product)

        # organize grid obj
        grid_obj = {'grid_cells_product': grid_product_cells, 'grid_gpis_product': grid_product_gpis,
                    'grid_reference_product': grid_reference_product,
                    'grid_reference_domain': grid_reference_domain}

        # info end method
        alg_logger.info(' ----> Organize static datasets ... DONE')

        return grid_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
