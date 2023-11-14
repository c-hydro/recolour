"""
Class Features

Name:          drv_map_data
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_io_geo import read_file_raster
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver map
class DrvMap:

    # method to initialize class
    def __init__(self, alg_settings,
                 tag_section_flags='flags',
                 tag_section_grid_reference='grid_reference',
                 tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_log='log'):

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_parameters = alg_settings[tag_section_params]
        self.alg_datasets = alg_settings[tag_section_datasets]['static']

        self.alg_grid_reference = self.alg_datasets[tag_section_grid_reference]

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'

        self.folder_name_grid_reference = self.alg_grid_reference[self.tag_folder_name]
        self.file_name_grid_reference = self.alg_grid_reference[self.tag_file_name]
        self.file_path_grid_reference = os.path.join(self.folder_name_grid_reference, self.file_name_grid_reference)

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize static datasets ... ')

        # read reference grid
        grid_reference_domain = read_file_raster(self.file_path_grid_reference)

        # organize grid obj
        grid_obj = {'grid_reference_domain': grid_reference_domain}

        # info end method
        logging.info(' ---> Organize static datasets ... DONE')

        return grid_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
