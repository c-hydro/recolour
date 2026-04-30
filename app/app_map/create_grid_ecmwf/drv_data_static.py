"""
Class Features

Name:          drv_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260428'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_io_geo import read_grid_data, read_grid_netcdf,  check_grid_data
from lib_info_args import geo_coord_name_x, geo_coord_name_y, geo_dim_name_x, geo_dim_name_y

# set debug level
logging.getLogger("findlibs").setLevel(logging.WARNING)
logging.getLogger("gribapi").setLevel(logging.WARNING)
logging.getLogger("eccodes").setLevel(logging.WARNING)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_settings,
                 tag_section_flags='flags',
                 tag_section_grid='grid_reference',
                 tag_section_datasets='datasets',
                 tag_section_log='log'):

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_datasets = alg_settings[tag_section_datasets]['static']

        self.alg_grid = self.alg_datasets[tag_section_grid]

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_format = 'format'

        self.folder_name_grid = self.alg_grid[self.tag_folder_name]
        self.file_name_grid = self.alg_grid[self.tag_file_name]
        self.file_path_grid = os.path.join(self.folder_name_grid, self.file_name_grid)

        self.format = None
        if self.tag_format in list(self.alg_grid.keys()):
            self.format = self.alg_grid[self.tag_format]
        else:
            self.format = "tiff"

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize static datasets ... ')

        # read grid data
        if self.format in ['tiff', 'tif', 'asc', 'txt']:

            # read data ascii or tiff
            grid_data, grid_attrs = read_grid_data(
                self.file_path_grid,
                coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y,
                dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y)

        elif self.format == "netcdf":

            # read data netcdf (create grid reference application for ecmwf datasets)
            grid_data, grid_attrs = read_grid_netcdf(
                 self.file_path_grid,
                var_name="mask", output_format="data_array", output_dtype="float32",
            )

        else:

            logging.error(' ===> Unknown format: grid data must be in tiff, netcdf')
            raise NotImplementedError('Case not implemented yet')

        # check grid data
        check_grid_data(self.file_path_grid, grid_data, grid_attrs, data_mandatory=True)

        # organize grid obj
        grid_obj = {'grid_data': grid_data, 'grid_attrs': grid_attrs}

        # info end method
        logging.info(' ---> Organize static datasets ... DONE')

        return grid_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
