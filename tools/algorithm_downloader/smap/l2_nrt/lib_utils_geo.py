"""
Library Features:

Name:          lib_utils_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20250410'
Version:       '1.1.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pyresample
import numpy as np

from lib_data_io_geo import read_grid_data
from lib_info_args import logger_name

# logger stream
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to set geographical information
def set_geo_info(file_name: str):

    if os.path.exists(file_name):

        grid_data, grid_attrs = read_grid_data(file_name)
        geo_proj_obj, geo_geotrans_obj = grid_attrs['crs'], grid_attrs['transform']
        geo_geotrans = (geo_geotrans_obj.xoff, geo_geotrans_obj.a, 0, geo_geotrans_obj.yoff, 0, geo_geotrans_obj.e)
        geo_proj_wkt = geo_proj_obj.wkt

        grid_lons_1d, grid_lats_1d = grid_data['longitude'].values, np.flipud(grid_data['latitude'].values)

        # define grid mesh
        grid_lons_2d, grid_lats_2d = np.meshgrid(grid_lons_1d, grid_lats_1d)
        # define grid data 2d
        grid_sm_2d = np.zeros(shape=(grid_lons_2d.shape[0], grid_lats_2d.shape[1]))
        grid_sm_2d[:, :] = np.nan

        # grid obj using pyresample definition
        grid_obj = pyresample.geometry.GridDefinition(lats=grid_lats_2d, lons=grid_lons_2d)

    else:
        logger_stream.error(' ===> Geographical file ' + file_name + ' not found')
        raise IOError('Geographical file location or name is wrong')

    return grid_obj, geo_proj_wkt, geo_geotrans, grid_sm_2d
# ----------------------------------------------------------------------------------------------------------------------