"""
Library Features:

Name:          lib_data_io_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210225'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import rasterio

import xarray as xr
import pandas as pd
import numpy as np

from copy import deepcopy
from rasterio.transform import Affine
from rasterio.crs import CRS
from osgeo import gdal, gdalconst

from lib_info_args import logger_name
from lib_info_args import proj_epsg, time_format_datasets

# set logger level
logging.getLogger('rasterio').setLevel(logging.WARNING)
# set logger obj
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize file tiff
def organize_file_tiff(obj_data, obj_time=None, obj_variable_in=None, obj_variable_out=None,
                       obj_transform=None, obj_proj=None,
                       var_attr_description='description', var_attr_time='time',
                       var_name_geo_x='longitude', var_name_geo_y='latitude',
                       ):

    if isinstance(obj_time, pd.Timestamp):
        string_time = obj_time.strftime(time_format_datasets)
    elif isinstance(obj_time, str):
        string_time = deepcopy(obj_time)
    else:
        alg_logger.error(' ===> Time obj format is not supported')
        raise NotImplemented('Case not implemented yet')

    if obj_variable_in is None:
        obj_variable_in = []
    if obj_variable_out is None:
        obj_variable_out = []

    var_height, var_width = None, None
    var_geo_x_west, var_geo_x_east, var_geo_y_south, var_geo_y_north = None, None, None, None
    var_geo_transform, var_geo_proj = None, None
    var_data_list, var_metadata_list = [], []

    # iterate over variable(s)
    for var_key_obj, var_key_file in zip(obj_variable_in, obj_variable_out):

        # check variable availability
        if var_key_obj in list(obj_data.keys()):
            var_obj = obj_data[var_key_obj]

            # check variable data
            if var_obj is not None:

                if (var_height is None) or (var_width is None):
                    var_height, var_width = var_obj.shape

                if obj_transform is None:
                    if (var_geo_x_west is None) or (var_geo_x_east is None):
                        var_geo_x = var_obj[var_name_geo_x].values
                        var_geo_x_west = np.min(np.min(var_geo_x))
                        var_geo_x_east = np.max(np.max(var_geo_x))
                    if (var_geo_y_south is None) or (var_geo_y_north is None):
                        var_geo_y = var_obj[var_name_geo_y].values
                        var_geo_y_south = np.min(np.min(var_geo_y))
                        var_geo_y_north = np.max(np.max(var_geo_y))
                    if var_geo_transform is None:
                        # TO DO: fix the 1/2 pixel of resolution in x and y ... using resolution/2
                        var_geo_transform = rasterio.transform.from_bounds(
                            var_geo_x_west, var_geo_y_south, var_geo_x_east, var_geo_y_north,
                            var_width, var_height)
                else:
                    var_geo_transform = deepcopy(obj_transform)

                if obj_proj is None:
                    var_geo_proj = deepcopy(proj_epsg)
                else:
                    var_geo_proj = deepcopy(obj_proj)

                if not isinstance(var_geo_proj, str):
                    var_geo_proj = var_geo_proj.to_string()

                if isinstance(var_obj, xr.DataArray):
                    var_data_step = var_obj.values
                else:
                    var_data_step = deepcopy(var_obj)
                var_metadata_step = {var_attr_description: var_key_file, var_attr_time: string_time}

                var_data_list.append(var_data_step)
                var_metadata_list.append(var_metadata_step)

            else:
                alg_logger.warning(' ===> Variable "' + var_key_obj + '" is defined by NoneType in the data object')

        else:
            alg_logger.warning(' ===> Variable "' + var_key_obj + '" is not available in the data object')

    var_attributes = {
        "file_wide": var_width, 'file_high': var_height,
        'file_transform': var_geo_transform, 'file_proj': var_geo_proj,
        'file_metadata': var_metadata_list
    }

    return var_data_list, var_attributes

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to write file tiff
def write_file_tiff(file_name, file_data, file_wide, file_high, file_transform, file_proj,
                    file_metadata=None, file_format=gdalconst.GDT_Float32):

    if not isinstance(file_data, list):
        file_data = [file_data]

    if file_metadata is None:
        file_metadata = {'description_field': 'data'}
    if not isinstance(file_metadata, list):
        file_metadata = [file_metadata] * file_data.__len__()

    if isinstance(file_transform, Affine):
        file_transform = file_transform.to_gdal()

    file_crs = rasterio.crs.CRS.from_string(file_proj)
    file_wkt = file_crs.to_wkt()

    file_n = file_data.__len__()
    dset_handle = gdal.GetDriverByName('GTiff').Create(file_name, file_wide, file_high, file_n, file_format,
                                                       options=['COMPRESS=DEFLATE'])
    dset_handle.SetGeoTransform(file_transform)
    dset_handle.SetProjection(file_wkt)

    for file_id, (file_data_step, file_metadata_step) in enumerate(zip(file_data, file_metadata)):
        dset_handle.GetRasterBand(file_id + 1).WriteArray(file_data_step)
        dset_handle.GetRasterBand(file_id + 1).SetMetadata(file_metadata_step)
    del dset_handle
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to read file tiff
def read_file_tiff(file_name):

    file_handle = rasterio.open(file_name)
    file_bounds, file_res = file_handle.bounds, file_handle.res
    file_proj = file_handle.crs.wkt
    file_geotrans = file_handle.transform

    file_tags = file_handle.tags()
    file_bands = file_handle.count
    file_metadata = file_handle.profile

    if file_bands == 1:
        file_data = file_handle.read(1)
    elif file_bands > 1:
        file_data = []
        for band_id in range(0, file_bands):
            file_data_tmp = file_handle.read(band_id + 1)
            file_data.append(file_data_tmp)
    else:
        alg_logger.error(' ===> File multi-band are not supported')
        raise NotImplementedError('File multi-band not implemented yet')

    return file_data, file_proj, file_geotrans
# ----------------------------------------------------------------------------------------------------------------------
