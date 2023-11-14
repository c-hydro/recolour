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

import pandas as pd
import numpy as np

from copy import deepcopy
from rasterio.transform import Affine
from rasterio.crs import CRS
from osgeo import gdal, gdalconst, osr

from lib_info_args import logger_name
from lib_info_args import proj_epsg, time_format_datasets

# set logger level
logging.getLogger('rasterio').setLevel(logging.WARNING)
# set logger obj
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize file tiff
def organize_file_tiff(obj_data, obj_time=None, obj_variable=None, obj_transform=None, obj_proj=None,
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

    var_list_file, var_list_obj = list(obj_variable.keys()), list(obj_variable.values())

    var_height, var_width = None, None
    var_geo_x_west, var_geo_x_east, var_geo_y_south, var_geo_y_north = None, None, None, None
    var_geo_transform, var_geo_proj = None, None
    var_data_list, var_metadata_list = [], []

    # iterate over variable(s)
    for var_key_obj, var_key_file in zip(var_list_obj, var_list_file):

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

                var_data_step = var_obj.values
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
# method to write file tiff
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
# method to read file tiff
def read_file_tiff(file_name, file_band_default='band_{}', file_proj='EPSG:4326'):

    src = gdal.Open(file_name)
    data_array = src.ReadAsArray()
    width, height, count = src.RasterXSize, src.RasterYSize, src.RasterCount
    geo_transform = src.GetGeoTransform()
    proj = src.GetProjectionRef()

    bound_left = geo_transform[0]
    bound_bottom = geo_transform[3] + width * geo_transform[4] + height * geo_transform[5]
    bound_right = geo_transform[0] + width * geo_transform[1] + height * geo_transform[2]
    bound_top = geo_transform[3]
    origin_x, origin_y = geo_transform[0], geo_transform[3]
    res_x, res_y = np.abs(geo_transform[1]), np.abs(geo_transform[5])

    srs = osr.SpatialReference(wkt=proj)
    if srs is None:
        crs = CRS.from_string(file_proj)
    else:
        crs = '{:}:{:}'.format(srs.GetAttrValue('AUTHORITY', 0), srs.GetAttrValue('AUTHORITY', 1))

    center_right = bound_right - (res_x / 2)
    center_left = bound_left + (res_x / 2)
    center_top = bound_top - (res_y / 2)
    center_bottom = bound_bottom + (res_y / 2)

    lon = np.arange(center_left, center_right + np.abs(res_x / 2), np.abs(res_x), float)
    lat = np.flip(np.arange(center_bottom, center_top + np.abs(res_y / 2), np.abs(res_y), float), axis=0)
    data_lons, data_lats = np.meshgrid(lon, lat)

    lat_flag, lat_upper, lat_lower = False, data_lats[0, 0], data_lats[-1, 0]
    if lat_lower > lat_upper:
        data_lats = np.flipud(data_lats)
        lat_flag = True

    common_attrs = {'transform': geo_transform, 'crs': crs,
                    'bbox': [bound_left, bound_bottom, bound_right, bound_top],
                    'bb_left': bound_left, 'bb_right': bound_right,
                    'bb_top': bound_top, 'bb_bottom': bound_bottom,
                    'res_lon': res_x, 'res_lat': res_y,
                    'width': width, 'height': height}

    data_obj, data_attrs = None, None
    if count == 1:

        band_name_select = file_band_default.format(str(1))

        if data_array.ndim == 2:
            data_tmp = deepcopy(data_array)
        elif data_array.ndim == 3:
            data_tmp = data_array[0, :, :]
        else:
            alg_logger.error(' ===> Dimensions are not expected for the datasets')
            raise NotImplemented('Case not implemented yet')

        if lat_flag:
            data_tmp = np.flipud(data_tmp)

        band_obj = src.GetRasterBand(1)
        metadata_obj = band_obj.GetMetadata()

        if data_obj is None:
            data_obj = {}
        if data_attrs is None:
            data_attrs = {}

        data_obj[band_name_select] = data_tmp
        data_attrs[band_name_select] = metadata_obj

    elif count > 1:
        for id in range(0, count):

            band_id = id + 1

            data_tmp = data_array[id, :, :]

            band_obj = src.GetRasterBand(band_id)
            metadata_obj = band_obj.GetMetadata()

            band_name_select = file_band_default.format(str(band_id))

            if data_obj is None:
                data_obj = {}
            if data_attrs is None:
                data_attrs = {}

            if lat_flag:
                data_tmp = np.flipud(data_tmp)

            data_obj[band_name_select] = data_tmp
            data_attrs[band_name_select] = metadata_obj

    else:
        alg_logger.error(' ===> File format is not supported')
        raise NotImplementedError('Unknown error related to the file format')

    if data_obj is None:
        alg_logger.error(' ===> Data object is defined by NoneType. All variable(s) are not included')
        raise RuntimeError('One or more variables must be included in the object to correctly run the algorithm')

    return data_obj, data_attrs, data_lons, data_lats, common_attrs
# ----------------------------------------------------------------------------------------------------------------------
