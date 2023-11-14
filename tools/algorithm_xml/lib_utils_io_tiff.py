"""
Library Features:

Name:          lib_utils_io_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import time
import numpy as np
import pandas as pd

from copy import deepcopy
from osgeo import gdal, gdalconst, osr
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file tiff
def read_file_tiff(file_name, file_band_default='band_{}', file_proj='EPSG:4326'):

    src = gdal.Open(file_name)
    data_array = src.ReadAsArray()
    width, height, count = src.RasterXSize, src.RasterYSize, src.RasterCount
    geo_transform = src.GetGeoTransform()
    proj = src.GetProjectionRef()

    datetime = src.GetMetadataItem("TIFFTAG_DATETIME")

    if datetime is None:
        datetime = pd.Timestamp(time.ctime(os.path.getmtime(file_name)))
    else:
        datetime = pd.Timestamp(datetime)

    bound_left = geo_transform[0]
    bound_bottom = geo_transform[3] + width * geo_transform[4] + height * geo_transform[5]
    bound_right = geo_transform[0] + width * geo_transform[1] + height * geo_transform[2]
    bound_top = geo_transform[3]
    origin_x, origin_y = geo_transform[0], geo_transform[3]
    res_x, res_y = np.abs(geo_transform[1]), np.abs(geo_transform[5])

    srs = osr.SpatialReference(wkt=proj)
    if srs is None:
        crs = '' # CRS.from_string(file_proj)
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
                    'width': width, 'height': height, 'datetime': datetime}

    data_obj, data_attrs = None, None
    if count == 1:

        band_name_select = file_band_default.format(str(1))

        if data_array.ndim == 2:
            data_tmp = deepcopy(data_array)
        elif data_array.ndim == 3:
            data_tmp = data_array[0, :, :]
        else:
            logging.error(' ===> Dimensions are not expected for the datasets')
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
        logging.error(' ===> File format is not supported')
        raise NotImplementedError('Unknown error related to the file format')

    if data_obj is None:
        logging.error(' ===> Data object is defined by NoneType. All variable(s) are not included')
        raise RuntimeError('One or more variables must be included in the object to correctly run the algorithm')

    return data_obj, data_attrs, data_lons, data_lats, common_attrs
# ----------------------------------------------------------------------------------------------------------------------
