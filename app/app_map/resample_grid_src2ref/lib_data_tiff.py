"""
Library Features:

Name:          lib_data_io_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240709'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import rasterio
import gdal
import numpy as np

from copy import deepcopy
from rasterio.transform import Affine
from osgeo import gdal, gdalconst

from lib_info_args import logger_name
from lib_info_args import proj_epsg as proj_default_epsg

# logging
logging.getLogger('rasterio').setLevel(logging.WARNING)
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to define file tiff metadata
def organize_file_tiff(file_data, file_geo_x, file_geo_y, file_geo_transform=None, file_geo_epsg=None):

    file_height, file_width = file_data.shape

    file_geo_x_west = np.min(np.min(file_geo_x))
    file_geo_x_east = np.max(np.max(file_geo_x))

    file_geo_y_south = np.min(np.min(file_geo_y))
    file_geo_y_north = np.max(np.max(file_geo_y))

    if file_geo_transform is None:
        # TO DO: fix the 1/2 pixel of resolution in x and y ... using resolution/2
        file_geo_transform = rasterio.transform.from_bounds(
            file_geo_x_west, file_geo_y_south, file_geo_x_east, file_geo_y_north,
            file_width, file_height)

    if file_geo_epsg is None:
        file_geo_epsg = deepcopy(proj_default_epsg)

    if not isinstance(file_geo_epsg, str):
        file_geo_epsg = file_geo_epsg.to_string()

    return file_height, file_width, file_geo_transform, file_geo_epsg
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to write file tiff
def write_file_tiff(file_name, file_data, file_wide, file_high, file_geotrans, file_proj,
                    file_metadata=None, file_format=gdalconst.GDT_Float32):

    if not isinstance(file_data, list):
        file_data = [file_data]

    if file_metadata is None:
        file_metadata = {'description_field': 'data'}
    if not isinstance(file_metadata, list):
        file_metadata = [file_metadata] * file_data.__len__()

    if isinstance(file_geotrans, Affine):
        file_geotrans = file_geotrans.to_gdal()

    file_n = file_data.__len__()
    dset_handle = gdal.GetDriverByName('GTiff').Create(file_name, file_wide, file_high, file_n, file_format,
                                                       options=['COMPRESS=DEFLATE'])
    dset_handle.SetGeoTransform(file_geotrans)
    dset_handle.SetProjection(file_proj)

    for file_id, (file_data_step, file_metadata_step) in enumerate(zip(file_data, file_metadata)):
        dset_handle.GetRasterBand(file_id + 1).WriteArray(file_data_step)
        dset_handle.GetRasterBand(file_id + 1).SetMetadata(file_metadata_step)
    del dset_handle
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to read file tiff
def read_file_tiff(file_name):

    file_handle = rasterio.open(file_name)

    file_proj = file_handle.crs.wkt
    file_bounds, file_res, file_transform = file_handle.bounds, file_handle.res, file_handle.transform
    file_tags = file_handle.tags()
    file_bands = file_handle.count
    file_metadata = file_handle.profile

    decimal_round_geo = 7

    center_right = file_bounds.right - (file_res[0] / 2)
    center_left = file_bounds.left + (file_res[0] / 2)
    center_top = file_bounds.top - (file_res[1] / 2)
    center_bottom = file_bounds.bottom + (file_res[1] / 2)

    if center_bottom > center_top:
        log_stream.warning(' ===> Coords "center_bottom": ' + str(center_bottom) + ' is greater than "center_top": '
                           + str(center_top) + '. Try to inverse the bottom and top coords. ')
        center_tmp = center_top
        center_top = center_bottom
        center_bottom = center_tmp

    lon = np.arange(center_left, center_right + np.abs(file_res[0] / 2), np.abs(file_res[0]), float)
    lat = np.flip(np.arange(center_bottom, center_top + np.abs(file_res[1] / 2), np.abs(file_res[1]), float), axis=0)
    file_lons, file_lats = np.meshgrid(lon, lat)

    lat_flipping = False
    lat_upper, lat_lower = file_lats[0, 0], file_lats[-1, 0]
    if lat_lower > lat_upper:
        file_lats = np.flipud(file_lats)
        lat_flipping = True

    min_lon_round = round(np.min(file_lons), decimal_round_geo)
    max_lon_round = round(np.max(file_lons), decimal_round_geo)
    min_lat_round = round(np.min(file_lats), decimal_round_geo)
    max_lat_round = round(np.max(file_lats), decimal_round_geo)

    center_right_round = round(center_right, decimal_round_geo)
    center_left_round = round(center_left, decimal_round_geo)
    center_bottom_round = round(center_bottom, decimal_round_geo)
    center_top_round = round(center_top, decimal_round_geo)

    assert min_lon_round == center_left_round
    assert max_lon_round == center_right_round
    assert min_lat_round == center_bottom_round
    assert max_lat_round == center_top_round

    file_height, file_width = None, None
    if file_bands == 1:
        file_data = file_handle.read(1)

        if lat_flipping:
            file_data = np.flipud(file_data)

        if file_height is None or file_width is None:
            file_height, file_width = file_data.shape

    elif file_bands > 1:
        file_data = []
        for band_id in range(0, file_bands):
            file_data_tmp = file_handle.read(band_id + 1)

            if lat_flipping:
                file_data_tmp = np.flipud(file_data_tmp)
            if file_height is None or file_width is None:
                file_height, file_width = file_data_tmp.shape

            file_data.append(file_data_tmp)
    else:
        log_stream.error(' ===> File multi-band are not supported')
        raise NotImplementedError('File multi-band not implemented yet')

    return file_data, file_lons, file_lats, file_height, file_width, file_proj, file_transform
# ----------------------------------------------------------------------------------------------------------------------
