# -------------------------------------------------------------------------------------
# Libraries
import logging
import rasterio

import pandas as pd
import numpy as np

from rasterio.transform import Affine
from rasterio.crs import CRS

from copy import deepcopy
from decimal import Decimal

from osgeo import gdal, gdalconst

logging.getLogger('rasterio').setLevel(logging.WARNING)

# Debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read file tiff
def read_file_tiff(file_name, file_scale_factor=1, file_format=gdalconst.GDT_Float32):

    file_handle = gdal.Open(file_name, file_format)
    file_metadata = file_handle.GetMetadata()
    file_tranform = file_handle.GetGeoTransform()
    file_proj = file_handle.GetProjection()

    file_width, file_height = file_handle.RasterXSize, file_handle.RasterYSize
    file_bands = file_handle.RasterCount

    file_values = np.zeros([file_height, file_width, file_bands])
    file_metadata = []
    if file_bands == 1:
        band_handle = file_handle.GetRasterBand(1)
        band_values = band_handle.ReadAsArray()

        band_values = band_values[:, :, np.newaxis]
        band_values = band_values / file_scale_factor

        file_values[:, :, 0] = band_values[:, :, 0]
    elif file_bands > 1:
        for band_id in range(0, file_bands):

            band_handle = file_handle.GetRasterBand(band_id + 1)

            band_values = band_handle.ReadAsArray()
            band_values = band_values / file_scale_factor

            band_metadata = band_handle.GetMetadata()

            file_values[:, :, band_id] = band_values
            file_metadata.append(band_metadata)

    else:
        logging.error(' ===> File multi-band are not supported')
        raise NotImplementedError('File multi-band not implemented yet')

    return file_values, file_metadata, file_proj, file_tranform
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to write file tiff
def write_file_tiff(file_name, file_data, file_time, file_wide, file_high, file_geotrans, file_epsg,
                     file_metadata=None, file_scale_factor=1,
                     file_format=gdalconst.GDT_Float32):

    file_shape = file_data.shape
    if not isinstance(file_data, list):
        file_data_list = []
        if file_data.shape.__len__() == 3:
            for id_step in range(0, file_shape[2]):
                file_data_list.append(file_data[:, :, id_step])
        elif file_data.shape.__len__() == 2:
            file_data_list.append(file_data)
        else:
            raise NotImplementedError('Write file tiff failed due to unsupported data format')
    else:
        file_data_list = file_data

    if file_metadata is None:
        file_metadata = {'description_field': 'data'}
    if not isinstance(file_metadata, list):
        if file_data.shape.__len__() == 3:
            file_metadata = [file_metadata] * file_shape[2]
        elif file_data.shape.__len__() == 2:
            file_metadata = [file_metadata]
        else:
            raise NotImplementedError('Write attributes tiff failed due to unsupported data format')

    if file_time is not None:
        for time_id, time_step in enumerate(file_time):
            dict_metadata = deepcopy(file_metadata[time_id])
            dict_metadata['time'] = pd.Timestamp(time_step).strftime('%Y-%m-%d %H:%M')
            file_metadata[time_id] = dict_metadata

    if isinstance(file_geotrans, Affine):
        file_geotrans = file_geotrans.to_gdal()

    if isinstance(file_epsg, int):
        file_crs = rasterio.crs.CRS.from_epsg(file_epsg)
    elif isinstance(file_epsg, str):
        file_crs = rasterio.crs.CRS.from_string(file_epsg)
    else:
        raise IOError('Geographical EPSG must be in string format "EPSG:4326" or integer format "4326".')

    file_wkt = file_crs.to_wkt()

    if file_data.shape.__len__() == 3:
        file_n = file_shape[2]
    elif file_data.shape.__len__() == 2:
        file_n = 1

    dset_handle = gdal.GetDriverByName('GTiff').Create(file_name, file_wide, file_high, file_n, file_format,
                                                       options=['COMPRESS=DEFLATE'])
    dset_handle.SetGeoTransform(file_geotrans)
    dset_handle.SetProjection(file_wkt)

    for file_id, (file_data_step, file_metadata_step) in enumerate(zip(file_data_list, file_metadata)):

        file_data_step = file_data_step / file_scale_factor

        dset_handle.GetRasterBand(file_id + 1).WriteArray(file_data_step)
        dset_handle.GetRasterBand(file_id + 1).SetMetadata(file_metadata_step)
    del dset_handle
# -------------------------------------------------------------------------------------
