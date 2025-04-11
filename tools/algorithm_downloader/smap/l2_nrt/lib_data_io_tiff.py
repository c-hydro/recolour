"""
Library Features:

Name:          lib_data_io_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20231110'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging

import os

import numpy as np
from osgeo import gdal, gdalconst

from lib_info_args import logger_name

# logger stream
logger_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read tiff file
def read_file_tiff(file_name_tiff):
    if os.path.exists(file_name_tiff):
        dset_tiff = gdal.Open(file_name_tiff, gdalconst.GA_ReadOnly)
        dset_proj = dset_tiff.GetProjection()
        dset_geotrans = dset_tiff.GetGeoTransform()
        dset_data = dset_tiff.ReadAsArray()
        dset_band = dset_tiff.GetRasterBand(1)
        dset_wide = dset_tiff.RasterXSize
        dset_high = dset_tiff.RasterYSize
    else:
        logger_stream.error(' ===> Tiff file ' + file_name_tiff + ' not found')
        raise IOError('Tiff file location or name is wrong')

    return dset_tiff, dset_proj, dset_geotrans, dset_data,\
        dset_wide, dset_high
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to write file tiff
def write_file_tiff(file_name, file_data, file_metadata, file_wide, file_high, file_geotrans, file_proj, flip=False):

    if not isinstance(file_data, list):
        file_data = [file_data]
    if flip: file_data = [np.flipud(data) for data in file_data]

    if file_metadata is None:
        file_metadata = {'description': 'data'}
    if not isinstance(file_metadata, list):
        file_metadata = [file_metadata] * file_data.__len__()

    layer_n = file_data.__len__()
    dset_handler = gdal.GetDriverByName('GTiff').Create(file_name, file_wide, file_high, layer_n,
                                                        gdalconst.GDT_Float32, options=['COMPRESS=DEFLATE'])
    dset_handler.SetGeoTransform(file_geotrans)
    dset_handler.SetProjection(file_proj)

    for layer_id, (layer_data, layer_metadata) in enumerate(zip(file_data, file_metadata)):
        dset_handler.GetRasterBand(layer_id + 1).SetNoDataValue(-9999.)
        dset_handler.GetRasterBand(layer_id + 1).WriteArray(layer_data)
        dset_handler.GetRasterBand(layer_id + 1).SetMetadata(layer_metadata)
    del dset_handler
# -------------------------------------------------------------------------------------