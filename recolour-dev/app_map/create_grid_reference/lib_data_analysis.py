"""
Library Features:

Name:          lib_data_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import shutil

import numpy as np
import xarray as xr

import gdal
from osgeo import gdal, gdalconst

from lib_data_io_tiff import write_file_tiff
from lib_utils_generic import create_filename_tmp, make_folder
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create tmp file name
def create_file_name(folder_name, file_name=None, file_prefix='tmp_geo_', file_suffix='.tiff'):

    if folder_name is not None:
        make_folder(folder_name)
    else:
        logging.error(' ===> Temporary folder is defined by NoneType')
        raise RuntimeError('Temporary folder must be defined according to the system environment')

    if (file_name is None) or (file_name == ''):
        file_path = create_filename_tmp(prefix=file_prefix, suffix=file_suffix, folder=folder_name)
    else:
        file_path = os.path.join(folder_name, file_name)

    return file_path
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to dump tmp file data
def dump_file_tiff(file_name, grid_data, grid_attrs,
                   grid_no_data=-9999, grid_format='integer'):

    if not isinstance(grid_data, list):
        grid_data = [grid_data]

    grid_ws = []
    for tmp_data in grid_data:
        if isinstance(tmp_data, xr.DataArray):
            grid_values = tmp_data.values
        else:
            logging.error(' ===> Data type is not supported')
            raise NotImplemented('Case not implemented yet')

        if grid_format == 'float':
            grid_values[grid_values == grid_no_data] = np.nan

        grid_ws.append(grid_values)

    if grid_format == 'integer':
        file_format = gdalconst.GDT_Int32
    elif grid_format == 'float':
        file_format = gdalconst.GDT_Float32
    else:
        logging.error(' ===> Grid format "' + grid_format + '" is not supported')
        raise NotImplemented('Case not implemented yet')

    write_file_tiff(file_name, file_data=grid_ws, file_format=file_format, **grid_attrs)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to execute mask
def compute_grid_mask(grid_da, grid_mask_min=0, grid_mask_max=None, grid_mask_data=1, grid_mask_format='integer',
                      grid_mask_no_data=-9999.0, grid_mask_active=True):
    if grid_mask_active:
        if grid_mask_min is not None:
            grid_da = xr.where(grid_da >= grid_mask_min, grid_mask_data, grid_mask_no_data)
        if grid_mask_max is not None:
            grid_da = xr.where(grid_da <= grid_mask_max, grid_mask_data, grid_mask_no_data)

        if grid_mask_format == 'integer':
            if grid_mask_no_data == np.nan:
                logging.error(' ===> Grid no_data equal to NaN is not supported in integer datasets')
                raise RuntimeError('Change the no_data value or the grid format in the settings file ')

        if grid_mask_format is not None:
            if grid_mask_format == 'integer':
                grid_da = grid_da.astype(int)
            elif grid_mask_format == 'float':
                grid_da = grid_da.astype(float)
            else:
                logging.error(' ===> Mask format "' + grid_mask_format + '" is not supported')
                raise NotImplemented('Case not implemented yet')
    else:
        logging.warning(' ===> Grid masking is not activated')
    return grid_da

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to warp grid data
def warp_grid_data(file_name_src, file_name_dst,
                   grid_warp_method='near', grid_warp_resolution=None,
                   grid_warp_format='GTiff', grid_warp_active=True, grid_warp_copy=True):

    if grid_warp_active:
        if grid_warp_resolution is not None:
            grid_warp_options = {
                'format': grid_warp_format,
                'xRes': grid_warp_resolution, 'yRes': grid_warp_resolution,
                'resampleAlg': grid_warp_method}

            file_ds = gdal.Warp(file_name_dst, file_name_src, **grid_warp_options)
            del file_ds
        else:
            logging.error(' ===> Warping grid is not possible due to the resolution defined by NoneType')
            raise RuntimeError('Set the resolution in the settings file')
    else:
        logging.warning(' ===> Grid warping is not activated')
        if grid_warp_copy:
            shutil.copy(file_name_src, file_name_dst)

# ----------------------------------------------------------------------------------------------------------------------
