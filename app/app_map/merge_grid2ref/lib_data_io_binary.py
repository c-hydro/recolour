"""
Library Features:

Name:          lib_data_io_binary
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210408'
Version:       '1.0.0'
"""
#######################################################################################
# Library
import logging
import os
import struct

from copy import deepcopy

import numpy as np
import pandas as pd
import xarray as xr

from lib_info_args import logger_name
from lib_info_args import proj_epsg, proj_wkt

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
import matplotlib.pylab as plt
#######################################################################################


# --------------------------------------------------------------------------------
# Method to read 2d variable in binary format (saved as 1d integer array)
def read_data_binary(file_name, var_geo_x, var_geo_y, var_geo_attrs=None,
                     var_format='i', var_scale_factor=10,
                     var_name=None, var_time=None, var_geo_1d=True, var_time_freq='H', var_time_steps_expected=1,
                     coord_name_geo_x='west_east', coord_name_geo_y='south_north', coord_name_time='time',
                     dim_name_geo_x='west_east', dim_name_geo_y='south_north', dim_name_time='time',
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_geo_y, dim_name_geo_x, dim_name_time]

    if os.path.exists(file_name):

        if isinstance(var_name, str):
            pass
        elif isinstance(var_name, list) and var_name.__len__() == 1:
            var_name = var_name[0]
        elif var_name is None:
            var_name = 'data'
        else:
            log_stream.error(' ===> The arguments "var_name" must be a string or a list with length equal to 1')
            raise NotImplemented('Case not implemented yet')

        if isinstance(var_scale_factor, list) and var_name.__len__() == 1:
            var_scale_factor = var_scale_factor[0]
        else:
            log_stream.error(' ===> The arguments "var_scale_factor" must be a scalar or list with length equal to 1')
            raise NotImplemented('Case not implemented yet')

        # Shape values 1d
        rows = var_geo_y.shape[0]
        cols = var_geo_x.shape[0]
        geo_n = rows * cols

        # Open and read binary file [OLD]
        file_handle = open(file_name, 'rb')
        data_format = var_format * geo_n
        data_stream = file_handle.read(-1)
        var_data_1d = struct.unpack(data_format, data_stream)
        file_handle.close()

        file_attrs = {}

        var_data_1d = np.asarray(var_data_1d, dtype=np.float32)
        var_data_1d = np.float32(var_data_1d / var_scale_factor)
        var_n = var_data_1d.shape[0]

        var_time_steps_cmp = int(var_n / geo_n)
        var_data_3d = np.reshape(var_data_1d, (rows, cols, var_time_steps_cmp), order='F')

        if var_geo_1d:
            var_geo_x_2d, var_geo_y_2d = np.meshgrid(var_geo_x, var_geo_y)
        else:
            var_geo_x_2d = var_geo_x
            var_geo_y_2d = var_geo_y

        geo_y_upper = var_geo_y_2d[0, 0]
        geo_y_lower = var_geo_y_2d[-1, 0]
        if geo_y_lower > geo_y_upper:
            var_geo_y_2d = np.flipud(var_geo_y_2d)

        var_dims = var_data_3d.shape
        var_high = var_dims[0]
        var_wide = var_dims[1]

        if var_geo_attrs is not None:

            nodata_value = var_geo_attrs.get('nodata_value', -9999)
            xll_corner = var_geo_attrs.get('xllcorner', None)
            yll_corner = var_geo_attrs.get('yllcorner', None)
            proj = var_geo_attrs.get('proj', proj_epsg)
            transform = var_geo_attrs.get('transform', proj_wkt)
            cellsize = var_geo_attrs.get('cellsize', None)

            geo_attrs = {'nrows': var_geo_y_2d.shape[0], 'ncols': var_geo_x_2d.shape[1],
                         'nodata_value': nodata_value,
                         'xllcorner': xll_corner,
                         'yllcorner': yll_corner, 'cellsize': abs(cellsize),
                         'proj': proj, 'transform': transform}
        else:
            geo_attrs = {}

        if var_time_steps_cmp == var_time_steps_expected:

            var_data = np.zeros(shape=[var_geo_x_2d.shape[0], var_geo_y_2d.shape[1], var_time_steps_cmp])
            var_data[:, :, :] = np.nan
            for step in np.arange(0, var_time_steps_cmp, 1):
                var_data_step = var_data_3d[:, :, step]
                var_data[:, :, step] = var_data_step

        elif (var_time_steps_cmp == 1) and (var_time_steps_cmp < var_time_steps_expected):

            log_stream.warning(' ===> File ' + file_name +
                               ' steps expected [' + str(var_time_steps_expected) +
                               '] and found [' + str(var_time_steps_cmp) + '] are different!')

            var_data = np.zeros(shape=[var_geo_x_2d.shape[0], var_geo_y_2d.shape[1], var_time_steps_expected])
            var_data[:, :, :] = np.nan
            for step in np.arange(0, var_time_steps_expected, 1):
                var_data_step = var_data_3d[:, :, 0]
                var_data[:, :, step] = var_data_step

            var_time_steps_cmp = deepcopy(var_time_steps_expected)

        else:
            log_stream.error(' ===> File ' + file_name + ' format are not expected!')
            raise NotImplemented('Case not implemented yet')

    else:
        log_stream.warning(' ===> File ' + file_name + ' not available in loaded datasets!')
        var_data = None

    if var_data is not None:

        if var_time is not None:

            if isinstance(var_time, pd.Timestamp):
                var_data_time = pd.DatetimeIndex([var_time])
            elif isinstance(var_time, pd.DatetimeIndex):
                var_data_time = deepcopy(var_time)
            else:
                log_stream.error(' ===> Time format is not allowed. Expected Timestamp or Datetimeindex')
                raise NotImplemented('Case not implemented yet')

            var_dset = xr.Dataset(coords={coord_name_time: ([dim_name_time], var_data_time)})
            var_dset.coords[coord_name_time] = var_dset.coords[coord_name_time].astype('datetime64[ns]')

            var_da = xr.DataArray(var_data, name=var_name, dims=dims_order,
                                  coords={coord_name_time: ([dim_name_time], var_data_time),
                                          coord_name_geo_x: ([dim_name_geo_x], var_geo_x_2d[0, :]),
                                          coord_name_geo_y: ([dim_name_geo_y], var_geo_y_2d[:, 0])})

            if file_attrs and geo_attrs:
                obj_attrs = {**file_attrs, **geo_attrs}
            elif (not file_attrs) and geo_attrs:
                obj_attrs = deepcopy(geo_attrs)
            elif file_attrs and (not geo_attrs):
                obj_attrs = deepcopy(file_attrs)
            else:
                obj_attrs = None

            if obj_attrs is not None:
                var_dset.attrs = obj_attrs

            var_dset[var_name] = var_da

        elif var_time is None:

            var_dset = xr.Dataset()
            var_da = xr.DataArray(var_data, name=var_name, dims=dims_order,
                                  coords={coord_name_geo_x: ([dim_name_geo_x], var_geo_x_2d[0, :]),
                                          coord_name_geo_y: ([dim_name_geo_y], var_geo_y_2d[:, 0])})

            if file_attrs and geo_attrs:
                obj_attrs = {**file_attrs, **geo_attrs}
            elif (not file_attrs) and geo_attrs:
                obj_attrs = deepcopy(geo_attrs)
            elif file_attrs and (not geo_attrs):
                obj_attrs = deepcopy(file_attrs)
            else:
                obj_attrs = None

            if obj_attrs is not None:
                var_dset.attrs = obj_attrs

            var_dset[var_name] = var_da

        else:
            log_stream.error(' ===> Error in creating time information for dataset object')
            raise RuntimeError('Unknown error in creating dataset. Check the procedure.')

    else:
        log_stream.warning(' ===> All filenames in the selected period are not available')
        var_dset = None

    return var_dset
# -------------------------------------------------------------------------------------
