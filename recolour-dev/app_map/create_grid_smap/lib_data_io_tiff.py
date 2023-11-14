"""
Library Features:

Name:          lib_data_io_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230915'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import rasterio

import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy
from rasterio.transform import Affine
from rasterio.crs import CRS
from osgeo import gdal, gdalconst

from lib_utils_io import create_darray_2d
from lib_info_args import proj_epsg, time_format_datasets

# logging
logging.getLogger('rasterio').setLevel(logging.WARNING)
# debugging
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize file tiff
def organize_file_tiff(obj_variable, obj_time=None, obj_transform=None, obj_proj=None,
                       var_attr_description='description', var_attr_time='time',
                       var_name_geo_x='longitude', var_name_geo_y='latitude',
                       ):

    if isinstance(obj_time, pd.Timestamp):
        string_time = obj_time.strftime(time_format_datasets)
    elif isinstance(obj_time, str):
        string_time = deepcopy(obj_time)
    else:
        logging.error(' ===> Time obj format is not supported')
        raise NotImplemented('Case not implemented yet')

    var_height, var_width = None, None
    var_geo_x_west, var_geo_x_east, var_geo_y_south, var_geo_y_north = None, None, None, None
    var_geo_transform, var_geo_proj = None, None
    var_data_list, var_metadata_list = [], []
    for var_key, var_obj in obj_variable.items():

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

        var_metadata_step = {var_attr_description: var_key, var_attr_time: string_time}
        var_metadata_list.append(var_metadata_step)

        var_data_step = var_obj.values
        var_data_list.append(var_data_step)

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
def read_file_tiff(file_name, file_time=None, output_format='data_array',
                   var_name_tmpl=None, var_dtype=None,
                   var_limit_min=None, var_limit_max=None, var_proj='EPSG:4326',
                   coord_name_time='time', coord_name_x='longitude', coord_name_y='latitude',
                   dim_name_time='time', dim_name_x='longitude', dim_name_y='latitude'):

    # fill undefined argument(s)
    if var_name_tmpl is None:
        var_name_tmpl = ['var_name_{var_id}']
    if not isinstance(var_name_tmpl, list):
        var_name_tmpl = [var_name_tmpl]
    if var_dtype is None:
        var_dtype = ['float32']
    if not isinstance(var_dtype, list):
        var_dtype = [var_dtype]
    if var_limit_min is None:
        var_limit_min = [None]
    if not isinstance(var_limit_min, list):
        var_limit_min = [var_limit_min]
    if var_limit_max is None:
        var_limit_max = [None]
    if not isinstance(var_limit_max, list):
        var_limit_max = [var_limit_max]

    var_time = None
    if file_time is not None:
        if isinstance(file_time, str):
            var_time = pd.DatetimeIndex([pd.Timestamp(file_time)])
        elif isinstance(file_time, pd.DatetimeIndex):
            var_time = deepcopy(file_time)
        elif isinstance(file_time, pd.Timestamp):
            var_time = pd.DatetimeIndex([file_time])
        else:
            logging.error(' ===> Time obj format is not supported')
            raise NotImplemented('Case not implemented yet')

    # get handle, data and geographical info
    dset = rasterio.open(file_name)
    bounds, res, transform = dset.bounds, dset.res, dset.transform
    tags, bands, metadata = dset.tags, dset.count, dset.profile
    proj = dset.crs.wkt

    if proj is None:
        crs = CRS.from_string(var_proj)
    else:
        crs = dset.crs

    # define geographical variable(s)
    decimal_round_geo = 7
    center_right = bounds.right - (res[0] / 2)
    center_left = bounds.left + (res[0] / 2)
    center_top = bounds.top - (res[1] / 2)
    center_bottom = bounds.bottom + (res[1] / 2)

    if center_bottom > center_top:
        logging.warning(' ===> Coords "center_bottom": ' + str(center_bottom) + ' is greater than "center_top": '
                        + str(center_top) + '. Try to inverse the bottom and top coords. ')
        center_tmp = center_top
        center_top = center_bottom
        center_bottom = center_tmp

    lon = np.arange(center_left, center_right + np.abs(res[0] / 2), np.abs(res[0]), float)
    lat = np.flip(np.arange(center_bottom, center_top + np.abs(res[1] / 2), np.abs(res[1]), float), axis=0)
    lons, lats = np.meshgrid(lon, lat)

    lat_flag, lat_upper, lat_lower = False, lats[0, 0], lats[-1, 0]
    if lat_lower > lat_upper:
        lats = np.flipud(lats)
        lat_flag = True

    min_lon_round = round(np.min(lons), decimal_round_geo)
    max_lon_round = round(np.max(lons), decimal_round_geo)
    min_lat_round = round(np.min(lats), decimal_round_geo)
    max_lat_round = round(np.max(lats), decimal_round_geo)

    center_right_round = round(center_right, decimal_round_geo)
    center_left_round = round(center_left, decimal_round_geo)
    center_bottom_round = round(center_bottom, decimal_round_geo)
    center_top_round = round(center_top, decimal_round_geo)

    assert min_lon_round == center_left_round
    assert max_lon_round == center_right_round
    assert min_lat_round == center_bottom_round
    assert max_lat_round == center_top_round

    var_geo_x_1d, var_geo_y_1d = lons[0, :], lats[:, 0]

    # define data variable(s)
    if bands == 1:
        var_data_list = [dset.read(1)]
    elif bands > 1:
        var_data_list, var_name_list = [], []
        for id in range(0, bands):
            data_tmp = dset.read(id + 1)
            name_tmp = var_name_tmpl.format(var_id=id)
            var_data_list.append(data_tmp)
            var_name_list.append(name_tmp)
    else:
        logging.error(' ===> Bands "' + str(bands) + '" are not supported')
        raise NotImplementedError('Case not implemented yet')

    # filter and organize data
    var_dtype_list = var_dtype * bands
    var_limit_min_list = var_limit_min * bands
    var_limit_max_list = var_limit_max * bands
    var_name_list = var_name_tmpl * bands

    var_data_obj = {}
    for var_id, (var_name_step, var_data_step, var_dtype_step, var_limit_min_step, var_limit_max_step) in enumerate(
                                                        zip(var_name_list, var_data_list, var_dtype_list,
                                                            var_limit_min_list, var_limit_max_list)):
        if var_dtype_step == 'float32':
            var_values_step = np.float32(var_data_step)

            if var_limit_min_step is not None:
                var_limit_min_step = np.float32(var_limit_min_step)
                var_values_step[var_values_step < var_limit_min_step] = np.nan
            if var_limit_max_step is not None:
                var_limit_max_step = np.float32(var_limit_max_step)
                var_values_step[var_values_step > var_limit_max_step] = np.nan

            if lat_flag:
                var_values_step = np.flipud(var_values_step)

            # '''
            plt.figure(); plt.imshow(var_values_step); plt.colorbar();
            plt.show()
            # '''
        else:
            logging.error(' ===> Data type is not allowed.')
            raise NotImplementedError('Case not implemented yet')

        var_data_obj[var_name_step.format(var_id=var_id)] = var_values_step

    # organize data attributes
    data_attrs = {'transform': transform, 'crs': crs,
                  'bbox': [bounds.left, bounds.bottom, bounds.right, bounds.top],
                  'bb_left': bounds.left, 'bb_right': bounds.right,
                  'bb_top': bounds.top, 'bb_bottom': bounds.bottom,
                  'res_lon': res[0], 'res_lat': res[1]}

    # organize data values
    if output_format == 'dictionary':

        # dictionary format
        if file_time is None:
            data_obj = {coord_name_x: var_geo_x_1d, coord_name_y: var_geo_y_1d, **var_data_obj}
        else:
            data_obj = {coord_name_time: var_time,
                        coord_name_x: var_geo_x_1d, coord_name_y: var_geo_y_1d,
                        **var_data_obj}

    elif output_format == 'data_array':

        # data array format
        data_obj = {}
        for var_name_step, var_data_step in var_data_obj.items():

            var_da_step = create_darray_2d(
                var_data_step, var_geo_x_1d, var_geo_y_1d, time=var_time,
                coord_name_time=coord_name_time, coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                dim_name_time=dim_name_time, dim_name_x=dim_name_x, dim_name_y=dim_name_y)

            var_da_step.attrs = data_attrs

            data_obj[var_name_step] = var_da_step

    elif output_format == 'dataset':

        # dataset format
        data_obj = None
        for var_name, var_data_step in var_data_obj.items():

            var_da_step = create_darray_2d(
                var_data_step, var_geo_x_1d, lats[:, 0], time=var_time,
                coord_name_time=coord_name_time, coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                dim_name_time=dim_name_time, dim_name_x=dim_name_x, dim_name_y=dim_name_y)

            var_da_step.attrs = data_attrs

            if data_obj is None:

                if var_time is None:
                    data_obj = xr.Dataset(
                        coords={coord_name_x: ([dim_name_x], var_geo_x_1d),
                                coord_name_y: ([dim_name_y], var_geo_y_1d)})
                else:
                    data_obj = xr.Dataset(
                        coords={coord_name_x: ([dim_name_x], var_geo_x_1d),
                                coord_name_y: ([dim_name_y], var_geo_y_1d),
                                coord_name_time: ([dim_name_time], var_time)})

            data_obj[var_name] = var_da_step.copy()

    else:

        logging.error(' ===> Format for datasets output of "' + file_name + '" is not supported')
        raise NotImplementedError('Case not implemented yet')

    return data_obj, data_attrs

# ----------------------------------------------------------------------------------------------------------------------
