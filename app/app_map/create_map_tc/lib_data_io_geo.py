"""
Library Features:

Name:          lib_data_io_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230307'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# Library
import logging
import os

import numpy as np
import rasterio
from rasterio.crs import CRS

from lib_fx_methods_common import interp_data
from lib_info_args import logger_name
from lib_utils_io import create_darray_2d

# set logger level
logging.getLogger('rasterio').setLevel(logging.WARNING)
# set logger obj
alg_logger = logging.getLogger(logger_name)

# debug
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remap grid data
def remap_grid_data(da_in,
                    lons_in, lats_in, lons_out, lats_out, attr_out=None,
                    interp_method='nearest', interp_no_data=-9999.0,
                    coord_name_x='longitude', coord_name_y='latitude',
                    dim_name_x='longitude', dim_name_y='latitude',
                    ):

    # interpolate input data array
    if attr_out is None:
        attr_out = {}

    values_in = da_in.values

    values_out = interp_data(lons_in, lats_in, values_in, lons_out, lats_out,
                             nodata=interp_no_data, interp_method=interp_method)

    ''' debug
    plt.figure(); plt.imshow(values_in); plt.colorbar()
    plt.figure(); plt.imshow(values_out); plt.colorbar()
    plt.show()
    '''

    # create output data array
    da_out = create_darray_2d(
        values_out, lons_out[0, :], lats_out[:, 0],
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    da_out.attrs = attr_out

    return da_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to check grid data
def check_grid_data(data_file_name, data_obj=None, data_attrs=None, data_mandatory=True):
    if (data_obj is None) and (data_attrs is None):
        if data_mandatory:
            alg_logger.error(' ===> Grid file "' + data_file_name + '" is not available')
            raise FileNotFoundError('File is mandatory to correctly run the algorithm')
        else:
            alg_logger.warning(' ===> Grid file "' + data_file_name + '" is not available')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read grid data
def read_grid_data(file_name, output_format='data_array', output_dtype='float32',
                   var_idx=1, var_limit_min=None, var_limit_max=None, var_proj='EPSG:4326',
                   coord_name_x='longitude', coord_name_y='latitude',
                   dim_name_x='longitude', dim_name_y='latitude',
                   binary_mask=False):

    try:
        dset = rasterio.open(file_name)
        bounds, res, transform = dset.bounds, dset.res, dset.transform
        data = dset.read()

        if dset.crs is None:
            crs = CRS.from_string(var_proj)
        else:
            crs = dset.crs

        if output_dtype == 'float32':
            values = np.float32(data[var_idx, :, :])
        else:
            alg_logger.error(' ===> Data type is not allowed.')
            raise NotImplementedError('Case not implemented yet')

        if var_limit_min is not None:
            var_limit_min = np.float32(var_limit_min)
            values[values < var_limit_min] = np.nan
        if var_limit_max is not None:
            var_limit_max = np.float32(var_limit_max)
            values[values > var_limit_max] = np.nan

        if binary_mask:
            values[(values >= 0) & (values <= 1)] = 1
            values[values < 0] = 0
            values[values > 1] = 0
            values[np.isnan(values)] = 0

        decimal_round_geo = 7
        # 6.5199,47.2175
        center_right = bounds.right - (res[0] / 2)
        center_left = bounds.left + (res[0] / 2)
        center_top = bounds.top - (res[1] / 2)
        center_bottom = bounds.bottom + (res[1] / 2)

        if center_bottom > center_top:
            alg_logger.warning(' ===> Coords "center_bottom": ' + str(center_bottom) + ' is greater than "center_top": '
                            + str(center_top) + '. Try to inverse the bottom and top coords. ')
            center_tmp = center_top
            center_top = center_bottom
            center_bottom = center_tmp

        lon = np.arange(center_left, center_right + np.abs(res[0] / 2), np.abs(res[0]), float)
        lat = np.flip(np.arange(center_bottom, center_top + np.abs(res[1] / 2), np.abs(res[1]), float), axis=0)
        lons, lats = np.meshgrid(lon, lat)

        lat_upper = lats[0, 0]
        lat_lower = lats[-1, 0]
        if lat_lower > lat_upper:
            lats = np.flipud(lats)
            values = np.flipud(values)

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

        data_attrs = {'transform': transform, 'crs': crs,
                      'bbox': [bounds.left, bounds.bottom, bounds.right, bounds.top],
                      'bb_left': bounds.left, 'bb_right': bounds.right,
                      'bb_top': bounds.top, 'bb_bottom': bounds.bottom,
                      'res_lon': res[0], 'res_lat': res[1]}

        if output_format == 'dictionary':

            data_var = {'values': values, 'longitude': lons[0, :], 'latitude': lats[:, 0]}
            data_obj = {**data_var, **data_attrs}

        elif output_format == 'data_array':

            data_obj = create_darray_2d(
                values, lons[0, :], lats[:, 0],
                coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                dim_name_x=dim_name_x, dim_name_y=dim_name_y)

            data_obj.attrs = data_attrs
        else:

            alg_logger.error(' ===> File static "' + file_name + '" output format not allowed')
            raise NotImplementedError('Case not implemented yet')

    except IOError as io_error:

        data_obj, data_attrs = None, None
        lons, lats = None, None
        alg_logger.warning(' ===> File static grid was not correctly open with error "' + str(io_error) + '"')
        logging.warning(' ===> Filename "' + os.path.split(file_name)[1] + '"')

    return data_obj, data_attrs, lons, lats
# ----------------------------------------------------------------------------------------------------------------------
