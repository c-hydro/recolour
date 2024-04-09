"""
Library Features:

Name:          lib_data_io
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
import rasterio
from rasterio.crs import CRS
from osgeo import gdal, gdalconst

from lib_utils_generic import create_darray_2d
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read geographical file
def read_file_geo(file_name_geo):
    """Function to read tif domain mask"""

    if os.path.exists(file_name_geo):
        if file_name_geo.endswith('tif') or file_name_geo.endswith('.tiff'):
            dset_geo = gdal.Open(file_name_geo, gdalconst.GA_ReadOnly)
            geo_proj = dset_geo.GetProjection()
            geo_geotrans = dset_geo.GetGeoTransform()
            geo_data = dset_geo.ReadAsArray()
            geo_band = dset_geo.GetRasterBand(1)
            geo_wide = dset_geo.RasterXSize
            geo_high = dset_geo.RasterYSize
            geo_min_x = geo_geotrans[0]
            geo_max_y = geo_geotrans[3]
            geo_max_x = geo_min_x + geo_geotrans[1] * dset_geo.RasterXSize
            geo_min_y = geo_max_y + geo_geotrans[5] * dset_geo.RasterYSize
            geo_mask_grid = np.where(geo_data == 1, 1, -9999) # in smap_grid no domain is defined, it is blank
            geo_mask_idx = np.where(geo_mask_grid.ravel() == -9999)
        else:
            logging.error(' ===> Geographical file ' + file_name_geo + ' format unknown')
            raise NotImplementedError('File type reader not implemented yet')
    else:
        logging.error(' ===> Geographical file ' + file_name_geo + ' not found')
        raise IOError('Geographical file location or name is wrong')

    return geo_proj, geo_geotrans, geo_data, \
           geo_wide, geo_high, geo_min_x, geo_max_y, geo_max_x, geo_min_y, geo_mask_idx
# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read grid data
def read_grid_data(file_name, output_format='data_array', output_dtype='float32',
                   var_limit_min=None, var_limit_max=None, var_proj='EPSG:4326',
                   coord_name_x='longitude', coord_name_y='latitude',
                   dim_name_x='longitude', dim_name_y='latitude'):

    try:
        dset = rasterio.open(file_name)
        bounds, res, transform = dset.bounds, dset.res, dset.transform
        data = dset.read()

        if dset.crs is None:
            crs = CRS.from_string(var_proj)
        else:
            crs = dset.crs

        if output_dtype == 'float32':
            values = np.float32(data[0, :, :])
        else:
            logging.error(' ===> Data type is not allowed.')
            raise NotImplementedError('Case not implemented yet')

        if var_limit_min is not None:
            var_limit_min = np.float32(var_limit_min)
            values[values < var_limit_min] = np.nan
        if var_limit_max is not None:
            var_limit_max = np.float32(var_limit_max)
            values[values > var_limit_max] = np.nan

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

            logging.error(' ===> File static "' + file_name + '" output format not allowed')
            raise NotImplementedError('Case not implemented yet')

    except IOError as io_error:

        data_obj, data_attrs = None, None
        logging.warning(' ===> File static grid was not correctly open with error "' + str(io_error) + '"')
        logging.warning(' ===> Filename "' + os.path.split(file_name)[1] + '"')

    return data_obj, data_attrs
# ----------------------------------------------------------------------------------------------------------------------
