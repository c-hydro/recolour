# lib_utils_io.py
"""
I/O routines for point data analysis:
- read_data_points: read point CSV
- write_geotiff: write interpolated grid to GeoTIFF
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json
import pandas as pd
import numpy as np
import rasterio

from rasterio.transform import Affine
from rasterio.transform import from_origin
from rasterio.crs import CRS

from osgeo import gdal, gdalconst

from lib_utils_generic import create_darray
from lib_info_args import logger_name

# logging
logging.getLogger('rasterio').setLevel(logging.WARNING)

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read a JSON configuration file
def read_json(json_path: str) -> dict:
    with open(json_path) as cf:
        cfg = json.load(cf)
    return cfg
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# method to read data grid
def read_data_grid(file_name, output_format='data_array', output_dtype='float32',
                   var_limit_min=None, var_limit_max=None, var_proj='EPSG:4326',
                   coord_name_x='longitude', coord_name_y='latitude',
                   dim_name_x='longitude', dim_name_y='latitude'):

    try:
        dset = rasterio.open(file_name)
        bounds = dset.bounds
        res = dset.res
        transform = dset.transform
        data = dset.read()

        if dset.crs is None:
            crs = CRS.from_string(var_proj)
        else:
            crs = dset.crs

        if output_dtype == 'float32':
            values = np.float32(data[0, :, :])
        else:
            alg_logger.error(' ===> Data type is not allowed.')
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
            alg_logger.warning(
                ' ===> Coords "center_bottom": ' + str(center_bottom) + ' is greater than "center_top": '
                + str(center_top) + '. Try to inverse the bottom and top coords. ')
            center_tmp = center_top
            center_top = center_bottom
            center_bottom = center_tmp

            values = np.flipud(values)

        ''' debug 
        plt.figure()
        plt.imshow(values)
        plt.colorbar()
        plt.show()
        '''

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

            data_obj = create_darray(
                values, lons[0, :], lats[:, 0],
                coord_name_x=coord_name_x, coord_name_y=coord_name_y,
                dim_name_x=dim_name_x, dim_name_y=dim_name_y)

            data_obj.attrs = data_attrs
        else:

            alg_logger.error(' ===> File static "' + file_name + '" output format not allowed')
            raise NotImplementedError('Case not implemented yet')

    except IOError as io_error:
        alg_logger.error(' ===> File static grid "' + file_name + '" was not correctly open with error "' + str(io_error) + '"')
        raise IOError(io_error)

    return data_obj
# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read point data from a semicolon-delimited CSV file
def read_data_points(csv_path: str, csv_var: str ='sm', csv_delimiter: str = ';') -> tuple:
    """
    Reads point data from a semicolon-delimited CSV with 'lon','lat','sm' columns.

    Parameters
    ----------
    csv_path : str
        Path to the input CSV file.

    Returns
    -------
    lons, lats, values : np.ndarray
        Arrays of longitude, latitude, and soil moisture values.
    """

    if not os.path.exists(csv_path):
        alg_logger.error(f" ===> CSV file not found: {csv_path}")
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df_complete = pd.read_csv(csv_path, sep=csv_delimiter)
    df_filter = df_complete[[csv_var, 'lon', 'lat']]

    return df_complete['lon'].values, df_complete['lat'].values, df_complete[csv_var].values, df_filter
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Method to write file tiff
def write_data_grid(file_name, file_data, file_wide, file_high, file_transform, file_proj,
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
