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
import xarray as xr
import rasterio
from rasterio.crs import CRS

from lib_utils_io import create_darray_2d

# logging
logging.getLogger('rasterio').setLevel(logging.WARNING)

# debug
try:
    import matplotlib.pylab as plt
except ImportError:
    pass
# ----------------------------------------------------------------------------------------------------------------------##


# ----------------------------------------------------------------------------------------------------------------------
# method to check grid data
def check_grid_data(data_file_name, data_obj=None, data_attrs=None, data_mandatory=True):
    if (data_obj is None) and (data_attrs is None):
        if data_mandatory:
            logging.error(' ===> Grid file "' + data_file_name + '" is not available')
            raise FileNotFoundError('File is mandatory to correctly run the algorithm')
        else:
            logging.warning(' ===> Grid file "' + data_file_name + '" is not available')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read grid data (ascii or tiff format)
def read_grid_netcdf(
    file_name,
    output_format="data_array",
    output_dtype="float32",
    var_name="mask",
    var_limit_min=None,
    var_limit_max=None,
    var_proj="EPSG:4326",
    coord_name_x="longitude",
    coord_name_y="latitude",
    dim_name_x="longitude",
    dim_name_y="latitude",
):
    """
    Read H26 NetCDF grid file with:

        longitude(lon)  -> 1D
        latitude(lat)   -> 1D
        mask(lat, lon)  -> 2D

    Returns:
        data_obj, data_attrs
    """

    try:
        dset = xr.open_dataset(file_name)

        if var_name not in dset:
            raise KeyError(f'Variable "{var_name}" not found in {file_name}')

        # Get coordinate names
        if coord_name_x in dset:
            lon = dset[coord_name_x].values
        elif "lon" in dset:
            lon = dset["lon"].values
        else:
            raise KeyError("Longitude coordinate not found")

        if coord_name_y in dset:
            lat = dset[coord_name_y].values
        elif "lat" in dset:
            lat = dset["lat"].values
        else:
            raise KeyError("Latitude coordinate not found")

        values = dset[var_name].values

        # Remove singleton dimensions if needed
        values = np.squeeze(values)

        if output_dtype == "float32":
            values = values.astype(np.float32)
        elif output_dtype == "int8":
            values = values.astype(np.int8)
        else:
            logging.error(" ===> Data type is not allowed.")
            raise NotImplementedError("Case not implemented yet")

        # Apply limits
        if var_limit_min is not None:
            var_limit_min = np.float32(var_limit_min)
            values = values.astype(np.float32)
            values[values < var_limit_min] = np.nan

        if var_limit_max is not None:
            var_limit_max = np.float32(var_limit_max)
            values = values.astype(np.float32)
            values[values > var_limit_max] = np.nan

        # Ensure latitude is north -> south
        if lat[0] < lat[-1]:
            lat = np.flip(lat)
            values = np.flipud(values)

        # Ensure longitude is west -> east
        if lon[0] > lon[-1]:
            lon = np.flip(lon)
            values = np.fliplr(values)

        res_lon = float(np.abs(lon[1] - lon[0])) if lon.size > 1 else np.nan
        res_lat = float(np.abs(lat[1] - lat[0])) if lat.size > 1 else np.nan

        bb_left = float(np.min(lon) - res_lon / 2)
        bb_right = float(np.max(lon) + res_lon / 2)
        bb_bottom = float(np.min(lat) - res_lat / 2)
        bb_top = float(np.max(lat) + res_lat / 2)

        data_attrs = {
            "crs": var_proj,
            "bbox": [bb_left, bb_bottom, bb_right, bb_top],
            "bb_left": bb_left,
            "bb_right": bb_right,
            "bb_top": bb_top,
            "bb_bottom": bb_bottom,
            "res_lon": res_lon,
            "res_lat": res_lat,
            "source_file": os.path.basename(file_name),
            "source_variable": var_name,
        }

        if output_format == "dictionary":

            data_obj = {
                "values": values,
                "longitude": lon,
                "latitude": lat,
                **data_attrs,
            }

        elif output_format == "data_array":

            data_obj = create_darray_2d(
                values,
                lon,
                lat,
                coord_name_x=coord_name_x,
                coord_name_y=coord_name_y,
                dim_name_x=dim_name_x,
                dim_name_y=dim_name_y,
            )

            data_obj.attrs = data_attrs

        else:
            logging.error(f' ===> File static "{file_name}" output format not allowed')
            raise NotImplementedError("Case not implemented yet")

        dset.close()

    except IOError as io_error:
        data_obj, data_attrs = None, None
        logging.warning(
            f' ===> File static grid was not correctly open with error "{io_error}"'
        )
        logging.warning(f' ===> Filename "{os.path.split(file_name)[1]}"')

    return data_obj, data_attrs
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read grid data (ascii or tiff format)
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
        # 6.5199,47.2175
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
