
import h5py

import xarray as xr
import numpy as np
import rasterio
from rasterio.crs import CRS

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


# ----------------------------------------------------------------------------------------------------------------------
# method to create a data array
def create_darray_2d(data, geo_x, geo_y, geo_1d=True, time=None, name=None,
                     coord_name_x='west_east', coord_name_y='south_north', coord_name_time='time',
                     dim_name_x='west_east', dim_name_y='south_north', dim_name_time='time',
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x]
    if time is not None:
        dims_order = [dim_name_y, dim_name_x, dim_name_time]

    if geo_1d:
        if geo_x.shape.__len__() == 2:
            geo_x = geo_x[0, :]
        if geo_y.shape.__len__() == 2:
            geo_y = geo_y[:, 0]

        if time is None:
            data_da = xr.DataArray(data,
                                   dims=dims_order,
                                   coords={coord_name_x: (dim_name_x, geo_x),
                                           coord_name_y: (dim_name_y, geo_y)})
        elif isinstance(time, pd.DatetimeIndex):

            if data.shape.__len__() == 2:
                data = np.expand_dims(data, axis=-1)

            data_da = xr.DataArray(data,
                                   dims=dims_order,
                                   coords={coord_name_x: (dim_name_x, geo_x),
                                           coord_name_y: (dim_name_y, geo_y),
                                           coord_name_time: (dim_name_time, time)})
        else:
            logging.error(' ===> Time obj is in wrong format')
            raise IOError('Variable time format not valid')

    else:
        logging.error(' ===> Longitude and Latitude must be 1d')
        raise IOError('Variable shape is not valid')

    if name is not None:
        data_da.name = name

    return data_da
# ----------------------------------------------------------------------------------------------------------------------

def read_smap_l2(file_name):

    with h5py.File(file_name, mode='r') as f:

        name = '/Soil_Moisture_Retrieval_Data/soil_moisture'
        data = f[name][:]
        units = f[name].attrs['units']
        units = units.decode('ascii', 'replace')
        long_name = f[name].attrs['long_name']
        long_name = long_name.decode('ascii', 'replace')
        _FillValue = f[name].attrs['_FillValue']
        valid_max = f[name].attrs['valid_max']
        valid_min = f[name].attrs['valid_min']
        invalid = np.logical_or(data > valid_max, data < valid_min)
        invalid = np.logical_or(invalid, data == _FillValue)
        data[invalid] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

        # Get the geolocation data
        latitude = f['/Soil_Moisture_Retrieval_Data/latitude'][:]
        longitude = f['/Soil_Moisture_Retrieval_Data/longitude'][:]

        col = f['/Soil_Moisture_Retrieval_Data/EASE_column_index'][:]
        row = f['/Soil_Moisture_Retrieval_Data/EASE_row_index'][:]

        return data, longitude, latitude, col, row