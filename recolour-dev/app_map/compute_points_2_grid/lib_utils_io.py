"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import pandas as pd
import xarray as xr
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


# ----------------------------------------------------------------------------------------------------------------------
# method to read file obj
def filter_dframe_by_vars(dframe_obj_in: pd.DataFrame,
                          dframe_parameter: str = None,
                          dframe_variables_data: list = None, dframe_variables_grid: list = None):
    if dframe_variables_data is None:
        dframe_variables_data = ['gpi', 'lon', 'lat', 'obs']
    if dframe_variables_grid is None:
        dframe_variables_grid = ['committed_area', 'land_flag']

    dframe_variables_in = dframe_variables_data + dframe_variables_grid
    dframe_variables_out = []
    for variable_name in dframe_variables_in:
        if variable_name in list(dframe_obj_in.columns):
            dframe_variables_out.append(variable_name)
        else:
            logging.warning(' ===> Variable "' + variable_name +
                            ' is not available in the DataFrame. Errors could be occurred')

    dframe_obj_out = dframe_obj_in[dframe_variables_out]

    if dframe_parameter not in dframe_obj_out.columns:
        logging.error(' ===> Variable "' + dframe_parameter + '" must be in the DataFrame')
        raise RuntimeError('Check your DataFrame to correctly set the columns names or data')
    return dframe_obj_out
# ----------------------------------------------------------------------------------------------------------------------


