"""
Library Features:

Name:          lib_utils_io_nc
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import netCDF4
import time
import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy

from lib_info_args import logger_name
from lib_info_args import proj_epsg, time_format_datasets

# set logger level
logging.getLogger('rasterio').setLevel(logging.WARNING)
# set logger obj
alg_logger = logging.getLogger(logger_name)

# default netcdf encoded attributes
attrs_encoded = ["_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping"]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize nc file
def organize_file_nc(obj_variable, obj_time=None,
                     var_name_time='time', var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_time='time', coord_name_x='longitude', coord_name_y='latitude',
                     dim_name_time='time', dim_name_x='longitude', dim_name_y='latitude'):

    var_data_time = None
    if obj_time is not None:
        if isinstance(obj_time, str):
            var_data_time = pd.DatetimeIndex([pd.Timestamp(obj_time)])
        elif isinstance(obj_time, pd.DatetimeIndex):
            var_data_time = deepcopy(obj_time)
        elif isinstance(obj_time, pd.Timestamp):
            var_data_time = pd.DatetimeIndex([obj_time])
        else:
            alg_logger.error(' ===> Time obj format is not supported')
            raise NotImplemented('Case not implemented yet')

    # iterate over variable(s)
    variable_dset = None
    for variable_name, variable_da in obj_variable.items():

        if variable_dset is None:

            var_geo_x_1d = variable_da[var_name_geo_x].values
            var_geo_y_1d = variable_da[var_name_geo_y].values

            if obj_time is None:
                variable_dset = xr.Dataset(
                    coords={coord_name_x: ([dim_name_x], var_geo_x_1d),
                            coord_name_y: ([dim_name_y], var_geo_y_1d)})
            else:
                variable_dset = xr.Dataset(
                    coords={coord_name_x: ([dim_name_x], var_geo_x_1d),
                            coord_name_y: ([dim_name_y], var_geo_y_1d),
                            coord_name_time: ([dim_name_time], var_data_time)})

        variable_dset[variable_name] = variable_da.copy()

    return variable_dset
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read nc file
def read_file_nc(file_name, file_variables_selected=None,
                 var_name_geo_x='lon', var_name_geo_y='lat'):

    if file_variables_selected is None:
        file_variables_selected = ['var40', 'var41', 'var42', 'var43',]

    # check file availability
    if os.path.exists(file_name):

        # open file
        file_dset_tmp = xr.open_dataset(file_name)

        # adjust coords
        file_geo_x_max = np.nanmax(file_dset_tmp.coords[var_name_geo_x])
        if file_geo_x_max > 180:
            file_dset_tmp.coords[var_name_geo_x] = (file_dset_tmp.coords[var_name_geo_x] + 180) % 360 - 180
            file_dset_tmp = file_dset_tmp.sortby(file_dset_tmp[var_name_geo_x])

        # select variable(s)
        file_variables_found = list(file_dset_tmp.variables)

        data_obj, data_attrs = {}, {}
        for var_name in file_variables_selected:
            if var_name in file_variables_found:
                file_values_tmp = np.squeeze(file_dset_tmp[var_name].values)
                file_attrs_tmp = file_dset_tmp[var_name].attrs
                data_obj[var_name] = file_values_tmp
                data_attrs[var_name] = file_attrs_tmp
            else:
                alg_logger.warning(' ===> Variable "' + var_name + '" not found in the filename "' + file_name + '"')

        data_lons = file_dset_tmp[var_name_geo_x].values
        data_lats = file_dset_tmp[var_name_geo_y].values

        common_attrs = file_dset_tmp.attrs

    else:
        alg_logger.warning(' ===> File "' + file_name + '" not found')
        data_obj, data_attrs = None, None
        data_lons, data_lats = None, None
        common_attrs = None

    return data_obj, data_attrs, data_lons, data_lats, common_attrs

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write nc file
def write_file_nc(file_name, dset_data,
                  dset_mode='w', dset_engine='netcdf4', dset_compression=0, dset_format='NETCDF4',
                  dim_key_time='time', no_data=-9999.0):

    dset_encoded = dict(zlib=True, complevel=dset_compression)

    dset_encoding = {}
    for var_name in dset_data.data_vars:

        if isinstance(var_name, bytes):
            tmp_name = var_name.decode("utf-8")
            dset_data.rename({var_name: tmp_name})
            var_name = deepcopy(tmp_name)

        var_data = dset_data[var_name]
        if len(var_data.dims) > 0:
            dset_encoding[var_name] = deepcopy(dset_encoded)

        var_attrs = dset_data[var_name].attrs
        if var_attrs:
            for attr_key, attr_value in var_attrs.items():
                if attr_key in attrs_encoded:

                    dset_encoding[var_name][attr_key] = {}

                    if isinstance(attr_value, list):
                        attr_string = [str(value) for value in attr_value]
                        attr_value = ','.join(attr_string)

                    dset_encoding[var_name][attr_key] = attr_value

        if '_FillValue' not in list(dset_encoding[var_name].keys()):
            dset_encoding[var_name]['_FillValue'] = no_data

    if dim_key_time in list(dset_data.coords):
        dset_encoding[dim_key_time] = {'calendar': 'gregorian'}

    dset_data.to_netcdf(path=file_name, format=dset_format, mode=dset_mode,
                        engine=dset_engine, encoding=dset_encoding)
# ----------------------------------------------------------------------------------------------------------------------
