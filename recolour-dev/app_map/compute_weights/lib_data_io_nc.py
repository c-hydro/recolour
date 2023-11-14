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

# default netcdf encoded attributes
attrs_encoded = ["_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping"]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize nc file
def organize_file_nc(obj_variable,
                     var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_x='longitude', coord_name_y='latitude',
                     dim_name_x='longitude', dim_name_y='latitude'):

    # iterate over variable(s)
    variable_dset = None
    for variable_name, variable_da in obj_variable.items():

        if variable_dset is None:

            var_geo_x_1d = variable_da[var_name_geo_x].values
            var_geo_y_1d = variable_da[var_name_geo_y].values

            variable_dset = xr.Dataset(
                coords={coord_name_x: ([dim_name_x], var_geo_x_1d),
                        coord_name_y: ([dim_name_y], var_geo_y_1d)})
                # ,coord_name_time: ([dim_name_time], var_data_time)})

        variable_dset[variable_name] = variable_da.copy()

    return variable_dset
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read nc file
def read_file_nc(file_name, file_variables_selected=None):

    # check file availability
    if os.path.exists(file_name):

        # open file
        file_dset_tmp = xr.open_dataset(file_name)
        # select variable(s)
        if file_variables_selected is not None:
            file_variables_found = list(file_dset_tmp.variables)

            file_variables_matched = []
            for variable_name in file_variables_selected:
                if variable_name in file_variables_found:
                    file_variables_matched.append(variable_name)
                else:
                    logging.warning(' ===> Variable "' + variable_name + '" not found in the filename "' +
                                    file_name + '"')
            # filter variable(s) and keep the selected dataset
            file_dset_def = file_dset_tmp[file_variables_matched]
        else:
            # copy the whole dataset
            file_dset_def = deepcopy(file_dset_tmp)

    else:
        logging.warning(' ===> File variables "' + file_name + '" not found')
        file_dset_def = None

    return file_dset_def

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
