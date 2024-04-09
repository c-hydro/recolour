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
import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy

from lib_info_args import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)

# default netcdf encoded attributes
attrs_encoded = ["_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping"]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize nc file
def organize_file_nc(obj_data_in, obj_time=None, obj_variable=None, obj_cell=None):

    if obj_variable is None:
        obj_variable = {
            'var40': "soil_moisture",
            "time": "time",
            "lon": "longitude", "lat": "latitude", "location_id": "location_id"
        }

    var_name_in_list, var_name_out_list = list(obj_variable.keys()), list(obj_variable.values())

    if obj_data_in is not None:

        obj_variable_check, alg_variable_missed = {}, None
        for var_name_in, var_name_out in obj_variable.items():
            if var_name_in in obj_data_in.variables:
                obj_variable_check[var_name_in] = var_name_out
            else:
                alg_logger.warning(' ===> File variable "' + var_name_in + '" is not available in the file')
                if alg_variable_missed is None:
                    alg_variable_missed = []
                alg_variable_missed.append(var_name_in)

        # check availability of all variable(s)
        if alg_variable_missed is None:

            if isinstance(obj_data_in, xr.Dataset):
                obj_data_tmp = obj_data_in.rename(obj_variable)
            else:
                alg_logger.error(' ===> Obj data format is not supported')
                raise NotImplemented('Case not implemented yet')

            # filter dataset by variable(s)
            obj_data_tmp = obj_data_tmp[var_name_out_list]

            if obj_time is not None:
                if 'time' in list(obj_data_in.variables):
                    time_stamp_max = pd.DatetimeIndex(obj_data_tmp['time'].values).max()

                else:
                    alg_logger.error(' ===> Time data is not available')
                    raise RuntimeError('Time data is needed by the algorithm to properly run')

                obj_data_out = obj_data_tmp.where(obj_data_tmp['time'] == time_stamp_max, drop=True)
            else:
                obj_data_out = deepcopy(obj_data_tmp)

            if obj_cell is not None:
                if 'locations' in list(obj_data_out.dims):
                    obj_n = obj_data_out.dims['locations']
                    obj_arr = [obj_cell] * obj_n

                    obj_da = xr.DataArray(obj_arr,
                                          dims=["locations"],
                                          coords={'locations': ('locations', range(obj_n))})

                    obj_data_out['cell'] = obj_da

                elif 'dim' in list(obj_data_out.dims):
                    obj_n = obj_data_out.dims['dim']
                    obj_arr = [obj_cell] * obj_n

                    obj_da = xr.DataArray(obj_arr,
                                          dims=["dim"],
                                          coords={'dim': ('dim', range(obj_n))})
                    obj_data_out['cell'] = obj_da
                else:
                    alg_logger.error(' ===> File dimensions must be "dim" or "locations" to be used by the algorithm')
                    raise NotImplemented('Case not implemented yet')

        else:
            string_variable_missed = ','.join(alg_variable_missed)
            alg_logger.warning(' ===> Object data does not have all variable(s). Some/All variable(s) "' +
                               string_variable_missed + '" is/are missed')
            obj_data_out = None

    else:
        alg_logger.warning(' ===> Object data is defined by NoneType. Data are not available')
        obj_data_out = None

    return obj_data_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read nc file
def read_file_nc(file_name):

    # check file availability
    if os.path.exists(file_name):
        # open file
        file_dset_def = xr.open_dataset(file_name)
    else:
        alg_logger.warning(' ===> File variables "' + file_name + '" not found')
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
