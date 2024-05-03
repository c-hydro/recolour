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
from lib_utils_grid import read_grid_file, find_grid_idx_to_data

# set logger
alg_logger = logging.getLogger(logger_name)

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

        variable_dset[variable_name] = variable_da

    return variable_dset
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


# ----------------------------------------------------------------------------------------------------------------------
# method to read file cell in netcdf format
def read_file_cell(file_name, selected_variables=None, time_variable='time', row_size_variable='row_size'):

    variable_workspace, variable_attrs = None, None
    with netCDF4.Dataset(file_name) as file_handle:

        file_variables = file_handle.variables
        for variable_tmp in selected_variables:
            if variable_tmp not in list(file_variables):
                if variable_tmp != row_size_variable:
                    alg_logger.error(
                        ' ===> Variable "' + variable_tmp + '" is not found in the cell file "' + file_name + '"')
                    raise RuntimeError('All the selected variables must be in the selected cell file')
                else:
                    alg_logger.warning(
                        ' ===> Variable "' + variable_tmp + '" is not found in the cell file "' +
                        file_name + '". In case of regular grid, this variable will be defined by the algorithm')

        for variable_id, variable_name in enumerate(file_variables):

            if selected_variables is not None:
                if variable_name in selected_variables:
                    variable_collections = file_variables[variable_name]
                    variable_data_tmp = variable_collections[:]
                else:
                    variable_data_tmp = None
            else:
                variable_collections = file_variables[variable_name]
                variable_data_tmp = variable_collections[:]

            if variable_data_tmp is not None:

                if variable_attrs is None:
                    variable_attrs = {}

                if isinstance(variable_data_tmp, np.ma.MaskedArray):
                    variable_data_def = variable_data_tmp.data
                else:
                    variable_data_def = deepcopy(variable_data_tmp)

                if variable_name == time_variable:
                    variable_data_def = netCDF4.num2date(
                        variable_data_def, units=variable_collections.units,
                        only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                    variable_attrs['time_units'] = variable_collections.units

                if variable_workspace is None:
                    variable_workspace = {}
                variable_workspace[variable_name] = variable_data_def
            else:
                # logging.warning(' ===> Variable "' + variable_name + '" is not saved in the workspace file')
                pass

    return variable_workspace, variable_attrs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file collection in netcdf format
def read_file_collection(file_name_data,
                         file_name_grid='TUW_WARP5_grid_info_2_3.nc', file_obj_grid=None,
                         var_dest_list=None, var_registry_list=None, var_extra_list=None,
                         variable_time='time',
                         variable_idx='location', variable_geo_x='longitude', variable_geo_y='latitude',
                         variable_type='ASCAT',
                         variable_no_data=-9999.0, variable_min_value=0, variable_max_value=100, variable_scale_factor=1,
                         variable_land_filter=True, variable_committed_filter=True):

    # organize variable(s) settings in list format
    if not isinstance(variable_no_data, list):
        variable_no_data = [variable_no_data]
    if not isinstance(variable_min_value, list):
        variable_min_value = [variable_min_value]
    if not isinstance(variable_max_value, list):
        variable_max_value = [variable_max_value]
    if not isinstance(variable_scale_factor, list):
        variable_scale_factor = [variable_scale_factor]

    # get reference grid
    if file_obj_grid is None:
        if os.path.exists(file_name_grid):
            grid_dframe = read_grid_file(file_name_grid)
        else:
            logging.error(' ===> File grid "' + file_name_grid + '" is not found!')
            raise FileNotFoundError('File must be available to correctly run the algorithm')
    else:
        grid_dframe = deepcopy(file_obj_grid)

    if var_dest_list is None:
        var_dest_list = ['sm_values']
    if var_registry_list is None:
        var_registry_list = ['longitude', 'latitude', 'location']
    if var_extra_list is None:
        var_extra_list = []

    variable_list = var_dest_list + var_registry_list + var_extra_list

    if 'row_size' in variable_list:
        variable_list.remove('row_size')

    with netCDF4.Dataset(file_name_data) as file_handle:

        variable_collection, attr_collection = {}, {}
        variable_data_idx, variable_data_geo_x, variable_data_geo_y = None, None, None
        for variable_id, variable_name in enumerate(variable_list):

            variable_tmp = file_handle.variables[variable_name]

            attrs_data = {}
            if hasattr(variable_tmp, 'units'):
                variable_units = variable_tmp.units
            else:
                variable_units = None

            variable_data = variable_tmp[:]

            if isinstance(variable_data, np.ma.MaskedArray):
                variable_data = variable_data.data

            if variable_idx == variable_name:
                variable_data_idx = np.array((deepcopy(variable_data)), dtype=int)
            elif variable_geo_x == variable_name:
                variable_data_geo_x = deepcopy(variable_data)
                variable_name_geo_x = deepcopy(variable_name)
            elif variable_geo_y == variable_name:
                variable_data_geo_y = deepcopy(variable_data)
                variable_name_geo_y = deepcopy(variable_name)
            elif variable_time in variable_name:

                variable_tmp = deepcopy(variable_data)
                variable_tmp = netCDF4.num2date(
                    variable_tmp, units=variable_units,
                    only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                variable_data = pd.DatetimeIndex(variable_tmp)

            attrs_data['time_units'] = variable_units

            attr_collection[variable_name] = attrs_data

            # apply variable settings
            if variable_name in var_dest_list:

                no_data, scale_factor = variable_no_data[variable_id], variable_scale_factor[variable_id]
                min_value, max_value = variable_min_value[variable_id], variable_max_value[variable_id]

                if no_data is not None:
                    variable_data[variable_data == no_data] = np.nan
                if min_value is not None:
                    variable_data[variable_data < min_value] = np.nan
                if max_value is not None:
                    variable_data[variable_data > max_value] = np.nan
                if scale_factor is not None:
                    variable_data = variable_data * scale_factor

            variable_collection[variable_name] = variable_data

    collections_dframe = pd.DataFrame(variable_collection, index=variable_data_idx)

    # check variable for the different type(s)
    if variable_type == 'ASCAT':

        # case ascat
        gpi_grid = grid_dframe['gpi'].values
        car_grid = grid_dframe['committed_area'].values
        land_grid = grid_dframe['land_flag'].values
        gpi_variable = collections_dframe.index.values

        idx_variable = find_grid_idx_to_data(gpi_grid, gpi_variable)

        car_data = car_grid[idx_variable]
        collections_dframe['committed_area'] = car_data
        land_data = land_grid[idx_variable]
        collections_dframe['land_flag'] = land_data

    else:

        # case default
        collections_dframe['committed_area'] = 1
        collections_dframe['land_flag'] = 1
        collections_dframe['longitude'] = variable_data_geo_x
        collections_dframe['latitude'] = variable_data_geo_y

    collections_dframe = collections_dframe.rename(
        columns={'longitude': variable_name_geo_x, 'latitude': variable_name_geo_y})

    if variable_land_filter:
        collections_dframe = collections_dframe[collections_dframe['land_flag'] == 1]
    if variable_committed_filter:
        collections_dframe = collections_dframe[collections_dframe['committed_area'] == 1]

    collections_dframe.attrs = attr_collection

    return collections_dframe

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write file collection in netcdf format
def write_file_collection(file_name, file_obj, file_tag_location='location'):

    # get data dimension
    file_n = file_obj[file_tag_location].__len__()
    file_var_list = list(file_obj.columns)

    file_attrs = file_obj.attrs
    if 'time_units' in list(file_attrs.keys()):
        file_time_units = file_attrs['time_units']
    else:
        file_time_units = 'days since 1858-11-17 00:00:00'

    # open and init file
    file_handle = netCDF4.Dataset(file_name, 'w')
    file_handle.createDimension('data', file_n)
    # add attr(s)
    file_handle.file_date = 'Created ' + time.ctime(time.time())

    # iterate over datasets
    for var_name in file_var_list:
        var_values = file_obj[var_name].values

        if isinstance(var_values[0], str):
            var_data = np.array(var_values, dtype=object)
            var_handle = file_handle.createVariable(varname=var_name, dimensions=('data',), datatype='str')
        elif isinstance(var_values[0], (float, np.floating)):
            var_data = var_values
            var_handle = file_handle.createVariable(varname=var_name, dimensions=('data',), datatype='f4')
        elif isinstance(var_values[0], (int, np.integer)):
            var_data = var_values
            var_handle = file_handle.createVariable(varname=var_name, dimensions=('data',), datatype='i4')
        elif isinstance(var_values[0], np.datetime64):
            var_tmp = pd.DatetimeIndex(var_values).to_pydatetime()
            var_data = netCDF4.date2num(var_tmp, units=file_time_units, calendar='gregorian')
            var_handle = file_handle.createVariable(varname=var_name, dimensions=('data',), datatype='f8')
            var_handle.units = file_time_units
        else:
            logging.error(' ===> Datasets format is not expected')
            raise IOError('Dump datasets failed due to the datasets format')

        var_handle[:] = var_data

    file_handle.close()

# ----------------------------------------------------------------------------------------------------------------------
