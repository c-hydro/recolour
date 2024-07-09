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
from lib_utils_io import create_darray
from lib_utils_grid import read_grid_file, find_grid_idx_to_data

# set logger
alg_logger = logging.getLogger(logger_name)

# default netcdf cell expected attributes
attrs_expected = ['long_name', 'standard_name', 'units', 'valid_range', 'name', 'sample_dimension']
# default netcdf file encoded attributes
attrs_encoded = ["_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping"]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize nc file
def organize_file_nc(obj_data, obj_time=None,
                     obj_geo_x=None, obj_geo_y=None,
                     obj_variable_in=None, obj_variable_out=None,
                     var_name_time='time', var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_time='time', coord_name_x='longitude', coord_name_y='latitude',
                     dim_name_time='time', dim_name_x='longitude', dim_name_y='latitude'):

    # organize time information
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

    # organize geo information
    var_geo_x_1d = np.unique(obj_geo_x.flatten())
    var_geo_y_1d = np.unique(obj_geo_y.flatten())

    # iterate over variable(s)
    variable_dset = xr.Dataset()

    for variable_name_in, variable_name_out in zip(obj_variable_in, obj_variable_out):

        if variable_name_in in list(obj_data.keys()):
            variable_data = obj_data[variable_name_in]

            variable_da = create_darray(
                variable_data, obj_geo_x, obj_geo_y,
                geo_1d=True, time=var_data_time,
                coord_name_x=coord_name_x, coord_name_y=coord_name_y, coord_name_time=coord_name_time,
                dim_name_x=dim_name_x, dim_name_y=dim_name_y, dim_name_time=dim_name_time)

            variable_dset[variable_name_out] = variable_da.copy()

        else:
            alg_logger.warning(' ===> Variable "' + variable_name_in + '" not found in the dataset')

    return variable_dset
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read nc file
def read_file_nc(file_name, file_variables=None):

    # check file availability
    if os.path.exists(file_name):

        # open file
        file_dset = xr.open_dataset(file_name)

        # organize file obj
        file_obj = {}
        for var_name_out, var_name_in in file_variables.items():
            if var_name_in in list(file_dset.data_vars):
                file_obj[var_name_out] = file_dset[var_name_in].values
            else:
                alg_logger.warning(' ===> Variable "' + var_name_in + '" not found in the dataset')
                file_obj[var_name_out] = None
    else:
        alg_logger.warning(' ===> File name "' + file_name + '" not found')
        file_obj = None

    return file_obj

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write nc file
def write_file_nc(file_name, dset_data,
                  dset_mode='w', dset_engine='netcdf4', dset_compression=9, dset_format='NETCDF4',
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
def read_file_cell(file_name, selected_variables=None, time_variable='time'):

    # open file and get variable(s)
    variable_workspace, attrs_workspace = None, None
    with netCDF4.Dataset(file_name) as file_handle:

        # check variables
        file_variables = file_handle.variables
        for variable_tmp in selected_variables:
            if variable_tmp not in list(file_variables):
                alg_logger.error(' ===> Variable "' + variable_tmp + '" is not found in the cell file "' + file_name + '"')
                raise RuntimeError('All the selected variables must be in the selected cell file')

        # iterate over variable(s)
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

                if isinstance(variable_data_tmp, np.ma.MaskedArray):
                    variable_data_def = variable_data_tmp.data
                else:
                    variable_data_def = deepcopy(variable_data_tmp)

                variable_attrs = {}
                for attr_name in attrs_expected:
                    if hasattr(variable_collections, attr_name):
                        variable_attrs[attr_name] = getattr(variable_collections, attr_name)

                if variable_name == time_variable:
                    variable_data_def = netCDF4.num2date(
                        variable_data_def, units=variable_collections.units,
                        only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                    variable_attrs['time_units'] = variable_collections.units

                if variable_workspace is None:
                    variable_workspace = {}
                variable_workspace[variable_name] = variable_data_def

                if attrs_workspace is None:
                    attrs_workspace = {}
                attrs_workspace[variable_name] = variable_attrs

            else:
                # alg_logger.warning(' ===> Variable "' + variable_name + '" is not saved in the workspace file')
                pass

    return variable_workspace, attrs_workspace
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file collection in netcdf format
def read_file_collection_OLD(file_name_data,
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
            alg_logger.error(' ===> File grid "' + file_name_grid + '" is not found!')
            raise FileNotFoundError('File must be available to correctly run the algorithm')
    else:
        grid_dframe = deepcopy(file_obj_grid)

    if var_dest_list is None:
        var_dest_list = ['sm_values', 'sm_noise', 'flag_corr', 'flag_proc']
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
# method to write file cell in netcdf format
def write_file_cell(file_name,
                    data_obj, registry_obj,
                    data_attrs, registry_attrs,
                    file_tag_location='location_id', file_tag_row_size='row_size',
                    file_tag_time='time',
                    file_feature_type='timeSeries'):

    # split file name parts
    file_folder, file_id = os.path.split(file_name)

    # get data dimension
    if file_tag_location in list(registry_obj.columns):
        loc_n = registry_obj[file_tag_location].__len__()
        data_n = data_obj[file_tag_time].__len__()
        mode = 'location'
    elif file_tag_time in list(data_obj.keys()):
        data_n = data_obj[file_tag_time].__len__()
        loc_n = registry_obj[file_tag_row_size].__len__()
        mode = 'row_size'
    else:
        logging.error(' ===> Data object does not contain location or row_size information')
        raise IOError('Check your data object')

    # define variable list
    if isinstance(data_obj, pd.DataFrame):
        file_var_list_data = list(data_obj.columns)
    elif isinstance(data_obj, xr.Dataset):
        file_var_list_data = list(data_obj.data_vars) + list(data_obj.dims)
        if file_tag_location in file_var_list_data:
            file_var_list_data.remove(file_tag_location)
    else:
        logging.error(' ===> Datasets format is not expected')
        raise IOError('Dump datasets failed due to the datasets format')
    # define registry list
    file_var_list_registry = list(registry_obj.columns)

    # search time attributes in data and registry attributes objects
    common_attrs = {**data_attrs, **registry_attrs}
    if 'time' in list(common_attrs.keys()):
        if 'time_units' in list(common_attrs[file_tag_time].keys()):
            file_time_units = common_attrs[file_tag_time]['time_units']
        else:
            file_time_units = 'days since 1858-11-17 00:00:00'
    else:
        alg_logger.error(' ===> Time information is not available in the data and registry attributes')
        raise IOError('Check your data object')

    # open and init file (location or row_size mode)
    file_handle = netCDF4.Dataset(file_name, 'w')
    if mode == 'location':
        dim_locations = file_handle.createDimension('locations', loc_n)
        dim_obs = file_handle.createDimension('obs', data_n)
    elif mode == 'row_size':
        dim_locations = file_handle.createDimension('locations', loc_n)
        dim_obs = file_handle.createDimension('obs', data_n)
    else:
        logging.error(' ===> Datasets format is not expected')
        raise IOError('Dump datasets failed due to the datasets format')

    # add common attributes defined for cell file
    file_handle.featureType = file_feature_type
    file_handle.file_date = 'Created ' + time.ctime(time.time())
    file_handle.id = file_id

    # iterate over variable data list
    for file_var_name in file_var_list_data:

        # get values and attributes
        var_values = data_obj[file_var_name].values
        var_name = data_obj[file_var_name].name

        var_attrs = {}
        if file_var_name in list(data_attrs.keys()):
            var_attrs = data_attrs[file_var_name]
        var_attrs['name'] = var_name

        if var_values.ndim == 1:
            var_check = var_values[0]
        elif var_values.ndim == 2:
            var_check = var_values[0, 0]
        else:
            logging.error(' ===> Values dimensions is not expected')
            raise IOError('Dump datasets failed due to the values dimensions')

        if isinstance(var_check, str):

            var_data = np.array(var_values, dtype=object)
            if var_data.ndim == 1:
                var_handle = file_handle.createVariable(
                    varname=file_var_name, dimensions=(dim_obs,), datatype='str')
            elif var_data.ndim == 2:
                var_handle = file_handle.createVariable(
                    varname=file_var_name, dimensions=(dim_locations, dim_obs), datatype='str')
            else:
                logging.error(' ===> Values dimensions is not expected')
                raise IOError('Dump datasets failed due to the values dimensions')

        elif isinstance(var_check, (float, np.floating)):

            var_data = var_values
            if var_data.ndim == 1:
                var_handle = file_handle.createVariable(
                    varname=file_var_name, dimensions=(dim_obs,), datatype='f4')
            elif var_data.ndim == 2:
                var_handle = file_handle.createVariable(
                    varname=file_var_name, dimensions=(dim_locations, dim_obs), datatype='f4')
            else:
                logging.error(' ===> Values dimensions is not expected')
                raise IOError('Dump datasets failed due to the values dimensions')

        elif isinstance(var_check, (int, np.integer)):

            var_data = var_values
            if var_data.ndim == 1:
                var_handle = file_handle.createVariable(
                    varname=file_var_name, dimensions=(dim_obs,), datatype='i4')
            elif var_data.ndim == 2:
                var_handle = file_handle.createVariable(
                    varname=file_var_name, dimensions=(dim_locations, dim_obs), datatype='i4')
            else:
                logging.error(' ===> Values dimensions is not expected')
                raise IOError('Dump datasets failed due to the values dimensions')

        elif isinstance(var_check, np.datetime64):

            var_tmp = pd.DatetimeIndex(var_values).to_pydatetime()
            var_data = netCDF4.date2num(var_tmp, units=file_time_units, calendar='gregorian')
            var_handle = file_handle.createVariable(varname=var_name, dimensions=(dim_obs,), datatype='f8')
            var_handle.units = file_time_units

        else:
            logging.error(' ===> Datasets format is not expected')
            raise IOError('Dump datasets failed due to the datasets format')

        if var_data.ndim == 1:
            var_handle[:] = var_data
        elif var_data.ndim == 2:
            var_handle[:, :] = var_data
        else:
            logging.error(' ===> Values dimensions is not expected')
            raise IOError('Dump datasets failed due to the values dimensions')

        if var_attrs is not None:
            var_handle.setncatts(var_attrs)

    # iterate over variable registry list
    for reg_name in file_var_list_registry:

        # get values and attributes
        reg_values = registry_obj[reg_name].values
        reg_attrs = None
        if reg_name in list(registry_attrs.keys()):
            reg_attrs = registry_attrs[reg_name]

        if isinstance(reg_values[0], str):
            reg_data = np.array(reg_values, dtype=object)
            reg_handle = file_handle.createVariable(varname=reg_name, dimensions=(dim_locations,), datatype='str')
        elif isinstance(reg_values[0], (float, np.floating)):
            reg_data = reg_values
            reg_handle = file_handle.createVariable(varname=reg_name, dimensions=(dim_locations,), datatype='f4')
        elif isinstance(reg_values[0], (int, np.integer)):
            reg_data = reg_values
            reg_handle = file_handle.createVariable(varname=reg_name, dimensions=(dim_locations,), datatype='i4')
        elif isinstance(reg_values[0], np.datetime64):
            reg_tmp = pd.DatetimeIndex(reg_values).to_pydatetime()
            reg_data = netCDF4.date2num(reg_tmp, units=file_time_units, calendar='gregorian')
            reg_handle = file_handle.createVariable(varname=reg_name, dimensions=(dim_locations,), datatype='f8')
            reg_handle.units = file_time_units
        else:
            logging.error(' ===> Datasets format is not expected')
            raise IOError('Dump datasets failed due to the datasets format')

        reg_handle[:] = reg_data

        if reg_attrs is not None:
            reg_handle.setncatts(reg_attrs)

    file_handle.close()

# ----------------------------------------------------------------------------------------------------------------------
