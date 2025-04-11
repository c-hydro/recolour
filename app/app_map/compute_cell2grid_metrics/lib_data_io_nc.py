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

from lib_utils_grid import read_grid_file, find_grid_idx_to_data

# default netcdf encoded attributes
attrs_encoded = ["_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping"]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize nc file
def organize_file_nc(obj_variable, obj_time=None, obj_variable_in=None, obj_variable_out=None,
                     var_name_geo_x='longitude', var_name_geo_y='latitude',
                     coord_name_x='longitude', coord_name_y='latitude',
                     dim_name_x='longitude', dim_name_y='latitude'):

    if obj_variable_in is None:
        obj_variable_in = []
    if obj_variable_out is None:
        obj_variable_out = []

    # iterate over variable(s)
    variable_dset = None
    for variable_name_in, variable_da in obj_variable.items():

        if variable_name_in in obj_variable_in:
            variable_id = obj_variable_in.index(variable_name_in)
            variable_name_out = obj_variable_out[variable_id]
        else:
            variable_name_out = variable_name_in

        if variable_dset is None:

            var_geo_x_1d = variable_da[var_name_geo_x].values
            var_geo_y_1d = variable_da[var_name_geo_y].values

            variable_dset = xr.Dataset(
                coords={coord_name_x: ([dim_name_x], var_geo_x_1d),
                        coord_name_y: ([dim_name_y], var_geo_y_1d)})
                # ,coord_name_time: ([dim_name_time], var_data_time)})

        variable_dset[variable_name_out] = variable_da

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
def read_file_cell(file_name, selected_variables=None):

    variable_workspace = None
    with netCDF4.Dataset(file_name) as file_handle:

        file_variables = file_handle.variables
        for variable_tmp in selected_variables:
            if variable_tmp not in list(file_variables):
                logging.error(' ===> Variable "' + variable_tmp + '" is not found in the cell file "' + file_name + '"')
                raise RuntimeError('All the selected variables must be in the selected cell file')

        for variable_id, variable_name in enumerate(file_variables):

            if selected_variables is not None:
                if variable_name in selected_variables:
                    variable_tmp = file_variables[variable_name][:]
                else:
                    variable_tmp = None
                    # logging.warning(' ===> Variable "' + variable_name + '" not selected ih the cell datasets.')
            else:
                variable_tmp = file_variables[variable_name][:]

            if variable_tmp is not None:
                if isinstance(variable_tmp, np.ma.MaskedArray):
                    variable_data = variable_tmp.data
                else:
                    variable_data = deepcopy(variable_tmp)

                if variable_workspace is None:
                    variable_workspace = {}
                variable_workspace[variable_name] = variable_data
            else:
                # logging.warning(' ===> Variable "' + variable_name + '" is not saved in the workspace file')
                pass

    return variable_workspace
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file collection in netcdf format
def read_file_collection(file_name_data,
                         file_name_grid='TUW_WARP5_grid_info_2_3.nc', file_obj_grid=None,
                         variable_list=None,
                         variable_idx='gpi', variable_geo_x='lon', variable_geo_y='lat',
                         variable_type='ASCAT', variable_land_filter=True):

    # get reference grid
    if file_obj_grid is None:
        if os.path.exists(file_name_grid):
            grid_dframe = read_grid_file(file_name_grid)
        else:
            logging.error(' ===> File grid "' + file_name_grid + '" is not found!')
            raise FileNotFoundError('File must be available to correctly run the algorithm')
    else:
        grid_dframe = deepcopy(file_obj_grid)

    if variable_list is None:
        variable_list = ['gpi', 'lon', 'lat', 'xyz_x_snr', 'obs', 'xy_pr', 'xz_pr']

    with netCDF4.Dataset(file_name_data) as file_handle:

        variable_collection = {}
        variable_data_idx, variable_data_geo_x, variable_data_geo_y = None, None, None
        for variable_name in variable_list:

            variable_tmp = file_handle.variables[variable_name]
            variable_data = variable_tmp[:]

            if isinstance(variable_data, np.ma.MaskedArray):
                variable_data = variable_data.data

            if variable_name == variable_idx:
                variable_data_idx = deepcopy(variable_data)
            elif variable_name == variable_geo_x:
                variable_data_geo_x = deepcopy(variable_data)
                variable_name_geo_x = deepcopy(variable_name)
            elif variable_name == variable_geo_y:
                variable_data_geo_y = deepcopy(variable_data)
                variable_name_geo_y = deepcopy(variable_name)

            variable_collection[variable_name] = variable_data

    variable_dframe = pd.DataFrame(variable_collection, index=variable_data_idx)
    variable_dframe.replace(-999999, np.nan, inplace=True)

    # check variable for the different type(s)
    if variable_type == 'ASCAT':

        gpi_grid = grid_dframe['gpi'].values
        car_grid = grid_dframe['committed_area'].values
        land_grid = grid_dframe['land_flag'].values
        gpi_variable = variable_dframe.index.values

        idx_variable = find_grid_idx_to_data(gpi_grid, gpi_variable)

        car_data = car_grid[idx_variable]
        variable_dframe['committed_area'] = car_data
        land_data = land_grid[idx_variable]
        variable_dframe['land_flag'] = land_data

    elif variable_type == 'RZSM':
        variable_dframe['committed_area'] = 1
        variable_dframe['land_flag'] = 1
        variable_dframe['lon'] = variable_data_geo_x
        variable_dframe['lat'] = variable_data_geo_y
    else:
        logging.error(' ===> Datasets type is not expected ("ASCAT" or "RZSM") ')
        raise NotImplemented('Case not implemented yet')

    variable_dframe = variable_dframe.rename(columns={'lon': variable_name_geo_x, 'lat': variable_name_geo_y})

    if variable_land_filter:
        variable_dframe = variable_dframe[variable_dframe['land_flag'] == 1]

    return variable_dframe

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write file collection in netcdf format
def write_file_collection(file_name, file_data, file_tag_location='gpi'):

    # get data dimension
    file_n = file_data[file_tag_location].__len__()

    # open and init file
    file_handle = netCDF4.Dataset(file_name, 'w')
    file_handle.createDimension('data', file_n)
    # add attr(s)
    file_handle.file_date = 'Created ' + time.ctime(time.time())

    # iterate over datasets
    for file_key, file_dict in file_data.items():

        file_values = list(file_dict.data)
        if isinstance(file_values[0], str):
            file_data = np.array(file_values, dtype=object)
            file_var = file_handle.createVariable(varname=file_key, dimensions=('data',), datatype='str')
        elif isinstance(file_values[0], (int, float, np.integer, np.floating)):
            file_data = file_values
            file_var = file_handle.createVariable(varname=file_key, dimensions=('data',), datatype='f4')
        else:
            logging.error(' ===> Datasets format is not expected')
            raise IOError('Dump datasets failed due to the datasets format')

        file_var[:] = file_data

    file_handle.close()

# ----------------------------------------------------------------------------------------------------------------------
