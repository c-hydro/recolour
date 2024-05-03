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

from copy import deepcopy

from lib_utils_grid import read_grid_file, find_grid_idx_to_data
# ----------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to read file cell in netcdf format
def read_file_cell(file_name, file_variables=None):

    variable_workspace = None
    with netCDF4.Dataset(file_name) as file_handle:

        for variable_id, variable_name in enumerate(file_handle.variables):

            if file_variables is not None:
                if variable_name in file_variables:
                    variable_tmp = file_handle.variables[variable_name][:]
                else:
                    variable_tmp = None

            else:
                variable_tmp = file_handle.variables[variable_name][:]

            if variable_tmp is not None:
                if isinstance(variable_tmp, np.ma.MaskedArray):
                    variable_data = variable_tmp.data
                else:
                    variable_data = deepcopy(variable_tmp)

                if variable_workspace is None:
                    variable_workspace = {}
                variable_workspace[variable_name] = variable_data

    return variable_workspace
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to get variable name
def get_variable_name_OLD(variable_list_in, variable_list_out, variable_name='gpi'):

    variable_idx = variable_list_out.index(variable_name)
    variable_name = variable_list_in[variable_idx]

    return variable_name
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to read file collection in netcdf format
def read_file_collection(file_name_data, file_name_grid,
                         variable_list=None,
                         variable_idx='gpi', variable_geo_x='lon', variable_geo_y='lat',
                         variable_cell='cell',
                         variable_type='ASCAT', variable_land_filter=True):

    # get reference grid
    if os.path.exists(file_name_grid):
        grid_dframe = read_grid_file(file_name_grid)
    else:
        logging.error(' ===> File grid "' + file_name_grid + '" is not found!')
        raise FileNotFoundError('File must be available to correctly run the algorithm')

    if variable_list is None:
        variable_list = ['gpi', 'lon', 'lat', 'xyz_x_snr', 'obs', 'xy_pr', 'xz_pr', 'yz_pr']
    # add extra variable(s)
    variable_list = variable_list + ['cell']

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

        # add extra variable(s)
        if variable_cell in list(file_handle.variables):
            variable_tmp = np.array(file_handle.variables[variable_cell][:], dtype=int)
            variable_collection[variable_cell] = variable_tmp

    variable_dframe = pd.DataFrame(variable_collection, index=variable_data_idx)
    variable_dframe.replace(-999999, np.nan, inplace=True)

    # check variable name (in case of redefinition)
    if variable_type == 'RZSM':
        logging.warning(' ===> Variable type "RZSM" redefined by "ECMWF" tag. '
                        'Please, check the variable type in the configuration file!')
        variable_type = 'ECMWF'

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

    elif variable_type == 'ECMWF' or variable_type == 'HMC':

        gpi_grid = grid_dframe['gpi'].values
        car_grid = grid_dframe['committed_area'].values
        land_grid = grid_dframe['land_flag'].values
        gpi_variable = variable_dframe.index.values

        idx_variable = find_grid_idx_to_data(gpi_grid, gpi_variable)

        try:
            car_data = car_grid[idx_variable]
            variable_dframe['committed_area'] = car_data
            land_data = land_grid[idx_variable]
            variable_dframe['land_flag'] = land_data
        except BaseException as base_error:
            logging.warning(' ===> Variable "committed_area" and "land_flag" not found in grid file. '
                            'The variable(s) will be set to 1')
            variable_dframe['committed_area'] = 1
            variable_dframe['land_flag'] = 1

        variable_dframe['lon'] = variable_data_geo_x
        variable_dframe['lat'] = variable_data_geo_y
    else:
        logging.error(' ===> Datasets type is not expected ("ASCAT" or "RZSM" or "HMC") ')
        raise NotImplemented('Case not implemented yet')

    variable_dframe = variable_dframe.rename(columns={'lon': variable_name_geo_x, 'lat': variable_name_geo_y})

    if variable_land_filter:
        variable_dframe = variable_dframe[variable_dframe['land_flag'] == 1]

    return variable_dframe

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
