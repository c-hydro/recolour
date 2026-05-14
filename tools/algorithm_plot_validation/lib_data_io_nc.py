"""
Library Features:

Name:          lib_utils_io_nc
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260512'
Version:       '1.1.0
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import netCDF4
import time
import numpy as np
import pandas as pd
from pandas.api.types import is_integer_dtype, is_float_dtype, is_string_dtype, is_object_dtype

from copy import deepcopy

from lib_utils_grid import read_grid_file, find_grid_idx_to_data
# ----------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# method to read file cell in netcdf format
def read_file_cell(file_name, expected_variables=None):

    # initialize variable workspace
    variable_workspace = None

    # check expected variables
    if expected_variables is None or not expected_variables:
        logging.error(
            f' ===> Expected variables for selecting in cell files are not defined. Exit.')
        raise SystemExit

    # open the cell file
    with netCDF4.Dataset(file_name) as file_handle:

        # define file variables
        file_variables = list(file_handle.variables)

        # analyze found or missing variables
        variables_found = [var for var in expected_variables if var in file_variables]
        variables_missing = [var for var in expected_variables if var not in file_variables]

        # print info
        logging.info(f' ::: File variables: {file_variables}')
        logging.info(f' ::: Expected variables: {expected_variables}')
        logging.info(f' ::: Variables found: {variables_found}')
        if variables_missing:
            logging.warning(f' ::: Variables missing: {variables_missing}')
        else:
            logging.info(' ::: All expected variables were found')

        # check file variables
        if file_variables is None or not file_variables:
            logging.warning(
                f' ===> Expected variables for selecting in cell files are not defined. Return NoneType.')
            return None

        # iterate over expected variable(s)
        for variable_id, variable_name in enumerate(expected_variables):

            # initialize workspace
            if variable_workspace is None:
                variable_workspace = {}

            # check if expected variables is in the file variable list
            if variable_name in file_variables:

                try:
                    variable_tmp = file_handle.variables[variable_name][:]
                except Exception as exc:
                    logging.warning(
                        f' ===> Reading variable "{variable_name}": {exc}. File not correctly saved. Return NoneType')
                    variable_tmp = None
            else:
                logging.warning(
                    f' ===> Variable name "{variable_name}" not in {file_variables}. Return NoneType')
                variable_tmp = None

            # check if variables is defined or not
            if variable_tmp is not None:
                # get values from masked arrays (if defined by mask)
                if isinstance(variable_tmp, np.ma.MaskedArray):
                    variable_data = variable_tmp.data
                else:
                    variable_data = deepcopy(variable_tmp)

                # saved variable
                variable_workspace[variable_name] = variable_data
            else:
                # saved nonetype if variable are not available
                variable_workspace[variable_name] = None

    # check workspace content
    variables_none = [
        variable_name
        for variable_name, variable_data in variable_workspace.items()
        if variable_data is None
    ]

    if variables_none:
        logging.warning(
            f' ===> One or more variables are NoneType: {variables_none}. '
            f'Set variable workspace to NoneType.')

        variable_workspace = None

    return variable_workspace
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
    file_n = len(file_data[file_tag_location])

    # open and init file
    file_handle = netCDF4.Dataset(file_name, 'w')
    file_handle.createDimension('data', file_n)

    # add attr(s)
    file_handle.file_date = 'Created ' + time.ctime(time.time())

    try:

        # iterate over dataframe columns
        for file_key in file_data.columns:

            file_series = file_data[file_key]
            file_dtype = file_series.dtype

            # string / object variables
            if is_string_dtype(file_dtype) or is_object_dtype(file_dtype):

                file_values = file_series.astype(str).values

                file_var = file_handle.createVariable(
                    varname=file_key,
                    dimensions=('data',),
                    datatype=str
                )

                # VLEN strings must be assigned element by element
                for i, value in enumerate(file_values):
                    file_var[i] = value

            # integer variables
            elif is_integer_dtype(file_dtype):

                file_values = file_series.values.astype(np.int32)

                file_var = file_handle.createVariable(
                    varname=file_key,
                    dimensions=('data',),
                    datatype='i4'
                )

                file_var[:] = file_values

            # float variables
            elif is_float_dtype(file_dtype):

                file_values = file_series.values.astype(np.float32)

                file_var = file_handle.createVariable(
                    varname=file_key,
                    dimensions=('data',),
                    datatype='f4'
                )

                file_var[:] = file_values

            # unsupported variables
            else:
                logging.error(
                    f' ===> Datasets format for "{file_key}" '
                    f'is not expected: {file_dtype}'
                )
                raise IOError(
                    'Dump datasets failed due to the datasets format'
                )

    finally:
        file_handle.close()

# -----------------------------------------------------------------------------
