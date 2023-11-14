"""
Library Features:

Name:          lib_fx_utils
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import xarray as xr

from lib_utils_io import create_darray_2d
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select variable from dataset
def select_variable(dset_obj: xr.Dataset,
                    variable_name: str, variable_mandatory: bool = True, variable_no_data=None):

    if variable_name in list(dset_obj.variables):
        variable_data = dset_obj[variable_name].values
    else:
        if variable_mandatory:
            logging.error(' ===> Variable "' + variable_name + '" is not available in the datasets')
            raise RuntimeError('Variable is needed by the algorithm to correctly run')
        else:
            logging.warning(' ===> Variable "' + variable_name + '" is not available in the datasets')
            variable_data = variable_no_data
    return variable_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize data
def organize_data(dset_obj: xr.Dataset,
                  var_name_geo_x: str = 'longitude', var_name_geo_y: str = 'latitude',
                  var_name_k1: str = 'xyz_x_err_var', var_name_k2: str = 'xyz_y_err_var',
                  var_name_k3: str = 'xyz_z_err_var'):

    # get variables data
    var_data_geo_x = select_variable(dset_obj, var_name_geo_x, variable_mandatory=True)
    var_data_geo_y = select_variable(dset_obj, var_name_geo_y, variable_mandatory=True)
    var_data_k1 = select_variable(dset_obj, var_name_k1, variable_mandatory=True)
    var_data_k2 = select_variable(dset_obj, var_name_k2, variable_mandatory=True)
    var_data_k3 = select_variable(dset_obj, var_name_k3, variable_mandatory=False, variable_no_data=None)
    # organize variables obj
    var_obj = {var_name_geo_x: var_data_geo_x, var_name_geo_y: var_data_geo_y,
               var_name_k1: var_data_k1, var_name_k2: var_data_k2,
               var_name_k3: var_data_k3}

    return var_obj

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert data from arrays to data arrays
def convert_data(data_obj_in, var_name_geo_x='longitude', var_name_geo_y='latitude',
                 coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude'):

    # initialize output obj
    data_obj_out = {}
    # check geo x variable
    if var_name_geo_x in list(data_obj_in.keys()):
        var_data_geo_x = data_obj_in[var_name_geo_x]
        data_obj_in.pop(var_name_geo_x)
    else:
        logging.error(' ===> Geographical variable "' + var_name_geo_x + '" is not available in the datasets')
        raise RuntimeError('Geographical variable is needed by the method')
    # check geo y variable
    if var_name_geo_y in list(data_obj_in.keys()):
        var_data_geo_y = data_obj_in[var_name_geo_y]
        data_obj_in.pop(var_name_geo_y)
    else:
        logging.error(' ===> Geographical variable "' + var_name_geo_y + '" is not available in the datasets')
        raise RuntimeError('Geographical variable is needed by the method')

    # iterate over data variable(s)
    for var_name, var_data in data_obj_in.items():
        # method to create data array
        var_da = create_darray_2d(
            var_data, var_data_geo_x, var_data_geo_y, name=var_name,
            coord_name_x=coord_name_x, coord_name_y=coord_name_y,
            dim_name_x=dim_name_x, dim_name_y=dim_name_y)
        # save in the output obj
        data_obj_out[var_name] = var_da

    return data_obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to mask data
def mask_data(da_obj_in, da_reference, mask_value_min=0, mask_value_max=None, mask_no_data=np.nan,
              var_name_data='variable', var_name_geo_x='longitude', var_name_geo_y='latitude',
              coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude',
              ):

    data_values = da_obj_in.values
    geo_x_values = da_obj_in[var_name_geo_x].values
    geo_y_values = da_obj_in[var_name_geo_y].values
    mask_values = da_reference.values

    if mask_value_min is not None:
        data_values[mask_values < mask_value_min] = mask_no_data
    if mask_value_max is not None:
        data_values[mask_values > mask_value_max] = mask_no_data

    # method to create data array
    da_obj_out = create_darray_2d(
        data_values, geo_x_values, geo_y_values, name=var_name_data,
        coord_name_x=coord_name_x, coord_name_y=coord_name_y,
        dim_name_x=dim_name_x, dim_name_y=dim_name_y)

    return da_obj_out

# ----------------------------------------------------------------------------------------------------------------------
