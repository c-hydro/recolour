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

from copy import deepcopy

from lib_utils_io import create_darray_2d

import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select variable from dataset
def select_variable(dset_obj: xr.Dataset,
                    variable_name: str, variable_mandatory: bool = True, variable_no_data=None):

    if variable_name in list(dset_obj.variables):
        variable_data = dset_obj[variable_name].values
        variable_data = np.squeeze(variable_data)
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
                  var_name_geo_x: str = 'lon', var_name_geo_y: str = 'lat',
                  var_name_k1: str = 'var40', var_name_k2: str = 'var41',
                  var_name_k3: str = 'var42', var_name_k4: str = 'var43',
                  var_name_k5: str = 'var200', var_name_upd=None):

    # get variables data
    var_data_geo_x = select_variable(dset_obj, var_name_geo_x, variable_mandatory=True)
    var_data_geo_y = select_variable(dset_obj, var_name_geo_y, variable_mandatory=True)
    var_data_k1 = select_variable(dset_obj, var_name_k1, variable_mandatory=True)
    var_data_k2 = select_variable(dset_obj, var_name_k2, variable_mandatory=True)
    var_data_k3 = select_variable(dset_obj, var_name_k3, variable_mandatory=True)
    var_data_k4 = select_variable(dset_obj, var_name_k4, variable_mandatory=True)
    var_data_k5 = select_variable(dset_obj, var_name_k5, variable_mandatory=False, variable_no_data=None)

    # organize variables obj
    var_obj_tmp = {var_name_geo_x: var_data_geo_x, var_name_geo_y: var_data_geo_y,
                   var_name_k1: var_data_k1, var_name_k2: var_data_k2,
                   var_name_k3: var_data_k3, var_name_k4: var_data_k4,
                   var_name_k5: var_data_k5}

    # update variable names (if requested)
    if var_name_upd is not None:
        var_obj_def = {}
        for var_key, var_data in var_obj_tmp.items():
            if var_key in list(var_name_upd.keys()):
                var_upd = var_name_upd[var_key]
                var_obj_def[var_upd] = var_data
            else:
                var_obj_def[var_key] = var_data
    else:
        var_obj_def = deepcopy(var_obj_tmp)

    return var_obj_def

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

        ''' debug
        plt.figure(); plt.imshow(var_data); plt.colorbar(); plt.show()
        '''

    return data_obj_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to check data
def check_data(obj_data, thr_data='all'):
    var_check_list = []
    for var_key, var_da in obj_data.items():
        if var_da is None:
            var_check_step = False
        else:
            var_check_step = True
        var_check_list.append(var_check_step)

    if thr_data == 'all':
        flag_data = all(var_check_list)
    elif thr_data == 'any':
        flag_data = any(var_check_list)
    else:
        logging.error(' ===> Check flag "' + thr_data + '" is not supported')
        raise NotImplemented('Case not implemented yet')
    return flag_data
# ----------------------------------------------------------------------------------------------------------------------
