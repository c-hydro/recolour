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
import pandas as pd
import numpy as np
import xarray as xr

from copy import deepcopy

from lib_data_fx import compute_rzsm
from lib_utils_obj import create_darray_2d

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)

# debugging
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert dataset to rzsm
def convert_dataset_to_rzsm(data_obj_in,
                            var_name_k1_in='var_data_k1', var_name_k1_out='var_rzsm_k1',
                            var_name_k2_in='var_data_k2', var_name_k2_out='var_rzsm_k2',
                            var_name_k3_in='var_data_k3', var_name_k3_out='var_rzsm_k3',
                            var_name_geo_x='longitude', var_name_geo_y='latitude',
                            coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude',
                            dim_name_y='latitude', **kwargs):

    var_list_in = [var_name_k1_in, var_name_k2_in, var_name_k3_in]
    var_list_out = [var_name_k1_out, var_name_k2_out, var_name_k3_out]

    # get datasets
    var_obj_in, var_obj_geo_x, var_obj_geo_y = {}, None, None
    for var_name_in in var_list_in:
        var_data_k = data_obj_in[var_name_in]
        if isinstance(var_data_k, xr.DataArray):

            var_data_values = var_data_k.values
            var_data_values = np.squeeze(var_data_values)

            if var_obj_geo_x is None:
                var_obj_geo_x = var_data_k[var_name_geo_x].values
            if var_obj_geo_y is None:
                var_obj_geo_y = var_data_k[var_name_geo_y].values
        else:
            log_stream.error(' ===> Datasets object is not supported')
            raise NotImplemented('Case not implemented yet')

        var_obj_in[var_name_in] = var_data_values

    # compute datasets
    var_rzsm_k1, var_rzsm_k2, var_rzsm_k3 = compute_rzsm(**var_obj_in)
    var_rzsm_out = [var_rzsm_k1, var_rzsm_k2, var_rzsm_k3]

    # organize datasets
    data_obj_out = {}
    for var_name_out, var_data_out in zip(var_list_out, var_rzsm_out):

        # method to create data array
        var_da = create_darray_2d(
            var_data_out, var_obj_geo_x, var_obj_geo_y,
            coord_name_x=coord_name_x, coord_name_y=coord_name_y,
            dim_name_x=dim_name_x, dim_name_y=dim_name_y)

        # save in the output obj
        data_obj_out[var_name_out] = var_da

    return data_obj_out
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
# method to convert obj from variable(s) to point(s)
def convert_obj_vars_to_point(file_obj_by_vars, no_data=-9999.0):

    # iterate over variable(s) and point(s)
    file_obj_by_points = {}
    for var_name, var_data in file_obj_by_vars.items():
        point_list = list(var_data.columns)
        for point_name in point_list:

            point_data = var_data[point_name]
            point_data = point_data.fillna(no_data)
            point_data.name = var_name

            if point_name not in list(file_obj_by_points.keys()):
                file_obj_by_points[point_name] = {}
                file_obj_by_points[point_name] = point_data.to_frame()
            else:
                tmp_data = file_obj_by_points[point_name]
                point_data = pd.concat([tmp_data, point_data], axis=1)  # .reset_index()
                file_obj_by_points[point_name] = point_data

    return file_obj_by_points

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
