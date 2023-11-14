"""
Class Features

Name:          drv_map_fx
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import inspect

import lib_fx_methods as fx_methods

from lib_utils_io import create_darray_2d
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class driver fx
class DrvFx:

    # method to initialize class
    def __init__(self, fx_datasets, fx_name, fx_vars_in, fx_vars_out,
                 fx_var_geo_x='longitude', fx_var_geo_y='latitude'):

        self.fx_datasets = fx_datasets
        self.fx_name = fx_name
        self.fx_vars_in, self.fx_vars_out = fx_vars_in, fx_vars_out
        self.fx_var_geo_x, self.fx_var_geo_y = fx_var_geo_x, fx_var_geo_y

    # method to organize fx source
    def organize_fx_src(self):

        # info start method
        logging.info(' ----> Organize fx source... ')

        # method to map fx collection data
        fx_collections_data = map_fx_kwargs(self.fx_datasets, self.fx_vars_in)
        # method to select fx method
        fx_handle = select_fx_method(self.fx_name, fx_methods)
        # method to select fx kwargs data
        fx_kwargs_data = select_fx_kwargs(fx_handle, fx_collections_data)
        # method to select fx kwargs geo
        fx_kwargs_geo = select_fx_geo(self.fx_datasets, fx_var_geo_x=self.fx_var_geo_x, fx_var_geo_y=self.fx_var_geo_y)

        # info end method
        logging.info(' ----> Organize fx source ... DONE')

        return fx_handle, fx_kwargs_data, fx_kwargs_geo

    # method to execute fx
    @staticmethod
    def execute_fx(fx_handle, fx_kwargs):

        # info start method
        logging.info(' ----> Execute fx ... ')

        fx_result = fx_handle(**fx_kwargs)

        # info end method
        logging.info(' ----> Execute fx ... DONE')

        return fx_result

    # method to organize fx result
    def organize_fx_dst(self, fx_result_data, fx_result_geo):

        # info start method
        logging.info(' ----> Organize fx destination ... ')

        fx_data_results = select_fx_result(
            fx_result_data, fx_result_geo, self.fx_vars_out)

        logging.info(' ----> Organize fx destination ... DONE')

        return fx_data_results

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select fx results
def select_fx_result(fx_result_data, fx_result_geo, fx_variables):
    fx_result_output = {}
    for fx_data, (fx_var_key, fx_var_name) in zip(fx_result_data, fx_variables.items()):
        fx_result_output[fx_var_name] = fx_data
    fx_result_output = {**fx_result_output, **fx_result_geo}
    return fx_result_output
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select fx geo
def select_fx_geo(fx_datasets_in, fx_var_geo_x='longitude', fx_var_geo_y='latitude'):
    # initialize output object
    fx_datasets_geo = {}
    # check var geo x
    if fx_var_geo_x in list(fx_datasets_in.keys()):
        fx_datasets_geo[fx_var_geo_x] = fx_datasets_in[fx_var_geo_x]
    else:
        logging.error(' ===> Geographical variable fx "' + fx_var_geo_x + '" is not available in the datasets')
        raise RuntimeError('Geographical variable is needed by the fx to correctly run')
    # check var geo y
    if fx_var_geo_y in list(fx_datasets_in.keys()):
        fx_datasets_geo[fx_var_geo_y] = fx_datasets_in[fx_var_geo_y]
    else:
        logging.error(' ===> Geographical variable fx "' + fx_var_geo_y + '" is not available in the datasets')
        raise RuntimeError('Geographical variable is needed by the fx to correctly run')
    return fx_datasets_geo
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select fx args
def select_fx_kwargs(fx_handle, fx_args):
    fx_sign_list = list(inspect.signature(fx_handle).parameters.keys())
    fx_kwargs = {}
    for fx_sign_name in fx_sign_list:
        if fx_sign_name != 'kwargs':
            if fx_sign_name in list(fx_args.keys()):
                fx_kwargs[fx_sign_name] = fx_args[fx_sign_name]
            else:
                logging.error(' ===> Function argument "' + fx_sign_name + '" is not available in the function signature')
                raise RuntimeError('Argument is needed in the signature by the method')
    return fx_kwargs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select fx
def select_fx_method(fx_name, fx_methods):
    if hasattr(fx_methods, fx_name):
        fx_handle = getattr(fx_methods, fx_name)
    else:
        logging.error(' ===> Function "' + fx_name + '" is not available')
        raise RuntimeError('Function must be included in the methods repository')
    return fx_handle
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to map fx variable(s)
def map_fx_kwargs(fx_datasets_in, fx_variables_data):
    # initialize output object
    fx_datasets_out = {}
    # iterate over data variable(s)
    for var_name_fx, var_name_dataset in fx_variables_data.items():
        if var_name_dataset in list(fx_datasets_in.keys()):
            fx_datasets_out[var_name_fx] = fx_datasets_in[var_name_dataset]
        else:
            logging.error(' ===> Dataset variable fx "' + var_name_dataset + '" is not available in the datasets')
            raise RuntimeError('Dataset variable is needed by the fx to correctly run')
    return fx_datasets_out
# ----------------------------------------------------------------------------------------------------------------------
