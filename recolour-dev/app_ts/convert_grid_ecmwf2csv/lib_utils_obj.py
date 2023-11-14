"""
Library Features:

Name:          lib_utils_obj
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""


# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import unicodedata
import re
import numpy as np
import pandas as pd
import xarray as xr
from copy import deepcopy

from lib_info_args import logger_name
from lib_info_args import (geo_coord_name_x, geo_coord_name_y, time_coord_name,
                           geo_dim_name_x, geo_dim_name_y, time_dim_name)

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to sanitize string
def sanitize_string(string_name):

    string_name = string_name.lower()
    string_name = re.sub(r"['.,-]", "", string_name)
    string_name = string_name.replace(' ', '')
    string_name = unicodedata.normalize('NFD', string_name).encode('ascii', 'ignore')
    string_name = string_name.decode("utf-8")

    return string_name
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to pad list collections
def pad_list(list_ref, list_in, default_value=None):
    list_ref_size = list_ref.__len__()
    list_in_size = list_in.__len__()

    list_out = deepcopy(list_in)
    if list_in_size < list_ref_size:
        list_out.extend([default_value] * (list_ref_size - len(list_in)))
        # alg_logger.warning(' ===> List is less than the reference size')
    elif list_in_size > list_ref_size:
        list_out = list_in[:list_ref_size]
        # alg_logger.warning(' ===> List is greater than the reference size')
    else:
        pass

    return list_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create a dictionary from list(s)
def create_dict_from_list(list_keys, list_values):

    if not isinstance(list_values, list):
        list_values = [list_values]

    list_values = pad_list(list_keys, list_values, default_value=list_values[0])
    obj_dict = {k: v for k, v in zip(list_keys, list_values)}

    return obj_dict
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert dataset to data array or data array collections
def convert_dataset_to_data_array(dset_obj):

    var_list_dset = list(dset_obj.data_vars)

    collections_obj = None
    if var_list_dset.__len__() == 1:
        collections_obj = dset_obj[var_list_dset[0]]
    else:
        for var_name_dset in var_list_dset:
            if collections_obj is None:
                collections_obj = {}
            da_obj = dset_obj[var_name_dset]
            collections_obj[var_name_dset] = da_obj

    return collections_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter dataset
def filter_dataset(dset_obj, dset_vars_filter=None):

    if dset_vars_filter is not None:

        var_list_out, var_list_in = list(dset_vars_filter.keys()), list(dset_vars_filter.values())
        var_list_dset = list(dset_obj.data_vars)

        dset_vars_map = None
        for var_name_dset in var_list_dset:

            if var_name_dset in var_list_in:
                var_idx_in = var_list_in.index(var_name_dset)
                var_name_out = var_list_out[var_idx_in]

                if dset_vars_map is None:
                    dset_vars_map = {}
                dset_vars_map[var_name_dset] = var_name_out

        dset_map_in = None
        if dset_vars_map is not None:
            dset_map_in = list(dset_vars_map.keys())

        if dset_map_in is not None:
            dset_select = dset_obj[dset_map_in]
            dset_select = dset_select.rename(dset_vars_map)
        else:
            dset_select = deepcopy(dset_obj)

    else:
        dset_select = deepcopy(dset_obj)

    return dset_select
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create dataset
def create_dataset(data_obj, data_geo_x, data_geo_y, data_time=None, data_attrs=None, common_attrs=None,
                   coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y, coord_name_time=time_coord_name,
                   dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y, dim_name_time=time_dim_name,
                   dims_order=None
                   ):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x]
    if data_time is not None:
        dims_order = [dim_name_y, dim_name_x, dim_name_time]

    if isinstance(data_time, pd.Timestamp):
        data_time = pd.DatetimeIndex([data_time])
    elif isinstance(data_time, pd.DatetimeIndex):
        pass
    else:
        log_stream.error(' ===> Time obj is in wrong format')
        raise IOError('Variable time format not valid')

    if data_geo_x.ndim == 2 and data_geo_y.ndim == 2:
        arr_geo_x, arr_geo_y = data_geo_x[0, :], data_geo_y[:, 0]
    elif data_geo_x.ndim == 1 and data_geo_y.ndim == 1:
        arr_geo_x, arr_geo_y = deepcopy(data_geo_x), deepcopy(data_geo_y)
    else:
        log_stream.error(' ===> Dimensions of geographical objects are not supported')
        raise NotImplemented('Case not implemented yet')

    var_dset = xr.Dataset(coords={coord_name_time: ([dim_name_time], data_time)})
    var_dset.coords[coord_name_time] = var_dset.coords[coord_name_time].astype('datetime64[ns]')
    if common_attrs is not None:
        var_dset.attrs = common_attrs

    for var_name, var_data in data_obj.items():

        var_attrs = None
        if data_attrs is not None:
            if var_name in list(data_attrs.keys()):
                var_attrs = data_attrs[var_name]

        if data_time is None:

            var_da = xr.DataArray(
                var_data, name=var_name,
                dims=dims_order,
                coords={coord_name_x: (dim_name_x, arr_geo_x),
                        coord_name_y: (dim_name_y, arr_geo_y)})

        elif isinstance(data_time, pd.DatetimeIndex):

            if var_data.shape.__len__() == 2:
                var_data = np.expand_dims(var_data, axis=-1)

            var_da = xr.DataArray(
                var_data, name=var_name,
                dims=dims_order,
                coords={coord_name_x: (dim_name_x, arr_geo_x),
                        coord_name_y: (dim_name_y, arr_geo_y),
                        coord_name_time: (dim_name_time, data_time)})
            if var_attrs is not None:
                var_da.attrs = var_attrs
        else:
            log_stream.error(' ===> Time obj is in wrong format')
            raise IOError('Variable time format not valid')

        var_dset[var_name] = var_da

    return var_dset
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create a data array
def create_darray_2d(data, geo_x, geo_y, geo_1d=True, time=None,
                     coord_name_x=geo_coord_name_x, coord_name_y=geo_coord_name_y, coord_name_time=time_coord_name,
                     dim_name_x=geo_dim_name_x, dim_name_y=geo_dim_name_y, dim_name_time=time_dim_name,
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x]
    if time is not None:
        dims_order = [dim_name_y, dim_name_x, dim_name_time]

    if geo_1d:
        if geo_x.shape.__len__() == 2:
            geo_x = geo_x[0, :]
        if geo_y.shape.__len__() == 2:
            geo_y = geo_y[:, 0]

        if time is None:
            data_da = xr.DataArray(data,
                                   dims=dims_order,
                                   coords={coord_name_x: (dim_name_x, geo_x),
                                           coord_name_y: (dim_name_y, geo_y)})
        elif isinstance(time, pd.DatetimeIndex):

            if data.shape.__len__() == 2:
                data = np.expand_dims(data, axis=-1)

            data_da = xr.DataArray(data,
                                   dims=dims_order,
                                   coords={coord_name_x: (dim_name_x, geo_x),
                                           coord_name_y: (dim_name_y, geo_y),
                                           coord_name_time: (dim_name_time, time)})
        else:
            log_stream.error(' ===> Time obj is in wrong format')
            raise IOError('Variable time format not valid')

    else:
        log_stream.error(' ===> Longitude and Latitude must be 1d')
        raise IOError('Variable shape is not valid')

    return data_da
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to map variable(s) in a dataframe object
def map_vars_dict(var_obj_in, var_map=None):

    if var_map is not None:
        # iterate over variable(s)
        var_obj_out = {}
        for var_name_in, var_dframe_tmp in var_obj_in.items():

            if var_name_in in list(var_map.keys()):
                var_idx = list(var_map.keys()).index(var_name_in)
                var_name_out = list(var_map.values())[var_idx]
                var_obj_out[var_name_out] = var_dframe_tmp
            else:
                log_stream.warning(' ===> Variable "' + var_name_in + '" not included in the dataframe obj')
    else:
        var_obj_out = deepcopy(var_obj_in)

    return var_obj_out
# ----------------------------------------------------------------------------------------------------------------------
