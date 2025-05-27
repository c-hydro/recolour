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
# method to remap dataframe object
def remap_dframe(obj_df_in, obj_map_in=None):

    index_name_in = obj_df_in.index.name

    if index_name_in is None:
        log_stream.warning(' ===> Index name is defined by NoneType. Try to remap to "time"')
        index_name_in = 'time'

    if index_name_in in list(obj_map_in):
        index_name_out = obj_map_in[index_name_in]
        obj_df_in.index.name = index_name_out
    else:
        log_stream.warning(' ===> Index "' + index_name_in + '" not in map object. Fields keep the standard name')

    obj_map_out = {}
    for field_name_in, field_name_out in obj_map_in.items():
        if field_name_in != index_name_in:
            if field_name_in in list(obj_df_in.columns):
                obj_map_out[field_name_in] = field_name_out
            else:
                log_stream.warning(' ===> Column "' + field_name_in + '" not in the dataframe object.')

    obj_df_out = obj_df_in.rename(columns=obj_map_out)

    return obj_df_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remap dict object
def remap_dict(obj_dict_in, obj_map=None):

    obj_dict_out = {}
    for field_name_in, field_name_out in obj_map.items():
        if field_name_in in list(obj_dict_in.keys()):
            obj_dict_out[field_name_out] = obj_dict_in[field_name_in]
        else:
            obj_dict_out[field_name_in] = obj_dict_in[field_name_in]

    return obj_dict_out

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
# method to find common elements between list
def find_element_in_list(list_1, list_2, occurrence='first'):

    if not isinstance(list_1, list):
        list_1 = [list_1]
    if not isinstance(list_2, list):
        list_2 = [list_2]

    select_element = None
    for list_el in list_1:
        if list_el in list_2:
            if occurrence == 'first':
                return list_el
            elif occurrence == 'all':
                if select_element is None:
                    select_element = []
                select_element.append(list_el)
    if select_element is None:
        log_stream.warning(' ===> All elements are not available in the expected list')
        select_element = None

    return select_element
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
