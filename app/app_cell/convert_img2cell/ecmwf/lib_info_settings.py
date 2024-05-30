"""
Library Features:

Name:          lib_info_settings
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230522'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import json
import logging
from copy import deepcopy
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get data by tag
def get_data_by_tag(data_settings, data_tag='log', data_default=None):
    if data_tag in list(data_settings.keys()):
        data_values = data_settings[data_tag]
    else:
        data_values = deepcopy(data_default)
    return data_values
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get data settings
def get_data_settings(file_name):
    if os.path.exists(file_name):
        with open(file_name) as file_handle:
            data_settings = json.load(file_handle)
    else:
        logging.error(' ===> Error in reading settings file "' + file_name + '"')
        raise IOError('File not found')
    return data_settings
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to parse data settings
def parse_data_settings(data_settings):

    # get product
    product_name = data_settings['product']['name']
    product_bbox = data_settings['product']['bbox']
    if 'image_buffer' in list(data_settings['product'].keys()):
        product_image_buffer = str(data_settings['product']['image_buffer'])
    else:
        product_image_buffer = str(4)
    if 'orientation' in list(data_settings['product'].keys()):
        product_orientation = str(data_settings['product']['orientation'])
    else:
        logging.error(' ===> Orientation not defined in settings. Impossible to set the grid orientation')
        raise NotImplementedError(
            'Set "orientation" parameter in settings file ["south-north-west-east", "north-south-west-east"]')
    # get flags
    reset_static = data_settings['flags']['reset_static']
    reset_ts = data_settings['flags']['reset_dynamic']
    # get info
    time_start = data_settings['time']['time_start']
    time_end = data_settings['time']['time_end']
    time_now = data_settings['time']['time_now']
    time_run = data_settings['time']['time_run']

    path_grid = data_settings['data']['path_grid']
    path_ts = data_settings['data']['path_ts']
    if 'path_stack' in list(data_settings['data'].keys()):
        path_stack = data_settings['data']['path_stack']
    else:
        path_stack = deepcopy(path_ts)

    grid_path = os.path.join(data_settings['grid']['folder_name'], data_settings['grid']['file_name'])

    file_name_src = data_settings['template']['file']['file_name_source']
    sub_path_src = data_settings['template']['time']['sub_path_source']
    datetime_src = data_settings['template']['time']['datetime_source']

    file_name_dst = data_settings['template']['file']['file_name_destination']
    sub_path_dst = data_settings['template']['time']['sub_path_destination']
    datetime_dst = data_settings['template']['time']['datetime_destination']

    if product_bbox is not None:
        geo_bbox = ''
        for item_bbox in product_bbox:
            geo_bbox += str(item_bbox)
            geo_bbox += ','
        geo_bbox = geo_bbox.strip(' ')
        geo_bbox = geo_bbox.strip(',')
    else:
        geo_bbox = ''

    # organize info
    product_args, geo_args, im_buffer_args, orient_args, path_args, time_args, tmpl_args_src, tmpl_args_dst = (
        [product_name], [geo_bbox], [product_image_buffer], [product_orientation],
        [path_grid, path_ts, path_stack], [time_start, time_end, time_run],
        [file_name_src, datetime_src, sub_path_src],
        [file_name_dst, datetime_dst, sub_path_dst],
    )
    flags_args = [str(reset_static), str(reset_ts)]
    grid_args = [grid_path]
    variables_args = data_settings['variables']

    # app settings
    app_settings = []
    app_settings.extend(product_args)
    app_settings.extend(geo_args)
    app_settings.extend(im_buffer_args)
    app_settings.extend(orient_args)
    app_settings.extend(flags_args)
    app_settings.extend(path_args)
    app_settings.extend(time_args)
    app_settings.extend(grid_args)
    app_settings.extend(tmpl_args_src)
    app_settings.extend(tmpl_args_dst)
    app_settings.extend(variables_args)

    return app_settings

# ----------------------------------------------------------------------------------------------------------------------
