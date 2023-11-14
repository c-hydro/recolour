"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import time
import re

import numpy as np
import pandas as pd

from copy import deepcopy
from datetime import datetime
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get geographical bounds
def get_geographical_bounds(lons, lats, decimal_precision=2):

    lon_west, lon_east = np.min(lons), np.max(lons)
    lat_south, lat_north = np.min(lats), np.max(lats)

    lon_west = str(round(lon_west, decimal_precision))
    lon_east = str(round(lon_east, decimal_precision))
    lat_south = str(round(lat_south, decimal_precision))
    lat_north = str(round(lat_north, decimal_precision))

    return lon_west, lon_east, lat_south, lat_north
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get time from filename
def extract_time_from_filename(file_name,
                               time_format='%Y%m%d%H%M', time_pattern='\d{4}\d{2}\d{2}\d{2}', time_sep='-'):
    time_match = re.search(pattern=time_pattern, string=file_name)

    if time_match is None:
        logging.error(' ===> Time pattern "' + time_pattern + '" did not match to the filename "' + file_name + '"')
        raise RuntimeError('Check the time pattern and the filename to pass the correct information')

    time_group = time_match.group()
    time_parts = time_group.split(time_sep)

    if time_parts.__len__() == 1:
        time_obj_start = datetime.strptime(time_parts[0], time_format)
        time_obj_end = deepcopy(time_obj_start)
    elif time_parts.__len__() == 2:
        time_obj_start = datetime.strptime(time_parts[0], time_format)
        time_obj_end = datetime.strptime(time_parts[1], time_format)
    else:
        logging.error(' ===> Time parts format is not supported')
        raise NotImplemented('Case not implemented yet')

    time_stamp_start = pd.Timestamp(time_obj_start)
    time_stamp_end = pd.Timestamp(time_obj_end)

    return time_stamp_start, time_stamp_end
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extrac part(s) from file path
def extract_parts_from_filename(file_path):
    folder_name, file_name = os.path.split(file_path)
    return folder_name, file_name
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get time xml
def get_time_xml(time_format="%Y-%m-%dT%H:%M:%S.%fZ"):
    time_string = datetime.utcnow().strftime(time_format)

    if 'Z' in time_format:
        time_string = time_string[:-3] + 'Z'

    return time_string
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get time product
def get_time_product(file_name, time_format="%Y-%m-%dT%H:%M:%S"):

    if not os.path.exists(file_name):
        logging.error(' ===> File "' + file_name + '" does not exist.')
        raise FileNotFoundError('File is mandatory')

    time_string_modif = time.ctime(os.path.getmtime(file_name))
    time_string_creation = time.ctime(os.path.getctime(file_name))

    time_stamp_modif = pd.Timestamp(time_string_modif)
    time_stamp_creation = pd.Timestamp(time_string_creation)

    time_string_modif = time_stamp_modif.strftime(time_format)
    time_string_creation = time_stamp_creation.strftime(time_format)

    if 'Z' in time_format:
        time_string_modif = time_string_modif[:-3] + 'Z'
    if 'Z' in time_format:
        time_string_creation = time_string_creation[:-3] + 'Z'

    return time_string_creation

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert time stamp to string
def convert_time_stamp_to_string(time_stamp, time_format="%Y-%m-%dT%H:%M:%S"):
    time_string = time_stamp.strftime(time_format)
    return time_string
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select the filename
def select_filename(file_name_selected, file_name_expected, flag_expected=True, flag_mandatory=True):

    if flag_expected:
        if file_name_selected != file_name_expected:
            if flag_mandatory:
                logging.error(' ===> File selected "' + file_name_selected + '" is not in the expected format "' +
                              file_name_expected + '"')
                raise RuntimeError('File selected must be equal to file expected')
            else:

                logging.warning(' ===> File selected "' + file_name_selected + '" is not in the expected format "' +
                                file_name_expected + '"')
                file_path_defined = deepcopy(file_name_expected)
        else:
            file_path_defined = deepcopy(file_name_selected)
    else:
        file_path_defined = deepcopy(file_name_selected)

    return file_path_defined
# ----------------------------------------------------------------------------------------------------------------------
