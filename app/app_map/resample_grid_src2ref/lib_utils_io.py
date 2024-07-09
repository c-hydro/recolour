"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240709'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_utils_generic import fill_tags2string

from lib_info_args import logger_name

# logging
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to define file name
def get_file_info(file_settings, tag_folder_name='folder_name', tag_file_name='file_name'):

    # get folder name
    if tag_folder_name in list(file_settings.keys()):
        folder_name = file_settings[tag_folder_name]
    else:
        alg_logger.error(' ===> Folder name is not defined in file settings')
        raise IOError('Folder name must be defined')
    # get file name
    if tag_file_name in list(file_settings.keys()):
        file_name = file_settings[tag_file_name]
    else:
        alg_logger.error(' ===> File name is not defined in file settings')
        raise IOError('File name must be defined')
    # define file path
    file_path = os.path.join(folder_name, file_name)

    return file_path
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to set time info
def set_file_info(file_name_template, template_time, template_tags=None):

    # fill template values
    template_values = {}
    if template_tags is not None:
        for template_tag, template_value in template_tags.items():
            template_values[template_tag] = template_time

    # define file name
    file_name_def = fill_tags2string(file_name_template, template_tags, template_values)

    # create file folder
    folder_name_def, _ = os.path.split(file_name_def)
    os.makedirs(folder_name_def, exist_ok=True)

    return file_name_def

# ----------------------------------------------------------------------------------------------------------------------
