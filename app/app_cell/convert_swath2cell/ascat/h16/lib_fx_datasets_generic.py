"""
Library Features:

Name:          lib_fx_datasets_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
from copy import deepcopy
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert sub path to root path
def convert_sub_path_2_root_path(sub_path_str, sub_path_delimiter='{'):
    sub_path_root = None
    if sub_path_str is not None:
        if sub_path_delimiter in sub_path_str:
            sub_path_root, _ = sub_path_str.split(sub_path_delimiter)
        else:
            sub_path_root = deepcopy(sub_path_str)
    return sub_path_root
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert template time string to time dictionary
def convert_sub_path_str_2_dict(sub_path_template_str='/%Y/%m/%d/%H/%M', sub_path_template_delimiter='/'):

    sub_path_template_dict = None
    if sub_path_template_str is not None:
        sub_path_template_list = sub_path_template_str.split(sub_path_template_delimiter)

        sub_path_template_dict = {}
        for sub_path_template_step in sub_path_template_list:

            sub_path_template_def = None
            if sub_path_template_step == '%Y' or sub_path_template_step == '%y':
                sub_path_template_def = {'years': '{year}'}
            elif sub_path_template_step == '%m':
                sub_path_template_def = {'months': '{month}'}
            elif sub_path_template_step == '%d':
                sub_path_template_def = {'days': '{day}'}
            elif sub_path_template_step == '%H':
                sub_path_template_def = {'hours': '{hour}'}
            elif sub_path_template_step == '%M':
                sub_path_template_def = {'minutes': '{minute}'}
            elif sub_path_template_step == '':
                pass
            else:
                logging.error(' ===> Sub path template "' + str(sub_path_template_step) + '" is not supported')
                raise IOError('Check the sub path template to correctly create the template obj')

            if sub_path_template_def is not None:
                sub_path_template_dict.update(sub_path_template_def)

    return sub_path_template_dict

# ----------------------------------------------------------------------------------------------------------------------

