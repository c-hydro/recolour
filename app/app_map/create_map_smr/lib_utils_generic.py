"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import shutil
from operator import is_not
from functools import partial
from datetime import datetime
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to check dictionary keys in other dictionaries
def check_key_between_dicts(dict_ref, dict_other_list):
    for key_ref, value_ref in dict_ref.items():
        for dict_other_step in dict_other_list:
            if key_ref not in list(dict_other_step.keys()):
                logging.error(' ===> Key ' + key_ref + ' not found in dictionary ' + dict_other_step)
                raise NotImplementedError('Check your dictionary keys')
    return True
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remove key with null value
def remove_key_with_null_value(dict_original):
    dict_filtered = {k: v for k, v in dict_original.items() if v is not None}
    dict_original.clear()
    dict_original.update(dict_filtered)
    return dict_original
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remove duplicates from list
def remove_duplicates_list(list_in):
    list_tmp = list(set(list_in))
    list_out = list(filter(partial(is_not, None), list_tmp))
    return list_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to combine variables in dictionary
def combine_vars2dict(var_pivot='name', var_dict=None):

    if var_pivot in list(var_dict.keys()):
        var_list = var_dict[var_pivot]
        var_dict.pop(var_pivot)
    else:
        logging.error(' ===> Variable pivot is not available in dictionary')
        raise NotImplementedError('Check your dictionary keys')

    var_parameters = {}
    for var_id, var_key in enumerate(var_list):

        if var_key not in list(var_parameters.keys()):
            var_parameters[var_key] = {}

        for var_name, var_values in var_dict.items():
            var_value = var_values[var_id]
            var_parameters[var_key][var_name] = var_value

    return var_parameters

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to pad list collections
def pad_list(list_ref, list_collections_in, default_value=None):
    list_ref_size = list_ref.__len__()
    list_collection_out = []
    for list_tmp in list_collections_in:
        list_tmp_size = list_tmp.__len__()

        if list_tmp_size < list_ref_size:
            list_tmp.extend([default_value] * (list_ref_size - len(list_tmp)))
            logging.warning(' ===> List is less than the reference size')
        elif list_tmp_size > list_ref_size:
            list_tmp = list_tmp[:list_ref_size]
            logging.warning(' ===> List is greater than the reference size')
        else:
            pass

        list_collection_out.append(list_tmp)

    return list_collection_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to add format(s) string (path or filename)
def fill_tags2string(string_raw, tags_format=None, tags_filling=None, tags_template='[TMPL_TAG_{:}]'):

    apply_tags = False
    if string_raw is not None:
        for tag in list(tags_format.keys()):
            if tag in string_raw:
                apply_tags = True
                break

    if apply_tags:

        string_filled = None
        tag_dictionary = {}
        for tag_id, (tag_key, tag_value) in enumerate(tags_format.items()):
            tag_key_tmp = '{' + tag_key + '}'
            if tag_value is not None:

                tag_id = tags_template.format(tag_id)
                tag_dictionary[tag_id] = {'key': None, 'value': None}

                if tag_key_tmp in string_raw:
                    tag_dictionary[tag_id] = {'key': tag_key, 'value': tag_value}
                    string_filled = string_raw.replace(tag_key_tmp, tag_id)
                    string_raw = string_filled
                else:
                    tag_dictionary[tag_id] = {'key': tag_key, 'value': None}

        dim_max = 1
        for tags_filling_values_tmp in tags_filling.values():
            if isinstance(tags_filling_values_tmp, list):
                dim_tmp = tags_filling_values_tmp.__len__()
                if dim_tmp > dim_max:
                    dim_max = dim_tmp

        string_filled_list = [string_filled] * dim_max

        string_filled_def = []
        for string_id, string_filled_step in enumerate(string_filled_list):

            for tag_dict_template, tag_dict_fields in tag_dictionary.items():
                tag_dict_key = tag_dict_fields['key']
                tag_dict_value = tag_dict_fields['value']

                if tag_dict_template in string_filled_step:
                    if tag_dict_value is not None:

                        if tag_dict_key in list(tags_filling.keys()):

                            value_filling_obj = tags_filling[tag_dict_key]

                            if isinstance(value_filling_obj, list):
                                value_filling = value_filling_obj[string_id]
                            else:
                                value_filling = value_filling_obj

                            tag_dict_key = '{' + tag_dict_key + '}'
                            string_filled_step = string_filled_step.replace(tag_dict_template, tag_dict_key)

                            if isinstance(value_filling, datetime):
                                tag_dict_value = value_filling.strftime(tag_dict_value)
                            elif isinstance(value_filling, (float, int)):
                                tag_dict_value = tag_dict_key.format(value_filling)
                            else:
                                tag_dict_value = value_filling

                            string_filled_step = string_filled_step.replace(tag_dict_key, tag_dict_value)

                        else:

                            # reverse the tag if not filled
                            tag_dict_key = '{' + tag_dict_key + '}'
                            string_filled_step = string_filled_step.replace(tag_dict_template, tag_dict_key)

            string_filled_def.append(string_filled_step)

        if dim_max == 1:
            string_filled_out = string_filled_def[0].replace('//', '/')
        else:
            string_filled_out = []
            for string_filled_tmp in string_filled_def:
                string_filled_out.append(string_filled_tmp.replace('//', '/'))

        return string_filled_out
    else:
        return string_raw
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to make folder
def make_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to reset folder
def reset_folder(folder_name, folder_reset=False, folder_make=True):
    if folder_reset:
        if os.path.exists(folder_name):
            shutil.rmtree(folder_name)
    if folder_make:
        os.makedirs(folder_name, exist_ok=True)
# ----------------------------------------------------------------------------------------------------------------------

