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

