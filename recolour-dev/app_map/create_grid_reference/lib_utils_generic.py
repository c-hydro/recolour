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
import tempfile
import shutil
# ----------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create a tmp name
def create_filename_tmp(prefix='tmp_', suffix='.tiff', folder=None):

    if folder is None:
        folder = '/tmp'

    with tempfile.NamedTemporaryFile(dir=folder, prefix=prefix, suffix=suffix, delete=False) as tmp:
        temp_file_name = tmp.name
    return temp_file_name
# -------------------------------------------------------------------------------------


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

