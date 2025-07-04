"""
Library Features:

Name:          lib_utils_gzip
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210408'
Version:       '1.0.0'
"""
#################################################################################
# Library
import logging
import gzip
import os
import random
import string

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
#################################################################################


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

# --------------------------------------------------------------------------------
# Method to unzip file
def unzip_filename(file_name_zip, file_name_unzip):

    if os.path.exists(file_name_unzip): os.remove(file_name_unzip)

    file_handle_zip = gzip.GzipFile(file_name_zip, "rb")
    file_handle_unzip = open(file_name_unzip, "wb")

    file_data_unzip = file_handle_zip.read()
    file_handle_unzip.write(file_data_unzip)

    file_handle_zip.close()
    file_handle_unzip.close()
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Method to zip file
def zip_filename(file_name_unzip, file_name_zip):

    file_handle_unzip = open(file_name_unzip, 'rb')
    file_handle_zip = gzip.open(file_name_zip, 'wb')

    file_handle_zip.writelines(file_handle_unzip)

    file_handle_zip.close()
    file_handle_unzip.close()
# --------------------------------------------------------------------------------
