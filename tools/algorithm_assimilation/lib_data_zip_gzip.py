"""
Library Features:

Name:          lib_data_zip_gzip
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20200401'
Version:       '3.0.0'
"""
#################################################################################
# Library
import logging
import gzip
#################################################################################


# --------------------------------------------------------------------------------
# Method to unzip file
def unzip_filename(file_name_zip, file_name_unzip):

    file_handle_zip = gzip.GzipFile(file_name_zip, "rb")
    file_handle_unzip = open(file_name_unzip, "wb")

    file_data_unzip = file_handle_zip.read()
    file_handle_unzip.write(file_data_unzip)

    file_handle_zip.close()
    file_handle_unzip.close()

# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Method to zip file
def zip_filename(file_name_unzip, file_name_zip, file_compression_level=9):

    file_handle_unzip = open(file_name_unzip, 'rb')
    with gzip.GzipFile(filename=file_name_zip, mode='w', compresslevel=file_compression_level) as file_handle_zip:
        file_handle_zip.writelines(file_handle_unzip)
    file_handle_unzip.close()
# --------------------------------------------------------------------------------
