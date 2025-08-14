"""
Class Features

Name:          drv_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20250813'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_io_csv import read_csv
from lib_utils_io import compose_paths

from lib_utils_info import logger_name

# set logger stream
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# geo class
class DrvGeo:

    def __init__(self, folder_name, file_name,
                 file_delimiter=';', file_index='tag', file_format='csv', labels=None):

        self.folder_name = folder_name
        self.file_name = file_name
        self.file_path = os.path.join(folder_name, file_name)

        self.file_format = file_format
        self.file_delimiter = file_delimiter
        self.file_index = file_index

        self.file_index = file_index
        self.labels = labels if labels is not None else [""]

    # method to organize data
    def organize_data(self, time_step=None):

        # info start method
        logger_stream.info(' ----> Organize geographical datasets ... ')

        # get file path based on setting information
        file_path = compose_paths(self.file_path, path_time=time_step, path_labels=self.labels)

        # Check if the file exists and read it
        if os.path.exists(file_path):
            if self.file_format == 'csv':
                file_df = read_csv(file_path, file_delimiter=self.file_delimiter)
            else:
                logger_stream.error(' ===> Unsupported file format: ' + self.file_format)
                raise ValueError(f"Only csv format is supported, but got {self.file_format}")
        else:
            logger_stream.error(' ===> File not found: ' + file_path)
            raise FileNotFoundError(f"Geographical file is mandatory")

        # Create a dict of DataFrames keyed by 'tag'
        file_dict = {tag: group for tag, group in file_df.groupby(self.file_index)}

        # info end method
        logger_stream.info(' ----> Organize geographical datasets ... DONE')

        return file_dict

# ----------------------------------------------------------------------------------------------------------------------