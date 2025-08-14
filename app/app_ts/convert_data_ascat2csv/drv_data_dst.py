"""
Class Features

Name:          drv_data_dst
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20250813'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd

from copy import deepcopy

from lib_io_csv import write_csv
from lib_utils_io import compose_paths, merge_by_data, check_dict_keys, filter_dataframe
from lib_utils_time import get_latest_feasible_date, define_time_reference

from lib_utils_info import logger_name

# set logger
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# data class
class DrvData:

    def __init__(self, folder_name, file_name,
                 file_variable='ssm_{type}', file_index='time', file_geo_x='longitude', file_geo_y='latitude',
                 file_delimiter=';', file_format='csv',
                 time_start=None, time_end=None, frequency='D',
                 mapping=None):

        self.folder_name = folder_name
        self.file_name = file_name
        self.file_path = os.path.join(folder_name, file_name)

        self.file_delimiter = file_delimiter
        self.file_format = file_format

        self.file_variable = file_variable
        self.file_index = file_index
        self.file_geo_x = file_geo_x
        self.file_geo_y = file_geo_y

        self.time_start = time_start
        self.time_end = time_end
        self.frequency = frequency

        self.time_match = get_latest_feasible_date(self.time_start, self.time_end, frequency=self.frequency)
        self.time_str_reference = define_time_reference(
            time_stamp=self.time_match, time_frequency=self.frequency, time_format='%Y%m')
        self.time_stamp_reference = pd.to_datetime(self.time_str_reference, format='%Y%m')

        self.mapping = mapping

    # method to organize data
    def organize_data(self, workspace_csv: dict, workspace_tiff: dict, no_data: (float, int) =-9999):

        # info start method
        logger_stream.info(' ----> Organize dst dynamic datasets [' + self.file_format + '] ... ')

        # method to merge dataframes of CSV and TIFF files by time index
        workspace_merged = merge_by_data(workspace_csv, workspace_tiff, no_data=no_data)
        # method to check keys in the workspace(s)
        workspace_checks = check_dict_keys([workspace_merged, workspace_csv, workspace_tiff])

        # check workspace consistency
        if not workspace_checks:
            logger_stream.warning(" ===> Keys are not available in all workspaces. Some time-series could be skipped.")

        # info end method
        logger_stream.info(' ----> Organize dst dynamic datasets [' + self.file_format + '] ... DONE')

        return workspace_merged

    # method to dump data
    def dump_data(self, workspace_df: dict, registry_df: pd.DataFrame) -> dict:

        # info start method
        logger_stream.info(' ----> Dump dst dynamic datasets [' + self.file_format + '] ... ')

        # define file path(s)
        file_path_raw = compose_paths(self.file_path, path_time=self.time_stamp_reference)

        # iterate over the DataFrames in the workspace
        for key, df_raw in workspace_df.items():
            # check if the DataFrame is empty
            if df_raw.empty:
                logger_stream.warning(f" ===> Empty DataFrame for key: {key}")
                continue

            # transform the DataFrame to the expected format
            df_def = filter_dataframe(
                df_raw, rename_map=self.mapping['rename'], filter_map=self.mapping['filter'])

            # define file path
            file_path_def = deepcopy(file_path_raw.format(point_name=key))
            # create folder if it does not exist
            folder_name_def, file_name_def = os.path.split(file_path_def)
            os.makedirs(folder_name_def, exist_ok=True)

            # method to write data to CSV
            write_csv(file_path_def, df_def,
                      file_delimiter=self.file_delimiter, keep_header=True, index=True,
                      float_precision=2, encoding='utf-8')

        # info end method
        logger_stream.info(' ----> Dump dst dynamic datasets [' + self.file_format + '] ... DONE')
# ----------------------------------------------------------------------------------------------------------------------
