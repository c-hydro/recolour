# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd

from copy import deepcopy

from lib_io_csv import write_csv
from lib_utils_io import compose_paths, adapt_dataframe_to_range
from lib_utils_time import get_latest_feasible_date, define_time_reference

from lib_utils_info import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# data class
class DriverData:

    def __init__(self, folder_name, file_name,
                 file_variable='ssm_{type}', file_index='time', file_geo_x='longitude', file_geo_y='latitude',
                 file_delimiter=';', file_format='csv',
                 time_start=None, time_end=None, frequency='D'):

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

    # method to organize data
    def organize_data(self, workspace_csv: dict, workspace_tiff: dict, no_data: (float, int) =-9999):

        workspace_merged = {}
        for (csv_key, csv_df), (tiff_key, tiff_df) in zip(workspace_csv.items(), workspace_tiff.items()):

            # check if the keys match
            if csv_key != tiff_key:
                alg_logger.error(f" ===> Keys do not match: {csv_key} != {tiff_key}")
                raise ValueError(f"Keys must have the same name, but got: {csv_key} and {tiff_key}")
            else:
                merged_key = csv_key = tiff_key

            # check if the DataFrame is empty
            if csv_df.empty or tiff_df.empty:
                alg_logger.warning(f" ===> Empty DataFrame for key: {csv_key}")
                continue

            # merge by time index (outer join to keep all dates)
            merged_df = pd.merge(csv_df, tiff_df, left_index=True, right_index=True, how='outer')
            # fill NaN values (e.g., with 0)
            merged_df = merged_df.fillna(no_data)

            # add data to workspace
            workspace_merged[merged_key] = merged_df

        return workspace_merged

    # method to dump data
    def dump_data(self, workspace_df: dict, registry_df: pd.DataFrame ) -> dict:

        # define file path(s)
        file_path_raw = compose_paths(self.file_path, path_time=self.time_stamp_reference)

        # iterate over the DataFrames in the workspace
        for key, df_raw in workspace_df.items():
            # check if the DataFrame is empty
            if df_raw.empty:
                alg_logger.warning(f" ===> Empty DataFrame for key: {key}")
                continue

            # transform the DataFrame to the expected format
            df_def = df_raw[['ssm_mean', 'ssm_filtered_mean']].rename(
                columns={'ssm_mean': 'ssm_ts', 'ssm_filtered_mean': 'ssm_grid'}
            )

            # define file path
            file_path_def = deepcopy(file_path_raw.format(point_name=key))
            # create folder if it does not exist
            folder_name_def, file_name_def = os.path.split(file_path_def)
            os.makedirs(folder_name_def, exist_ok=True)

            # method to write data to CSV
            write_csv(file_path_def, df_def,
                      file_delimiter=self.file_delimiter, keep_header=True, index=True,
                      float_precision=2, encoding='utf-8')

# ----------------------------------------------------------------------------------------------------------------------
