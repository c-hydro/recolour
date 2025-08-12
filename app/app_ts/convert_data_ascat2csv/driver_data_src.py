# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd
import xarray as xr

from lib_io_csv import read_csv, check_csv
from lib_io_tiff import read_tiff
from lib_utils_io import compose_paths, merge_by_rows, adapt_dataframe_to_range, merge_by_time
from lib_utils_decoretors import iterate_time_steps
from lib_utils_analysis import (search_values_points,  search_values_maps,
                                aggregate_values_ts_by_frequency, aggregate_values_maps_by_frequency)
from lib_utils_time import compute_time_by_labels

from lib_utils_info import logger_name

# set logger
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# data class
class DriverData:

    def __init__(self, folder_name, file_name,
                 file_variable='ssm', file_index='time', file_geo_x='longitude', file_geo_y='latitude',
                 file_delimiter=';', file_format='csv',
                 search_radius_km=12.5,  # Default search radius in kilometers
                 time_start=None, time_end=None, frequency='D', time_window=None):

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

        if time_window is None:
            self.labels = [""]
        else:
            self.labels = [time_window[k]["label"] for k in sorted(time_window.keys())]

        self.search_radius_m = search_radius_km * 1000  # Convert kilometers to meters

    # method to read data
    @iterate_time_steps
    def read_data(self, time_step=None):

        # define file path(s)
        file_path = compose_paths(self.file_path, path_time=time_step, path_labels=self.labels)
        # define file time(s)
        file_times = compute_time_by_labels(time_step, labels=self.labels, ref_point='end')

        # check file format
        if self.file_format == 'csv':
            # read data csv
            file_data = read_csv(
                file_path, file_delimiter=self.file_delimiter, )

            # check and prepare data
            file_data = check_csv(file_data)

            # merge data
            file_obj = merge_by_rows(obj=file_data, ignore_index=False)

        elif self.file_format == 'tiff' or self.file_format == 'tif':

            # read data tiff
            file_data = read_tiff(file_path, file_times, variable='ssm_filtered')
            # merge data
            file_obj = merge_by_time(file_data[0], file_data[1])

        else:
            raise ValueError(f"Unsupported file format: {self.file_format}")

        return file_obj

    # method to organize data maps
    def organize_data_maps(self, file_das: (xr.DataArray, list), registry_df: (dict, pd.DataFrame) ) -> pd.DataFrame:

        workspace_da_global = None
        for da in file_das:
            if workspace_da_global is None:
                workspace_da_global = da
            else:
                workspace_da_global = merge_by_time(workspace_da_global, da)

        workspace_df_selected = search_values_maps(
            obj=registry_df, source_da=workspace_da_global,
            radius_m=self.search_radius_m,  # Default radius in meters
            target_lat_col="latitude",
            target_lon_col="longitude",
            value_col=self.file_variable,
            lat_name="latitude",
            lon_name="longitude",
            time_name='time')

        workspace_df_aggregated = aggregate_values_maps_by_frequency(
            obj=workspace_df_selected,
            time_col='time', value_col=self.file_variable, freq=self.frequency, min_frac=0.75)

        return workspace_df_aggregated

    # method to organize data time-series
    def organize_data_ts(self, file_dfs: (pd.DataFrame, list) = None, registry_df: (dict, pd.DataFrame) = None) -> pd.DataFrame:

        # concatenate dfs to a single df
        workspace_df_global = None
        for df in file_dfs:
            if workspace_df_global is None:
                workspace_df_global = df
            else:
                workspace_df_global = pd.concat([workspace_df_global, df], ignore_index=False)

        # select values using a radius search
        workspace_df_selected = search_values_points(
            obj=registry_df, source_df=workspace_df_global,
            radius_m=self.search_radius_m,  # Default radius in meters
            target_lat_col="latitude",
            target_lon_col="longitude",
            source_lat_col="latitude",
            source_lon_col="longitude")

        workspace_df_aggregated = aggregate_values_ts_by_frequency(
            workspace_df_selected, time_col='time', value_cols=self.file_variable,  freq=self.frequency)

        return workspace_df_aggregated

    # method to sync data
    def sync_data(self, workspace_df_aggregated: pd.DataFrame ) -> dict:

        workspace_df_adapted = adapt_dataframe_to_range(
            workspace_df_aggregated, time_start=self.time_start, time_end=self.time_end, freq=self.frequency)

        return workspace_df_adapted

# ----------------------------------------------------------------------------------------------------------------------