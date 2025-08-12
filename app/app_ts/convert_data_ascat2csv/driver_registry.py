# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os

from lib_io_csv import read_csv
from lib_utils_io import compose_paths
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# registry class
class DriverRegistry:

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

    # method to read data
    def read_data(self, time_step=None):

        file_path = compose_paths(self.file_path, path_time=time_step, path_labels=self.labels)

        if os.path.exists(file_path):
            if self.file_format == 'csv':
                file_df = read_csv(file_path, file_delimiter=self.file_delimiter)
            else:
                raise ValueError(f"Unsupported file format: {self.file_format}")
        else:
            raise FileNotFoundError(f"File not found: {file_path}")

        # Create a dict of DataFrames keyed by 'tag'
        file_dict = {tag: group for tag, group in file_df.groupby(self.file_index)}

        return file_dict

# ----------------------------------------------------------------------------------------------------------------------