# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import sys
import pandas as pd

from lib_utils_decoretors import iterate_file_list, iterate_dict, iterate_items
from lib_utils_info import logger_name

# set logger
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to read a CSV file into a DataFrame
@iterate_file_list
def read_csv(file_name: str, file_delimiter=',') -> pd.DataFrame:

    # check if the file exists
    if os.path.exists(file_name):

        # Read the CSV into a DataFrame with a custom delimiter
        file_df = pd.read_csv(
            file_name,
            delimiter=file_delimiter,  # Change to your delimiter (e.g., ',', '\t', '|')
            encoding="utf-8",  # Optional: ensures proper text encoding
        )
        file_df.columns = file_df.columns.str.strip().str.replace(" ", "")

    else:
        logger_stream.warning(f" ===> File data csv not found: {file_name}")
        file_df = None

    return file_df
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to check and prepare a DataFrame from a CSV file
@iterate_items(iter_types=(list, tuple), strict_zip=True)
def check_csv(file_df = None, file_variable='ssm',
             file_geo_x='latitude', file_geo_y='latitude', file_index='time'):

    # check if the DataFrame is provided
    if file_df is None:
        return None

    # check required columns
    required = [file_variable, file_geo_x, file_geo_y]
    missing = [c for c in required if c not in file_df.columns]
    if missing:
        logger_stream.warning(f" ===> Missing required columns: {missing}. Available: {list(file_df.columns)}")

    # check if the index is already set
    if ((file_df.index.name == file_index) or
            (hasattr(file_df.index, "names") and file_index in (file_df.index.names or []))):
        file_df = file_df.reset_index()
    # check if the index is a column
    if file_index not in file_df.columns:
        logger_stream.error(f" ===> Missing required column: {file_index}")
        raise KeyError(f"'{file_index}' not found as column or index.")

    # convert time column to datetime if not already
    file_df[file_index] = pd.to_datetime(file_df[file_index], errors="raise")
    # set time as the index
    file_df = file_df.set_index(file_index)
    # sort by time
    file_df = file_df.sort_index()

    return file_df

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to write a DataFrame to a CSV file
def write_csv(
    file_name: str,
    file_data: pd.DataFrame,
    file_delimiter: str = ",",
    keep_header: bool = True,
    encoding: str = "utf-8",
    index: bool = False,
    float_precision: int = 2
) -> None:
    """
    Writes a pandas DataFrame to a CSV file with customizable options.

    Parameters:
        file_name (str): Path and name of the CSV file to create.
        file_data (pd.DataFrame): DataFrame containing the data.
        file_delimiter (str): Delimiter used in CSV (default: ',').
        keep_header (bool): If False, the column names will not be written.
        encoding (str): File encoding (default: 'utf-8').
        index (bool): Whether to include DataFrame's index in the file (default: False).

    Returns:
        None
    """
    if not isinstance(file_data, pd.DataFrame):
        logger_stream.error(" ===> The file_data must be a pandas DataFrame.")
        raise NotImplementedError("Case not implemented yet")

    if file_data.empty:
        logger_stream.warning(" ===> The file_data cannot be empty.")
        return None

    try:
        file_data.to_csv(
            file_name,
            sep=file_delimiter,
            header=keep_header,
            encoding=encoding,
            index=index,
            float_format=f"%.{float_precision}f"
        )

    except Exception as e:

        logger_stream.error(f" ===> CSV file '{file_name}' written failed.")
        raise IOError(f"An error occurred while writing the file: {e}")
# ----------------------------------------------------------------------------------------------------------------------
