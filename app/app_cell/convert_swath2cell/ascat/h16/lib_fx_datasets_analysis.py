"""
Library Features:

Name:          lib_fx_datasets_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import numpy as np
import pandas as pd

from copy import deepcopy
from datetime import timedelta
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to build points dataframe
def build_points(data: dict) -> pd.DataFrame:
    """
    Given a dict mapping key → {arr, lon, lat, gpi_global, gpi_local},
    builds and returns a single DataFrame with all blocks concatenated,
    tagging each row with its source key.
    """

    # info algorithm start
    logging.info(' ----> Build point(s) ... ')

    #`iterate over the data dict, creating a DataFrame for each block
    dfs = []
    for key, block in data.items():

        # info point block start
        logging.info(' -----> Point block "' + str(key) + '" ... ' )

        # 1) DataFrame from structured array
        df = pd.DataFrame(block['data'])

        # 2) Add masked-array columns, filling masks with NaN
        df['lon'] = block['lon'].filled(np.nan)
        df['lat'] = block['lat'].filled(np.nan)
        df['gpi_global'] = block['gpi_global'].filled(np.nan)

        # 3) Add plain ndarray
        df['gpi_local'] = block['gpi_local']

        # 4) (Optional) convert flags to bool
        if 'corr_flag' in df: df['corr_flag'] = df['corr_flag'].astype(bool)
        if 'proc_flag' in df: df['proc_flag'] = df['proc_flag'].astype(bool)

        # 5) Tag with the id
        df['id'] = key

        df['datetime_n'] = pd.to_datetime(df['jd'], unit='D', origin='1858-11-17')
        df['datetime_str'] = df['datetime_n'].dt.strftime('%Y-%m-%d %H:%M:%S')

        dfs.append(df)

        # info point block end
        logging.info(' -----> Point block "' + str(key) + '" ... DONE')

    # info algorithm end
    logging.info(' ----> Build point(s) ... DONE')

    # 6) Concatenate all
    return pd.concat(dfs, ignore_index=True)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to save points dataframe
def save_points(df: pd.DataFrame, file_path: str, file_name: str, file_update=True) -> None:
    """
    Saves the DataFrame to a CSV file at `filepath`, using semicolon as delimiter,
    including the full header, without writing an extra index column.
    """

    formats = {
        'sm': '{:.3f}',
        'sm_noise': '{:.3f}',
        'dir': '{:.0f}',            # integer
        'gpi_global': '{:.0f}',     # integer
        # note: 'jd' is *not* here, so it stays as a numeric column
    }

    # info algorithm start
    logging.info(' ----> Save point(s) ... ')

    # format the DataFrame columns
    df = format_columns(df, formats)

    # make sure you have a datetime column
    df['datetime_n'] = pd.to_datetime(df['datetime_n'])

    # extract the date for grouping
    df['date'] = df['datetime_n'].dt.date

    # iterate over date(s)
    for single_date, day_group in df.groupby('date'):

        # info date start
        logging.info(' -----> Time "' + str(single_date) + '" ... ')

        # get hour
        hours = day_group['datetime_n'].dt.hour

        # day window: 00:00 ≤ h < 12:00
        df_day = day_group[(hours >= 0) & (hours < 12)]
        logging.info(' :: DAY POINTS ')
        logging.info(' :: Data: ' + str(len(df_day)))
        if not df_day.empty:

            single_date_day = deepcopy(single_date)

            sub_path_ts = single_date_day.strftime('%Y/%m/%d')
            day_file = single_date_day.strftime('%Y%m%d')

            out_dir = file_path.format(sub_path_ts=sub_path_ts)
            os.makedirs(out_dir, exist_ok=True)

            out_name = file_name.format(cell_day_time=day_file, cell_day_hour='1200')
            out_path = os.path.join(out_dir, out_name)

            if file_update and os.path.exists(out_path):
                logging.info(' :: File "' + out_path + '" already exists, removing it before writing new data')
                os.remove(out_path)

            # drop column 'datetime_n' in place
            df_day.drop('datetime_n', axis=1, inplace=True)
            # rename column in place
            df_day.rename(columns={'datetime_str': 'time'}, inplace=True)

            df_day.to_csv(out_path, sep=';', index=False, header=True)
            logging.info(' :: Wrote to "' + out_path + '"')
            logging.info(' :: ')

        # night window: 12:00 ≤ h or h < 00:00
        df_night = day_group[(hours >= 12) | (hours < 0)]

        logging.info(' :: NIGHT POINTS ')
        logging.info(' :: Data: ' + str(len(df_night)))
        if not df_night.empty:

            # ensure the date is incremented by one day for night data
            single_date_night = deepcopy(single_date + timedelta(days=1))

            sub_path_ts = single_date_night.strftime('%Y/%m/%d')
            day_file = single_date_night.strftime('%Y%m%d')

            out_dir = file_path.format(sub_path_ts=sub_path_ts)
            os.makedirs(out_dir, exist_ok=True)

            out_name = file_name.format(cell_day_time=day_file, cell_day_hour='0000')
            out_path = os.path.join(out_dir, out_name)

            if file_update and os.path.exists(out_path):
                logging.info(' :: File "' + out_path + '" already exists, removing it before writing new data')
                os.remove(out_path)

            # drop column 'datetime_n' in place
            df_night.drop('datetime_n', axis=1, inplace=True)
            # rename column in place
            df_night.rename(columns={'datetime_str': 'time'}, inplace=True)

            df_night.to_csv(out_path, sep=';', index=False, header=True)
            logging.info(' :: Wrote to "' + out_path + '"')

        # info date end
        logging.info(' -----> Time "' + str(single_date) + '" ... DONE')

    # info algorithm end
    logging.info(' ----> Save point(s) ... DONE')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to format columns in a DataFrame
def format_columns(df: pd.DataFrame, fmt_map: dict) -> pd.DataFrame:
    """
    Returns a copy of `df` where each column in fmt_map
    has been converted to strings via its format specifier.
    Columns *not* in fmt_map are left alone.
    """
    df2 = df.copy()
    for col, fmt in fmt_map.items():
        if col in df2.columns:
            df2[col] = df2[col].map(lambda x: fmt.format(x) if pd.notna(x) else '')
    return df2
# ----------------------------------------------------------------------------------------------------------------------
