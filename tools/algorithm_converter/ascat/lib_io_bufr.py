# libraries
import logging
import datetime
import geopandas as gpd
import numpy as np
import pandas as pd

from eccodes import (
    codes_bufr_new_from_file, codes_set, codes_release,
    codes_get_array, CodesInternalError
)

from lib_utils_time import datetime_to_jd

# method to decode BUFR files and extract soil moisture data
def decode_bufr(filepath):
    lats, lons = [], []
    ssm, flag_proc, flag_corr, quality = [], [], [], []
    dates, times, jd = [], [], []
    with open(filepath, 'rb') as f:
        while True:
            try:
                bufr = codes_bufr_new_from_file(f)
                if bufr is None:
                    break
                codes_set(bufr, "unpack", 1)
                try:
                    lat = codes_get_array(bufr, 'latitude').astype(np.float32)
                    lon = codes_get_array(bufr, 'longitude').astype(np.float32)
                    lat = np.clip(lat, -90.0, 90.0)
                    lon = np.clip(lon, -180.0, 180.0)

                    v_ssm = codes_get_array(bufr, 'surfaceSoilMoisture')
                    v_proc = codes_get_array(bufr, 'soilMoistureProcessingFlag')
                    v_corr = codes_get_array(bufr, 'soilMoistureCorrectionFlag')
                    v_qual = codes_get_array(bufr, 'soilMoistureQuality')

                    v_second = codes_get_array(bufr, 'second').astype(int)
                    v_minute = codes_get_array(bufr, 'minute').astype(int)
                    v_hour = codes_get_array(bufr, 'hour').astype(int)
                    v_day = codes_get_array(bufr, 'day').astype(int)
                    v_month = codes_get_array(bufr, 'month').astype(int)
                    v_year = codes_get_array(bufr, 'year').astype(int)

                    # Stack the scalars with the seconds array
                    v_dates = np.vstack((np.full_like(v_second, v_year),  # Year (same for all)
                                         np.full_like(v_second, v_month),  # Month (same for all)
                                         np.full_like(v_second, v_day),  # Day (same for all)
                                         np.full_like(v_second, v_hour),  # Hour (same for all)
                                         np.full_like(v_second, v_minute),  # Minute (same for all)
                                         v_second)).T

                    # Check if the lengths of v_dates and v_ssm match
                    if len(v_dates) != len(v_ssm):
                        # Log the mismatch and display the current lengths for tracking
                        logging.warning(
                            " ===> Mismatch in lengths of v_dates ({} entries) and v_ssm ({} entries). Trying to adapt dates to ssm.".format(
                                len(v_dates), len(v_ssm)))

                        # Case 1: v_dates has only 1 entry (single timestamp for all ssm values)
                        if len(v_dates) == 1:
                            logging.warning(
                                "  ::: v_dates has only 1 entry. Repeating v_dates for all {} entries in v_ssm.".format(
                                    v_ssm.shape[0]))
                            v_dates = np.repeat(v_dates, v_ssm.shape[0], axis=0)

                        # Case 2: v_dates is shorter than v_ssm (need to repeat dates)
                        elif len(v_dates) < len(v_ssm):
                            logging.warning(
                                " ::: v_dates is shorter than v_ssm. Tiling v_dates to match {} entries.".format(
                                v_ssm.shape[0]))
                            # Tile v_dates to match the number of ssm observations, repeating if necessary
                            v_dates = np.tile(v_dates, (v_ssm.shape[0] // len(v_dates) + 1, 1))[:v_ssm.shape[0]]

                        # Case 3: v_dates is longer than v_ssm (truncate v_dates)
                        elif len(v_dates) > len(v_ssm):
                            logging.warning(
                                " ::: v_dates is longer than v_ssm. Truncating v_dates to match {} entries.".format(
                                    v_ssm.shape[0]))
                            v_dates = v_dates[:v_ssm.shape[0]]

                    # Final check to ensure v_dates and v_ssm have the same length
                    if len(v_dates) == len(v_ssm):
                        logging.warning(" ::: v_dates successfully adjusted to match v_ssm with {} entries.".format(len(v_ssm)))
                        logging.warning(" ::: v_dates adjusted successfully to match v_ssm.")
                    else:
                        # If lengths still don't match, log the error and skip the BUFR message
                        logging.error(" ::: Failed to adjust v_dates to match v_ssm. Skipping this BUFR message.")
                        raise NotImplementedError(
                            'Case not implemented: v_dates and v_ssm lengths do not match after adjustment.')

                    # Construct the datetime strings in 'YYYY-MM-DD HH:MM:SS' format
                    v_time_strings = ['{}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}'.format(
                        v_dates[i, 0], v_dates[i, 1], v_dates[i, 2], v_dates[i, 3], v_dates[i, 4], v_dates[i, 5])
                        for i in range(len(v_dates))]

                    # Convert to pandas DateTimeIndex
                    v_time_index = pd.to_datetime(v_time_strings)

                    # Convert v_dates to Julian Dates
                    v_jd = np.array([datetime_to_jd(datetime.datetime(int(v_dates[i, 0]),
                                                                      int(v_dates[i, 1]),
                                                                      int(v_dates[i, 2]),
                                                                      int(v_dates[i, 3]),
                                                                      int(v_dates[i, 4]),
                                                                      int(v_dates[i, 5])))
                                     for i in range(len(v_dates))])

                    if all(len(x) == len(lat) for x in [lon, v_ssm, v_proc, v_corr, v_qual]):
                        lats.extend(lat)
                        lons.extend(lon)
                        ssm.extend(v_ssm)
                        flag_proc.extend(v_proc)
                        flag_corr.extend(v_corr)
                        quality.extend(v_qual)
                        dates.extend(v_dates)
                        times.extend(v_time_strings)
                        jd.extend(v_jd)
                except CodesInternalError:
                    pass
                codes_release(bufr)
            except Exception as e:
                logging.error(f"Error decoding BUFR: {e}")
                break
    workspace = {
        "longitude": np.array(lons),
        "latitude": np.array(lats),
        "ssm": np.array(ssm),
        "flag_processing": np.array(flag_proc),
        "flag_corrections": np.array(flag_corr),
        "quality": np.array(quality),
        "jd": np.array(jd),
        "time": pd.to_datetime(times),
        "date": np.array(dates)
    }

    return workspace
