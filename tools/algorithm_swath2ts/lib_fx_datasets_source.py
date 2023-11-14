"""
Library Features:

Name:          lib_fx_datasets_source
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from datetime import datetime
from datetime import timedelta
from copy import deepcopy

import numpy as np
import xarray as xr
import pandas as pd

from ascat.eumetsat.level2 import AscatL2File
from ascat.read_native.cdr import AscatGriddedNcTs
from ascat.file_handling import ChronFiles

from pygeobase.utils import split_daterange_in_intervals
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to handle eps bufr file
class AscatEpsBufrFileList(ChronFiles):

    """
    Class reading ASCAT Eps files.
    """

    def __init__(self, root_path, sat_name, product_id,
                 subfolder_template=None, filename_template=None,
                 filetime_template='%Y%m%d%H%M%S', filetime_locs=[16, 30],
                 chunk_minutes=100, var_ref='sm', var_dim='obs',
                 compression=True
                 ):

        """
        Initialize.
        """

        sat_lut = {'a': 2, 'b': 1, 'c': 3, '?': '?'}

        self.sat_id = sat_lut[sat_name.lower()]
        self.product_id = product_id

        self.filetime_locs = filetime_locs
        self.filetime_template = filetime_template
        self.chunk_minutes = chunk_minutes

        self.sf_tmpl_fmt = subfolder_template

        self.var_ref = var_ref
        self.var_dim = var_dim

        self.var_to_rename = {'sat_track_azi': 'dir'}
        self.var_to_keep = ['sm', 'sm_noise', 'lon', 'lat', 'dir', 'corr_flag', 'proc_flag', 'jd']

        # Example file: ASCA_SMR_02_M03_20190531233300Z_20190601011458Z_N_C_20190601011516Z.nat.gz
        if filename_template is None:
            if compression:
                filename_template = 'ASCA_{product_id}_02_M0{sat_id}_{date}Z_*_*_*_*.nat.gz'
            else:
                filename_template = 'ASCA_{product_id}_02_M0{sat_id}_{date}Z_*_*_*_*.nat'

        super().__init__(root_path, AscatL2File,
                         filename_template, sf_templ=subfolder_template)

    def tstamps_for_daterange(self, dt_start, dt_end, dt_buffer=timedelta(minutes=30), dt_delta=timedelta(days=1)):

        file_names, file_timestamps = [], []

        dt_buffer_totalseconds = dt_buffer.total_seconds()
        dt_buffer = timedelta(seconds=dt_buffer_totalseconds)

        dt_start = dt_start - dt_buffer
        dt_end = dt_end + dt_buffer

        ts_start, ts_end = pd.Timestamp(dt_start), pd.Timestamp(dt_end)
        dt_delta_ts_totalseconds = (ts_end - ts_start).total_seconds()
        dt_delta_arg_totalseconds = dt_delta.total_seconds()

        if dt_delta_ts_totalseconds < dt_delta_arg_totalseconds:
            dt_delta = timedelta(seconds=dt_delta_ts_totalseconds)
        else:
            dt_delta = timedelta(seconds=dt_delta_arg_totalseconds)

        dt_range = np.arange(dt_start, dt_end, dt_delta)
        for dt_cur in dt_range.astype(datetime):
            files, dates = self.search_date(dt_cur, return_date=True)
            for f, dt in zip(files, dates):
                if (f not in file_names) and (dt >= dt_start and dt <= dt_end):
                    file_names.append(f)
                    file_timestamps.append(dt)

        file_time_intervals = split_daterange_in_intervals(file_timestamps[0], file_timestamps[-1],
                                                           self.chunk_minutes)

        return file_timestamps, file_time_intervals

    def search_period(self, dt_start, dt_end, dt_delta=timedelta(days=1)):
        filenames = []
        for dt_cur in np.arange(dt_start, dt_end, dt_delta).astype(datetime):
            files, dates = self.search_date(dt_cur, return_date=True)
            for f, dt in zip(files, dates):
                if f not in filenames and dt >= dt_start and dt <= dt_end:
                    filenames.append(f)
        return filenames

    def read_period(self, dt_start, dt_end, dt_delta=timedelta(days=1),
                    dt_buffer=timedelta(days=1), **kwargs):

        filenames = self.search_period(dt_start-dt_buffer, dt_end, dt_delta)
        data = []
        for filename in filenames:
            self._open(filename)
            d = self.fid.read(None, None, **kwargs)
            if d is not None:
                data.append(d)

        if data:
            data = self._merge_data(data)

        return data

    def resample_image(self, array, index, distance, weights):

        total_weights = np.nansum(weights, axis=1)

        resOrbit = {}
        # resample backscatter
        for name in array.dtype.names:
            if name != 'time':
                if name in ['corr_flag', 'proc_flag']:
                    # The flags are resampled by taking the minimum flag This works
                    # since any totally valid observation has the flag 0 and
                    # overrides the flagged observations. This is true in cases
                    # where the data was set to NaN by the flag as well as when the
                    # data was set to 0 or 100. The last image element is the one
                    # standing for NaN so we fill it with all flags filled to not
                    # interfere with the minimum.
                    array[name][-1] = 255
                    bits = np.unpackbits(array[name].reshape(
                        (-1, 1)).astype(np.uint8), axis=1)
                    resampled_bits = np.min(bits[index, :], axis=1)
                    resOrbit[name] = np.packbits(resampled_bits)
                else:
                    resOrbit[name] = np.nansum(
                        array[name][index] * weights, axis=1) / total_weights

        return resOrbit

    def flush(self):
        pass

    def close(self):
        pass


    def _fmt(self, timestamp):
        """
        Definition of filename and subfolder format.

        Parameters
        ----------
        timestamp : datetime
            Time stamp.

        Returns
        -------
        fn_fmt : dict
            Filename format.
        sf_fmt : dict
            Subfolder format.
        """
        fn_read_fmt = {'date': timestamp.strftime(self.filetime_template),
                       'sat_id': self.sat_id, 'product_id': self.product_id.upper()}
        fn_write_fmt = None
        sf_read_fmt = None
        sf_write_fmt = sf_read_fmt

        return fn_read_fmt, sf_read_fmt, fn_write_fmt, sf_write_fmt

    def _parse_date(self, filename):
        """
        Parse date from filename.

        Parameters
        ----------
        filename : str
            Filename.

        Returns
        -------
        date : datetime
            Parsed date.
        """
        return datetime.strptime(os.path.basename(filename)[self.filetime_locs[0]:self.filetime_locs[1]],
                                 '%Y%m%d%H%M%S')

    def _merge_data(self, data):
        """
        Merge data.

        Parameters
        ----------
        data : list
            List of array.

        Returns
        -------
        data : numpy.ndarray
            Data.
        """

        if not isinstance(data, list):
            raise IOError('merge data accepted only list')

        data_merged = None
        if isinstance(data[0], tuple):

            if isinstance(data[0][0], xr.Dataset):
                for element in data:

                    data_step = element[0]
                    data_valid = data_step.where(data_step[self.var_ref] >= 0, drop=True)
                    # idx_valid = list(np.where(np.isfinite(data_step['sm']))[0])

                    if data_valid[self.var_dim].size > 0:

                        if data_merged is None:
                            data_merged = deepcopy(data_valid)
                        else:
                            data_merged = xr.concat([data_merged, data_valid], dim=self.var_dim)
            else:
                metadata = [element[1] for element in data]
                data_merged = np.hstack([element[0] for element in data])
                data_merged = (data_merged, metadata)
        else:
            data_merged = np.hstack(data)

        if isinstance(data_merged, xr.Dataset):
            data_merged['jd'] = pd.DatetimeIndex(data_merged['time'].values).to_julian_date().values

        if isinstance(data_merged, xr.Dataset) and (self.var_to_rename is not None):
            data_merged = data_merged.rename_vars(self.var_to_rename)

        if isinstance(data_merged, xr.Dataset) and (self.var_to_keep is not None):
            data_merged = data_merged[self.var_to_keep]

        return data_merged
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to handle nrt bufr file
class AscatNrtBufrFileList(ChronFiles):

    """
    Class reading ASCAT NRT BUFR files.
    """

    def __init__(self, root_path, product_id='*', filename_template=None,
                 subfolder_template=None,
                 filetime_template='%Y%m%d_%H%M%S', filetime_locs=[4, 19],
                 chunk_minutes=50, var_ref='sm', var_dim='obs'
                 ):
        """
        Initialize.
        """

        self.product_id = product_id
        self.sf_tmpl_fmt = subfolder_template

        self.filetime_locs = filetime_locs
        self.filetime_template = filetime_template
        self.chunk_minutes = chunk_minutes

        self.var_ref = var_ref
        self.var_dim = var_dim

        self.var_to_rename = {'sat_track_azi': 'dir'}
        self.var_to_keep = ['sm', 'sm_noise', 'lon', 'lat', 'dir', 'corr_flag', 'proc_flag', 'jd']

        if filename_template is None:
            filename_template = '{product_id}_{date}*.buf'

        super().__init__(root_path, AscatL2File, filename_template,
                         sf_templ=subfolder_template)

    def tstamps_for_daterange(self, dt_start, dt_end, dt_delta=timedelta(days=1)):

        file_names, file_timestamps = [], []

        ts_start, ts_end = pd.Timestamp(dt_start), pd.Timestamp(dt_end)
        ts_totalseconds = (ts_end - ts_start).total_seconds()
        dt_totalseconds = dt_delta.total_seconds()

        if ts_totalseconds < dt_totalseconds:
            dt_delta = timedelta(seconds=ts_totalseconds)
        else:
            dt_delta = timedelta(seconds=dt_totalseconds)

        dt_range = np.arange(dt_start, dt_end + dt_delta, dt_delta)

        if ts_end < dt_range[-1]:
            dt_range[-1] = ts_end

        file_time_stamps, file_time_intervals = None, None
        for dt_cur in dt_range.astype(datetime):
            files, dates = self.search_date(dt_cur, return_date=True)
            for f, dt in zip(files, dates):
                if (files not in file_names) and (dt >= dt_start and dt <= dt_end):

                    if file_time_stamps is None:
                        file_time_stamps = []

                    file_names.append(f)
                    file_time_stamps.append(dt)

        if file_time_stamps is not None:
            file_time_intervals = split_daterange_in_intervals(file_time_stamps[0], file_time_stamps[-1],
                                                               self.chunk_minutes)
        else:
            logging.warning(' ===> File list selected is empty. Check your time information')

        return file_time_stamps, file_time_intervals

    def search_period(self, dt_start, dt_end, dt_delta=timedelta(days=1)):
        filenames = []
        for dt_cur in np.arange(dt_start, dt_end, dt_delta).astype(datetime):
            files, dates = self.search_date(dt_cur, return_date=True)
            for f, dt in zip(files, dates):
                if f not in filenames and dt >= dt_start and dt <= dt_end:
                    filenames.append(f)
        return filenames

    def resample_image(self, array, index, distance, weights):

        total_weights = np.nansum(weights, axis=1)

        resOrbit = {}
        # resample backscatter
        for name in array.dtype.names:
            if name != 'time':
                if name in ['corr_flag', 'proc_flag']:
                    # The flags are resampled by taking the minimum flag This works
                    # since any totally valid observation has the flag 0 and
                    # overrides the flagged observations. This is true in cases
                    # where the data was set to NaN by the flag as well as when the
                    # data was set to 0 or 100. The last image element is the one
                    # standing for NaN so we fill it with all flags filled to not
                    # interfere with the minimum.
                    array[name][-1] = 255
                    bits = np.unpackbits(array[name].reshape(
                        (-1, 1)).astype(np.uint8), axis=1)
                    resampled_bits = np.min(bits[index, :], axis=1)
                    resOrbit[name] = np.packbits(resampled_bits)
                else:
                    resOrbit[name] = np.nansum(
                        array[name][index] * weights, axis=1) / total_weights

        return resOrbit

    def flush(self):
        pass

    def close(self):
        pass

    def _fmt(self, timestamp):
        """
        Definition of filename and subfolder format.

        Parameters
        ----------
        timestamp : datetime
            Time stamp.

        Returns
        -------
        fn_fmt : dict
            Filename format.
        sf_fmt : dict
            Subfolder format.
        """
        fn_read_fmt = {'date': timestamp.strftime(self.filetime_template),
                       'product_id': self.product_id}

        sf_tmp_fmt = {'years': {'year': timestamp.strftime('%Y')},
                      'months': {'month': timestamp.strftime('%m')},
                      'days': {'day': timestamp.strftime('%d')},
                      'hours': {'hour': timestamp.strftime('%H')}}

        sf_step_fmt = None
        if self.sf_tmpl_fmt is not None:
            for tmpl_key, tmpl_value in self.sf_tmpl_fmt.items():
                if tmpl_key in list(sf_tmp_fmt.keys()):
                    if sf_step_fmt is None:
                        sf_step_fmt = {}
                    sf_step_fmt[tmpl_key] = sf_tmp_fmt[tmpl_key]

        sf_read_fmt = sf_step_fmt
        fn_write_fmt = None
        sf_write_fmt = None

        return fn_read_fmt, sf_read_fmt, fn_write_fmt, sf_write_fmt

    def _parse_date(self, filename):
        """
        Parse date from filename.

        Parameters
        ----------
        filename : str
            Filename.

        Returns
        -------
        date : datetime
            Parsed date.
        """

        filedate = os.path.basename(filename)[self.filetime_locs[0]:self.filetime_locs[1]].replace('_', '')

        return datetime.strptime(filedate, '%Y%m%d%H%M%S')

    def _merge_data(self, data):
        """
        Merge data.

        Parameters
        ----------
        data : list
            List of array.

        Returns
        -------
        data : numpy.ndarray
            Data.
        """

        if not isinstance(data, list):
            raise IOError('merge data accepted only list')

        if isinstance(data[0], tuple):

            data_merged = None
            for obj_step in data:

                data_step = obj_step[0]
                data_valid = data_step.where(np.isfinite(data_step[self.var_ref]), drop=True)
                # idx_valid = list(np.where(np.isfinite(data_step['sm']))[0])

                if data_valid[self.var_dim].size > 0:
                    if data_merged is None:
                        data_merged = deepcopy(data_valid)
                    else:
                        data_merged = xr.concat([data_merged, data_valid], dim=self.var_dim)
        else:
            data_merged = np.hstack(data)

        if isinstance(data_merged, xr.Dataset):
            data_merged['jd'] = pd.DatetimeIndex(data_merged['time'].values).to_julian_date().values

        if isinstance(data_merged, xr.Dataset) and (self.var_to_rename is not None):
            data_merged = data_merged.rename_vars(self.var_to_rename)

        if isinstance(data_merged, xr.Dataset) and (self.var_to_keep is not None):
            data_merged = data_merged[self.var_to_keep]

        return data_merged
    # ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
