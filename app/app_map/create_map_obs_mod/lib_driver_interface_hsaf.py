"""
Library Features:

Name:          lib_driver_interface_hsaf
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240209'
Version:       '2.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import warnings
import pdbufr

from datetime import datetime
from datetime import timedelta
from copy import deepcopy
from cadati.cal_date import cal2dt

import numpy as np
import xarray as xr
import pandas as pd

try:
    import pygrib
except ImportError:
    warnings.warn('pygrib can not be imported GRIB files (H14) can not be read.')

from ascat.eumetsat.level2 import AscatL2File
from ascat.utils import get_toi_subset, get_roi_subset
from ascat.file_handling import ChronFiles

from pygeobase.utils import split_daterange_in_intervals

from lib_info_args import logger_name

# from ascat.utils import tmp_unzip
# from ascat.file_handling import ChronFiles
# from ascat.read_native.bufr import BUFRReader
# from ascat.read_native.bufr import AscatL2BufrFile
# from ascat.read_native.cdr import AscatGriddedNcTs
# from pygeobase.io_base import IntervalReadingMixin, ImageBase

# defined nan values
bufr_nan = 1.7e+38
float32_nan = -999999.
uint8_nan = np.iinfo(np.uint8).max
uint16_nan = np.iinfo(np.uint16).max
int32_nan = np.iinfo(np.int32).max

nan_val_dict = {
    np.float32: float32_nan,
    np.uint8: uint8_nan,
    np.uint16: uint16_nan,
    np.int32: int32_nan
}

# logging
log_stream = logging.getLogger(logger_name)
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

    def search_period(self, dt_start, dt_end, dt_delta=timedelta(days=1),
                      search_date_str="%Y%m%d*", date_str="%Y%m%d", date_field="date",
                      ):

        filenames = []
        for dt_cur in np.arange(dt_start, dt_end, dt_delta).astype(datetime):
            files, dates = self.search_date(dt_cur, return_date=True)
            for f, dt in zip(files, dates):
                if f not in filenames and dt_start <= dt <= dt_end:
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

    def _parse_date(self, filename, date_str, date_field):
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

class AscatNrtBufrWrapper(AscatL2File):

    msg_name_lookup = {
            'satelliteIdentifier': "Satellite Identifier",
            'orbitNumber': "Orbit Number",
            'beamIdentifier': "f_Beam Identifier",
            'directionOfMotionOfMovingObservingPlatform': 'Direction Of Motion Of Moving Observing Platform',
            'crossTrackCellNumber': 'Cross-Track Cell Number',
            'latitude': "latitude",
            'longitude': "longitude",
            'year': "year",
            'month': "month",
            'day': "day",
            'hour': "hour",
            'minute': "minute",
            'second': "second",
            'estimatedErrorInSurfaceSoilMoisture': "Estimated Error In Surface Soil Moisture",
            'surfaceSoilMoisture': "Surface Soil Moisture (Ms)",
            'meanSurfaceSoilMoisture': "Mean Surface Soil Moisture",
            'soilMoistureCorrectionFlag': "Soil Moisture Correction Flag",
            'soilMoistureProcessingFlag': "Soil Moisture Processing Flag",
            'soilMoistureQuality': "Soil Moisture Quality",
            'snowCover': "Snow Cover",
            'frozenLandSurfaceFraction': "Frozen Land Surface Fraction",
            'topographicComplexity': "Topographic Complexity",
        }

    def read(self, toi=None, roi=None, generic=True, to_xarray=False):

        data, metadata = self.select(generic=generic, to_xarray=to_xarray)

        if data is not None:
            if toi:
                data = get_toi_subset(data, toi)

            if roi:
                data = get_roi_subset(data, roi)

        return data, metadata

    def select(self, generic=False, to_xarray=False):

        # select datasets columns
        columns_name = tuple(list(self.msg_name_lookup.keys()))

        # get dataset
        try:
            df = pdbufr.read_bufr(self.filename, columns=columns_name)
        except Exception as e:
            log_stream.warning(' ===> Errors in reading BUFR file: "' + str(e) + '". Data are not available.')
            data, metadata = None, None
            return data, metadata
        # check if dataset is not empty
        if df.empty:
            log_stream.warning(' ===> Empty dataframe. Data are not available.')
            data, metadata = None, None
            return data, metadata

        # rename dataset columns
        col_rename = {}
        for i, col in enumerate(df.columns.to_list()):
            name = self.msg_name_lookup[col]
            if name is not None:
                col_rename[col] = name
        data = df.rename(columns=col_rename)

        # check if all selected data are available
        var_na = []
        for var in list(self.msg_name_lookup.values()):
            if var not in data.columns:
                var_na.append(var)
        if var_na:
            var_str = ', '.join(var_na)
            log_stream.warning(' ===> Variable(s) "' + var_str + '" not available. Data are not complete.')
            data, metadata = None, None
            return data, metadata

        # adapt, convert and organize dataset
        data["lat"] = df["latitude"].values.astype(np.float32)
        data["lon"] = df["longitude"].values.astype(np.float32)

        year = df["year"].values.astype(int)
        month = df["month"].values.astype(int)
        day = df["day"].values.astype(int)
        hour = df["hour"].values.astype(int)
        minute = df["minute"].values.astype(int)
        seconds = df["second"].values.astype(int)
        milliseconds = np.zeros(seconds.size)
        cal_dates = np.vstack(
            (year, month, day, hour, minute, seconds, milliseconds)).T

        data['time'] = cal2dt(cal_dates)
        data = data.to_records(index=False)
        data = {name:data[name] for name in data.dtype.names}

        # control to the soil moisture datasets
        if not np.all(data["Mean Surface Soil Moisture"]) is None: # add to avoid error for all data equal to None
            data["Mean Surface Soil Moisture"] *= 100.
        elif np.all(data["Mean Surface Soil Moisture"]) is None:
            log_stream.warning(' ===> Soil Moisture data are not available')
            data, metadata = None, None
            return data, metadata
        else:
            log_stream.error(' ===> Soil moisture data format is not expected')
            raise NotImplemented('Case not implemented yet')

        # control to the snow cover datasets
        if np.all(data["Snow Cover"]) is None:
            data["Snow Cover"] = np.zeros(data["Snow Cover"].size, dtype=np.uint8)
        # control to the topographic complexity datasets
        if np.all(data["Topographic Complexity"] ) is None:
            data["Topographic Complexity"] = np.full(data["Topographic Complexity"].size, np.nan)

        # organize metadata
        metadata = {
            'platform_id': data['Satellite Identifier'][0].astype(int),
            'orbit_start': np.uint32(data['Orbit Number'][0]),
            'filename': os.path.basename(self.filename)}

        # add/rename/remove fields according to generic format
        if generic:
            data = conv_bufrl2_generic(data, metadata)

        # convert dict to xarray.Dataset or numpy.ndarray
        if to_xarray:
            for k in data.keys():
                if len(data[k].shape) == 1:
                    dim = ['obs']
                elif len(data[k].shape) == 2:
                    dim = ['obs', 'beam']

                data[k] = (dim, data[k])

            coords = {}
            coords_fields = ['lon', 'lat', 'time']
            for cf in coords_fields:
                coords[cf] = data.pop(cf)

            data = xr.Dataset(data, coords=coords, attrs=metadata)
        else:
            # collect dtype info
            dtype = []
            # fill_value = []

            for var_name in data.keys():

                if len(data[var_name].shape) == 1:
                    dtype.append((var_name, data[var_name].dtype.str))
                    # fill_value.append(data[var_name].fill_value)

                elif len(data[var_name].shape) > 1:
                    dtype.append((var_name, data[var_name].dtype.str,
                                  data[var_name].shape[1:]))
                    # fill_value.append(data[var_name].shape[1] *
                    #                   [data[var_name].fill_value])

            ds = np.ma.empty(data['time'].size, dtype=np.dtype(dtype))
            # fill_value_arr = np.array((*fill_value, ), dtype=np.dtype(dtype))

            for k, v in data.items():
                ds[k] = v

            # ds.fill_value = fill_value_arr
            data = ds

        return data, metadata

def conv_bufrl2_generic(data, metadata):
    """
    Rename and convert data types of dataset.

    Spacecraft_id vs sat_id encoding

    BUFR encoding - Spacecraft_id
    - 1 ERS 1
    - 2 ERS 2
    - 3 Metop-1 (Metop-B)
    - 4 Metop-2 (Metop-A)
    - 5 Metop-3 (Metop-C)

    Internal encoding - sat_id
    - 1 ERS 1
    - 2 ERS 2
    - 3 Metop-2 (Metop-A)
    - 4 Metop-1 (Metop-B)
    - 5 Metop-3 (Metop-C)

    Parameters
    ----------
    data: dict of numpy.ndarray
        Original dataset.
    metadata: dict
        Metadata.

    Returns
    -------
    data: dict of numpy.ndarray
        Converted dataset.
    """
    skip_fields = ['Satellite Identifier']

    gen_fields_beam = {
        #'Radar Incidence Angle': ('inc', np.float32),
        #'Backscatter': ('sig', np.float32),
        #'Antenna Beam Azimuth': ('azi', np.float32),
        #'ASCAT Sigma-0 Usability': ('f_usable', np.uint8),
        'Beam Identifier': ('beam_num', np.uint8),
        #'Radiometric Resolution (Noise Value)': ('kp_noise', np.float32),
        #'ASCAT KP Estimate Quality': ('kp', np.float32),
        #'ASCAT Land Fraction': ('f_land', np.float32)
    }

    gen_fields_lut = {
        'Orbit Number': ('abs_orbit_nr', np.int32),
        'Cross-Track Cell Number': ('node_num', np.uint8),
        'Direction Of Motion Of Moving Observing Platform': ('sat_track_azi', np.float32),
        'Surface Soil Moisture (Ms)': ('sm', np.float32),
        'Estimated Error In Surface Soil Moisture': ('sm_noise', np.float32),
        #'Backscatter': ('sig40', np.float32),
        #'Estimated Error In Sigma0 At 40 Deg Incidence Angle': ('sig40_noise', np.float32),
        #'Slope At 40 Deg Incidence Angle': ('slope40', np.float32),
        #'Estimated Error In Slope At 40 Deg Incidence Angle': ('slope40_noise', np.float32),
        #'Soil Moisture Sensitivity': ('sm_sens', np.float32),
        #'Dry Backscatter': ('dry_sig40', np.float32),
        #'Wet Backscatter': ('wet_sig40', np.float32),
        'Mean Surface Soil Moisture': ('sm_mean', np.float32),
        # 'Rain Fall Detection': ('rf', np.float32),
        'Soil Moisture Correction Flag': ('corr_flag', np.uint8),
        'Soil Moisture Processing Flag': ('proc_flag', np.uint8),
        'Soil Moisture Quality': ('agg_flag', np.uint8),
        'Snow Cover': ('snow_prob', np.uint8),
        'Frozen Land Surface Fraction': ('frozen_prob', np.uint8),
        #'Inundation And Wetland Fraction': ('wetland', np.uint8),
        'Topographic Complexity': ('topo', np.uint8)
    }

    for var_name in skip_fields:
        if var_name in data:
            data.pop(var_name)

    for var_name, (new_name, new_dtype) in gen_fields_lut.items():
        mask = (data[var_name] == bufr_nan) |  (np.isnan(data[var_name]))
        data[var_name][mask] = nan_val_dict[new_dtype]

        data[new_name] = np.ma.array(data.pop(var_name).astype(new_dtype), mask=mask)
        data[new_name].fill_value = nan_val_dict[new_dtype]

    '''
    for var_name, (new_name, new_dtype) in gen_fields_beam.items():

        f = ['{}_{}'.format(b, var_name) for b in ['f', 'm', 'a']]

        mask = np.vstack((data[f[0]] == bufr_nan, data[f[1]] == bufr_nan,
                          data[f[2]] == bufr_nan)).T

        data[new_name] = np.ma.vstack((data.pop(f[0]), data.pop(f[1]),
                                       data.pop(f[2]))).T.astype(new_dtype)

        data[new_name].mask = mask
        data[new_name][mask] = nan_val_dict[new_dtype]

        data[new_name].fill_value = nan_val_dict[new_dtype]
    '''

    if data['node_num'].max() == 82:
        data['swath_indicator'] = np.ma.array(1 * (data['node_num'] > 41),
                                              dtype=np.uint8,
                                              mask=data['node_num'] > 82)
    elif data['node_num'].max() == 42:
        data['swath_indicator'] = np.ma.array(1 * (data['node_num'] > 21),
                                              dtype=np.uint8,
                                              mask=data['node_num'] > 42)
    else:
        raise ValueError('Cross-track cell number size unknown')

    n_lines = data['lat'].shape[0] / data['node_num'].max()
    line_num = np.arange(n_lines).repeat(data['node_num'].max())
    data['line_num'] = np.ma.array(line_num,
                                   dtype=np.uint16,
                                   mask=np.zeros_like(line_num),
                                   fill_value=uint16_nan)

    sat_id = np.ma.array([0, 0, 0, 4, 3, 5], dtype=np.uint8)
    data['sat_id'] = np.ma.zeros(data['time'].size,
                                 dtype=np.uint8) + sat_id[int(
                                     metadata['platform_id'])]
    data['sat_id'].mask = np.zeros(data['time'].size)
    data['sat_id'].fill_value = uint8_nan

    # compute ascending/descending direction
    data['as_des_pass'] = np.ma.array(data['sat_track_azi'] < 270,
                                      dtype=np.uint8,
                                      mask=np.zeros(data['time'].size),
                                      fill_value=uint8_nan)

    mask = data['lat'] == bufr_nan
    data['lat'] = np.ma.array(data['lat'], mask=mask, fill_value=float32_nan)

    mask = data['lon'] == bufr_nan
    data['lon'] = np.ma.array(data['lon'], mask=mask, fill_value=float32_nan)

    data['time'] = np.ma.array(data['time'], mask=mask, fill_value=0)

    return data


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

        super().__init__(root_path, AscatNrtBufrWrapper, filename_template,
                         sf_templ=subfolder_template)
        '''
        super().__init__(root_path, AscatL2File, filename_template,
                         sf_templ=subfolder_template)
        '''
    def tstamps_for_daterange(self, dt_start, dt_end, dt_delta=timedelta(days=1)):

        file_names, file_timestamps = [], []

        if not isinstance(dt_start, pd.Timestamp):
            ts_start = pd.Timestamp(dt_start)
        else:
            ts_start = deepcopy(dt_start)
        if not isinstance(dt_end, pd.Timestamp):
            ts_end = pd.Timestamp(dt_end)
        else:
            ts_end = deepcopy(dt_end)

        ts_totalseconds = (ts_end - ts_start).total_seconds()
        dt_totalseconds = dt_delta.total_seconds()

        if ts_totalseconds < dt_totalseconds:
            dt_delta = timedelta(seconds=ts_totalseconds)
        else:
            dt_delta = timedelta(seconds=dt_totalseconds)

        dt_range = np.arange(ts_start, ts_end + dt_delta, dt_delta)
        for dt_cur in dt_range.astype(datetime):
            files, dates = self.search_date(
                dt_cur,
                return_date=True, search_date_str="%Y%m%d*", date_str="%Y%m%d", date_field="date")
            for f, dt in zip(files, dates):

                ts_step = pd.to_datetime(dt)

                if (files not in file_names) and (ts_start <= ts_step <= ts_end):
                    file_names.append(f)
                    file_timestamps.append(dt)

        if file_timestamps:
            file_time_intervals = split_daterange_in_intervals(file_timestamps[0], file_timestamps[-1],
                                                               self.chunk_minutes)
        else:
            log_stream.warning(' ===> No files found in the selected period')
            file_time_intervals = []

        return file_timestamps, file_time_intervals

    def search_period(self, dt_start, dt_end, dt_delta=timedelta(days=1),
                      search_date_str="%Y%m%d*", date_str="%Y%m%d", date_field="date"):

        if not isinstance(dt_start, pd.Timestamp):
            ts_start = pd.Timestamp(dt_start)
        else:
            ts_start = deepcopy(dt_start)
        if not isinstance(dt_end, pd.Timestamp):
            ts_end = pd.Timestamp(dt_end)
        else:
            ts_end = deepcopy(dt_end)

        filenames = []
        for dt_cur in np.arange(ts_start, ts_end, dt_delta).astype(datetime):
            files, dates = self.search_date(dt_cur, return_date=True)
            for f, dt in zip(files, dates):

                ts_step = pd.to_datetime(dt)

                if (f not in filenames) and (ts_start <= ts_step <= ts_end):
                    filenames.append(f)
        return filenames

    def read_period(self, dt_start, dt_end, dt_delta=timedelta(days=1),
                    dt_buffer=timedelta(days=1), **kwargs):

        filenames = self.search_period(dt_start-dt_buffer, dt_end, dt_delta)
        data = []
        for filename in filenames:

            log_stream.info(' ------> Read file: {}'.format(filename) + ' ... ')

            self._open(filename)
            d = self.fid.read(None, None, **kwargs)
            if d[0] is not None and d[1] is not None:
                data.append(d)
                log_stream.info(' ------> Read file: {}'.format(filename) + ' ... DONE')
            else:
                log_stream.info(' ------> Read file: {}'.format(filename) + ' ... SKIPPED. DATA NOT FOUND.')

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

    def _parse_date(self, filename, date_str, date_field):
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
        objdate = datetime.strptime(filedate, '%Y%m%d%H%M%S')

        return objdate

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

                if not isinstance(data_step, xr.Dataset):
                    log_stream.error(' ===> Data not in xarray dataset format')
                    raise RuntimeError('Set "to_xarray" to True in the read method')

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
