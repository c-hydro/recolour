"""
Library Features:

Name:          lib_resampler_base
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240209'
Version:       2.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import numpy as np

from datetime import timedelta
from pygeogrids.netcdf import load_grid

from lib_info_args import logger_name
from lib_utils_generic import fill_tags2string
from lib_data_io_generic import delete_file_cell

from lib_resempler_apps import OrbitResampler
from lib_driver_interface_obs import AscatDriverIndexed #AscatDriverH16, AscatDriverH101,
from lib_driver_interface_hsaf import AscatNrtBufrFileList

# suppress warnings
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to set ascat swath data (h16, h101, h102, h103)
class AscatResamplerConfigure:

    # -------------------------------------------------------------------------------------
    # Method to initialize class
    def __init__(self, time_run,
                 ascat_product,
                 ascat_path_swath,
                 ascat_path_raw,
                 ascat_path_ts,
                 ascat_filename_raw='{time_start}_{time_end}.workspace',
                 ascat_filename_ts='{cell}.nc',
                 grid_path='/',
                 grid_filename='TUW_WARP5_grid_info_2_3.nc',
                 ascat_month_path_str='%Y/%m/%d/%H/',
                 writing_mode='w',
                 mask=True,
                 tmpl_tags=None,
                 tmpl_values=None,
                 time_start=None, time_end=None,
                 ascat_updating_ts=False, ascat_updating_raw=True):

        # Set time information
        self.time_run = time_run

        if time_start is None:
            time_start = self.time_run
        if time_end is None:
            time_end = self.time_run

        self.time_start = time_start
        self.time_end = time_end

        # Set tag(s) and value(s) template
        if tmpl_tags is None:
            tmpl_tags = {'datetime': '%Y%m%d%H%M', 'sub_path_time': '%Y/%m/%d', 'cell': '%04d',
                         'time_start': 'from_%Y%m%d%H%M', 'time_end': 'to_%Y%m%d%H%M'}
        self.tmpl_tags = tmpl_tags
        if tmpl_values is None:
            tmpl_values = {'datetime':  None, 'sub_path_time': None, 'cell': None,
                           'time_start': self.time_start, 'time_end': self.time_end}
        self.tmpl_values = tmpl_values

        # Select product driver(s)
        ascat_io_lut = {
            'h16': AscatNrtBufrFileList,
            'h101': None,
            'h102': None,
            'h103': AscatNrtBufrFileList}

        # Compute target grid
        self.target_grid = load_grid(os.path.join(grid_path, grid_filename))
        self.target_cells = np.unique(self.target_grid.activearrcell)

        # Initialize product swath class
        ascat_io_class = ascat_io_lut[ascat_product]

        if ascat_product == 'h16' or ascat_product == 'h103':

            subfolder_template = {'years': '{year}', 'months': '{month}', 'days': '{day}', 'hours': '{hour}'}
            filename_template = '{product_id}_{date}*.buf'

            product_id = ascat_product

            ''' debug h16
            self.time_start, self.time_end = '2023-12-07 07:00:00', '2023-12-07 13:00:00' # test case
            self.time_start, self.time_end = '2023-05-24 09:00:00', '2023-05-24 13:00:00' # error in empty dataframe
            self.time_start, self.time_end = '2023-05-25 10:00:00', '2023-05-25 12:00:00' # error in wrong message length (gribapi)
            print('DEBUG DATETIME --> DELETE ME -- lib_resampler_base.py ll: 87')
            '''

            self.ascat_io_swath = ascat_io_class(
                root_path=ascat_path_swath, product_id=product_id,
                filename_template=filename_template,
                subfolder_template=subfolder_template)

            self.ascat_time_stamps, self.ascat_time_intervals = np.array(
                self.ascat_io_swath.tstamps_for_daterange(self.time_start, self.time_end, frequency, rounding))
            ascat_file_names = self.ascat_io_swath.search_period(self.time_start, self.time_end)

            self.dt_delta = timedelta(minutes=3)
            self.dt_buffer = timedelta(days=0)

        else:
            log_stream.error(' ===> Product "' + ascat_product + '" is not expected')
            raise NotImplemented('Case not implemented yet')

        # Define ts path and filename
        tmpl_values['sub_path_time'] = self.time_run
        tmpl_values['datetime'] = self.time_run
        self.ascat_path_ts = ascat_path_ts
        self.ascat_filename_ts = ascat_filename_ts

        ascat_path_raw = fill_tags2string(ascat_path_raw, self.tmpl_tags, tmpl_values)
        #ascat_filename_raw = fill_tags2string(ascat_filename_raw, self.tmpl_tags, tmpl_values)
        if not os.path.exists(ascat_path_raw):
            os.makedirs(ascat_path_raw)
        self.ascat_path_raw = ascat_path_raw
        self.ascat_filename_raw = ascat_filename_raw

        # Initialize product indexed class
        self.ascat_io_ts = AscatDriverIndexed(
            ascat_path_ts,
            grid=self.target_grid,
            mode=writing_mode,
            ioclass_kws={'time_units': "days since 1858-11-17 00:00:00"})

        # Define time stamps
        #self.time_stamps = np.array(self.ascat_io_swath.tstamps_for_daterange(time_start, time_end))
        # Set mask flag
        self.mask = mask
        # Set updating flag
        self.ascat_updating_raw = ascat_updating_raw
        self.ascat_updating_ts = ascat_updating_ts

        if self.ascat_updating_ts:
            delete_file_cell(self.ascat_path_ts, filename_ts=self.ascat_filename_ts,
                             cells=self.target_cells)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to resample data from swath to ts format
    def resampler(self, write_n_resampled=2000, write_orbit_buffer=True):

        # initialize resampler
        resampler = AscatResamplerInterface(
            self.ascat_io_swath, self.ascat_io_ts,
            spatial_res=25000,
            dt=15,
            wfunc='hamming',
            write_orbit_buffer=write_orbit_buffer,
            mask=self.mask,
            dt_buffer=self.dt_buffer,dt_delta=self.dt_delta,
            upd_raw=self.ascat_updating_raw, folder_raw=self.ascat_path_raw, fileraw_tmpl=self.ascat_filename_raw,
            filename_tmpl=self.ascat_filename_ts, to_xarray=True)
        # run resampler
        resampler.resample(self.ascat_time_intervals, write_n_resampled=write_n_resampled)

    # -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to resample ascat swath data (h16, h101, h102, h103)
class AscatResamplerInterface(OrbitResampler):
    """
    Resample ASCAT Orbits to time series while also performing SSF
    retrieval.

    Parameters
    ----------
    orbit_io : Subclass of MultiTemporalImageBase
        Orbit file reader.
    resam_usort_io : Subclass of GriddedTsBase
        Re-sampled data writer, has to support the write_cell method.
    spatial_res: int, optional
        spatial resolution in meters
    wfunc: string, optional
        one of 'hamming', 'boxcar'
    min_num_nn: int, optional
        Minimum number of nearest neighbors to be found to use the data
        from the swath/orbit
    dt: int, optional
        Temporal distance in minutes between two files so that they
        should be considered temporal neighbors
    """

    def __init__(self, orbit_io, resam_io, spatial_res=25000.,
                 wfunc='hamming', min_num_nn=3, dt=15,
                 write_orbit_buffer=False,
                 to_xarray=False,
                 dt_delta=timedelta(minutes=3), dt_buffer=timedelta(days=0),
                 mask=False, filename_tmpl='{:}.nc',
                 upd_raw=True, folder_raw=None, fileraw_tmpl='{:}.workspace'):

        super(AscatResamplerInterface, self).__init__(
            orbit_io, resam_io,
            spatial_res=spatial_res,
            wfunc=wfunc,
            min_num_nn=min_num_nn,
            dt=dt,
            write_orbit_buffer=write_orbit_buffer,
            to_xarray=to_xarray,
            dt_delta=dt_delta, dt_buffer=dt_buffer,
            mask=mask,
            upd_raw=upd_raw, folder_raw=folder_raw, fileraw_tmpl=fileraw_tmpl,
            filename_tmpl=filename_tmpl)


    def resample(self, time_stamps, write_n_resampled=14 * 2 * 365):
        """
        Run re-sampling for given time_stamps.
        Data is written if available memory becomes low
        or if write_n_resampled files have been resampled.

        Parameters
        ----------
        time_stamps : numpy.ndarray
            Orbit time stamp information.
        write_n_resampled: int, optional
            Write data if more than n timestamps have been resampled.
            The default is one year of ASCAT data at the moment.
        proc_param : dict
            Processing parameters.

        Returns
        -------
        resampled_timestamps: list
            list of resampled timestamps that actually contained
            data on the target grid
        """
        super(AscatResamplerInterface, self).resample(
            time_stamps,
            write_n_resampled=write_n_resampled)

    def resample_orbit(self, timestamp):
        """
        Resample orbit.
        """

        gpis, orbit = super(AscatResamplerInterface,
                            self).resample_orbit(timestamp)

        return gpis, orbit

    def write_resampled_data(self, gpi_data, orbit_data, orbit_workspace):
        """
        write the data after name conversion.

        Parameters
        ----------
        gpi_data : numpy.recarray
            Grid point information. A grid point indices
            for each element in the orbit data array.
        orbit_data : numpy.recarray
            record array containing the resampled data at the grid point
            indices of gpi_data
        """

        '''
        orbit_data['Soil Moisture Processing Flag'] = (
            translate_processing_flag(orbit_data['Soil Moisture Processing Flag']))

        fields_to_rename = {
            'Surface Soil Moisture (Ms)': 'sm',
            'Estimated Error In Surface Soil Moisture': 'sm_noise',
            'Soil Moisture Correction Flag': 'corr_flag',
            'Soil Moisture Processing Flag': 'proc_flag',
            'Direction Of Motion Of Moving Observing Platform': 'dir'}
        orbit_data = rename_fields(orbit_data, fields_to_rename)

        fields_to_keep = ['sm', 'sm_noise', 'dir',
                          'corr_flag', 'proc_flag', 'jd']
                          
        fields_to_drop = set(orbit_data.dtype.names).difference(fields_to_keep)
        orbit_data = drop_fields(orbit_data, fields_to_drop)
        '''

        # update processing flag
        orbit_data['proc_flag'] = (translate_processing_flag(orbit_data['proc_flag']))

        # convert to modified julian date
        orbit_data['jd'] = orbit_data['jd'] - 2400000.5
        # convert direction to ascending/descending flag
        orbit_data['dir'] = orbit_data['dir'] < 270

        orbit_workspace['jd'] = orbit_workspace['jd'] - 2400000.5
        orbit_workspace['dir'] = orbit_workspace['dir'] < 270

        super(AscatResamplerInterface, self).write_resampled_data(
            gpi_data, orbit_data, orbit_workspace=orbit_workspace)

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to shift bits
def translate_processing_flag(proc_flag_orbit):
    """
    Shift bits by one position since orbit processing
    flag contain a flag at bit 0 that does not exist
    in the time series format.
    """
    proc_flag_ts = np.right_shift(proc_flag_orbit.astype(np.uint8), 1)
    return proc_flag_ts
# -------------------------------------------------------------------------------------
