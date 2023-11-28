"""
Class Features

Name:          lib_fx_orbit_resampler_ascat
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230807'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
from numpy.lib.recfunctions import append_fields, rename_fields, drop_fields
from datetime import timedelta

from lib_fx_orbit_resampler_base import OrbitResamplerBase
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class orbit resampler ascat
class OrbitResamplerAscat(OrbitResamplerBase):
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
                 write_orbit_buffer=False, to_xarray=False,
                 dt_delta=timedelta(minutes=3), dt_buffer=timedelta(days=0)):

        super(OrbitResamplerAscat, self).__init__(
            orbit_io, resam_io,
            spatial_res=spatial_res,
            wfunc=wfunc,
            min_num_nn=min_num_nn,
            dt=dt,
            write_orbit_buffer=write_orbit_buffer, to_xarray=to_xarray,
            dt_delta=dt_delta, dt_buffer=dt_buffer)

    # method to resample data
    def resample(self, time_stamps, write_n_resampled=14 * 2 * 365, **kwargs):
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
        super(OrbitResamplerAscat, self).resample(
            time_stamps,
            write_n_resampled=write_n_resampled, **kwargs)

    def resample_orbit(self, timestamp):
        """
        Resample orbit.
        """

        gpis, orbit = super(OrbitResamplerAscat,
                            self).resample_orbit(timestamp)

        return gpis, orbit

    def write_resampled_data(self, gpi_data, orbit_data,
                             init_file_cell=True, use_file_cell=True, name_file_cell='cell_{:}.data',
                             name_file_ts='{cell_n}.nc'):

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
        if 'sm' not in list(orbit_data.dtype.names):
            orbit_data['Soil Moisture Processing Flag'] = (
                translate_processing_flag(
                    orbit_data['Soil Moisture Processing Flag']))

            fields_to_rename = {
                'Surface Soil Moisture (Ms)': 'sm',
                'Estimated Error In Surface Soil Moisture': 'sm_noise',
                'Soil Moisture Correction Flag': 'corr_flag',
                'Soil Moisture Processing Flag': 'proc_flag',
                'Direction Of Motion Of Moving Observing Platform': 'dir'}
            orbit_data = rename_fields(orbit_data, fields_to_rename)

        # convert to modified julian date
        orbit_data['jd'] = orbit_data['jd'] - 2400000.5

        # convert direction to ascending/descending flag
        if 'dir' in list(orbit_data.dtype.names):
            orbit_data['dir'] = orbit_data['dir'] < 270
        elif 'sat_track_azi' in list(orbit_data.dtype.names):
            orbit_data['sat_track_azi'] = orbit_data['sat_track_azi'] < 270
        else:
            logging.error(' ===> Satellite track direction name is not supported')
            raise RuntimeError('Only "dir" or "sat_track_azi" names are supported by the algorithm')

        fields_to_keep = ['sm', 'sm_noise', 'dir', 'corr_flag', 'proc_flag', 'jd']
        fields_to_drop = set(orbit_data.dtype.names).difference(fields_to_keep)
        orbit_data = drop_fields(orbit_data, fields_to_drop)

        super(OrbitResamplerAscat, self).write_resampled_data(
            gpi_data, orbit_data,
            init_file_cell=init_file_cell, use_file_cell=use_file_cell, name_file_cell=name_file_cell,
            name_file_ts=name_file_ts)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to translate asdat processing flag
def translate_processing_flag(proc_flag_orbit):
    """
    Shift bits by one position since orbit processing
    flag contain a flag at bit 0 that does not exist
    in the time series format.
    """
    proc_flag_ts = np.right_shift(proc_flag_orbit.astype(np.uint8), 1)
    return proc_flag_ts
# ----------------------------------------------------------------------------------------------------------------------
