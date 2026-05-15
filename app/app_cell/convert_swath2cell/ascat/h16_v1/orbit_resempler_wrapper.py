import numpy as np
from numpy.lib.recfunctions import append_fields, rename_fields, drop_fields

from orbit_resempler_base import OrbitResampler

class ASCATOrbitResampler(OrbitResampler):
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
                 write_orbit_buffer=False):

        super(ASCATOrbitResampler, self).__init__(
            orbit_io, resam_io,
            spatial_res=spatial_res,
            wfunc=wfunc,
            min_num_nn=min_num_nn,
            dt=dt,
            write_orbit_buffer=write_orbit_buffer)

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
        super(ASCATOrbitResampler, self).resample(
            time_stamps,
            write_n_resampled=write_n_resampled)

    def resample_orbit(self, timestamp):
        """
        Resample orbit.
        """

        gpis, orbit = super(ASCATOrbitResampler,
                            self).resample_orbit(timestamp)

        return gpis, orbit

    def write_resampled_data(self, gpi_data, orbit_data):
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

        fields_to_keep = ['sm', 'sm_noise', 'dir',
                          'corr_flag', 'proc_flag', 'jd']
        fields_to_drop = set(orbit_data.dtype.names).difference(fields_to_keep)
        orbit_data = drop_fields(orbit_data, fields_to_drop)

        # convert to modified julian date
        orbit_data['jd'] = orbit_data['jd'] - 2400000.5
        # convert direction to ascending/descending flag
        orbit_data['dir'] = orbit_data['dir'] < 270
        super(ASCATOrbitResampler, self).write_resampled_data(
            gpi_data, orbit_data)

def translate_processing_flag(proc_flag_orbit):
    """
    Shift bits by one position since orbit processing
    flag contain a flag at bit 0 that does not exist
    in the time series format.
    """
    proc_flag_ts = np.right_shift(proc_flag_orbit.astype(np.uint8), 1)
    return proc_flag_ts
