"""
Class Features

Name:          lib_fx_orbit_resampler_base
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230807'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import sys
import warnings

import pandas as pd
import progressbar

from functools import partial
from datetime import timedelta, datetime

import numpy as np
import xarray as xr

from pygeobase.io_base import MultiTemporalImageBase, GriddedTsBase, Image
from pykdtree.kdtree import KDTree
from memon import MemoryMonitor

import ascat

from lib_data_io_pickle import write_file_obj, read_file_obj
from lib_fx_datasets_source import AscatNrtBufrFileList, AscatEpsBufrFileList

# suppress warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to resample orbit base
class OrbitResamplerBase(object):
    """
    Resample Orbit/Swath data onto a target grid defined by the resam_io
    object. The resam_io object has to implement a write_cell method
    and must also work on a :py:class:`pygeogrids.grids.CellGrid` grid object.

    Parameters
    ----------
    orbit_io : Subclass of MultiTemporalImageBase
        Orbit file reader.
    resam_io : Subclass of GriddedTsBase
        Re-sampled data writer, has to support the write_cell method.
    spatial_res: int, optional
        spatial resolution in meters
    wfunc: string, optional
        one of 'hamming', 'boxcar'
    min_num_nn: int, optional
        Minimum number of nearest neighbors to be found to use the data
        from the swath/orbit
    max_num_nn: int, optional
        Maximum number of nearest neighbors to search for. If not given
        then this is estimated from target grid.
    dt: int, optional
        Temporal distance in minutes between two files so that they
        should be considered temporal neighbors
    write_orbit_buffer: boolean, optional
        Specify if the orbit butter numpy array should be written to disk
        after the resampling is finished. Set to True if you want
        to continue the resampling in a later processing run.
    resample_first: boolean, optional
        if set to False the data from the first sub-orbit that is read
        is not resampled, except for the points that contribute to the
        next sub-orbit.
    """

    def __init__(self, orbit_io, resam_io, spatial_res=25000.,
                 wfunc='hamming', min_num_nn=3, max_num_nn=None, dt=15,
                 write_orbit_buffer=False, buffer_saved=False, buffer_init=True,
                 resample_first=True, to_xarray=False,
                 dt_delta=timedelta(minutes=3), dt_buffer=timedelta(days=0)):

        if isinstance(orbit_io, AscatNrtBufrFileList):
            pass
        elif isinstance(orbit_io, MultiTemporalImageBase):
            pass
        elif isinstance(orbit_io, AscatEpsBufrFileList):
            pass
        else:
            raise ValueError('Orbit data must be of type AscatNrtBufrFileList or MultiTemporalImageBase')

        if not isinstance(resam_io, GriddedTsBase):
            name = GriddedTsBase.__name__
            raise ValueError('Output data must be of type{:}'.format(name))

        self.orbit_io = orbit_io
        self.resam_io = resam_io

        # set maximum time delta of concurrent measurements
        self.dt = dt / 60. / 24.

        # set minimum number of nearest neighbours for one target grid point
        # to be considered for resampling
        self.min_num_nn = min_num_nn

        self.resample_first = resample_first
        self.count_instructed_resampled = 0

        # set flag to check if orbit buffer is saved to disk before closing
        # exiting object
        self.buffer_saved = False

        # set flag to indicate if orbit buffer needs to written
        self.write_orbit_buffer = write_orbit_buffer

        # get orbit buffer for append if available
        self.orbit_buffer_file = os.path.join(resam_io.path, "orbitBuffer.npz")
        if buffer_init:
            if os.path.exists(self.orbit_buffer_file):
                os.remove(self.orbit_buffer_file)

        self.buffer_keys = ['jds_prev', 'data_prev', 'nn_flag_prev',
                            'index_prev', 'distance_prev']

        self.orbit_buffer = {}
        if os.path.isfile(self.orbit_buffer_file):
            prev = np.load(self.orbit_buffer_file)
            for key in self.buffer_keys:
                self.orbit_buffer[key] = prev[key]
        else:
            for key in self.buffer_keys:
                self.orbit_buffer[key] = None

        # define window function and corresponding search radius
        self.window_function_name = wfunc
        self.spatial_res = spatial_res

        self.search_radius = get_window_radius(self.window_function_name, self.spatial_res / 2.)

        self.target_grid = resam_io.grid
        self.ref_geodatum = self.target_grid.geodatum
        self.xyz_target_grid = \
            self.ref_geodatum.toECEF(self.target_grid.activearrlon,
                                     self.target_grid.activearrlat)

        orig_window_function = {'hamming': hamming_window, 'boxcar': boxcar}[self.window_function_name]

        self.window_function = partial(
            orig_window_function, self.search_radius)

        self.xyz_target_grid = np.transpose(np.vstack(self.xyz_target_grid))

        if max_num_nn is None:
            # use first grid point to estimate number of neighbours within
            # search radius
            tmp_diff = self.xyz_target_grid[0, :] - self.xyz_target_grid
            tmp_dist = np.sqrt(
                tmp_diff[:, 0]**2 + tmp_diff[:, 1]**2 + tmp_diff[:, 2]**2)

            self.n_neighbours = 3 * np.sum(tmp_dist <= self.search_radius)
        else:
            self.n_neighbours = max_num_nn

        self.dt_delta = dt_delta
        self.dt_buffer = dt_buffer
        self.to_xarray = to_xarray

        self.timefield = None

    def __enter__(self):
        """
        Enter the runtime context related to the object.
        """
        return self

    def close(self, *args):
        """
        Close method.
        """
        return self.save_orbit_buffer()

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Exit the runtime context related to the object
        """
        self.save_orbit_buffer()
        return self.buffer_saved

    def __del__(self):
        """
        Distructor, save orbit files.
        """
        self.save_orbit_buffer()
        return self.buffer_saved

    def save_orbit_buffer(self):
        """
        Save orbit buffer on disk.
        """
        if self.buffer_saved is False and self.write_orbit_buffer is True:

            if os.path.isdir(os.path.dirname(self.orbit_buffer_file)) is False:
                os.makedirs(os.path.dirname(self.orbit_buffer_file))

            if os.path.isfile(self.orbit_buffer_file):
                os.remove(self.orbit_buffer_file)

            np.savez_compressed(self.orbit_buffer_file, **self.orbit_buffer)
            self.buffer_saved = True

    def nn_search(self, lon, lat, **kwargs):
        """
        Method to perform nearest neighbour search.

        Parameters
        ----------
        lon : list, np.array
            Longitudes of source grid.
        lat : list, np.array
            Latitudes of source grid.

        Returns
        -------
        nn_flag : np.array
            Flag indicating if a nearest neighbour was found for the given
            target grid locations.
        nn_num : int
            Number of neighbours found.
        index : np.array
            Indices of neighbours to source
        distance : np.array
            Distance to neighbors.
        """
        xyz = np.transpose(np.vstack(self.ref_geodatum.toECEF(lon, lat)))

        init_KDTree = KDTree(xyz)

        distance, index = \
            init_KDTree.query(self.xyz_target_grid, k=self.n_neighbours,
                              distance_upper_bound=self.search_radius,
                              **kwargs)

        # calculate flag and number of neighbours
        nn_flag = np.any(distance != np.inf, axis=1)
        nn_num = np.sum(distance != np.inf, axis=1)

        return nn_flag, nn_num, index, distance

    def resample_orbit(self, timestamp, to_array=False):
        """
        This method performs the spatial interpolation (resampling) of one
        orbit file to the earth fixed grid defined in the resam_io.

        Parameters
        ----------
        timestamp : datetime
            Orbit time information.
        """
        resampled_gpis = []
        resampled_data = []

        self.count_instructed_resampled = self.count_instructed_resampled + 1
        try:

            if float(ascat.__version__.split('.')[0]) <= 1 and sys.version_info[0] == 2:
                orbit = self.orbit_io.read(timestamp, mask=True)
            elif float(ascat.__version__.split('.')[0]) >= 2.0 and sys.version_info[0] == 3:
                orbit = self.orbit_io.read_period(timestamp[0], timestamp[1],
                                                  dt_buffer=self.dt_buffer,
                                                  dt_delta=self.dt_delta,
                                                  **{'to_xarray': self.to_xarray})
            else:
                raise RuntimeError('python and ascat libraries are not tested for this case')

            if orbit is None:
                raise Exception("No data for timestamp")
        except Exception as ex:
            print(ex)
            return np.array([]), np.array([])

        if isinstance(orbit, Image):
            data = self.rectDictToNdarray(orbit.data)
            self.timefield = orbit.timekey
            orbit_jd = orbit.data[orbit.timekey]
            lon = orbit.lon
            lat = orbit.lat
        elif isinstance(orbit, xr.Dataset):
            data = self.rectDsetToNdarray(orbit)
            # data_sm = orbit.sm.values
            orbit_jd = orbit.jd.values
            lon = orbit.lon.values
            lat = orbit.lat.values

            self.timefield = 'jd'

        if not orbit:
            return np.array([]), np.array([])

        # Check for empty messages file
        if not data.__len__() > 0:
            return np.array([]), np.array([])

        jds = np.hstack((orbit_jd, np.nan))

        # search for nearest neighbours
        nn_flag, nn_num, index, distance = self.nn_search(lon, lat)

        # check if the data is self overlapping
        jd_diff = np.nanmax(jds[index], axis=1) - np.nanmin(jds[index], axis=1)
        n_overlapping_points = np.nansum(jd_diff[np.isfinite(jd_diff)] > self.dt)

        if n_overlapping_points > 0:
            warnings.warn("{} points overlap in one sub-orbit".format(n_overlapping_points),
                          UserWarning)
        # resample previous orbit by considering boundary to
        # new orbit part.
        # Check if previous data is available and nearest
        # neighbors were found.
        if self.orbit_buffer['data_prev'] is not None and \
                np.sum(self.orbit_buffer['nn_flag_prev']) > 0.:

            # check which grid points where found for both orbits parts
            pos = (self.orbit_buffer['nn_flag_prev'] & nn_flag)

            # check if previous sub-orbit is the chronological
            # previous orbit
            # calculate mean temporal difference of duplicates
            jd_diff = np.abs(
                np.nanmean(jds[index[pos, :]], axis=1) -
                np.nanmean(self.orbit_buffer['jds_prev'][
                    self.orbit_buffer['index_prev'][pos, :]],
                    axis=1))

            # check if positions are also temporal neighbors (15 minutes)
            include_gpis = pos.copy()
            include_gpis[include_gpis] = (jd_diff <= self.dt)

            # exclude positions which are not temporal neighbors
            exclude_gpis = pos.copy()
            exclude_gpis[exclude_gpis] = (jd_diff > self.dt)

            # attach data from grid points found in this sub orbit and the
            # previous sub orbit to the data from the previous sub orbit
            # Add a last element for indexing when no neighbor was found.
            data_tmp = self.orbit_buffer['data_prev']
            data_res = np.hstack((data_tmp, data, data[0]))

            # prepare nearest neighbour index array for resampling previous
            # index and index for new points have to be harmonized. the last
            # index of each array always indicates no neighbor found.
            index_tmp = index.copy()
            # shift the index of the new data by the size of the previous data
            # so that they still fit together after hstack
            shift_index = self.orbit_buffer['data_prev'].shape[0]
            index_tmp += shift_index
            index_res = np.hstack((self.orbit_buffer['index_prev'],
                                   index_tmp))
            # now values that were zero in the new index have the same
            # values as the NaN values in the previous orbit
            prevNaN = np.where(self.orbit_buffer['index_prev'] == shift_index)
            # since the old data is in the beginning of the index_res after hstack
            # we can use the index in prevNaN unchanged.
            # set the previous NaN values to the new last possible indices
            index_res[prevNaN] = data_res.shape[0] - 1

            # prepare distance array of nearest neighbour for resampling
            distance_tmp = distance.copy()
            distance_tmp[exclude_gpis, :] = np.inf
            distance_res = np.hstack((self.orbit_buffer['distance_prev'],
                                      distance_tmp))

            # only resample points of the last orbit including the new
            # found neighbors in this orbit part.
            gpis_to_resample = self.orbit_buffer['nn_flag_prev']
            # if the it is the second orbit and we should discard the first
            # one, save for the overlap then we select only the overlap gpis
            if self.count_instructed_resampled == 2 and not self.resample_first:
                gpis_to_resample = include_gpis

            index_res = index_res[gpis_to_resample, :]
            distance_res = distance_res[gpis_to_resample, :]

            # check if number of nearest neighbours is greater than
            # defined minimum number
            num_nn = (np.sum(np.isfinite(distance_res), axis=1) >=
                      self.min_num_nn)
            if np.sum(num_nn) > 0.:
                index_res = index_res[num_nn, :]
                distance_res = distance_res[num_nn, :]

                resampled_gpis.append(self.resam_io.grid.activegpis[
                    gpis_to_resample][num_nn])

                # calculated weights from distances
                weights_res, _ = self.window_function(distance_res)

                # resampled Data
                resdata = self.orbit_io.resample_image(data_res,
                                                       index_res,
                                                       distance_res,
                                                       weights_res)
                if type(resdata) == dict:
                    resdata = self.rectDictToNdarray(resdata)
                resampled_data.append(resdata)

            # remove already resampled grid points at the orbit split
            # boundary from further processing
            distance[include_gpis, :] = np.inf
            nn_flag[np.sum(np.isfinite(distance), axis=1) == 0] = False

        self.orbit_buffer['jds_prev'] = jds.copy()
        self.orbit_buffer['data_prev'] = data.copy()
        self.orbit_buffer['nn_flag_prev'] = nn_flag.copy()
        self.orbit_buffer['index_prev'] = index.copy()
        self.orbit_buffer['distance_prev'] = distance.copy()

        if len(resampled_data) > 0:
            return np.hstack(resampled_gpis), np.hstack(resampled_data)
        else:
            return np.array([]), np.array([])

    # method to resample datasets
    def resample(self, time_stamps,
                 write_n_resampled=14 * 2 * 365,
                 use_memon=True,
                 init_file_ws=False, use_file_ws=True,
                 name_file_ws='ws_{datetime_workspace_start}_{datetime_workspace_end}.workspace',
                 init_file_cell=False, use_file_cell=True,
                 name_file_cell='cell_{cell_n}.data',
                 name_file_ts='{cell_n}.nc',
                 time_format_file='%Y%m%d%H%M', time_format_print='%Y%m%d %H:%M',
                 **kwargs):

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
        use_memon: boolean, optional
            If True then use the MemoryMonitor class to decide
            when to write to disk. This can fail if the machine is shared
            in which case it should be disabled.

        Returns
        -------
        resampled_timestamps: list
            list of resampled timestamps that actually contained
            data on the target grid
        """

        # initialize variable(s)
        orbit_data, gpi_data = [], []
        resample_times_list, all_resampled_times = [], []
        # active memory monitoring
        mem_mon = MemoryMonitor()
        mem_mon.start()

        # get time information
        time_period_start, time_period_end = time_stamps[0][0], time_stamps[-1][1]
        time_file_start, time_file_end = time_period_start.strftime(time_format_file), \
            time_period_end.strftime(time_format_file)

        # define workspace file
        name_file_ws = name_file_ws.format(
            datetime_workspace_start=time_file_start, datetime_workspace_end=time_file_end)

        # condition(s) to restore or delete tmp file(s)
        if use_file_ws:
            if init_file_ws:
                if os.path.exists(name_file_ws):
                    os.remove(name_file_ws)
        else:
            if os.path.exists(name_file_ws):
                os.remove(name_file_ws)

        # condition to active resampling
        if not os.path.exists(name_file_ws):

            # iterate over time stamp
            for t in time_stamps:

                # time information
                dt_start, dt_end = t[0], t[1]
                ts_start, ts_end = dt_start.strftime(time_format_print), dt_end.strftime(time_format_print)
                # time info start
                logging.info(' ----> Resample data from "' + ts_start + '" to "' + ts_end + '" ...')

                try:
                    # start memory recording
                    mem_mon.start_recording()

                    # get gpis and orbit data
                    gpis, orbit = self.resample_orbit(t)
                    if gpis.size > 0 and orbit.size > 0:
                        gpi_data.append(gpis)
                        orbit_data.append(orbit)
                        logging.info(' ----> ... Data found ...')
                    else:
                        logging.info(' ----> ... Data not found ...')

                    resample_times_list.append(t)

                except:
                    # stop memory recording (due to errors)
                    mem_mon.stop_recording()
                    mem_mon.stop()
                    logging.error(' ===> Resampling data failed at the "' + t + '" time step')
                    raise RuntimeError('Algorithm failed in resampling data. Check the datasets')

                # stop memory recording
                mem_mon.stop_recording()

                # manage memory available
                try:
                    if use_memon:
                        memory_full = not mem_mon.memory_available()
                    else:
                        memory_full = False
                except ValueError:
                    memory_full = False

                # dump data to disk if free memory is short
                if ((memory_full is True or
                     len(resample_times_list) > write_n_resampled) and
                        len(gpi_data) > 0):

                    logging.info(' -----> Write resampled data to disk (memory full) ... ')
                    gpi_data = np.hstack(gpi_data)
                    orbit_data = np.hstack(orbit_data)

                    self.write_resampled_data(
                        gpi_data, orbit_data,
                        init_file_cell=init_file_cell, use_file_cell=use_file_cell, name_file_cell=name_file_cell,
                        name_file_ts=name_file_ts)

                    all_resampled_times.extend(resample_times_list)
                    resample_times_list = []
                    gpi_data = []
                    orbit_data = []

                    logging.info(' ----> Write resampled data to disk (memory full) ... DONE')

                mem_mon.clear_recording_history()

                # time info end
                logging.info(' ----> Resample data from "' + ts_start + '" to "' + ts_end + '" ... DONE')

            # save ancillary workspace file
            logging.info(' ----> Save ancillary workspace ... ')
            if use_file_ws:
                data_ws = {'gpi': gpi_data, 'orbit': orbit_data, 'time': resample_times_list, 'timefield': 'jd'}
                write_file_obj(name_file_ws, data_ws)
                logging.info(' ----> Save ancillary workspace ... DONE')
            else:
                logging.info(' ----> Save ancillary workspace ... SKIPPED. Flag is not activated')
        else:

            # read ancillary workspace file (previously created)
            logging.info(' ----> Read ancillary workspace ... ')
            file_obj = read_file_obj(name_file_ws)

            if file_obj is None:
                logging.error(' ===> Workspace is defined by NoneType object')
                raise RuntimeError('Workspace obj must be defined by dictionary object with satellite information')

            gpi_data, orbit_data, \
                resample_times_list, self.timefield = file_obj['gpi'], file_obj['orbit'], \
                                                      file_obj['time'], file_obj['timefield']
            logging.info(' ----> Read ancillary workspace ... DONE')

        # write remaining data to disk
        logging.info(' ----> Write resampled data to disk ... ')
        if len(gpi_data) > 0:
            gpi_data = np.hstack(gpi_data)
            orbit_data = np.hstack(orbit_data)
            self.write_resampled_data(
                gpi_data, orbit_data,
                init_file_cell=init_file_cell, use_file_cell=use_file_cell, name_file_cell=name_file_cell,
                name_file_ts=name_file_ts)
            all_resampled_times.extend(resample_times_list)

            logging.info(' ----> Write resampled data to disk ... DONE')
        else:
            logging.info(' ----> Write resampled data to disk ... SKIPPED. No data found')

        mem_mon.stop()
        self.orbit_io.close()
        self.resam_io.close()
        self.save_orbit_buffer()

        return all_resampled_times

    def write_resampled_data(self, gpi_data, orbit_data,
                             init_file_cell=True, use_file_cell=True,
                             name_file_cell='{cell_n}.data',
                             name_file_ts='{cell_n}.nc'):
        """
        write the data.

        Parameters
        ----------
        gpi_data : numpy.recarray
            Grid point information. A grid point indices
            for each element in the orbit data array.
        orbit_data : numpy.recarray
            record array containing the resampled data at the grid point
            indices of gpi_data
        """
        cells = self.target_grid.gpi2cell(gpi_data)
        indSort = np.argsort(cells)
        gpi_data = gpi_data[indSort]
        orbit_data = orbit_data[indSort]

        cells = cells[indSort]
        uniq_cells, st_index, n_meas = np.unique(cells, return_index=True, return_counts=True)
        n_cells = uniq_cells.size
        end_index = st_index + n_meas

        gpi_act = self.resam_io.grid.activegpis
        gpi_idx = np.apply_along_axis(lambda f: gpi_act.searchsorted(f), 0, gpi_data)

        # iterate over cell(s)
        for i in np.arange(0, n_cells):

            logging.info(' ------> Dump cell: ' + str(uniq_cells[i]) + ' ... ')

            cell_i = uniq_cells[i]

            gpi_i_local = gpi_idx[st_index[i]:end_index[i]]
            orbit_i = orbit_data[st_index[i]:end_index[i]]

            lon_i_global = self.resam_io.grid.arrlon[gpi_i_local]
            lat_i_global = self.resam_io.grid.arrlat[gpi_i_local]
            gpi_i_global = self.resam_io.grid.gpis[gpi_i_local]

            if use_file_cell:
                name_file_step = name_file_cell.format(cell_n=cell_i)

                if init_file_cell:
                    if os.path.exists(name_file_step):
                        os.remove(name_file_step)

                cell_obj = {'gpi': gpi_i_local, 'orbit': orbit_i}
                write_file_obj(name_file_step, cell_obj)


            name_file_ts = name_file_ts.format(cell_n=cell_i)

            self.resam_io.write_cell(cell_i, gpi_i_local, orbit_i, self.timefield)

            logging.info(' ------> Dump cell: ' + str(uniq_cells[i]) + ' ... DONE')

    @staticmethod
    def rectDictToNdarray(dd):
        """
        Method to convert a rectangular shaped dictionary into a structured
        nd array.

        Parameters
        ----------
        dd : dict
            Rectangluar dictionary to convert.

        Returns
        -------
        data : numpy.array (structured)
            Converted dictionary as numpy array.
        """
        dtypeList = []
        for key in dd.keys():
            dtypeList.append((key, dd[key].dtype.type))

        if sys.version_info[0] == 2:
            data = np.empty(dd[dd.keys()[0]].size, dtype=np.dtype(dtypeList))
        elif sys.version_info[0] == 3:
            var_tmp = list(dd.keys())[0]
            data = np.empty(dd[var_tmp].size, dtype=np.dtype(dtypeList))
        else:
            raise RuntimeError('python version "' + str(sys.version_info[0]) + '" is not allowed')

        for key in dd.keys():
            data[key] = dd[key]

        return data

    @staticmethod
    def rectDsetToNdarray(dset):

        dtypeList = []
        for var in list(dset.variables):
            if var != 'time':
                #dtypeList.append((var, dset[var].dtype.type))
                dtypeList.append((var, np.dtype(np.float64)))

        n = dset[list(dset.variables)[0]].size
        data = np.empty(n, dtype=np.dtype(dtypeList))
        for var in list(dset.variables):
            if var != 'time':
                if dset[var].values.ndim == 1:
                    data[var] = dset[var].values
                else:
                    no_data = np.zeros((n,))
                    no_data[:] = -9999
                    data[var] = no_data
                    warnings.warn(' variable ' + var + ' skipped due to dims or type')

        return data

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get window radius
def get_window_radius(window, hp_radius):
    """
    Calculates the required radius of a window function in order to achieve
    the provided half power radius.

    Parameters
    ----------
    window : string
        Window function name.
        Current supported windows:
            - Hamming
            - Boxcar
    hp_radius : float32
        Half power radius. Radius of window function for weight
        equal to 0.5 (-3 dB). In the spatial domain this corresponds to
        half of the spatial resolution one would like to achieve with the
        given window.
    Returns
    -------
    r : float32
        Window radius needed to achieve the given half power radius

    """
    window = window.lower()
    hp_weight = 0.5
    if window == 'hamming':
        alpha = 0.54
        r = (np.pi * hp_radius) / np.arccos((hp_weight - alpha) / (1 - alpha))
    elif window == 'boxcar':
        r = hp_radius
    else:
        raise ValueError('Window name not supported.')

    return r
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set hamming window
def hamming_window(radius, distances):
    """
    Hamming window filter.

    Parameters
    ----------
    radius : float32
        Radius of the window.
    distances : numpy.ndarray
        Array with distances.

    Returns
    -------
    weights : numpy.ndarray
        Distance weights.
    tw : float32
        Sum of weigths.
    """
    alpha = 0.54
    weights = alpha + (1 - alpha) * np.cos(np.pi / radius * distances)

    return weights, np.sum(weights)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set boxcar window
def boxcar(radius, distance):
    """
    Boxcar filter

    Parameters
    ----------
    n : int
        Length.

    Returns
    -------
    weights : numpy.ndarray
        Distance weights.
    tw : float32
        Sum of weigths.
    """
    weights = np.zeros(distance.size)
    weights[distance <= radius] = 1.

    return weights, np.sum(weights)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get window weight(s)
def get_window_weights(window, radius, distance, norm=False):
    """
    Function returning weights for the provided window function

    Parameters
    ----------
    window : str
        Window function name
    radius : float
        Radius of the window.
    distance : numpy.ndarray
        Distance array
    norm : boolean
        If true, normalised weights will be returned.

    Returns
    -------
    weights : numpy.ndarray
        Weights according to distances and given window function

    """
    if window == 'hamming':
        weights, w_sum = hamming_window(radius, distance)
    elif window == 'boxcar':
        weights, w_sum = boxcar(radius, distance)
    else:
        raise ValueError('Window name not supported.')

    if norm is True:
        weights = weights / w_sum

    return weights
# ----------------------------------------------------------------------------------------------------------------------
