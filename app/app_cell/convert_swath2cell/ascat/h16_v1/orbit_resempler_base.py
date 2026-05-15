import os
import warnings
from functools import partial

import numpy as np
import psutil

from pygeobase.io_base import MultiTemporalImageBase, GriddedTsBase
from pykdtree.kdtree import KDTree


class MemoryMonitor(object):
    """
    Lightweight replacement for memon.MemoryMonitor using psutil.
    """

    def __init__(self, min_available_mb=500):
        self.min_available_mb = min_available_mb

    def start(self):
        pass

    def stop(self):
        pass

    def start_recording(self):
        pass

    def stop_recording(self):
        pass

    def clear_recording_history(self):
        pass

    def memory_available(self):
        mem = psutil.virtual_memory()
        available_mb = mem.available / 1024.0 / 1024.0
        return available_mb > self.min_available_mb


class OrbitResampler(object):
    """
    Resample Orbit/Swath data onto a target grid defined by the resam_io
    object. The resam_io object has to implement a write_cell method
    and must also work on a pygeogrids.grids.CellGrid grid object.
    """

    def __init__(self, orbit_io, resam_io, spatial_res=25000.,
                 wfunc='hamming', min_num_nn=3, max_num_nn=None, dt=15,
                 write_orbit_buffer=False,
                 resample_first=True):

        if not isinstance(orbit_io, MultiTemporalImageBase):
            name = MultiTemporalImageBase.__name__
            raise ValueError('Orbit data must be of type {:}'.format(name))

        if not isinstance(resam_io, GriddedTsBase):
            name = GriddedTsBase.__name__
            raise ValueError('Output data must be of type {:}'.format(name))

        self.orbit_io = orbit_io
        self.resam_io = resam_io

        self.dt = dt / 60. / 24.
        self.min_num_nn = min_num_nn

        self.resample_first = resample_first
        self.count_instructed_resampled = 0

        self.buffer_saved = False
        self.write_orbit_buffer = write_orbit_buffer

        self.orbit_buffer_file = os.path.join(resam_io.path, "orbitBuffer.npz")

        self.buffer_keys = [
            'jds_prev',
            'data_prev',
            'nn_flag_prev',
            'index_prev',
            'distance_prev'
        ]

        self.orbit_buffer = {}

        if os.path.isfile(self.orbit_buffer_file):
            try:
                prev = np.load(self.orbit_buffer_file, allow_pickle=True)

                if all(key in prev.files for key in self.buffer_keys):
                    for key in self.buffer_keys:
                        self.orbit_buffer[key] = prev[key]
                else:
                    warnings.warn(
                        "Orbit buffer file is empty or incompatible: {}. "
                        "Reinitializing buffer. Available keys: {}".format(
                            self.orbit_buffer_file,
                            prev.files
                        ),
                        UserWarning
                    )

                    for key in self.buffer_keys:
                        self.orbit_buffer[key] = None

            except Exception as exc:
                warnings.warn(
                    "Could not read orbit buffer file: {}. "
                    "Reinitializing buffer. Error: {}".format(
                        self.orbit_buffer_file,
                        exc
                    ),
                    UserWarning
                )

                for key in self.buffer_keys:
                    self.orbit_buffer[key] = None

        else:
            for key in self.buffer_keys:
                self.orbit_buffer[key] = None

        self.window_function_name = wfunc
        self.spatial_res = spatial_res

        self.search_radius = get_window_radius(
            self.window_function_name,
            self.spatial_res / 2.
        )

        self.target_grid = resam_io.grid
        self.ref_geodatum = self.target_grid.geodatum

        self.xyz_target_grid = self.ref_geodatum.toECEF(
            self.target_grid.activearrlon,
            self.target_grid.activearrlat
        )

        orig_window_function = {
            'hamming': hamming_window,
            'boxcar': boxcar
        }[self.window_function_name]

        self.window_function = partial(
            orig_window_function,
            self.search_radius
        )

        self.xyz_target_grid = np.transpose(np.vstack(self.xyz_target_grid))

        if max_num_nn is None:
            tmp_diff = self.xyz_target_grid[0, :] - self.xyz_target_grid

            tmp_dist = np.sqrt(
                tmp_diff[:, 0] ** 2 +
                tmp_diff[:, 1] ** 2 +
                tmp_diff[:, 2] ** 2
            )

            self.n_neighbours = 3 * np.sum(tmp_dist <= self.search_radius)

        else:
            self.n_neighbours = max_num_nn

    def __enter__(self):
        return self

    def close(self, *args):
        return self.save_orbit_buffer()

    def __exit__(self, exc_type, exc_value, traceback):
        self.save_orbit_buffer()
        return self.buffer_saved

    def __del__(self):
        self.save_orbit_buffer()
        return self.buffer_saved

    def save_orbit_buffer(self):
        if self.buffer_saved is False and self.write_orbit_buffer is True:

            if os.path.isdir(os.path.dirname(self.orbit_buffer_file)) is False:
                os.makedirs(os.path.dirname(self.orbit_buffer_file))

            if os.path.isfile(self.orbit_buffer_file):
                os.remove(self.orbit_buffer_file)

            np.savez_compressed(
                self.orbit_buffer_file,
                **self.orbit_buffer
            )

            self.buffer_saved = True

    def nn_search(self, lon, lat, **kwargs):
        xyz = np.transpose(np.vstack(self.ref_geodatum.toECEF(lon, lat)))

        init_KDTree = KDTree(xyz)

        distance, index = init_KDTree.query(
            self.xyz_target_grid,
            k=self.n_neighbours,
            distance_upper_bound=self.search_radius,
            **kwargs
        )

        nn_flag = np.any(distance != np.inf, axis=1)
        nn_num = np.sum(distance != np.inf, axis=1)

        return nn_flag, nn_num, index, distance

    def resample_orbit(self, timestamp):
        resampled_gpis = []
        resampled_data = []

        self.count_instructed_resampled = self.count_instructed_resampled + 1

        try:
            orbit = self.orbit_io.read(timestamp, mask=True)

            if orbit is None:
                raise Exception("No data for timestamp")

        except Exception as ex:
            print(ex)
            return np.array([]), np.array([])

        data = self.rectDictToNdarray(orbit.data)
        self.timefield = orbit.timekey
        orbit_jd = orbit.data[orbit.timekey]
        lon = orbit.lon
        lat = orbit.lat

        if not data.__len__() > 0:
            return np.array([]), np.array([])

        jds = np.hstack((orbit_jd, np.nan))

        nn_flag, nn_num, index, distance = self.nn_search(lon, lat)

        jd_diff = np.nanmax(jds[index], axis=1) - np.nanmin(jds[index], axis=1)

        n_overlapping_points = np.nansum(
            jd_diff[np.isfinite(jd_diff)] > self.dt
        )

        if n_overlapping_points > 0:
            warnings.warn(
                "{} points overlap in one sub-orbit".format(
                    n_overlapping_points
                ),
                UserWarning
            )

        if self.orbit_buffer['data_prev'] is not None and \
                np.sum(self.orbit_buffer['nn_flag_prev']) > 0.:

            pos = self.orbit_buffer['nn_flag_prev'] & nn_flag

            jd_diff = np.abs(
                np.nanmean(jds[index[pos, :]], axis=1) -
                np.nanmean(
                    self.orbit_buffer['jds_prev'][
                        self.orbit_buffer['index_prev'][pos, :]
                    ],
                    axis=1
                )
            )

            include_gpis = pos.copy()
            include_gpis[include_gpis] = jd_diff <= self.dt

            exclude_gpis = pos.copy()
            exclude_gpis[exclude_gpis] = jd_diff > self.dt

            data_res = np.hstack((
                self.orbit_buffer['data_prev'],
                data,
                data[0]
            ))

            index_tmp = index.copy()

            shift_index = self.orbit_buffer['data_prev'].shape[0]
            index_tmp += shift_index

            index_res = np.hstack((
                self.orbit_buffer['index_prev'],
                index_tmp
            ))

            prevNaN = np.where(
                self.orbit_buffer['index_prev'] == shift_index
            )

            index_res[prevNaN] = data_res.shape[0] - 1

            distance_tmp = distance.copy()
            distance_tmp[exclude_gpis, :] = np.inf

            distance_res = np.hstack((
                self.orbit_buffer['distance_prev'],
                distance_tmp
            ))

            gpis_to_resample = self.orbit_buffer['nn_flag_prev']

            if self.count_instructed_resampled == 2 and not self.resample_first:
                gpis_to_resample = include_gpis

            index_res = index_res[gpis_to_resample, :]
            distance_res = distance_res[gpis_to_resample, :]

            num_nn = (
                np.sum(np.isfinite(distance_res), axis=1) >= self.min_num_nn
            )

            if np.sum(num_nn) > 0.:
                index_res = index_res[num_nn, :]
                distance_res = distance_res[num_nn, :]

                resampled_gpis.append(
                    self.resam_io.grid.activegpis[gpis_to_resample][num_nn]
                )

                weights_res, _ = self.window_function(distance_res)

                resdata = self.orbit_io.resample_image(
                    data_res,
                    index_res,
                    distance_res,
                    weights_res
                )

                if type(resdata) == dict:
                    resdata = self.rectDictToNdarray(resdata)

                resampled_data.append(resdata)

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

    def resample(self, time_stamps,
                 write_n_resampled=14 * 2 * 365,
                 use_memon=True):

        orbit_data = []
        gpi_data = []
        resample_times_list = []
        all_resampled_times = []

        mem_mon = MemoryMonitor()
        mem_mon.start()

        for t in time_stamps:
            print(t)

            try:
                mem_mon.start_recording()

                gpis, orbit = self.resample_orbit(t)

                if gpis.size > 0 and orbit.size > 0:
                    gpi_data.append(gpis)
                    orbit_data.append(orbit)

                resample_times_list.append(t)

            except Exception:
                mem_mon.stop_recording()
                mem_mon.stop()
                raise

            mem_mon.stop_recording()

            try:
                if use_memon:
                    memory_full = not mem_mon.memory_available()
                else:
                    memory_full = False

            except ValueError:
                memory_full = False

            if ((memory_full is True or
                 len(resample_times_list) > write_n_resampled) and
                    len(gpi_data) > 0):

                gpi_data = np.hstack(gpi_data)
                orbit_data = np.hstack(orbit_data)

                self.write_resampled_data(gpi_data, orbit_data)

                all_resampled_times.extend(resample_times_list)

                resample_times_list = []
                gpi_data = []
                orbit_data = []

            mem_mon.clear_recording_history()

        if len(gpi_data) > 0:
            gpi_data = np.hstack(gpi_data)
            orbit_data = np.hstack(orbit_data)

            self.write_resampled_data(gpi_data, orbit_data)

            all_resampled_times.extend(resample_times_list)

        mem_mon.stop()

        self.orbit_io.close()
        self.resam_io.close()

        self.save_orbit_buffer()

        return all_resampled_times

    def write_resampled_data(self, gpi_data, orbit_data):
        cells = self.target_grid.gpi2cell(gpi_data)

        indSort = np.argsort(cells)

        gpi_data = gpi_data[indSort]
        orbit_data = orbit_data[indSort]

        cells = cells[indSort]

        uniq_cells, st_index, n_meas = np.unique(
            cells,
            return_index=True,
            return_counts=True
        )

        n_cells = uniq_cells.size
        end_index = st_index + n_meas

        gpi_act = self.resam_io.grid.activegpis

        gpi_idx = np.apply_along_axis(
            lambda f: gpi_act.searchsorted(f),
            0,
            gpi_data
        )

        for i in np.arange(0, n_cells):
            self.resam_io.write_cell(
                uniq_cells[i],
                gpi_idx[st_index[i]:end_index[i]],
                orbit_data[st_index[i]:end_index[i]],
                self.timefield
            )

    @staticmethod
    def rectDictToNdarray(dd):
        dtypeList = []

        for key in dd.keys():
            dtypeList.append((key, dd[key].dtype.type))

        first_key = list(dd.keys())[0]

        data = np.empty(
            dd[first_key].size,
            dtype=np.dtype(dtypeList)
        )

        for key in dd.keys():
            data[key] = dd[key]

        return data


def get_window_radius(window, hp_radius):
    window = window.lower()
    hp_weight = 0.5

    if window == 'hamming':
        alpha = 0.54

        r = (
            np.pi * hp_radius
        ) / np.arccos(
            (hp_weight - alpha) / (1 - alpha)
        )

    elif window == 'boxcar':
        r = hp_radius

    else:
        raise ValueError('Window name not supported.')

    return r


def hamming_window(radius, distances):
    alpha = 0.54

    weights = alpha + \
        (1 - alpha) * \
        np.cos(np.pi / radius * distances)

    return weights, np.sum(weights)


def boxcar(radius, distance):
    weights = np.zeros(distance.size)

    weights[distance <= radius] = 1.

    return weights, np.sum(weights)


def get_window_weights(window, radius, distance, norm=False):
    if window == 'hamming':
        weights, w_sum = hamming_window(radius, distance)

    elif window == 'boxcar':
        weights, w_sum = boxcar(radius, distance)

    else:
        raise ValueError('Window name not supported.')

    if norm is True:
        weights = weights / w_sum

    return weights