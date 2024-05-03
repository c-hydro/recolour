"""
Library Features:

Name:          lib_interface_ascat
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240502'
Version:       '1.5.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import os
import glob

import numpy as np
import pandas as pd

from pygeobase.io_base import GriddedTsBase
from ascat.read_native.cdr import StaticLayers, load_grid
from pynetcf.time_series import OrthoMultiTs, ContiguousRaggedTs

from ascat.read_native.cdr import AscatGriddedNcTs

float32_nan = -999999.0
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to set ascat gridded netcdf ts base
class GriddedNcTs(GriddedTsBase):

    def __init__(self, *args, **kwargs):

        self.parameters = None
        if 'parameters' in kwargs:
            self.parameters = kwargs.pop('parameters')

        self.offsets = None
        if 'offsets' in kwargs:
            self.offsets = kwargs.pop('offsets')

        self.scale_factors = None
        if 'scale_factors' in kwargs:
            self.scale_factors = kwargs.pop('scale_factors')

        self.dtypes = None
        if 'dtypes' in kwargs:
            self.dtypes = kwargs.pop('dtypes')

        self.autoscale = True
        if 'autoscale' in kwargs:
            self.autoscale = kwargs.pop('autoscale')

        self.automask = True
        if 'automask' in kwargs:
            self.automask = kwargs.pop('automask')

        # add options to read_bulk
        self.read_bulk = False
        if 'read_bulk' in kwargs:
            self.read_bulk = kwargs.pop('read_bulk')

        super(GriddedNcTs, self).__init__(*args, **kwargs)

        self.ioclass_kws.update(
            {'autoscale': self.autoscale,
             'automask': self.automask,
             'read_bulk': self.read_bulk
             }
        )

        self.dates = None

        if self.ioclass == OrthoMultiTs:
            self.read_dates = False
        else:
            self.read_dates = True

    def _open(self, gp):
        """
        Open file.

        Parameters
        ----------
        gp : int
            Grid point.

        Returns
        -------
        success : boolean
            Flag if opening the file was successful.
        """
        success = True
        cell = self.grid.gpi2cell(gp)
        filename = os.path.join(self.path,
                                '{:}.nc'.format(self.fn_format.format(cell)))

        if self.mode == 'r':
            if self.previous_cell != cell:
                self.close()

                try:
                    self.fid = self.ioclass(filename, mode=self.mode,
                                            **self.ioclass_kws)
                    self.previous_cell = cell
                except (IOError, RuntimeError):
                    success = False
                    self.fid = None
                    msg = "I/O error {:}".format(filename)
                    warnings.warn(msg, RuntimeWarning)

        if self.mode in ['w', 'a']:
            if self.previous_cell != cell:
                self.flush()
                self.close()

                try:
                    if self.mode == 'w':
                        if 'n_loc' not in self.ioclass_kws:
                            n_loc = self.grid.grid_points_for_cell(cell)[
                                0].size
                            self.ioclass_kws['n_loc'] = n_loc
                    self.fid = self.ioclass(filename, mode=self.mode,
                                            **self.ioclass_kws)
                    self.previous_cell = cell
                    self.ioclass_kws.pop('n_loc', None)
                except (IOError, RuntimeError):
                    success = False
                    self.fid = None
                    msg = "I/O error {:}".format(filename)
                    warnings.warn(msg, RuntimeWarning)

        return success

    def _read_gp(self, gpi, period=None, **kwargs):
        """
        Method reads data for given gpi, additional keyword arguments
        are passed to ioclass.read_ts

        Parameters
        ----------
        gp : int
            Grid point.
        period : list
            2 element array containing datetimes [start, end]

        Returns
        -------
        ts : pandas.DataFrame
            Time series data.
        """
        if self.mode in ['w', 'a']:
            raise IOError("trying to read file is in 'write/append' mode")

        if not self._open(gpi):
            return None

        if self.parameters is None:
            data = self.fid.read_all_ts(gpi, **kwargs)
        else:
            data = self.fid.read_ts(self.parameters, gpi, **kwargs)

        if self.dates is None or self.read_dates:
            if "dates_direct" in kwargs:
                self.dates = self.fid.read_time(gpi)
            else:
                self.dates = self.fid.read_dates(gpi)

        time = self.dates

        # remove time column from dataframe, only index should contain time
        try:
            data.pop('time')
        except KeyError:
            # if the time value is not found then do nothing
            pass

        ts = pd.DataFrame(data, index=time)

        if period is not None:
            ts = ts[period[0]:period[1]]

        if self.dtypes is not None:
            for dtype_column in self.dtypes:
                if dtype_column in ts.columns:
                    try:
                        ts[dtype_column] = ts[dtype_column].astype(
                            self.dtypes[dtype_column])
                    except ValueError:
                        raise ValueError(
                            "Dtype conversion did not work. Try turning off automatic masking.")

        if self.scale_factors is not None:
            for scale_column in self.scale_factors:
                if scale_column in ts.columns:
                    ts[scale_column] *= self.scale_factors[scale_column]

        if self.offsets is not None:
            for offset_column in self.offsets:
                if offset_column in ts.columns:
                    ts[offset_column] += self.offsets[offset_column]

        return ts
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to set ascat gridded netcdf ts contiguous
class GriddedNcContiguousRaggedTs(GriddedNcTs):

    def __init__(self, *args, **kwargs):
        kwargs['ioclass'] = ContiguousRaggedTs
        super(GriddedNcContiguousRaggedTs, self).__init__(*args, **kwargs)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to set ascat data record
class AscatGriddedNcTs(GriddedNcContiguousRaggedTs):

    """
    Class reading Metop ASCAT soil moisture Climate Data Record (CDR).

    Parameters
    ----------
    path : str
        Path to Climate Data Record (CDR) data set.
    fn_format : str
        Filename format string, typical '<prefix>_{:04d}'
    grid_filename : str
        Grid filename.
    static_layer_path : str, optional
        Path to static layer files (default: None).
    thresholds : dict, optional
        Thresholds for topographic complexity (default 50) and
        wetland fraction (default 50).

    Attributes
    ----------
    grid : pygeogrids.CellGrid
        Cell grid.
    thresholds : dict
        Thresholds for topographic complexity (default 50) and
        wetland fraction (default 50).
    slayer : str
        StaticLayer object
    """

    def __init__(self, path, fn_format, grid_filename, static_layer_path=None,
                 cache_static_layer=False, thresholds=None, **kwargs):

        grid = load_grid(grid_filename)

        self.thresholds = {'topo_complex': 50, 'wetland_frac': 50}

        if thresholds is not None:
            self.thresholds.update(thresholds)

        self.slayer = None

        if static_layer_path is not None:
            self.slayer = StaticLayers(static_layer_path,
                                       cache=cache_static_layer)

        super().__init__(path, grid, fn_format=fn_format,
                         **kwargs)

    def _read_gp(self, gpi, **kwargs):
        """
        Read time series for specific grid point.

        Parameters
        ----------
        gpi : int
            Grid point index.
        mask_ssf : boolean, optional
            Default False, if True only SSF values of 1 and 0 will be
            allowed, all others are removed
        mask_frozen_prob : int, optional
            If included in kwargs then all observations taken when
            frozen probability > mask_frozen_prob are removed from the data
            Default: no masking
        mask_snow_prob : int, optional
            If included in kwargs then all observations taken when
            snow probability > mask_snow_prob are removed from the data

        Returns
        -------
        ts : AscatTimeSeries
            Time series object.
        """
        absolute_sm = kwargs.pop('absolute_sm', None)
        mask_frozen_prob = kwargs.pop('mask_frozen_prob', None)
        mask_snow_prob = kwargs.pop('mask_snow_prob', None)
        mask_ssf = kwargs.pop('mask_ssf', None)

        data = super()._read_gp(gpi, **kwargs)
        data.attrs = {}
        data.attrs['gpi'] = gpi
        data.attrs['lon'], data.attrs['lat'] = self.grid.gpi2lonlat(gpi)
        data.attrs['cell'] = self.grid.gpi2cell(gpi)

        if self.slayer is not None:
            data.attrs['topo_complex'] = self.slayer.topo_wetland[gpi]['topo']
            data.attrs['wetland_frac'] = self.slayer.topo_wetland[gpi]['wetland']
            snow_prob = self.slayer.frozen_snow_prob[gpi]['snow_prob']
            frozen_prob = self.slayer.frozen_snow_prob[gpi]['frozen_prob']
            data.attrs['porosity_gldas'] = self.slayer.porosity[gpi]['por_gldas']
            data.attrs['porosity_hwsd'] = self.slayer.porosity[gpi]['por_hwsd']

            if data.attrs['porosity_gldas'] == float32_nan:
                data.attrs['porosity_gldas'] = np.nan

            if data.attrs['porosity_hwsd'] == float32_nan:
                data.attrs['porosity_hwsd'] = np.nan

            if data is not None:
                data['snow_prob'] = snow_prob[data.index.dayofyear - 1]
                data['frozen_prob'] = frozen_prob[data.index.dayofyear - 1]
        else:
            data.attrs['topo_complex'] = np.nan
            data.attrs['wetland_frac'] = np.nan
            data.attrs['porosity_gldas'] = np.nan
            data.attrs['porosity_hwsd'] = np.nan
            data['snow_prob'] = np.nan
            data['frozen_prob'] = np.nan

        if absolute_sm:
            # no error assumed for porosity values, i.e. variance = 0
            por_var = 0.

            data['abs_sm_gldas'] = data['sm'] / \
                100.0 * data.attrs['porosity_gldas']
            data['abs_sm_noise_gldas'] = np.sqrt(
                por_var * (data['sm'] / 100.0)**2 + data['sm_noise']**2 *
                (data.attrs['porosity_gldas'] / 100.0)**2)

            data['abs_sm_hwsd'] = data['sm'] / \
                100.0 * data.attrs['porosity_hwsd']
            data['abs_sm_noise_hwsd'] = np.sqrt(
                por_var * (data['sm'] / 100.0)**2 + data['sm_noise']**2 *
                (data.attrs['porosity_hwsd'] / 100.0)**2)
        else:
            data['abs_sm_gldas'] = np.nan
            data['abs_sm_noise_gldas'] = np.nan
            data['abs_sm_hwsd'] = np.nan
            data['abs_sm_noise_hwsd'] = np.nan

        if mask_ssf is not None:
            data = data[data['ssf'] < 2]

        if mask_frozen_prob is not None:
            data = data[data['frozen_prob'] < mask_frozen_prob]

        if mask_snow_prob is not None:
            data = data[data['snow_prob'] < mask_snow_prob]

        if (data.attrs['topo_complex'] is not None and
                data.attrs['topo_complex'] >= self.thresholds['topo_complex']):
            msg = "Topographic complexity >{:2d} ({:2d})".format(
                self.thresholds['topo_complex'], data.attrs['topo_complex'])
            warnings.warn(msg)

        if (data.attrs['wetland_frac'] is not None and
                data.attrs['wetland_frac'] >= self.thresholds['wetland_frac']):
            msg = "Wetland fraction >{:2d} ({:2d})".format(
                self.thresholds['wetland_frac'], data.attrs['wetland_frac'])
            warnings.warn(msg)

        return data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to interface ascat data record
class AscatCdr(AscatGriddedNcTs):

    """
    Class reading Metop ASCAT soil moisture Climate Data Record (CDR).

    Parameters
    ----------
    cdr_path : str
        Path to Climate Data Record (CDR) data set.
    grid_path : str
        Path to grid file.
    grid_filename : str
        Name of grid file.
    static_layer_path : str
        Path to static layer files.

    Attributes
    ----------
    grid : pygeogrids.CellGrid
        Cell grid.
    """

    def __init__(self, cdr_path, grid_path,
                 grid_filename='TUW_WARP5_grid_info_2_3.nc',
                 static_layer_path=None, read_bulk=False,  **kwargs):

        list_file = glob.glob(os.path.join(cdr_path, '*.nc'))
        if list_file.__len__() == 0:
            logging.error(' ===> ASCAT files not found in the folder "' + cdr_path + '"')
            raise FileNotFoundError('File must be available in the selected folder')

        first_file = list_file[0]
        version = os.path.basename(first_file).rsplit('_', 1)[0]
        fn_format = '{:}_{{:04d}}'.format(version)
        grid_filename = os.path.join(grid_path, grid_filename)

        kwargs = {'read_bulk': read_bulk}  # add options to read_bulk to speed up the reading of the file(s)

        super(AscatCdr, self).__init__(cdr_path, fn_format, grid_filename, static_layer_path, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------
