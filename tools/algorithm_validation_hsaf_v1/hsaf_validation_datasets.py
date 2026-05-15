
"""
This module implements datasets used for validation.
"""

# ----------------------------------------------------------------------------------------------------------------------
# Libraries
import os

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from pandas import read_pickle

try:
    from ascat.timeseries import AscatSsmCdr
except BaseException:
    from ascat.read_native.cdr import AscatSsmCdr
from gldas.interface import GLDASTs

from library_cci.interface import CCITs
from library_rzsm.interface import RZSMTs

#from library.cci.interface import CCITs
# from library.rzsm.interface import RZSMTs

# import pytesmo.temporal_matching as temp_match

from hsaf_validation_utils import get_bit
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap RZSM time series
class RZSM_Dataset(RZSMTs):

    saturation = None

    # Method to init class
    def __init__(self, **kwargs):

        # data record path
        if "dr_path" in kwargs:
            self.dr_path = kwargs.pop('dr_path')
        else:
            raise ValueError('data record path not defined')

        # initialise static layer paths
        if 'grid_path' in kwargs:
            self.grid_path = kwargs.pop('grid_path')
        else:
            raise ValueError('grid layer path not defined')

        if "temp_path" in kwargs:
            self.temp_path = kwargs.pop('temp_path')

            if not os.path.exists(self.temp_path):
                os.makedirs(self.temp_path)
        else:
            raise ValueError('temp path not defined')

        if "flag_clean_tmp" in kwargs:
            self.flag_clean_tmp = kwargs.pop('flag_clean_tmp')
        else:
            raise ValueError('flag_clean_tmp not defined')

        super(RZSM_Dataset, self).__init__(ts_path=self.dr_path, **kwargs)

    # Method to read time-series
    def read_ts(self, *args, **kwargs):

        try:
            # get rzsm gpi
            gpi_rzsm = args[0]

            temp_file = os.path.join(self.temp_path, 'rzsm_ts_gpi_' + str(gpi_rzsm) + '.pkl')

            if self.flag_clean_tmp:
                os.remove(temp_file)

            if not os.path.exists(temp_file):
                ts = super(RZSM_Dataset, self).read(*args)
                ts.to_pickle(temp_file)
            else:
                ts = read_pickle(temp_file)

            if ts.size == 0:
                print(' ----> WARNING: No data valid for RZSM dataset')
                ts = pd.DataFrame()

            else:
                print(' ----> Data valid for RZSM dataset')

                n_value = ts.shape[0]

                if n_value > 0:
                    time_start = ts.index[0].strftime('%Y%m%d')
                    time_end = ts.index[-1].strftime('%Y%m%d')
                else:
                    time_start, time_end = 'NA', 'NA'

                print(' ----> N Data valid for RZSM dataset: ' + str(n_value))
                print(' ----> TimePeriod Data valid for RZSM dataset: from ' + time_start + ' to ' + time_end)

        except Exception as exc:
            print(' ----> WARNING: RunTime error for RZSM dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap ASCAT nrt time series
class ASCAT_Dataset_NRT(AscatSsmCdr):

    por = None

    def __init__(self, *args, **kwargs):

        # initialise static layer paths
        if 'dr_path' in kwargs:
            self.dr_path = kwargs.pop('dr_path')
        else:
            raise ValueError('data record path not defined')

        # initialise static layer paths
        if 'grid_path' in kwargs:
            self.grid_path = kwargs.pop('grid_path')
        else:
            raise ValueError('grid layer path not defined')

        # initialise static layer paths
        if 'static_layer_path' in kwargs:
            self.static_layer_path = kwargs.pop('static_layer_path')
        else:
            raise ValueError('static layer path not defined')

        if "flag_clean_tmp" in kwargs:
            self.flag_clean_tmp = kwargs.pop('flag_clean_tmp')
        else:
            raise ValueError('flag_clean_tmp not defined')

        if 'file_porosity' in kwargs:
            self._por_file = kwargs.pop('file_porosity')
        else:
            self._por_file = 'porosity.nc'

        self._por_path = os.path.join(self.static_layer_path, self._por_file)
        self._read_porosity()

        #self.lsm = GLDAS_Dataset(self._lsm_path)
        #self.temp_match = BasicTemporalMatching(window=self.window)

        super(ASCAT_Dataset_NRT, self).__init__(cdr_path=self.dr_path,
                                            grid_path=self.grid_path,
                                            grid_filename='TUW_WARP5_grid_info_2_3.nc',
                                            static_layer_path=self.static_layer_path,
                                            **kwargs)

        #self.luts = self.grid.calc_lut(self.lsm.grid, max_dist=35000)

    def _read_porosity(self):
        """
        Read global porosity from NOAH GLDAS.
        """

        if self.por is None:

            ncFile = Dataset(self._por_path, mode='r')
            por_gpi = ncFile['location_id'][:]
            por = ncFile['por_gldas'][:]
            self.por = por[~por.mask]
            self.por_gpi = por_gpi[~por.mask]

            ncFile.close()

    def get_porosity(self, *args):
        """
        Read porosity for given location.

        Takes either 1 or 2 arguments and calls the correct function
        which is either reading the gpi directly or finding
        the nearest gpi from given lat/lon coordinates and then reading it.

        Returns
        -------
        por : float32
            Porosity.
        """
        if len(args) == 1:
            gpi = args[0]
        if len(args) == 2:
            gpi, _ = self.grid.find_nearest_gpi(args[0], args[1])
        if len(args) < 1 or len(args) > 2:
            raise ValueError('Wrong number of arguments.')

        ind = np.where(self.por_gpi == gpi)[0]

        if ind.size == 0:
            por = np.nan
        else:
            por = self.por[ind]

        return por

    def read_ts(self, *args, **kwargs): # SoilTMP0_10cm_inst
        """
        Method to read time series and mask the data and convert it volumetric
        soil moisture information by making use of the porosity layer.
        """
        try:
            # 2128135
            ts_obj = super(ASCAT_Dataset_NRT, self).read(*args, **kwargs)
            ts = ts_obj.data

            """
            gpi = args[0]
            lsm_gpi = self.luts[gpi]
            lsm_ts = self.lsm.read(lsm_gpi)
            match = temp_match.df_match(ts, lsm_ts, window=0.5)
            """

            # convert to absolute soil moisture
            porosity = self.get_porosity(*args)

            if porosity is not np.nan:
                ts['sm'] = ts['sm'] / 100. * porosity
                ts['sm_noise'] = ts['sm_noise'] / 100. * porosity
            else:
                print(' ----> WARNING: No porosity valid for ASCAT dataset')
                # logging.info(' -----> WARNING: No porosity valid for ASCAT dataset')
                ts = pd.DataFrame()

            if ts.size == 0:
                print(' ----> WARNING: No data valid for ASCAT dataset')
                # logging.info(' ----> WARNING: No data valid for ASCAT dataset')
                ts = pd.DataFrame()
            else:
                print(' ----> Data valid for ASCAT dataset')

                n_value = len(ts)
                print(' ----> N Data valid for ASCAT dataset: ' + str(n_value))

        except Exception as exc:
            print(' ----> WARNING: RunTime error for ASCAT dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap ASCAT data record
class ASCAT_Dataset_DR(AscatSsmCdr):

    por = None

    def __init__(self, *args, **kwargs):

        # initialise static layer paths
        if 'dr_path' in kwargs:
            self.dr_path = kwargs.pop('dr_path')
        else:
            raise ValueError('data record path not defined')

        # initialise static layer paths
        if 'grid_path' in kwargs:
            self.grid_path = kwargs.pop('grid_path')
        else:
            raise ValueError('grid layer path not defined')

        # initialise static layer paths
        if 'static_layer_path' in kwargs:
            self.static_layer_path = kwargs.pop('static_layer_path')
        else:
            raise ValueError('static layer path not defined')

        if 'file_porosity' in kwargs:
            self._por_file = kwargs.pop('file_porosity')
        else:
            self._por_file = 'porosity.nc'

        self._por_path = os.path.join(self.static_layer_path, self._por_file)

        self._read_porosity()

        super(ASCAT_Dataset_DR, self).__init__(cdr_path=self.dr_path,
                                            grid_path=self.grid_path,
                                            grid_filename='TUW_WARP5_grid_info_2_2.nc',
                                            static_layer_path=self.static_layer_path,
                                            **kwargs)

    def _read_porosity(self):
        """
        Read global porosity from NOAH GLDAS.
        """

        if self.por is None:

            ncFile = Dataset(self._por_path, mode='r')
            por_gpi = ncFile['location_id'][:]
            por = ncFile['por_gldas'][:]
            self.por = por[~por.mask]
            self.por_gpi = por_gpi[~por.mask]

            ncFile.close()

        # with Dataset(self._por_path, mode='r') as ncFile:
        #     por_gpi = ncFile['location_id'][:]
        #     por = ncFile['por_gldas'][:]
        #     self.por = por[~por.mask]
        #     self.por_gpi = por_gpi[~por.mask]
        #
        #     ncFile.close()

    def get_porosity(self, *args):
        """
        Read porosity for given location.

        Takes either 1 or 2 arguments and calls the correct function
        which is either reading the gpi directly or finding
        the nearest gpi from given lat/lon coordinates and then reading it.

        Returns
        -------
        por : float32
            Porosity.
        """
        if len(args) == 1:
            gpi = args[0]
        if len(args) == 2:
            gpi, _ = self.grid.find_nearest_gpi(args[0], args[1])
        if len(args) < 1 or len(args) > 2:
            raise ValueError('Wrong number of arguments.')

        ind = np.where(self.por_gpi == gpi)[0]

        if ind.size == 0:
            por = np.nan
        else:
            por = self.por[ind]

        return por

    def read_ts(self, *args, **kwargs):
        """
        Method to read time series and mask the data and convert it volumetric
        soil moisture information by making use of the porosity layer.
        """
        try:
            ts_obj = super(ASCAT_Dataset_DR, self).read(*args, **kwargs)
            ts = ts_obj.data

            bit_mask = ((get_bit(ts['proc_flag'], 1)) |
                        (get_bit(ts['proc_flag'], 2)) |
                        (get_bit(ts['proc_flag'], 3)))

            ts = ts[((ts['ssf'] == 0) | (ts['ssf'] == 1)) & (bit_mask == 0)]

            # convert to absolute soil moisture
            porosity = self.get_porosity(*args)
            
            if porosity is not np.nan:
                ts['sm'] = ts['sm'] / 100. * porosity
                ts['sm_noise'] = ts['sm_noise'] / 100. * porosity
            else:
                print(' ----> WARNING: No porosity valid for ASCAT dataset')
                ts = pd.DataFrame()

            if ts.size == 0:
                print(' ----> WARNING: No data valid for ASCAT dataset')
                ts = pd.DataFrame()
            else:
                print(' ----> Data valid for ASCAT dataset')

                n_value = ts.loc['2018':'2019'].shape
                print(' ----> N Data valid for ASCAT dataset: ' + str(n_value))

        except Exception as exc:
            print(' ----> WARNING: RunTime error for ASCAT dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap CCI time series
class CCI_Dataset(CCITs):

    # Method to init class
    def __init__(self, cci_data_folder, **kwargs):

        self.cci_data_folder = cci_data_folder

        if 'only_valid' in kwargs:
            self.only_valid = kwargs.pop('only_valid')
        else:
            self.only_valid = True

        if 'mask_sm_nan' in kwargs:
            self.mask_sm_nan = kwargs.pop('mask_sm_nan')
        else:
            self.mask_sm_nan = False

        if 'mask_invalid_flags' in kwargs:
            self.mask_invalid_flags = kwargs.pop('mask_invalid_flags')
        else:
            self.mask_invalid_flags = False

        if 'sm_nan' in kwargs:
            self.sm_nan = kwargs.pop('sm_nan')
        else:
            self.sm_nan = -9999.0

        if 'valid_flag' in kwargs:
            self.valid_flag = kwargs.pop('valid_flag')
        else:
            self.valid_flag = 0

        super(CCI_Dataset, self).__init__(ts_path=cci_data_folder, **kwargs)

    # Method to read time-series
    def read_ts(self, *args, **kwargs):

        try:
            ts = super(CCI_Dataset, self).read(*args)

            if self.only_valid:
                self.mask_sm_nan = True
                self.mask_invalid_flags = True
            if self.mask_sm_nan:
                ts = ts[ts['sm'] != self.sm_nan]
            if self.mask_invalid_flags:
                ts = ts[ts['flag'] == self.valid_flag]

            if ts.size == 0:
                print(' ----> WARNING: No data valid for CCI dataset')
                #logging.info(' ----> WARNING: No data valid for CCI dataset')
                ts = pd.DataFrame()
            else:
                print(' ----> Data valid for CCI dataset')

                n_value = len(ts)
                print(' ----> N Data valid for CCI dataset: ' + str(n_value))

        except Exception as exc:
            print(' ----> WARNING: RunTime error for CCI dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap GLDAS time series
class GLDAS_Dataset(GLDASTs):

    # Method to init class
    def __init__(self, gldas_data_folder, **kwargs):

        self.gldas_data_folder = gldas_data_folder

        if 'mask_snow' in kwargs:
            self.mask_snow = kwargs.pop('mask_snow')
        else:
            self.mask_snow = 0

        if 'mask_soil_temp' in kwargs:
            self.mask_soil_temp = kwargs.pop('mask_soil_temp')
        else:
            self.mask_soil_temp = 277.15    # kelvin degree in gldas data (= 4 celsius degree)

        super(GLDAS_Dataset, self).__init__(ts_path=gldas_data_folder, **kwargs)

    # Method to read time-series
    def read_ts(self, *args, **kwargs):

        # Call read time-series
        try:
            ts = super(GLDAS_Dataset, self).read(*args)

            # soil moisture [kg/m^2] converted into [m^3/m^3]
            d = 0.10  # thickness of soil layer in m
            ts['SoilMoi0_10cm_inst'] = ts['SoilMoi0_10cm_inst'] * 0.001 * 1 / d

            # Filter(s)
            ts = ts[ts['SWE_inst'] == self.mask_snow]
            
            try:
                ts = ts[ts['SoilTMP0_10cm_inst'] > self.mask_soil_temp]
            except:
                ts = ts[ts['Tair_f_inst'] > self.mask_soil_temp]

            if ts.size == 0:
                print(' ----> WARNING: No data valid for GLDAS dataset')
                ts = pd.DataFrame()
            else:
                print(' ----> Data valid for GLDAS dataset')

                n_value = ts.loc['2018':'2019'].shape
                print(' ----> N Data valid for GLDAS dataset: ' + str(n_value))

        except Exception as exc:
            print(' ----> WARNING: RunTime error for GLDAS dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------
