
"""
Library Features:

Name:          lib_datasets_ascat
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import tempfile

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from lib_utils_generic import read_obj, write_obj
from lib_interface_ascat import AscatCdr

# import pytesmo.temporal_matching as temp_match

from lib_utils_generic import get_bit
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap ASCAT nrt time series
class ASCAT_Dataset_NRT(AscatCdr):

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

        # set porosity file
        self._por_file = 'porosity.nc'
        if 'file_porosity' in kwargs:
            self._por_file = kwargs.pop('file_porosity')

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        super(ASCAT_Dataset_NRT, self).__init__(
            cdr_path=self.dr_path,
            grid_path=self.grid_path, grid_filename='TUW_WARP5_grid_info_2_3.nc',
            static_layer_path=self.static_layer_path,
            **kwargs)

    # method to set tmp info
    @staticmethod
    def _set_tmp_info(tmp_info):

        active_tmp = False
        if 'active_tmp' in list(tmp_info.keys()):
            active_tmp = tmp_info['active_tmp']
        path_tmp = tempfile.mkdtemp()
        if 'path_tmp' in list(tmp_info.keys()):
            path_tmp = tmp_info['path_tmp']
        temp_obj = tempfile.NamedTemporaryFile(prefix='ascat_', suffix='.dframe')
        _, file_tmp = os.path.split(temp_obj.name)
        if 'file_tmp' in list(tmp_info.keys()):
            file_tmp = tmp_info['file_tmp']
        clean_tmp = True
        if 'clean_tmp' in list(tmp_info.keys()):
            clean_tmp = tmp_info['clean_tmp']

        if active_tmp:
            os.makedirs(path_tmp, exist_ok=True)

        return active_tmp, os.path.join(path_tmp, file_tmp), clean_tmp

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

    def read(self, *args, **kwargs):
        """
        Method to read time series and mask the data and convert it volumetric
        soil moisture information by making use of the porosity layer.
        """
        try:

            active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp

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

                n_value = ts.loc['2018':'2019'].shape
                print(' ----> N Data valid for ASCAT dataset: ' + str(n_value))

        except Exception as exc:
            print(' ----> WARNING: RunTime error for ASCAT dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap ASCAT data record
class ASCAT_Dataset_DR(AscatCdr):

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

        # set porosity file
        self._por_file = 'porosity.nc'
        if 'file_porosity' in kwargs:
            self._por_file = kwargs.pop('file_porosity')

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        # set read bulk
        self.read_bulk = False
        if "bulk" in kwargs:
            self.read_bulk = kwargs.pop('bulk')

        # get porosity information
        self._por_path = os.path.join(self.static_layer_path, self._por_file)
        self._read_porosity()

        super(ASCAT_Dataset_DR, self).__init__(
            cdr_path=self.dr_path,
            grid_path=self.grid_path, grid_filename='TUW_WARP5_grid_info_2_2.nc',
            static_layer_path=self.static_layer_path,
            read_bulk=self.read_bulk,
            **kwargs)

    # method to set tmp info
    @staticmethod
    def _set_tmp_info(tmp_info):

        active_tmp = False
        if 'active_tmp' in list(tmp_info.keys()):
            active_tmp = tmp_info['active_tmp']
        path_tmp = tempfile.mkdtemp()
        if 'path_tmp' in list(tmp_info.keys()):
            path_tmp = tmp_info['path_tmp']
        temp_obj = tempfile.NamedTemporaryFile(prefix='ascat_', suffix='.dframe')
        _, file_tmp = os.path.split(temp_obj.name)
        if 'file_tmp' in list(tmp_info.keys()):
            file_tmp = tmp_info['file_tmp']
        clean_tmp = True
        if 'clean_tmp' in list(tmp_info.keys()):
            clean_tmp = tmp_info['clean_tmp']

        if active_tmp:
            os.makedirs(path_tmp, exist_ok=True)

        return active_tmp, os.path.join(path_tmp, file_tmp), clean_tmp

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

    # method to read ascat time-series
    def read(self, *args, **kwargs):

        # get tmp information
        active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp
        gpi = args[0]

        # info time-series
        logging.info(' ------> Read ASCAT time-series for GPI "' + str(gpi) + '" ... ')

        # organize dataframe
        try:

            file_gpi = file_tmp.format(gpi_n=str(gpi))
            if clean_tmp:
                if os.path.exists(file_gpi):
                    os.remove(file_gpi)

            if not os.path.exists(file_gpi):

                # info get data start
                logging.info(' -------> Get data ... ')
                ts_dframe = super(ASCAT_Dataset_DR, self).read(*args, **kwargs)

                if ts_dframe is None:
                    logging.warning(' ===> ASCAT Dataset is defined by NoneType')
                    logging.warning(' ===> ASCAT time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()

                    # info get data end (skipped)
                    logging.info(' -------> Get data ... SKIPPED. ASCAT time-series is not defined.')
                else:

                    # info get data end
                    logging.info(' -------> Get data ... DONE')

                    # info apply ssf start
                    logging.info(' -------> Apply ssf ...')
                    bit_mask = ((get_bit(ts_dframe['proc_flag'], 1)) |
                                (get_bit(ts_dframe['proc_flag'], 2)) |
                                (get_bit(ts_dframe['proc_flag'], 3)))

                    # check surface state flag (available in the data record published by the authors))
                    if 'ssf' in ts_dframe.columns:
                        ts_dframe = ts_dframe[((ts_dframe['ssf'] == 0) | (ts_dframe['ssf'] == 1)) & (bit_mask == 0)]
                        # info apply ssf end
                        logging.info(' -------> Apply ssf ... DONE')
                    else:
                        logging.warning(' ===> Surface state flag "ssf" not available in ASCAT dataset')
                        # info apply ssf end
                        logging.info(' -------> Apply ssf ... SKIPPED')

                    # info apply porosity start
                    logging.info(' -------> Apply porosity ...')

                    # convert to absolute soil moisture
                    porosity = self.get_porosity(*args)

                    # check porosity value
                    if porosity is not np.nan:
                        # update sm and sm_noise
                        ts_dframe['sm'] = ts_dframe['sm'] / 100. * porosity
                        ts_dframe['sm_noise'] = ts_dframe['sm_noise'] / 100. * porosity
                        # info apply porosity end
                        logging.info(' -------> Apply porosity ... DONE')
                    else:
                        logging.warning(' ===> No porosity valid for ASCAT dataset')
                        logging.warning(' ===> ASCAT time-series will be initialized by empty dataframe')

                        # info apply porosity end
                        logging.info(' -------> Apply porosity ... FAILED')
                        ts_dframe = pd.DataFrame()

                    # check data valid
                    logging.info(' -------> Check data ... ')
                    if ts_dframe.size == 0:
                        logging.warning(' ===> No data valid for ASCAT dataset')
                        logging.warning(' ===> ASCAT time-series will be initialized by empty dataframe')
                        ts_dframe = pd.DataFrame()
                        logging.info(' -------> Check data ... FAILED. No data available')
                    else:
                        valid_value = ts_dframe.shape[0]
                        start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                        end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                        logging.info(' --------> Data valid for ASCAT dataset (N: "' + str(valid_value) + '")')
                        logging.info(' --------> Time valid for ASCAT dataset from "' + start_index +
                                     '" to "' + end_index + '"')
                        logging.info(' -------> Check data ... DONE')

                # save time-series if flag is active
                if active_tmp:
                    write_obj(file_gpi, ts_dframe)

                # info time-series
                logging.info(' ------> Read ASCAT time-series for GPI "' + str(gpi) + '" ... DONE ')

            else:

                # read time-series previously save
                ts_dframe = read_obj(file_gpi)

                # check data valid
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for ASCAT dataset')
                    logging.warning(' ===> ASCAT time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' -------> Data valid for ASCAT dataset (N: "' + str(valid_value) + '")')
                    logging.info(' -------> Time valid for ASCAT dataset from "' + start_index +
                                 '" to "' + end_index + '"')

                # info time-series
                logging.info(' ------> Read ASCAT time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

        except BaseException as base_exc:
            logging.warning(' ===> Error in reading ASCAT time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> ASCAT time-series will be initialized by empty dataframe')
            ts_dframe = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read ASCAT time-series for GPI "' + str(gpi) + '" ... FAILED')

        # ts must be defined by dataframe
        if ts_dframe is None:
            logging.warning(' ===> ASCAT time-series is defined by NoneType. It will be format using DataFrame')
            ts_dframe = pd.DataFrame()

        return ts_dframe
# ----------------------------------------------------------------------------------------------------------------------
