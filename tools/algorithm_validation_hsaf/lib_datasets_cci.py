
"""
Library Features:

Name:          lib_datasets_cci
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240318'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import tempfile

import numpy as np
import pandas as pd

from lib_utils_generic import read_obj, write_obj
from lib_interface_cci import CCITs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to wrap CCI time series
class CCI_Dataset(CCITs):

    # method to init class
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

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        super(CCI_Dataset, self).__init__(ts_path=cci_data_folder, **kwargs)

    # method to set tmp info
    @staticmethod
    def _set_tmp_info(tmp_info):

        active_tmp = False
        if 'active_tmp' in list(tmp_info.keys()):
            active_tmp = tmp_info['active_tmp']
        path_tmp = tempfile.mkdtemp()
        if 'path_tmp' in list(tmp_info.keys()):
            path_tmp = tmp_info['path_tmp']
        temp_obj = tempfile.NamedTemporaryFile(prefix='cci_', suffix='.dframe')
        _, file_tmp = os.path.split(temp_obj.name)
        if 'file_tmp' in list(tmp_info.keys()):
            file_tmp = tmp_info['file_tmp']
        clean_tmp = True
        if 'clean_tmp' in list(tmp_info.keys()):
            clean_tmp = tmp_info['clean_tmp']

        if active_tmp:
            os.makedirs(path_tmp, exist_ok=True)

        return active_tmp, os.path.join(path_tmp, file_tmp), clean_tmp

    # method to read CCI time-series
    def read(self, *args, **kwargs):

        # get tmp information
        active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp
        gpi = args[0]

        # info time-series
        logging.info(' ------> Read CCI time-series for GPI "' + str(gpi) + '" ... ')

        # organize dataframe
        try:

            # define tmp file
            file_gpi = file_tmp.format(gpi_n=str(gpi))
            if clean_tmp:
                if os.path.exists(file_gpi):
                    os.remove(file_gpi)

            # check tmp file and read time-series
            if not os.path.exists(file_gpi):

                # read time-series
                logging.info(' -------> Get data ... ')
                ts_dframe = super(CCI_Dataset, self).read(*args)
                # check time-series extracted using gpi
                if ts_dframe is not None:
                    ts_dframe[ts_dframe['sm'] < 0] = np.nan
                    ts_dframe.dropna(inplace=True)
                    if ts_dframe.empty:
                        lon, lat = self.grid.gpi2lonlat(gpi)
                        logging.warning(' ===> ALL DATA ARE DEFINED BY NAN(S) :: POINT :: LON: ' +
                                        str(lon) + ' -- LAT: ' + str(lat) + ' -- GPI: ' + str(gpi))
                        logging.warning(' ===> All data for this gpi are defined by no_data time-series')
                        logging.info(' -------> Get data ... FAILED. Time-series will be initialized by empty dataframe')
                        ts_dframe = pd.DataFrame()
                    else:
                        logging.info(' -------> Get data ... DONE')
                else:
                    lon, lat = self.grid.gpi2lonlat(gpi)
                    logging.warning(' ===> ALL DATA ARE DEFINED BY NONETYPE :: POINT :: LON: ' +
                                    str(lon) + ' -- LAT: ' + str(lat) + ' -- GPI: ' + str(gpi))
                    logging.warning(' ===> All data for this gpi are defined by NoneType')
                    logging.info(' -------> Get data ... FAILED. time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()

                # testing for debugging
                # lon, lat = self.grid.gpi2lonlat(gpi)
                # ts_dframe_test = super(CCI_Dataset, self).read(lon, lat)

                # info apply filter(s) start
                logging.info(' -------> Apply filter(s) ...')
                if not ts_dframe.empty:
                    # filter datasets
                    if self.only_valid:
                        self.mask_sm_nan = True
                        self.mask_invalid_flags = True
                    if self.mask_sm_nan:
                        ts_dframe = ts_dframe[ts_dframe['sm'] != self.sm_nan]
                    if self.mask_invalid_flags:
                        ts_dframe = ts_dframe[ts_dframe['flag'] == self.valid_flag]
                    # info apply filter(s) end
                    logging.info(' -------> Apply filter(s) ... DONE')
                else:
                    # info apply filter(s) end
                    logging.info(' -------> Apply filter(s) ... SKIPPED. Dataframe is empty')

                # check data valid
                logging.info(' -------> Check data ... ')
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for CCI dataset')
                    logging.warning(' ===> CCI time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                    logging.info(' -------> Check data ... FAILED. No data available')
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' -------> Data valid for CCI dataset (N: "' + str(valid_value) + '")')
                    logging.info(' -------> Time valid for CCI dataset from "' + start_index +
                                 '" to "' + end_index + '"')
                    logging.info(' -------> Check data ... DONE')

                # save time-series if flag is active
                if active_tmp:
                    write_obj(file_gpi, ts_dframe)

                # info time-series
                logging.info(' ------> Read CCI time-series for GPI "' + str(gpi) + '" ... DONE ')

            else:

                # read time-series previously saved
                ts_dframe = read_obj(file_gpi)

                # check data valid
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for CCI dataset')
                    logging.warning(' ===> CCI time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' -------> Data valid for CCI dataset (N: "' + str(valid_value) + '")')
                    logging.info(' -------> Time valid for CCI dataset from "' + start_index +
                                 '" to "' + end_index + '"')

                # info time-series
                logging.info(' ------> Read CCI time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

        except BaseException as base_exc:

            logging.warning(' ===> Error in reading CCI time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> CCI time-series will be initialized by empty dataframe')
            ts_dframe = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read CCI time-series for GPI "' + str(gpi) + '" ... FAILED')

        return ts_dframe
# ----------------------------------------------------------------------------------------------------------------------
