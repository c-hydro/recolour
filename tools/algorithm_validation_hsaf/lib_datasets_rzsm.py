
"""
Library Features:

Name:          lib_datasets_rzsm
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import tempfile
import numpy as np
import pandas as pd

from lib_utils_generic import read_obj, write_obj
from lib_interface_rzsm import RZSMTs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to wrap RZSM time series
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

        # tag = 'reference', 'k1', 'k2'
        if 'dset_tag' in kwargs:
            self.tag = kwargs.pop('dset_tag')

        # deprecated
        # if 'sm_nan' in kwargs:
        #     self.sm_nan = kwargs.pop('sm_nan')
        # else:
        #     self.sm_nan = -9999.

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        super(RZSM_Dataset, self).__init__(ts_path=self.dr_path, grid_path=None, **kwargs)

    # method to set tmp info
    @staticmethod
    def _set_tmp_info(tmp_info):

        active_tmp = False
        if 'active_tmp' in list(tmp_info.keys()):
            active_tmp = tmp_info['active_tmp']
        path_tmp = tempfile.mkdtemp()
        if 'path_tmp' in list(tmp_info.keys()):
            path_tmp = tmp_info['path_tmp']
        temp_obj = tempfile.NamedTemporaryFile(prefix='rzsm_', suffix='.dframe')
        _, file_tmp = os.path.split(temp_obj.name)
        if 'file_tmp' in list(tmp_info.keys()):
            file_tmp = tmp_info['file_tmp']
        clean_tmp = True
        if 'clean_tmp' in list(tmp_info.keys()):
            clean_tmp = tmp_info['clean_tmp']

        if active_tmp:
            os.makedirs(path_tmp, exist_ok=True)

        return active_tmp, os.path.join(path_tmp, file_tmp), clean_tmp

    # Method to read time-series
    def read(self, *args, **kwargs):

        # get tmp information
        active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp
        gpi = args[0]

        # info time-series
        logging.info(' ------> Read ECMWF-RZSM time-series for GPI "' + str(gpi) + '" ... ')

        # organize dataframe
        try:

            file_gpi = file_tmp.format(gpi_n=str(gpi))
            if clean_tmp and os.path.exists(file_gpi):
                os.remove(file_gpi)

            if not os.path.exists(file_gpi):

                # read time-series
                logging.info(' -------> Get data ... ')
                ts_dframe = super(RZSM_Dataset, self).read(*args)

                # check time-series extracted using gpi
                if ts_dframe is not None:
                    ts_dframe[ts_dframe['var40'] < 0] = np.nan
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
                # ts_dframe_test = super(GLDAS_Dataset, self).read(lon, lat)

                # check data valid
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for ECMWF-RZSM dataset')
                    logging.warning(' ===> ECMWF-RZSM time-series will be initialized by empty dataframe')

                    ts_dframe = pd.DataFrame()
                    # info get data end (skipped)
                    logging.info(' -------> Get data ... SKIPPED. ECMWF-RZSM time-series is not defined.')
                else:

                    # info get data end
                    logging.info(' -------> Get data ... DONE')

                # info apply mask start
                logging.info(' -------> Apply mask ...')

                # RZSM is expressed as soil wetness index, which takes values in [0, 1]
                mask = (ts_dframe < 0) | (ts_dframe > 1)
                ts_dframe[mask] = np.nan
                # drop nan(s)
                ts_dframe = ts_dframe.dropna()

                # info apply mask end
                logging.info(' -------> Apply mask ... DONE')

                # check data valid
                logging.info(' -------> Check data ... ')
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for ECMWF-RZSM dataset')
                    logging.warning(' ===> ECMWF-RZSM time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                    logging.info(' -------> Check data ... FAILED. No data available')
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' --------> Data valid for ECMWF-RZSM dataset (N: "' + str(valid_value) + '")')
                    logging.info(' --------> Time valid for ECMWF-RZSM dataset from "' + start_index +
                                 '" to "' + end_index + '"')
                    logging.info(' -------> Check data ... DONE')

                # save time-series if flag is active
                if active_tmp:
                    write_obj(file_gpi, ts_dframe)

                # info time-series
                logging.info(' ------> Read ECMWF-RZSM time-series for GPI "' + str(gpi) + '" ... DONE ')

            else:

                # read time-series previously saved
                ts_dframe = read_obj(file_gpi)

                # check data valid
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for ECMWF-RZSM  dataset')
                    logging.warning(' ===> ECMWF-RZSM time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' -------> Data valid for ECMWF-RZSM dataset (N: "' + str(valid_value) + '")')
                    logging.info(' -------> Time valid for ECMWF-RZSM dataset from "' + start_index +
                                 '" to "' + end_index + '"')

                # info time-series
                logging.info(' ------> Read ECMWF-RZSM time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

        except BaseException as base_exc:

            logging.warning(' ===> Error in reading ECMWF-RZSM time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> ECMWF-RZSM time-series will be initialized by empty dataframe')
            ts_dframe = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read ECMWF-RZSM time-series for GPI "' + str(gpi) + '" ... FAILED')

        # ts must be defined by dataframe
        if ts_dframe is None:
            logging.warning(' ===> ECMWF-RZSM time-series is defined by NoneType. It will be format using DataFrame')
            ts_dframe = pd.DataFrame()

        return ts_dframe
# ----------------------------------------------------------------------------------------------------------------------
