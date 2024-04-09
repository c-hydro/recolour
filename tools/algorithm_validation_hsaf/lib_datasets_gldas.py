
"""
Library Features:

Name:          lib_datasets_gldas
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import tempfile

import pandas as pd

from lib_utils_generic import read_obj, write_obj
from lib_interface_gldas import GLDASTs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to wrap GLDAS time series
class GLDAS_Dataset(GLDASTs):

    # method to init class
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

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        super(GLDAS_Dataset, self).__init__(ts_path=gldas_data_folder, **kwargs)

    # method to set tmp info
    @staticmethod
    def _set_tmp_info(tmp_info):

        active_tmp = False
        if 'active_tmp' in list(tmp_info.keys()):
            active_tmp = tmp_info['active_tmp']
        path_tmp = tempfile.mkdtemp()
        if 'path_tmp' in list(tmp_info.keys()):
            path_tmp = tmp_info['path_tmp']
        temp_obj = tempfile.NamedTemporaryFile(prefix='gldas_', suffix='.dframe')
        _, file_tmp = os.path.split(temp_obj.name)
        if 'file_tmp' in list(tmp_info.keys()):
            file_tmp = tmp_info['file_tmp']
        clean_tmp = True
        if 'clean_tmp' in list(tmp_info.keys()):
            clean_tmp = tmp_info['clean_tmp']

        if active_tmp:
            os.makedirs(path_tmp, exist_ok=True)

        return active_tmp, os.path.join(path_tmp, file_tmp), clean_tmp

    # method to read GLDAS time-series
    def read(self, *args, **kwargs):

        # get tmp information
        active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp
        gpi = args[0]

        # info time-series
        logging.info(' ------> Read GLDAS time-series for GPI "' + str(gpi) + '" ... ')

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
                ts_dframe = super(GLDAS_Dataset, self).read(*args)

                # soil moisture [kg/m^2] converted into [m^3/m^3]
                d = 0.10  # thickness of soil layer in m
                ts_dframe['SoilMoi0_10cm_inst'] = ts_dframe['SoilMoi0_10cm_inst'] * 0.001 * 1 / d

                # Filter(s)
                ts_dframe = ts_dframe[ts_dframe['SWE_inst'] == self.mask_snow]

                if 'SoilTMP0_10cm_inst' in ts_dframe.columns:
                    ts_dframe = ts_dframe[ts_dframe['SoilTMP0_10cm_inst'] > self.mask_soil_temp]
                elif 'Tair_f_inst' in ts_dframe.columns:
                    ts_dframe = ts_dframe[ts_dframe['Tair_f_inst'] > self.mask_soil_temp]
                else:
                    logging.error(' ===> GLDAS temperature variable is not expected.')
                    raise NotImplemented('Case not implemented yet')

                # check data valid
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for GLDAS dataset')
                    logging.warning(' ===> GLDAS time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' -------> Data valid for GLDAS dataset (N: "' + str(valid_value) + '")')
                    logging.info(' -------> Time valid for GLDAS dataset from "' + start_index +
                                 '" to "' + end_index + '"')

                # save time-series if flag is active
                if active_tmp:
                    write_obj(file_gpi, ts_dframe)

                # info time-series
                logging.info(' ------> Read GLDAS time-series for GPI "' + str(gpi) + '" ... DONE ')

            else:

                # read time-series previously saved
                ts_dframe = read_obj(file_gpi)

                # check data valid
                if ts_dframe.size == 0:
                    logging.warning(' ===> No data valid for GLDAS dataset')
                    logging.warning(' ===> GLDAS time-series will be initialized by empty dataframe')
                    ts_dframe = pd.DataFrame()
                else:
                    valid_value = ts_dframe.shape[0]
                    start_index = ts_dframe.index[0].strftime('%Y-%m-%d %H:%M:%S')
                    end_index = ts_dframe.index[-1].strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(' -------> Data valid for GLDAS dataset (N: "' + str(valid_value) + '")')
                    logging.info(' -------> Time valid for GLDAS dataset from "' + start_index +
                                 '" to "' + end_index + '"')

                # info time-series
                logging.info(' ------> Read GLDAS time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

        except BaseException as base_exc:

            logging.warning(' ===> Error in reading GLDAS time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> GLDAS time-series will be initialized by empty dataframe')
            ts_dframe = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read GLDAS time-series for GPI "' + str(gpi) + '" ... FAILED')

        # ts must be defined by dataframe
        if ts_dframe is None:
            logging.warning(' ===> GLDAS time-series is defined by NoneType. It will be format using DataFrame')
            ts_dframe = pd.DataFrame()

        return ts_dframe
# ----------------------------------------------------------------------------------------------------------------------
