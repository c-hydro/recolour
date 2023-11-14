
"""
Library Features:

Name:           lib_datasets_smap
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org),
                Martina Natali (martina01.natali@edu.unife.it)
Date:           '20230719'
Version:        '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import tempfile
import pandas as pd
import numpy as np

from lib_utils_generic import read_obj, write_obj

from lib_interface_smap import SMAPTs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap SMAP time series
class SMAP_Dataset(SMAPTs):

    # Method to init class
    def __init__(self, smap_data_folder, **kwargs):

        self.smap_data_folder = smap_data_folder

        # deprecated
        # if 'sm_nan' in kwargs:
        #     self.sm_nan = kwargs.pop('sm_nan')
        # else:
        #     self.sm_nan = np.nan # -9999.

        # tag = 'reference', 'k1', 'k2'
        if 'dset_tag' in kwargs:
            self.tag = kwargs.pop('dset_tag')

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        super(SMAP_Dataset, self).__init__(ts_path=smap_data_folder, **kwargs)

    # method to set tmp info
    @staticmethod
    def _set_tmp_info(tmp_info):

        active_tmp = False
        if 'active_tmp' in list(tmp_info.keys()):
            active_tmp = tmp_info['active_tmp']
        path_tmp = tempfile.mkdtemp()
        if 'path_tmp' in list(tmp_info.keys()):
            path_tmp = tmp_info['path_tmp']
        temp_obj = tempfile.NamedTemporaryFile(prefix='smap_', suffix='.dframe')
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
        logging.info(' ------> Read SMAP time-series for GPI "' + str(gpi) + '" ... ')

        # organize dataframe
        try:

            file_gpi = file_tmp.format(gpi_n=str(gpi))
            if clean_tmp:
                if os.path.exists(file_gpi):
                    os.remove(file_gpi)

            if not os.path.exists(file_gpi):

                ts = super(SMAP_Dataset, self).read(*args)

                if ts.size == 0:
                    logging.warning(' ===> No data valid for SMAP dataset')
                    logging.warning(' ===> SMAP time-series will be initialized by empty dataframe')
                    ts = pd.DataFrame()

                # save time-series if flag is active
                if active_tmp:
                    write_obj(file_gpi, ts)

                # info time-series
                logging.info(' ------> Read SMAP time-series for GPI "' + str(gpi) + '" ... DONE ')

            else:

                # read time-series previously saved
                ts = read_obj(file_gpi)
                # info time-series
                logging.info(' ------> Read SMAP time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

        except BaseException as base_exc:

            logging.warning(' ===> Error in reading SMAP time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> SMAP time-series will be initialized by empty dataframe')
            ts = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read SMAP time-series for GPI "' + str(gpi) + '" ... FAILED')

        # ts must be defined by dataframe
        if ts is None:
            logging.warning(' ===> SMAP time-series is defined by NoneType. It will be format using DataFrame')
            ts = pd.DataFrame()

        # data cleaning
        # SMAP si expressed in [m3m-3], which take values in [0, 1]
        # so all values outside this range need to be put to nan
        mask = (ts < 0) | (ts > 1)
        ts[mask] = np.nan

        valid_value = np.sum(~np.isnan(ts.values))
        logging.info(' ------> Data valid for SMAP dataset (N: "' + str(valid_value) + '")')

        # to speed up validation, if reference dataset has less than 10 valid values
        # it is substituted with an empty dataframe to make validation skip the whole gpi
        # which is exactly what would happen in the same case when computing metrics
        if self.tag == 'reference' and valid_value < 10:
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------
