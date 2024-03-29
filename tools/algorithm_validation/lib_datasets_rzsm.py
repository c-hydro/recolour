
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

# import pytesmo.temporal_matching as temp_match
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

        super(RZSM_Dataset, self).__init__(
            ts_path=self.dr_path,
            grid_path=None,
            **kwargs
        )

    # --------------------------------------------------------------------------

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
    # --------------------------------------------------------------------------

    # Method to read time-series
    def read(self, *args, **kwargs):

        # get tmp information
        active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp
        gpi = args[0]

        # info time-series
        logging.info(' ------> Read RZSM time-series for GPI "' + str(gpi) + '" ... ')

        # organize dataframe
        try:

            file_gpi = file_tmp.format(gpi_n=str(gpi))
            if clean_tmp and os.path.exists(file_gpi):
                os.remove(file_gpi)

            if not os.path.exists(file_gpi):

                ts = super(RZSM_Dataset, self).read(*args)

                if ts.size == 0:
                    logging.warning(' ===> No data valid for RZSM dataset')
                    logging.warning(' ===> RZSM time-series will be initialized by empty dataframe')
                    ts = pd.DataFrame()


                # save time-series if flag is active
                if active_tmp:
                    write_obj(file_gpi, ts)

                # info time-series
                logging.info(' ------> Read RZSM time-series for GPI "' + str(gpi) + '" ... DONE ')

            else:

                # read time-series previously saved
                ts = read_obj(file_gpi)
                # info time-series
                logging.info(' ------> Read RZSM time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

        except BaseException as base_exc:

            logging.warning(' ===> Error in reading RZSM time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> RZSM time-series will be initialized by empty dataframe')
            ts = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read RZSM time-series for GPI "' + str(gpi) + '" ... FAILED')

            # ts must be defined by dataframe
        if ts is None:
            logging.warning(' ===> RZSM time-series is defined by NoneType. It will be format using DataFrame')
            ts = pd.DataFrame()

        # data cleaning
        # RZSM is expressed as soil wetness index, which takes values in [0, 1]
        # so all values outside this range need to be put to nan
        mask = (ts < 0) | (ts > 1)
        ts[mask] = np.nan

        valid_value = np.sum(~np.isnan(ts.values))
        logging.info(' ------> Data valid for RZSM dataset (N: "' + str(valid_value) + '")')

        # to speed up validation, if reference dataset has less than 10 valid values
        # it is substituted with an empty dataframe to make validation skip the whole gpi
        # which is exactly what would happen in the same case when computing metrics
        if self.tag == 'reference' and valid_value < 10:
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------
