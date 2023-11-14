
"""
Library Features:

Name:          lib_datasets_cci
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from pandas import read_pickle

from lib_interface_cci import CCITs

# import pytesmo.temporal_matching as temp_match

from lib_utils_generic import get_bit
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
    def read(self, *args, **kwargs):

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

                n_value = ts.loc['2018':'2019'].shape
                print(' ----> N Data valid for CCI dataset: ' + str(n_value))

        except Exception as exc:
            print(' ----> WARNING: RunTime error for CCI dataset -- ' + repr(exc))
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------
