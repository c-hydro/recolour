
"""
Library Features:

Name:          lib_datasets_hmc
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import tempfile
import pandas as pd
import numpy as np

from lib_utils_generic import read_obj, write_obj
from lib_interface_hmc import HMCTs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to wrap HMC time series
class HMC_Dataset(HMCTs):

    # Method to init class
    def __init__(self, hmc_data_folder, **kwargs):

        self.hmc_data_folder = hmc_data_folder

        # deprecated
        # if 'sm_nan' in kwargs:
        #     self.sm_nan = kwargs.pop('sm_nan')
        # else:
        #     self.sm_nan = -9999.0

        # tag = 'reference', 'k1', 'k2'
        if 'dset_tag' in kwargs:
            self.tag = kwargs.pop('dset_tag')

        # set tmp info
        tmp_info = {}
        if "tmp_info" in kwargs:
            tmp_info = kwargs.pop('tmp_info')
        self.active_tmp, self.file_tmp, self.clean_tmp = self._set_tmp_info(tmp_info)

        # set read bulk
        self.read_bulk = False
        if "bulk" in kwargs:
            self.read_bulk = kwargs.pop('bulk')

        # set stack part
        self.stack_obj = {}
        if "path_stack" in kwargs:
            self.path_stack = kwargs.pop('path_stack')
        if "file_stack" in kwargs:
            self.file_stack = kwargs.pop('file_stack')
        else:
            self.file_stack = '{:}.stack'
        self.cell_templ = '%04d'

        super(HMC_Dataset, self).__init__(
            ts_path=hmc_data_folder, read_bulk=self.read_bulk,
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
        temp_obj = tempfile.NamedTemporaryFile(prefix='hmc_', suffix='.dframe')
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

    # method to read time-series
    def read(self, *args, **kwargs):

        # get tmp information
        active_tmp, file_tmp, clean_tmp = self.active_tmp, self.file_tmp, self.clean_tmp
        # get stack information
        file_stack_tmp, path_stack = self.file_stack, self.path_stack
        gpi = args[0]
        cell_num = self.grid.gpi2cell(gpi)
        cell_string = self.cell_templ % cell_num

        # info time-series
        logging.info(' ------> Read HMC time-series for GPI "' + str(gpi) + '" ... ')
        logging.info(' ------> Bulk Activated: ' + str(self.read_bulk))

        # organize dataframe
        try:

            file_gpi = file_tmp.format(gpi_n=str(gpi))
            if clean_tmp and os.path.exists(file_gpi):
                os.remove(file_gpi)

            # stack information
            stack_value = -1
            file_stack = os.path.join(path_stack, file_stack_tmp.format(cell_string))
            if os.path.exists(file_stack):
                stack_dframe = read_obj(file_stack)
                stack_row = stack_dframe.loc[stack_dframe.index == gpi]
                if not stack_row.empty:
                    stack_value = stack_row.values[0]

            # debug65
            stack_value = stack_dframe.loc[stack_dframe.values == 1].values[0][0]
            gpi = stack_dframe.loc[stack_dframe.values == 1].index[0]

            if stack_value == 1:
                if not os.path.exists(file_gpi):

                    ts = super(HMC_Dataset, self).read(*args)

                    if ts.size == 0:
                        logging.warning(' ===> No data valid for HMC dataset')
                        logging.warning(' ===> HMC time-series will be initialized by empty dataframe')
                        ts = pd.DataFrame()

                    # save time-series if flag is active
                    if active_tmp:
                        write_obj(file_gpi, ts)

                    # info time-series
                    logging.info(' ------> Read HMC time-series for GPI "' + str(gpi) + '" ... DONE ')

                else:

                    # read time-series previously saved
                    ts = read_obj(file_gpi)
                    # info time-series
                    logging.info(' ------> Read HMC time-series for GPI "' + str(gpi) + '" ... PREVIOUSLY SAVED')

            else:
                logging.info(' ------> Read HMC time-series for GPI "' + str(gpi) +
                             '" ... SKIPPED. GPI is classified in the stack with having no data values for the analysis period')
                ts = pd.DataFrame()

        except BaseException as base_exc:

            logging.warning(' ===> Error in reading HMC time-series "' + repr(base_exc) + '"')
            logging.warning(' ===> HMC time-series will be initialized by empty dataframe')
            ts = pd.DataFrame()

            # info time-series
            logging.info(' ------> Read HMC time-series for GPI "' + str(gpi) + '" ... FAILED')

        # ts must be defined by dataframe
        if ts is None:
            logging.warning(' ===> HMC time-series is defined by NoneType. It will be format using DataFrame')
            ts = pd.DataFrame()

        # data cleaning
        # HMC is expressed in [m3m-3], which take values in [0, 1]
        # so all values outside this range need to be put to nan
        mask = (ts < 0) | (ts > 1)
        ts[mask] = np.nan

        valid_value = np.sum(~np.isnan(ts.values))
        logging.info(' ------> Data valid for HMC dataset (N: "' + str(valid_value) + '")')

        # to speed up validation, if reference dataset has less than 10 valid values
        # it is substituted with an empty dataframe to make validation skip the whole gpi
        # which is exactly what would happen in the same case when computing metrics
        if valid_value < 10:
            ts = pd.DataFrame()

        return ts
# ----------------------------------------------------------------------------------------------------------------------
