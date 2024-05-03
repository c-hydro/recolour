"""
Class Features

Name:           cpl_process_datasets
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org),
                Martina Natali (martina01.natali@edu.unife.it)
Date:           '20230719'
Version:        '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import tempfile

from copy import deepcopy

from lib_datasets_ascat import ASCAT_Dataset_DR, ASCAT_Dataset_NRT
from lib_datasets_gldas import GLDAS_Dataset
from lib_datasets_cci import CCI_Dataset
from lib_datasets_rzsm import RZSM_Dataset
from lib_datasets_hmc import HMC_Dataset
from lib_datasets_smap import SMAP_Dataset

from lib_utils_generic import make_folder
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# datasets types
datasets_type = {
    'ASCAT': ['data_record', 'nrt'],
    'GLDAS': ['data_record'], 'HMC': ['data_record'], 'CCI': ['data_record'],
    'RZSM': ['data_record', 'nrt'], 'ECMWF': ['data_record', 'nrt'], 'SMAP': ['data_record']
}
datasets_variable = {
    'ASCAT': ['sm'],
    'GLDAS': ['SoilMoi0_10cm_inst'], 'HMC': ['soil_moisture'], 'CCI': ['sm'],
    'RZSM': ['var40'], 'ECMWF': ['var40'], 'SMAP': ['soil_moisture']
}
datasets_name = ['ASCAT', 'RZSM', 'ECMWF', 'GLDAS', 'HMC', 'CCI', 'SMAP']
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class coupler process datasets
class CplDatasets:

    # method to initialize class
    def __init__(self, dset_collections_src, dset_collections_dst,
                 tag_dset_reference='ref', tag_dset_k1='k1', tag_dset_k2='k2'):

        self.dimensions = len(dset_collections_src)
        self.dset_collections_src = dset_collections_src
        self.dset_collections_dst = dset_collections_dst

        if self.dimensions == 2:
            self.dset_mode = ['ref', 'other']
            self.dset_tag = [tag_dset_reference, tag_dset_k1]
        elif self.dimensions == 3:
            self.dset_mode = ['ref', 'other', 'other']
            self.dset_tag = [tag_dset_reference, tag_dset_k1, tag_dset_k2]
        else:
            raise NotImplemented('Case not implemented yet.')

    # method to setup process datasets source
    def setup_datasets_src(self):

        # info source datasets start
        logging.info(' ----> Setup source datasets ... ')

        # iterate over source datasets
        dset_interfaces, dset_modes = {}, {}
        for dset_id, (dset_mode, dset_tag) in enumerate(zip(self.dset_mode, self.dset_tag)):

            # info datasets setup start
            logging.info(' -----> Datasets "' + dset_tag + '" ... ')

            # get dset obj
            dset_obj = self.dset_collections_src[dset_tag]

            # check name, type and variable
            dset_name = self.__check_key(dset_obj, dset_key='name')
            self.__check_field(dset_name, dset_expected_obj=datasets_name)
            dset_type = self.__check_key(dset_obj, dset_key='type')
            self.__check_field(dset_name, dset_type, dset_expected_obj=datasets_type)
            dset_variable = self.__check_key(dset_obj, dset_key='variable')
            self.__check_field(dset_name, dset_variable, dset_expected_obj=datasets_variable)
            dset_bulk = self.__check_key(dset_obj, dset_key='bulk', dset_mandatory=False, dset_default=None)

            # set tmp
            dset_tmp = self.set_dset_tmp(dset_obj)
            # set reader
            dset_reader = self.set_dset_reader(dset_obj, dset_name, dset_type, dset_tag,
                                               dset_tmp=dset_tmp, dset_bulk=dset_bulk)
            # set args and kwargs
            dset_args, dset_kwargs = self.set_dset_arguments(dset_obj)
            # set interface
            dset_interface = self.set_dset_interface(
                dset_reader, dset_variable, dset_mode,
                dset_args=dset_args, dset_kwargs=dset_kwargs)

            # collect interface(s) and mode(s)
            dset_interfaces[dset_tag] = dset_interface
            dset_modes[dset_tag] = dset_mode

            # info datasets setup end
            logging.info(' -----> Datasets "' + dset_tag + '" ... DONE')

        # info source datasets end
        logging.info(' ----> Setup source datasets ... DONE')

        return dset_interfaces, dset_modes

    # method to setup process datasets destination
    def setup_datasets_dst(self):

        # info destination datasets start
        logging.info(' ----> Setup destination datasets ... ')

        dset_obj = self.dset_collections_dst

        path_analysis = None
        if 'path_analysis' in list(dset_obj.keys()):
            path_analysis = dset_obj['path_analysis']

        if (path_analysis == "") or (path_analysis is None):
            path_analysis = tempfile.mkdtemp()
            logging.warning(' ===> Destination folder is null or defined by NoneType. Create temporary folder.'
                            ' Datasets will be saved in "' + path_analysis + '"')

        # make analysis folder
        make_folder(path_analysis)

        # info destination datasets end
        logging.info(' ----> Setup destination datasets ... DONE')

        return path_analysis

    # method to set datasets interface
    @staticmethod
    def set_dset_interface(dset_reader, dset_variable, dset_mode,
                           dset_args=None, dset_kwargs=None,
                           dset_grids_compatible=False,
                           dset_use_lut=True, dset_lut_max_dist=None):

        if isinstance(dset_variable, str):
            dset_variable = [dset_variable]
        if dset_args is None:
            dset_args = []
        if dset_kwargs is None:
            dset_kwargs = {}
        else:
            if 'lut_max_dist' in dset_kwargs:
                dset_lut_max_dist = dset_kwargs['lut_max_dist']
            else:
                dset_lut_max_dist = 35000

        if dset_mode == 'ref':

            dset_interface = {
                'class': dset_reader, 'columns': dset_variable,
                'type': dset_mode, 'args': dset_args, 'kwargs': dset_kwargs}

        elif dset_mode == 'other':

            dset_interface = {
                'class': dset_reader, 'columns': dset_variable,
                'type': dset_mode, 'args': dset_args, 'kwargs': dset_kwargs,
                'grids_compatible': dset_grids_compatible,
                'use_lut': dset_use_lut, 'lut_max_dist': dset_lut_max_dist}

        else:

            logging.error(' ===> Dataset mode "' + dset_mode + '" is not expected by the datasets interface')
            raise NotImplemented('Case not implemented yet')

        return dset_interface

    # method to set datasets tmp info
    @staticmethod
    def set_dset_tmp(dset_obj, tag_tmp='tmp'):

        if tag_tmp in list(dset_obj.keys()):
            tmp_kwargs = dset_obj[tag_tmp]
        else:
            logging.warning(' ===> The tmp "kwargs" is not defined in the dataset obj. Use empty dictionary.')
            tmp_kwargs = {}
        return tmp_kwargs

    # method to set datasets arguments
    @staticmethod
    def set_dset_arguments(dset_obj, tag_args='args', tag_kwargs='kwargs'):

        if tag_args in list(dset_obj.keys()):
            dset_args = dset_obj[tag_args]
        else:
            logging.warning(' ===> The obj "args" is not defined in the dataset obj. Use empty list.')
            dset_args = []
        if tag_kwargs in list(dset_obj.keys()):
            dset_kwargs = dset_obj[tag_kwargs]
        else:
            logging.warning(' ===> The obj "kwargs" is not defined in the dataset obj. Use empty dictionary.')
            dset_kwargs = {}

        return dset_args, dset_kwargs

    # method to set dataset reader
    @staticmethod
    def set_dset_reader(dset_obj, dset_name, dset_type, dset_tag, dset_tmp=None, dset_bulk=None):

        if dset_tmp is None:
            dset_tmp = {}

        path_ts, path_grid, path_static = dset_obj['path_ts'], dset_obj['path_static'], dset_obj['path_static']
        path_stack = None
        if 'path_stack' in list(dset_obj.keys()):
            path_stack = dset_obj['path_stack']

        if dset_name == 'ASCAT':
            if dset_type == 'data_record':
                dset_reader = ASCAT_Dataset_DR(
                    dr_path=path_ts, grid_path=path_grid,
                    static_layer_path=path_static,
                    tmp_info=dset_tmp, bulk=dset_bulk)
            elif dset_type == 'nrt':
                dset_reader = ASCAT_Dataset_NRT(
                    dr_path=path_ts, grid_path=path_grid,
                    static_layer_path=path_static,
                    tmp_info=dset_tmp, bulk=dset_bulk)
            else:
                logging.error(' ===> Dataset type "' + dset_type +
                              '" for the Dataset name "' + dset_name + '" is not expected by the datasets driver')
                raise NotImplemented('Case not implemented yet')

        elif (dset_name == 'RZSM') or (dset_name == 'ECMWF'):
            if dset_type == 'data_record' or dset_type == 'nrt':
                dset_reader = RZSM_Dataset(
                    dr_path=path_ts, grid_path=path_grid,
                    tmp_info=dset_tmp, dset_tag=dset_tag)
            else:
                logging.error(' ===> Dataset type "' + dset_type +
                              '" for the Dataset name "' + dset_name + '" is not expected by the datasets driver')
                raise NotImplemented('Case not implemented yet')

        elif dset_name == 'GLDAS':
            if dset_type == 'data_record':
                dset_reader = GLDAS_Dataset(path_ts, tmp_info=dset_tmp)
            else:
                logging.error(' ===> Dataset type "' + dset_type +
                              '" for the Dataset name "' + dset_name + '" is not expected by the datasets driver')
                raise NotImplemented('Case not implemented yet')

        elif dset_name == 'CCI':
            if dset_type == 'data_record':
                dset_reader = CCI_Dataset(path_ts, tmp_info=dset_tmp)
            else:
                logging.error(' ===> Dataset type "' + dset_type +
                              '" for the Dataset name "' + dset_name + '" is not expected by the datasets driver')
                raise NotImplemented('Case not implemented yet')

        elif dset_name == 'HMC':
            if dset_type == 'data_record':
                dset_reader = HMC_Dataset(path_ts, tmp_info=dset_tmp, dset_tag=dset_tag, path_stack=path_stack,
                                          bulk=dset_bulk)
            else:
                logging.error(' ===> Dataset type "' + dset_type +
                              '" for the Dataset name "' + dset_name + '" is not expected by the datasets driver')
                raise NotImplemented('Case not implemented yet')

        elif dset_name == 'SMAP':
            if dset_type == 'data_record':
                dset_reader = SMAP_Dataset(path_ts, tmp_info=dset_tmp, dset_tag=dset_tag)
            else:
                logging.error(' ===> Dataset type "' + dset_type +
                              '" for the Dataset name "' + dset_name + '" is not expected by the datasets driver')
                raise NotImplemented('Case not implemented yet')

        else:
            logging.error(' ===> Dataset name "' + dset_name + '" is not expected by the datasets driver')
            raise NotImplemented('Case not implemented yet')

        return dset_reader

    # method to check key availability
    @staticmethod
    def __check_key(dset_obj, dset_key='name', dset_mandatory=True, dset_default=None):
        if dset_key in list(dset_obj.keys()):
            dset_value = dset_obj[dset_key]
        else:
            if dset_mandatory:
                logging.error(' ===> Dataset key "' + dset_key + '" is not available in the datasets obj')
                raise RuntimeError('The key is needed by the algorithm. Please check your settings file')
            else:
                dset_value = deepcopy(dset_default)
        return dset_value

    # method to check field availability
    @staticmethod
    def __check_field(dset_field_1, dset_field_2=None, dset_expected_obj=None):
        if dset_expected_obj is not None:
            if isinstance(dset_expected_obj, list):
                if dset_field_1 not in dset_expected_obj:
                    logging.error(' ===> Dataset field "' + dset_field_1 +
                                  '" is not available in the datasets expected obj')
                    raise RuntimeError('The name is needed by the algorithm. Please check your datasets driver')
                else:
                    return True
            if isinstance(dset_expected_obj, dict):
                if dset_field_1 in list(dset_expected_obj.keys()):
                    dset_values = dset_expected_obj[dset_field_1]
                    if dset_field_2 is not None:
                        if dset_field_2 not in dset_values:
                            logging.error(' ===> Dataset value "' + dset_field_2 + '" is not in the expected values')
                            raise RuntimeError('The value is not expected by the algorithm')
                        else:
                            return True
                else:
                    logging.error(' ===> Dataset key "' + dset_field_1 + '" is not in the expected keys')
                    raise RuntimeError('The key is not expected by the algorithm')
        else:
            logging.error(' ===> Dataset expected obj is defined by NoneType')
            raise RuntimeError('The obj must be defined by list or dictionary')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
