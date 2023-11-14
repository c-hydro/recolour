"""
Class Features

Name:          drv_map_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_utils_generic import make_folder, reset_folder
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_nc import read_file_nc, organize_file_nc, write_file_nc

from lib_fx_utils import organize_data, convert_data, mask_data

from drv_map_fx import DrvFx
# debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver map
class DrvMap:

    # method to initialize class
    def __init__(self, alg_obj_static, alg_settings,
                 tag_section_flags='flags',
                 tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_time='time', tag_section_log='log'):

        self.alg_obj_static = alg_obj_static

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_parameters = alg_settings[tag_section_params]

        # self.alg_time = alg_settings[tag_section_time]

        self.alg_datasets_src = alg_settings[tag_section_datasets]['dynamic']['source']
        self.alg_datasets_anc_data = alg_settings[tag_section_datasets]['dynamic']['ancillary']['data']
        self.alg_datasets_anc_weights = alg_settings[tag_section_datasets]['dynamic']['ancillary']['weights']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_file_vars = 'variables'

        self.reset_datasets_anc_data = self.alg_flags['reset_ancillary_datasets_data']
        self.reset_datasets_anc_weights = self.alg_flags['reset_ancillary_datasets_weights']
        self.reset_datasets_dst = self.alg_flags['reset_destination_datasets']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_src = self.alg_datasets_src[self.tag_folder_name]
        self.file_name_src = self.alg_datasets_src[self.tag_file_name]
        self.variables_src = self.alg_datasets_src[self.tag_file_vars]
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)

        self.folder_name_anc_data = self.alg_datasets_anc_data[self.tag_folder_name]
        self.file_name_anc_data = self.alg_datasets_anc_data[self.tag_file_name]
        self.file_path_anc_data = os.path.join(self.folder_name_anc_data, self.file_name_anc_data)

        self.folder_name_anc_weights = self.alg_datasets_anc_weights[self.tag_folder_name]
        self.file_name_anc_weights = self.alg_datasets_anc_weights[self.tag_file_name]
        self.file_path_anc_weights = os.path.join(self.folder_name_anc_weights, self.file_name_anc_weights)

        self.folder_name_dst = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_datasets_dst[self.tag_file_name]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.grid_reference_domain = self.alg_obj_static['grid_reference_domain']

        self.weights_variables_in = self.alg_parameters['weights_variables_in']
        self.weights_variables_out = self.alg_parameters['weights_variables_out']
        self.weights_fx = self.alg_parameters['weights_fx']

    # method to check datasets (if clean or not)
    @staticmethod
    def _clean_ancillary_folders(dset_obj,
                                 dset_key_root=None, dset_key_sub=None, dset_clean=True):

        if (dset_key_root is not None) and (dset_key_sub is not None):
            for dset_key, dset_fields in dset_obj.items():
                dset_path = dset_fields[dset_key_root][dset_key_sub]
                reset_folder(dset_path, folder_reset=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is not None):
            dset_path = dset_obj[dset_key_sub]
            reset_folder(dset_path, folder_reset=dset_clean)
        elif (dset_key_root is not None) and (dset_key_sub is None):
            dset_path = dset_obj[dset_key_root]
            reset_folder(dset_path, folder_reset=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is None):
            pass
        else:
            pass

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize dynamic datasets ... ')

        # get file variables
        variables_src = self.variables_src
        # get file path
        file_path_src = self.file_path_src
        file_path_anc_data, file_path_anc_weights = self.file_path_anc_data, self.file_path_anc_weights
        file_path_dst = self.file_path_dst

        # clean ancillary and destination datasets (if ancillary flag(s) is activated)
        if self.reset_datasets_anc_data or self.reset_datasets_anc_weights:
            if os.path.exists(file_path_anc_data):
                os.remove(file_path_anc_data)
            if os.path.exists(file_path_anc_weights):
                os.remove(file_path_anc_weights)
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)
        # clean ancillary and destination datasets if are not available together
        if (not os.path.exists(file_path_anc_data)) or (not os.path.exists(file_path_dst)):
            if os.path.exists(file_path_anc_data):
                os.remove(file_path_anc_data)
            if os.path.exists(file_path_anc_weights):
                os.remove(file_path_anc_weights)
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)

        # check file point ancillary availability
        if not os.path.exists(file_path_anc_data):

            # get dataset source
            dset_src = read_file_nc(file_path_src, file_variables_selected=variables_src)

            # organize variable(s) source
            variables_src = organize_data(dset_src)

            # save variable(s) obj to ancillary file
            folder_name_anc_data, file_name_anc_data = os.path.split(file_path_anc_data)
            make_folder(folder_name_anc_data)
            write_file_obj(file_path_anc_data, variables_src)

            # info end method
            logging.info(' ---> Organize dynamic datasets ... DONE')

        else:
            # info end method
            logging.info(' ---> Organize dynamic datasets ... SKIPPED. Data points previously saved')

    # method to analyze data
    def analyze_data(self):

        # info start method
        logging.info(' ---> Analyze dynamic data ... ')

        # get file path
        file_path_anc_data, file_path_anc_weights = self.file_path_anc_data, self.file_path_anc_weights

        # get weights information
        fx_weights_variables_in = self.weights_variables_in
        fx_weights_variables_out = self.weights_variables_out
        fx_weights_name = self.weights_fx

        # get reference domain
        grid_reference_domain = self.grid_reference_domain

        # check file ancillary points availability
        if os.path.exists(file_path_anc_data):

            # check file ancillary grid availability
            if not os.path.exists(file_path_anc_weights):

                # get variable(s) obj
                fx_source_data = read_file_obj(file_path_anc_data)

                # initialize fx driver
                driver_fx = DrvFx(fx_source_data, fx_name=fx_weights_name,
                                  fx_vars_in=fx_weights_variables_in,
                                  fx_vars_out=fx_weights_variables_out)
                # organize fx data source
                fx_handle, fx_kwargs_data, fx_kwargs_geo = driver_fx.organize_fx_src()
                # execute fx method
                fx_result_data = driver_fx.execute_fx(fx_handle, fx_kwargs_data)
                # organize fx data destination
                fx_destination_data = driver_fx.organize_fx_dst(fx_result_data, fx_kwargs_geo)

                # method to convert data
                fx_darray_data = convert_data(fx_destination_data)

                # save weights collection in pickle format
                folder_name_anc_weights, file_name_anc_weights = os.path.split(file_path_anc_weights)
                make_folder(folder_name_anc_weights)
                write_file_obj(file_path_anc_weights, fx_darray_data)

                # info end method
                logging.info(' ---> Analyze dynamic data ... DONE')

            else:
                # info end method
                logging.info(' ---> Analyze dynamic data ... SKIPPED. Data grid previously saved')

        else:
            # info end method
            logging.info(' ---> Analyze dynamic data... FAILED. Data points not available')

    # method to save data
    def dump_data(self):

        # info start method
        logging.info(' ---> Dump dynamic datasets ... ')

        # get file path
        file_path_anc_weights, file_path_dst = self.file_path_anc_weights, self.file_path_dst

        # check file ancillary grid availability
        if os.path.exists(file_path_anc_weights):
            if not os.path.exists(file_path_dst):

                # method to get data obj
                variable_collection = read_file_obj(file_path_anc_weights)

                # method to organize dataset
                variable_dset = organize_file_nc(variable_collection)

                # method to write dataset
                folder_name_dst, file_name_dst = os.path.split(file_path_dst)
                make_folder(folder_name_dst)
                write_file_nc(file_path_dst, variable_dset)

                # info end method
                logging.info(' ---> Dump dynamic datasets ... DONE')

            else:
                # info end method
                logging.info(' ---> Dump dynamic datasets ... SKIPPED. Data previously saved')
        else:
            # info end method
            logging.info(' ---> Dump dynamic datasets ... FAILED. Data grid not available')


# -------------------------------------------------------------------------------------




