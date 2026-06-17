"""
Class Features

Name:          driver_data_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20241018'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

# from lib_data_io_ascii import read_parameters_ascii
from lib_data_io_csv import read_registry_csv, read_parameters_csv
from lib_data_io_pickle import read_obj, write_obj

from lib_utils_io import fill_string_with_time, fill_string_with_info
from lib_utils_generic import make_folder
from lib_utils_obj import join_dframe

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DriverData:

    # -------------------------------------------------------------------------------------
    # initialize class
    def __init__(self, time_reference, alg_datasets, alg_info, alg_template, alg_flags):

        # set time reference
        self.time_reference = time_reference

        # set algorithm information
        self.alg_flags = alg_flags
        self.alg_info = alg_info
        self.alg_datasets_src_registry = alg_datasets['source']['registry']
        self.alg_datasets_src_parameters = alg_datasets['source']['parameters']
        self.alg_datasets_dst = alg_datasets['destination']
        self.alg_template_time = alg_template['time']
        self.alg_template_datasets = alg_template['datasets']

        # reset flags
        self.reset_data_static = self.alg_flags['reset_data_static']
        self.reset_data_dynamic = self.alg_flags['reset_data_dynamic']

        # registry and datasets tag(s)
        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.filters_tag, self.delimiter_tag = 'filters', 'delimiter'
        self.fields_tag, self.format_tag = 'fields', 'format'

        # source registry object(s)
        self.folder_name_src_registry = self.alg_datasets_src_registry['folder_name']
        self.file_name_src_registry = self.alg_datasets_src_registry['file_name']
        self.format_registry = self.alg_datasets_src_registry[self.format_tag]
        self.fields_registry = self.alg_datasets_src_registry[self.fields_tag]
        self.filters_registry = self.alg_datasets_src_registry[self.filters_tag]
        if self.delimiter_tag in self.alg_datasets_src_registry.keys():
            self.delimiter_registry = self.alg_datasets_src_registry[self.delimiter_tag]
        else:
            self.delimiter_registry = ';'
        self.file_path_src_registry = os.path.join(self.folder_name_src_registry, self.file_name_src_registry)

        # source parameters object(s)
        self.folder_name_src_params = self.alg_datasets_src_parameters['folder_name']
        self.file_name_src_params = self.alg_datasets_src_parameters['file_name']
        self.format_params = self.alg_datasets_src_parameters[self.format_tag]
        self.fields_params = self.alg_datasets_src_parameters[self.fields_tag]
        self.filters_params = self.alg_datasets_src_parameters[self.filters_tag]
        if self.delimiter_tag in self.alg_datasets_src_registry.keys():
            self.delimiter_params = self.alg_datasets_src_parameters[self.delimiter_tag]
        else:
            self.delimiter_params = ';'
        self.file_path_src_params = os.path.join(self.folder_name_src_params, self.file_name_src_params)

        # destination object(s)
        self.folder_name_dst = self.alg_datasets_dst['folder_name']
        self.file_name_dst = self.alg_datasets_dst['file_name']
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get registry object
    def get_obj_registry(self, file_name, file_fields=None, file_filters=None):

        # info start method
        log_stream.info(' -----> Read file registry "' + file_name + '" ... ')

        # check file existence
        if not os.path.exists(file_name):
            log_stream.error(' ===> File "' + file_name + '" does not exist')
            raise IOError('File parameters must be available')

        # check file format
        if self.format_registry == 'csv':
            # get registry in csv format
            fields_obj = read_registry_csv(
                file_name, file_fields=file_fields, file_filters=file_filters, file_sep=self.delimiter_registry)
        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Read file registry "' + file_name + '" ... DONE')

        return fields_obj

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get parameters object
    def get_obj_parameters(self, file_name, file_fields=None, file_filters=None):

        # info start method
        log_stream.info(' -----> Read file parameters "' + file_name + '" ... ')

        # check file existence
        if not os.path.exists(file_name):
            log_stream.error(' ===> File "' + file_name + '" does not exist')
            raise IOError('File parameters must be available')

        # check file format
        if self.format_params == 'ascii':
            # get params in ascii format
            # fields_dframe = read_parameters_ascii(file_name)
            log_stream.error(' ===> File ascii format is not longer supported')
            raise NotImplemented('Case not longer supported by the library')
        elif self.format_params == 'csv':
            # get params in csv format
            fields_dframe = read_parameters_csv(
                file_name, file_fields=file_fields, file_filters=file_filters, file_sep=self.delimiter_params)
        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Read file parameters "' + file_name + '" ... DONE')

        return fields_dframe

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump collections object
    def dump_obj_collections(self, file_name, file_obj, file_format='pickle'):

        # info start method
        log_stream.info(' -----> Dump collections "' + file_name + '" ... ')

        # check file format
        if file_format == 'pickle':
            # dump dframe obj
            folder_name, _ = os.path.split(file_name)
            make_folder(folder_name)
            write_obj(file_name, file_obj)
        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Dump collections "' + file_name + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to read collections object
    def read_obj_collections(self, file_name, file_format='pickle'):
        # info start method
        log_stream.info(' -----> Read collections "' + file_name + '" ... ')

        # check file format
        if file_format == 'pickle':
            # dump dframe obj
            file_obj = read_obj(file_name)
        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Read collections "' + file_name + '" ... DONE')

        return file_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to define file string
    def __define_file_string(self, file_string_tmpl):
        file_string_def = fill_string_with_time(file_string_tmpl, self.time_reference, self.alg_template_time)
        file_string_def = fill_string_with_info(file_string_def, self.alg_info, self.alg_template_datasets)
        return file_string_def
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize data static object(s) ... ')

        # get time reference
        time_step_reference = self.time_reference

        # get path(s
        file_path_src_registry_tmpl = self.file_path_src_registry
        file_path_src_params_tmpl = self.file_path_src_params
        file_path_dst_tmpl = self.file_path_dst

        # get flag(s)
        reset_data_static = self.reset_data_static

        # method to fill the filename(s)
        file_path_src_registry_def = self.__define_file_string(file_path_src_registry_tmpl)
        file_path_src_params_def = self.__define_file_string(file_path_src_params_tmpl)
        file_path_dst_def = self.__define_file_string(file_path_dst_tmpl)

        # reset destination file if required
        if reset_data_static:
            if os.path.exists(file_path_dst_def):
                os.remove(file_path_dst_def)

        # check ancillary file availability
        if not os.path.exists(file_path_dst_def):

            # get registry obj
            obj_registry = self.get_obj_registry(
                file_path_src_registry_def, self.fields_registry, self.filters_registry)
            # get parameters obj
            obj_parameters = self.get_obj_parameters(
                file_path_src_params_def, self.fields_params, self.filters_params)

            # join registry and parameters obj
            obj_registry = join_dframe(obj_registry, obj_parameters, column_ref='tag', column_suffix='_tmp')

            # organize destination obj
            obj_collections = {'registry': obj_registry}

            # dump destination obj
            self.dump_obj_collections(file_path_dst_def, obj_collections)

            # method end info
            log_stream.info(' ----> Organize data static object(s) ... DONE')

        else:
            # read source obj
            obj_collections = self.read_obj_collections(file_path_dst_def)
            # method end info
            log_stream.info(' ----> Organize data static object(s) ... DONE. Datasets previously saved')

        return obj_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
