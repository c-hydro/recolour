"""
Class Features

Name:          driver_data_destination
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_io_csv import write_file_csv
from lib_utils_obj import map_vars_dframe
from lib_utils_system import make_folder
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DriverData:

    # -------------------------------------------------------------------------------------
    # initialize class
    def __init__(self, registry_dict, datasets_dict, ancillary_dict,
                 flags_dict=None, template_dict=None, params_dict=None, tmp_dict=None):

        self.registry_dict = registry_dict
        self.ancillary_dict = ancillary_dict
        self.datasets_dict = datasets_dict

        self.flags_dict = flags_dict
        self.template_dict = template_dict
        self.params_dict = params_dict
        self.tmp_dict = tmp_dict

        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'
        self.fields_tag = 'fields'
        self.no_data_tag = 'no_data'

        self.reset_source = flags_dict['reset_source']
        self.reset_destination = flags_dict['reset_destination']

        self.format_source = self.params_dict['format_source']
        self.format_destination = self.params_dict['format_destination']

        # ancillary object(s)
        folder_name_ancillary = ancillary_dict[self.folder_name_tag]
        file_name_ancillary = ancillary_dict[self.file_name_tag]
        self.file_path_ancillary = os.path.join(folder_name_ancillary, file_name_ancillary)

        # registry object(s)
        folder_name_registry = self.registry_dict[self.folder_name_tag]
        file_name_registry = self.registry_dict[self.file_name_tag]
        self.file_path_registry = os.path.join(folder_name_registry, file_name_registry)
        self.fields_registry = self.registry_dict[self.fields_tag]
        self.no_data_registry = self.registry_dict[self.no_data_tag]

        # datasets object(s)
        folder_name_datasets = self.datasets_dict[self.folder_name_tag]
        file_name_datasets = self.datasets_dict[self.file_name_tag]
        self.file_path_datasets = os.path.join(folder_name_datasets, file_name_datasets)
        self.fields_datasets = self.datasets_dict[self.fields_tag]
        self.no_data_datasets = self.datasets_dict[self.no_data_tag]

        # tmp object(s)
        self.folder_name_tmp = tmp_dict[self.folder_name_tag]
        self.file_name_tmp = tmp_dict[self.file_name_tag]

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump datasets object
    def dump_obj_datasets(self, file_path, file_dframe_raw, file_fields=None):

        # info start method
        log_stream.info(' -----> Dump file datasets "' + file_path + '" ... ')

        if file_path.endswith('csv'):

            # map datasets fields
            # file_dframe_map = map_vars_dframe(file_dframe_raw, file_fields)
            # fill nan(s) with no_data
            file_dframe_map = file_dframe_raw.fillna(self.no_data_datasets)

            # create datasets folder
            folder_name, file_name = os.path.split(file_path)
            make_folder(folder_name)

            # dump datasets dframe
            write_file_csv(
                file_path, file_dframe_map,
                dframe_sep=',', dframe_decimal='.', dframe_float_format='%.1f',
                dframe_index=True, dframe_header=True)

        else:
            log_stream.info(' -----> Dump file datasets "' + file_path + '" ... FAILED')
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Dump file datasets "' + file_path + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump registry object
    def dump_obj_registry(self, file_path, file_dframe_raw, file_fields=None):

        # info start method
        log_stream.info(' -----> Dump file registry "' + file_path + '" ... ')

        if file_path.endswith('csv'):

            # map registry fields
            file_dframe_map = map_vars_dframe(file_dframe_raw, file_fields)

            # fill nan(s) with no_data
            file_dframe_map = file_dframe_map.fillna(self.no_data_registry)

            # create registry folder
            folder_name, file_name = os.path.split(file_path)
            make_folder(folder_name)

            # dump registry dframe
            write_file_csv(
                file_path, file_dframe_map,
                dframe_sep=',', dframe_decimal='.', dframe_float_format='%.4f',
                dframe_index=False, dframe_header=True)

        else:
            log_stream.info(' -----> Dump file registry "' + file_path + '" ... FAILED')
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Dump file registry "' + file_path + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self, data_source_obj):

        # method start info
        log_stream.info(' ----> Organize destination object(s) ... ')

        # get path(s)
        file_path_registry, file_path_datasets = self.file_path_registry, self.file_path_datasets

        # get flag(s)
        reset_destination = self.reset_destination

        # reset destination file(s)
        if reset_destination:
            if os.path.exists(file_path_registry):
                os.remove(file_path_registry)
            if os.path.exists(file_path_datasets):
                os.remove(file_path_datasets)

        # check destination format
        if self.format_destination == 'csv_dr':

            # get source obj
            dframe_registry, dframe_datasets = data_source_obj['registry'], data_source_obj['datasets']

            # check availability of destination file(s)
            if (not os.path.exists(file_path_registry)) or (not os.path.exists(file_path_datasets)):

                # dump registry dframe
                self.dump_obj_registry(file_path_registry, dframe_registry, self.fields_registry)
                # dump datasets dframe
                self.dump_obj_datasets(file_path_datasets, dframe_datasets, self.fields_datasets)

                # method end info
                log_stream.info(' ----> Organize destination object(s) ... DONE')

            else:
                # method end info
                log_stream.info(' ----> Organize destination object(s) ... SKIPPED. '
                                'Destination file were previously dumped.')
        else:
            # method end info
            log_stream.info(' ----> Organize destination object(s) ... FAILED.')
            log_stream.error(' ===> Destination format "' + self.format_source + '" is not supported')
            raise NotImplemented('Case not implemented yet')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
