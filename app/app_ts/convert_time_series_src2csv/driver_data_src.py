"""
Class Features

Name:          driver_data_source
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231120'
Version:       '1.5.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_data_io_generic import combine_data_point_by_time, range_data_point
from lib_data_io_ascii import wrap_datasets_ascii
from lib_data_io_csv import wrap_registry_csv
from lib_data_io_mat import wrap_registry_mat, wrap_datasets_mat
from lib_data_io_pickle import read_obj, write_obj

from lib_utils_system import make_folder

from lib_info_args import logger_name, time_format_algorithm, time_format_datasets

# logging
log_stream = logging.getLogger(logger_name)

# debugging
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DriverData:

    # -------------------------------------------------------------------------------------
    # initialize class
    def __init__(self, time_reference, time_range, registry_dict, datasets_dict, ancillary_dict,
                 flags_dict=None, template_dict=None, tmp_dict=None):

        # set time reference
        self.time_reference = time_reference
        self.time_range = time_range
        # set registry, ancillary and datasets dictionary
        self.registry_dict = registry_dict
        self.ancillary_dict = ancillary_dict
        self.datasets_dict = datasets_dict
        # set flags, template and tmp dictionary
        self.flags_dict = flags_dict
        self.template_dict_time = template_dict['time']
        self.template_dict_datasets = template_dict['datasets']
        self.tmp_dict = tmp_dict

        # registry and datasets tag(s)
        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'
        self.fields_tag = 'fields'
        self.format_tag = 'format'

        self.time_start_tag = "time_start"
        self.time_end_tag = "time_end"
        self.time_frequency_tag = "time_frequency"
        self.time_rounding_tag = "time_rounding"

        # reset flags
        self.reset_source = flags_dict['reset_source']
        self.reset_destination = flags_dict['reset_destination']

        # registry object(s)
        folder_name_registry = self.registry_dict[self.folder_name_tag]
        file_name_registry = self.registry_dict[self.file_name_tag]
        self.format_registry = self.registry_dict[self.format_tag]
        self.fields_registry = self.registry_dict[self.fields_tag]
        self.file_path_registry = os.path.join(folder_name_registry, file_name_registry)

        # datasets object(s)
        folder_name_datasets = self.datasets_dict[self.folder_name_tag]
        file_name_datasets = self.datasets_dict[self.file_name_tag]
        self.format_datasets = self.datasets_dict[self.format_tag]
        self.fields_datasets = self.datasets_dict[self.fields_tag]

        self.time_start = None
        if self.time_start_tag in list(self.datasets_dict.keys()):
            self.time_start = self.datasets_dict[self.time_start_tag]
        self.time_end = None
        if self.time_end_tag in list(self.datasets_dict.keys()):
            self.time_end = self.datasets_dict[self.time_end_tag]
        self.time_rounding = None
        if self.time_rounding_tag in list(self.datasets_dict.keys()):
            self.time_rounding = self.datasets_dict[self.time_rounding_tag]
        self.time_frequency = None
        if self.time_frequency_tag in list(self.datasets_dict.keys()):
            self.time_frequency = self.datasets_dict[self.time_frequency_tag]

        # to use in realtime run
        if self.time_end is None:
            self.time_end = self.time_reference.strftime(time_format_algorithm)

        self.file_path_datasets = os.path.join(folder_name_datasets, file_name_datasets)

        # ancillary object(s)
        folder_name_ancillary = ancillary_dict[self.folder_name_tag]
        file_name_ancillary = ancillary_dict[self.file_name_tag]
        self.file_path_ancillary = os.path.join(folder_name_ancillary, file_name_ancillary)

        # tmp object(s)
        self.folder_name_tmp = tmp_dict[self.folder_name_tag]
        self.file_name_tmp = tmp_dict[self.file_name_tag]

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get datasets object
    def get_obj_datasets(self, file_path_generic, file_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' -----> Read file datasets ... ')

        # check file format
        if self.format_datasets == 'mat':

            # get datasets in mat format
            fields_dframe_obj, fields_time_start, fields_time_end = wrap_datasets_mat(
                file_path_generic,
                file_fields=file_fields,
                registry_fields=registry_fields, folder_tmp=self.folder_name_tmp)
            fields_time_file = self.time_reference

        elif self.format_datasets == 'ascii':

            # get datasets in ascii format
            fields_dframe_obj, fields_time_start, fields_time_end, fields_time_file = wrap_datasets_ascii(
                file_path_generic,
                time_reference=self.time_reference, time_start=self.time_start, time_end=self.time_end,
                time_range=self.time_range,
                file_fields=file_fields, registry_fields=registry_fields,
                template_time_tags=self.template_dict_time, template_datasets_tags=self.template_dict_datasets,
                time_rounding=self.time_rounding, time_frequency=self.time_frequency,
                time_format=time_format_datasets,
                file_sep=' ', file_decimal='.')

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + self.format_datasets + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        #  check dataframe object
        if fields_dframe_obj is not None:

            # method to range data point
            time_frequency_expected, time_start_expected, time_end_expected = range_data_point(
                fields_dframe_obj, time_run_reference=fields_time_file,
                time_start_reference=fields_time_start, time_end_reference=fields_time_end)

            # method to combine data point to the expected time range
            fields_dframe_obj = combine_data_point_by_time(
                fields_dframe_obj, registry_fields,
                time_start_expected=time_start_expected, time_end_expected=time_end_expected,
                time_frequency_expected=time_frequency_expected, time_reverse=True)

            time_file_expected = fields_time_file

            # info end method
            log_stream.info(' -----> Read file datasets ... DONE')

        else:

            # info end method
            fields_dframe_obj, time_start_expected, time_end_expected = None, None, None
            log_stream.info(' -----> Read file datasets ... SKIPPED. All dynamic datasets are not defined.')

        return fields_dframe_obj, time_start_expected, time_end_expected, time_file_expected
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get registry object
    def get_obj_registry(self, file_name, file_fields=None):

        # info start method
        log_stream.info(' -----> Read file registry "' + file_name + '" ... ')

        # check file format
        if self.format_registry == 'mat':
            # get registry in mat format
            fields_dframe = wrap_registry_mat(file_name, file_fields=file_fields, folder_tmp=self.folder_name_tmp)
        elif self.format_registry == 'csv':
            # get registry in csv format
            fields_dframe = wrap_registry_csv(file_name, file_fields=file_fields)
        else:
            # exit with error if file format is not supported
            log_stream.info(' -----> Read file registry "' + file_name + '" ... FAILED')
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Read file registry "' + file_name + '" ... DONE')

        return fields_dframe

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize source object(s) ... ')

        # get path(s)
        file_path_registry, file_path_datasets = self.file_path_registry, self.file_path_datasets
        file_path_ancillary = self.file_path_ancillary

        # get flag(s)
        reset_source = self.reset_source

        # reset ancillary file
        if reset_source:
            if os.path.exists(file_path_ancillary):
                os.remove(file_path_ancillary)

        # check ancillary file availability
        if not os.path.exists(file_path_ancillary):

            # get registry dframe
            dframe_registry = self.get_obj_registry(file_path_registry, self.fields_registry)
            # get datasets dframe
            dframe_datasets, dframe_time_start, dframe_time_end, dframe_time_file = self.get_obj_datasets(
                file_path_datasets, self.fields_datasets, dframe_registry)

            # check dframe datasets
            if dframe_datasets is not None:

                # organize dframe obj
                dframe_obj = {
                    'registry': dframe_registry, 'datasets': dframe_datasets,
                    'time_start': dframe_time_start, 'time_end': dframe_time_end,
                    'time_file': dframe_time_file,
                }

                # dump dframe obj
                folder_name_ancillary, file_name_ancillary = os.path.split(file_path_ancillary)
                make_folder(folder_name_ancillary)
                write_obj(file_path_ancillary, dframe_obj)

                # method end info
                log_stream.info(' ----> Organize source object(s) ... DONE')

            else:
                # info end method
                dframe_obj = None
                log_stream.info(' ----> Organize source object(s) ... SKIPPED. '
                                'All dynamic datasets are not defined.')

        else:
            # read source obj
            dframe_obj = read_obj(file_path_ancillary)

            # method end info
            log_stream.info(' ----> Organize source object(s) ... DONE. '
                            'Datasets are read from file "' + str(file_path_ancillary) + '"')

        return dframe_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
