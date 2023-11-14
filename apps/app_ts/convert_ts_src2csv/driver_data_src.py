"""
Class Features

Name:          driver_data_source
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import os

import pandas as pd

from lib_data_io_mat import read_file_mat, convert_var_mat
from hmc.lib_data_io_pickle import read_obj, write_obj

from lib_utils_obj import map_vars_dict, sanitize_string
from lib_utils_system import make_folder

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)

# debugging
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

        self.reset_source = flags_dict['reset_source']
        self.reset_destination = flags_dict['reset_destination']

        self.format_source = self.params_dict['format_source']
        self.format_destination = self.params_dict['format_destination']

        # registry object(s)
        folder_name_registry = self.registry_dict[self.folder_name_tag]
        file_name_registry = self.registry_dict[self.file_name_tag]
        self.fields_registry = self.registry_dict[self.fields_tag]
        self.file_path_registry = os.path.join(folder_name_registry, file_name_registry)

        # ancillary object(s)
        folder_name_ancillary = ancillary_dict[self.folder_name_tag]
        file_name_ancillary = ancillary_dict[self.file_name_tag]
        self.file_path_ancillary = os.path.join(folder_name_ancillary, file_name_ancillary)

        # datasets object(s)
        folder_name_datasets = self.datasets_dict[self.folder_name_tag]
        file_name_datasets = self.datasets_dict[self.file_name_tag]
        self.fields_datasets = self.datasets_dict[self.fields_tag]
        self.file_path_datasets = os.path.join(folder_name_datasets, file_name_datasets)

        # tmp object(s)
        self.folder_name_tmp = tmp_dict[self.folder_name_tag]
        self.file_name_tmp = tmp_dict[self.file_name_tag]

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get datasets object
    def get_obj_datasets(self, file_name, file_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' -----> Read file datasets "' + file_name + '" ... ')

        if file_name.endswith('mat'):

            # get fields data raw
            fields_data_raw = read_file_mat(
                file_name, file_mandatory=True,
                file_fields_excluded=['__header__', '__version__', '__globals__'])
            # map fields
            fields_data_map = map_vars_dict(fields_data_raw, file_fields)

            # parser time(s)
            var_data_time = convert_var_mat(
                fields_data_map['time'], var_type='time',
                reset_tmp=False, folder_tmp=self.folder_name_tmp, file_tmp='time.workspace')
            fields_data_map.pop('time')

            # parser name(s)
            var_data_name = convert_var_mat(
                fields_data_map['name'], var_type='other',
                reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='name.workspace')
            fields_data_map.pop('name')

            # parser description(s)
            var_data_description = convert_var_mat(
                fields_data_map['description'], var_type='other',
                reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='description.workspace')
            fields_data_map.pop('description')

            # parser tag
            if 'tag' not in list(fields_data_map.keys()):

                if 'tag' in list(registry_fields.columns):
                    var_data_tag = list(registry_fields['tag'].values)
                else:
                    var_data_tag = []
                    for string_name in var_data_name:
                        string_tag = sanitize_string(string_name)
                        var_data_tag.append(string_tag)
            else:
                var_data_tag = convert_var_mat(
                    fields_data_map['tag'], var_type='other',
                    reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='tag.workspace')
                fields_data_map.pop('tag')

            # create fields data update
            fields_data_upd = {
                'sm': fields_data_map['sm'].T,
                'time': var_data_time, 'tag': var_data_tag,
                'name': var_data_name, 'description': var_data_description
            }

            # create fields data merge
            fields_data_merge = {**fields_data_map, **fields_data_upd}

            # adapt data dict to dataframe
            data_values = fields_data_merge['sm']
            data_time = fields_data_merge['time']
            data_name, data_description = fields_data_merge['name'], fields_data_merge['description']
            data_tag = fields_data_merge['tag']
            # no data values
            data_values[data_values < -99] = np.nan

            data_obj, columns_obj = {}, []
            for var_id, (var_tag, var_description) in enumerate(zip(data_tag, data_description)):
                data_obj[var_tag] = data_values[:, var_id]
                columns_obj.append([(var_tag, var_description)])

            # create dframe
            fields_dframe = pd.DataFrame(data=data_obj, index=data_time)
            fields_dframe.index.name = 'time'

            # check dframe consistency
            start_time, end_time = data_time[0], data_time[-1]
            time_range_expected = pd.date_range(start=start_time, end=end_time, freq='H')
            expected_dframe = pd.DataFrame(index=time_range_expected)
            expected_dframe = expected_dframe.join(fields_dframe)

            if expected_dframe.__len__() != fields_dframe.__len__():
                log_stream.error(' ===> Time range consistency is not verified')
                raise RuntimeError('Check the time range and the time format to adjust the consistency')

            # search indexes for columns with values greater than 100
            fields_index = np.argwhere(fields_dframe.max().values > 100)[:, 0]
            # update dframe with scaling factor 10 for selected indexes
            fields_dframe.iloc[:, fields_index] = fields_dframe.iloc[:, fields_index] / 10
            # update columns for multiindex (not needed)
            # fields_dframe.columns = pd.MultiIndex.from_tuples(columns_obj)

        else:
            log_stream.info(' -----> Read file datasets "' + file_name + '" ... FAILED')
            log_stream.error(' ===> File format is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Read file datasets "' + file_name + '" ... DONE')

        '''
        # punto 71-72 buco dal 2023-06-06 23 al 2023-06-13 00
        plt.figure()
        fields_dframe['LaThuile']['2023-06-05':'2023-07-01'].plot()
        plt.show()
        '''

        return fields_dframe
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get registry object
    def get_obj_registry(self, file_name, file_fields=None):

        # info start method
        log_stream.info(' -----> Read file registry "' + file_name + '" ... ')

        if file_name.endswith('mat'):

            # get fields data raw
            fields_data_raw = read_file_mat(
                file_name, file_mandatory=True,
                file_fields_excluded=['__header__', '__version__', '__globals__', 'data', 'times'])
            # map fields
            fields_data_map = map_vars_dict(fields_data_raw, file_fields)

            # parser name(s)
            var_data_name = convert_var_mat(
                fields_data_map['name'], var_type='other',
                reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='names.workspace')
            fields_data_map.pop('name')

            # parser description(s)
            var_data_description = convert_var_mat(
                fields_data_map['description'], var_type='other',
                reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='description.workspace')
            fields_data_map.pop('description')

            # parser amm bnd 1
            var_data_amm_bnd_1 = convert_var_mat(
                fields_data_map['amm_level_1'], var_type='other',
                reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='amm_level_1.workspace')
            fields_data_map.pop('amm_level_1')

            # parser amm bnd 2
            var_data_amm_bnd_2 = convert_var_mat(
                fields_data_map['amm_level_2'], var_type='other',
                reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='amm_level_2.workspace')
            fields_data_map.pop('amm_level_2')

            # parser tag
            if 'tag' not in list(fields_data_map.keys()):
                var_data_tag = []
                for string_name in var_data_name:
                    string_tag = sanitize_string(string_name)
                    var_data_tag.append(string_tag)
            else:
                var_data_tag = convert_var_mat(
                    fields_data_map['tag'], var_type='other',
                    reset_tmp=True, folder_tmp=self.folder_name_tmp, file_tmp='tag.workspace')
                fields_data_map.pop('tag')

            # define valid fields
            fields_n = np.squeeze(fields_data_map['altitude']).shape[0]
            var_data_valid = np.ones(shape=[fields_n], dtype=int)

            # create fields data update
            fields_data_upd = {
                'altitude': np.squeeze(fields_data_map['altitude']),
                'longitude': np.squeeze(fields_data_map['longitude']),
                'latitude': np.squeeze(fields_data_map['latitude']),
                'code': np.squeeze(fields_data_map['code']),
                'name': var_data_name, 'description': var_data_description,
                'amm_level_1': var_data_amm_bnd_1, 'amm_level_2': var_data_amm_bnd_2,
                'valid': var_data_valid, 'tag': var_data_tag
            }

            # create fields data merge
            fields_data_merge = {**fields_data_map, **fields_data_upd}
            # create fields dframe
            fields_dframe = pd.DataFrame(data=fields_data_merge)
        else:
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
            if self.format_source == 'drops_dr':

                # get registry dframe
                dframe_registry = self.get_obj_registry(file_path_registry, self.fields_registry)
                # get datasets dframe
                dframe_datasets = self.get_obj_datasets(file_path_datasets, self.fields_datasets, dframe_registry)

                # organize dframe obj
                dframe_obj = {'registry': dframe_registry, 'datasets': dframe_datasets}

                # dump dframe obj
                folder_name_ancillary, file_name_ancillary = os.path.split(file_path_ancillary)
                make_folder(folder_name_ancillary)
                write_obj(file_path_ancillary, dframe_obj)

            else:
                log_stream.error(' ===> Source format "' + self.format_source + '" is not supported')
                raise NotImplemented('Case not implemented yet')

        else:
            # read source obj
            dframe_obj = read_obj(file_path_ancillary)

        # method end info
        log_stream.info(' ----> Organize source object(s) ... DONE')

        return dframe_obj

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
