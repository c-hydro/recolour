"""
Class Features

Name:          driver_data_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd
from copy import deepcopy

from lib_data_io_pickle import read_obj, write_obj
from lib_data_io_json import read_file_json, write_file_json
from lib_data_io_generic import (adjust_data_point,
                                 convert_data_to_vars, convert_vars_to_data, convert_vars_to_dict)
from lib_data_io_csv import read_file_csv, write_file_csv

from lib_utils_analysis import (apply_time_series_scaling, apply_time_series_filter, apply_time_series_metrics,
                                join_time_series)
from lib_utils_graph import configure_time_series_info, view_time_series
from lib_utils_obj import create_dict_from_list
from lib_utils_system import fill_tags2string, make_folder

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
    def __init__(self, time_run, time_range,
                 static_obj, source_dict, ancillary_dict, destination_dict=None,
                 flags_dict=None, template_dict=None, params_dict=None, methods_dict=None, tmp_dict=None):

        self.time_run = time_run
        self.time_range = time_range

        self.time_start = pd.DatetimeIndex([time_range[0], time_range[-1]]).min()
        self.time_end = pd.DatetimeIndex([time_range[0], time_range[-1]]).max()

        self.time_end, self.time_start = time_range[0], time_range[-1]
        self.static_obj = static_obj

        self.source_dict = source_dict
        self.ancillary_dict = ancillary_dict
        self.destination_dict = destination_dict

        self.flags_dict = flags_dict
        self.template_dict = template_dict
        self.params_dict = params_dict
        self.methods_dict = methods_dict
        self.tmp_dict = tmp_dict

        # data tag(s)
        self.file_name_tag, self.folder_name_tag, self.path_name_tag = 'file_name', 'folder_name', 'path_name'
        self.value_min_tag, self.value_max_tag, self.scale_factor_tag = 'value_min', 'value_max', 'scale_factor'
        self.value_nodata_tag = 'value_no_data'
        self.fields_tag, self.format_tag = 'fields', 'format'

        # figure tag(s)
        self.title_tag = 'title'
        self.label_axis_x_tag, self.label_axis_y_tag = 'label_axis_x', 'label_axis_y'
        self.legend_tag, self.style_tag = 'legend', "style"
        # other tags
        self.metrics_tag, self.group_data_tag, self.group_time_tag = 'metrics', 'groups_data', 'groups_time'

        # get reset flags
        self.reset_src = flags_dict['reset_dynamic_source']
        self.reset_dst = flags_dict['reset_dynamic_destination']

        # get static point
        self.point_obj = self.static_obj['point_obj']

        # get source data
        (folder_name_src_data, file_name_src_data, fields_src_data, format_src_data, scale_factor_src_data,
         value_min_src_data, value_max_src_data, value_no_data_src_data) = self.get_info_data(
            self.source_dict['data'])
        file_path_src_data = os.path.join(folder_name_src_data, file_name_src_data)
        # zip source data
        self.dset_obj_src_data = self.zip_info_data(
            folder_name_src_data, file_name_src_data, file_path_src_data,
            fields_src_data, format_src_data, scale_factor_src_data,
            value_min_src_data, value_max_src_data, value_no_data_src_data)

        # get source metrics
        (folder_name_src_met, file_name_src_met, fields_src_met, format_src_met, _, _, _, _) = self.get_info_data(
            self.source_dict['metrics'])
        file_path_src_met = os.path.join(folder_name_src_met, file_name_src_met)
        # zip source metrics
        self.dset_obj_src_met = self.zip_info_data(
            folder_name_src_met, file_name_src_met, file_path_src_met,
            fields_src_met, format_src_met)

        # get ancillary data
        (folder_name_anc, file_name_anc, _, _, _, _, _, _) = self.get_info_data(self.ancillary_dict)
        file_path_anc = os.path.join(folder_name_anc, file_name_anc)
        # zip ancillary data
        self.dset_obj_anc = self.zip_info_data(folder_name_anc, file_name_anc, file_path_anc)

        # get destination figure
        (folder_name_dst_fig, file_name_dst_fig, fields_dst_fig, format_dst_fig, scale_factor_dst_fig,
         _, _, _) = self.get_info_data(self.destination_dict)
        file_path_dst_fig = os.path.join(folder_name_dst_fig, file_name_dst_fig)
        # zip destination figure
        self.dset_obj_dst_file = self.zip_info_data(
            folder_name_dst_fig, file_name_dst_fig, file_path_dst_fig,
            fields_dst_fig, format_dst_fig)

        (title_dst_fig, label_axis_x_dst_fig, label_axis_y_dst_fig,
         fig_legend, fig_style) = self.get_info_figure(self.destination_dict)
        # zip tmp data
        self.dset_obj_dst_fig = self.zip_info_figure(
            title_dst_fig, label_axis_x_dst_fig, label_axis_y_dst_fig,
            fig_legend, fig_style)

        # get destination other
        self.metrics_settings = self.destination_dict[self.metrics_tag]
        self.groups_data_settings = self.destination_dict[self.group_data_tag]
        self.groups_time_settings = self.destination_dict[self.group_time_tag]

        if self.groups_data_settings is None:
            self.groups_data_settings = {'group_default': None}

        # get tmp data
        (folder_name_tmp, file_name_tmp, _, _, _, _, _, _) = self.get_info_data(self.tmp_dict)
        # zip tmp data
        self.dset_obj_tmp = self.zip_info_data(folder_name_tmp, file_name_tmp, None)

        # time-series name(s)
        self.ts_name_ref = 'ref'
        ts_name_tmp = list(fields_src_data.keys())

        if self.ts_name_ref not in ts_name_tmp:
            log_stream.error(' ===> Datasets reference tag "' + self.ts_name_ref +
                             '" must be in the source datasets')
            raise RuntimeError('Define the datasets reference tag in the source datasets')
        ts_name_other = []
        for ts_name_step in ts_name_tmp:
            if 'other' in ts_name_step:
                ts_name_other.append(ts_name_step)
        self.ts_name_other = sorted(ts_name_other)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to zip info data
    def zip_info_data(self, folder_name, file_name, path_name, fields=None, file_format=None, scale_factor=None,
                      value_min=None, value_max=None, value_no_data=None):
        info_obj = {self.folder_name_tag: folder_name, self.file_name_tag: file_name, self.path_name_tag: path_name,
                    self.fields_tag: fields, self.format_tag: file_format, self.scale_factor_tag: scale_factor,
                    self.value_min_tag: value_min, self.value_max_tag: value_max, self.value_nodata_tag: value_no_data}
        return info_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get info data
    def get_info_data(self, obj_data):

        folder_name = check_key_of_obj(self.folder_name_tag, obj_data, value_data_default=None)
        file_name = check_key_of_obj(self.file_name_tag, obj_data, value_data_default=None)
        fields = check_key_of_obj(self.fields_tag, obj_data, value_data_default={})
        scale_factor = check_key_of_obj(self.scale_factor_tag, obj_data, value_data_default=1)
        file_format = check_key_of_obj(self.format_tag, obj_data, value_data_default=None)
        v_min = check_key_of_obj(self.value_min_tag, obj_data, value_data_default=None)
        v_max = check_key_of_obj(self.value_max_tag, obj_data, value_data_default=None)
        v_no_data = check_key_of_obj(self.value_nodata_tag, obj_data, value_data_default=np.nan)

        if v_no_data is None:
            v_no_data = np.nan

        return folder_name, file_name, fields, file_format, scale_factor, v_min, v_max, v_no_data
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to zip info figure
    def zip_info_figure(self, fig_title, fig_label_axis_x, fig_label_axis_y, fig_legend, fig_style):
        info_obj = {self.title_tag: fig_title,
                    self.label_axis_x_tag: fig_label_axis_x, self.label_axis_y_tag: fig_label_axis_y,
                    self.legend_tag: fig_legend, self.style_tag: fig_style}
        return info_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get info data
    def get_info_figure(self, obj_data):

        fig_title = check_key_of_obj(self.title_tag, obj_data, value_data_default='title')
        fig_label_axis_x = check_key_of_obj(self.label_axis_x_tag, obj_data, value_data_default='x [-]')
        fig_label_axis_y = check_key_of_obj(self.label_axis_y_tag, obj_data, value_data_default='y [-]')
        fig_legend = check_key_of_obj(self.legend_tag, obj_data, value_data_default=None)
        fig_style = check_key_of_obj(self.style_tag, obj_data, value_data_default=None)

        return fig_title, fig_label_axis_x, fig_label_axis_y, fig_legend, fig_style

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to define file name
    def define_file_name(self, file_path_raw,
                         time_step=None, time_start=None, time_end=None,
                         point_name=None, point_code=None, point_tag=None):

        if time_step is None:
            time_step = deepcopy(time_end)

        template_time_generic = self.template_dict['time']
        template_datasets_generic = self.template_dict['datasets']

        generic_time_list = list(template_time_generic.keys())
        values_time_dict_generic = create_dict_from_list(generic_time_list, time_step)

        if time_end is not None:
            time_end_list = [key for key, value in values_time_dict_generic.items() if 'end' in key.lower()]
            values_dict_tmp = create_dict_from_list(time_end_list, time_end)
            values_time_dict_generic.update(**values_dict_tmp)

        if time_start is not None:
            time_start_list = [key for key, value in values_time_dict_generic.items() if 'start' in key.lower()]
            values_dict_tmp = create_dict_from_list(time_start_list, time_start)
            values_time_dict_generic.update(**values_dict_tmp)

        # generic_datasets_list = list(template_datasets_generic.keys())
        values_datasets_dict_generic = {}
        if point_name is not None:
            values_datasets_dict_generic['point_name'] = point_name
        if point_code is not None:
            values_datasets_dict_generic['point_code'] = str(point_code)
        if point_tag is not None:
            values_datasets_dict_generic['point_tag'] = point_tag

        file_path_def = fill_tags2string(file_path_raw, template_time_generic, values_time_dict_generic)[0]
        file_path_def = fill_tags2string(file_path_def, template_datasets_generic, values_datasets_dict_generic)[0]

        return file_path_def

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to view data
    def view_data(self, point_data_collections):

        # method start info
        log_stream.info(' ----> View time-series object(s) ... ')

        # get time
        time_run = self.time_run

        # get method(s) object
        point_methods = self.methods_dict
        # get static object
        point_registry = self.point_obj

        # get time-series name(s)
        ts_name_ref, ts_name_other = self.ts_name_ref, self.ts_name_other

        # get file settings
        file_path_dst_raw = self.dset_obj_dst_file[self.path_name_tag]
        file_format_dst = self.dset_obj_dst_file[self.format_tag]
        file_fields_dst = self.dset_obj_dst_file[self.fields_tag]
        # get figure settings
        fig_legend = self.dset_obj_dst_fig[self.legend_tag]
        fig_style = self.dset_obj_dst_fig[self.style_tag]
        fig_title = self.dset_obj_dst_fig[self.title_tag]
        fig_label_x_tag = self.dset_obj_dst_fig[self.label_axis_x_tag]
        fig_label_y_tag = self.dset_obj_dst_fig[self.label_axis_y_tag]

        # iterate over point(s)
        for point_fields in point_registry.to_dict(orient="records"):

            # debug
            # point_fields = point_registry.to_dict('records')[8] # valzemola

            # get point information
            point_tag, point_name, point_code = point_fields['tag'], point_fields['name'], point_fields['code']
            # point_geo_x, point_geo_y = point_fields['longitude'], point_fields['latitude']

            # point info start
            log_stream.info(' -----> Point "' + point_tag + '" ... ')

            # define file name
            file_path_dst_point = self.define_file_name(
                file_path_dst_raw,
                time_step=time_run, time_start=None, time_end=None,
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            # get point data and metrics
            point_data_raw = point_data_collections[point_tag]['data']
            point_metrics_raw = point_data_collections[point_tag]['metrics']
            # select point data and metrics
            point_data_selected, point_metrics_selected = configure_time_series_info(
                point_data_raw, point_metrics_raw,
                fields=self.dset_obj_dst_file['fields'], metrics=self.metrics_settings)

            # iterate over groups
            for group_tag, group_data in self.groups_data_settings.items():

                # group info start
                log_stream.info(' ------> Group "' + group_tag + '" ... ')

                # get point data and metrics by group
                if group_data is not None:
                    point_data_group = point_data_selected[group_data]

                    point_metrics_group = {}
                    for group_step in group_data:
                        point_metrics_group[group_step] = point_metrics_selected[group_step]

                else:
                    point_data_group = deepcopy(point_data_selected)
                    point_metrics_group = deepcopy(point_metrics_selected)

                # view point time-series
                view_time_series(file_path_dst_point,
                                 point_ts_data=point_data_group, point_metrics=point_metrics_group,
                                 point_registry=point_fields,
                                 file_fields=file_fields_dst,
                                 file_groups_name=group_tag, file_groups_time=self.groups_time_settings,
                                 file_metrics=self.metrics_settings,
                                 fig_title=fig_title,
                                 fig_legend=fig_legend, fig_style=fig_style,
                                 fig_label_axis_x=fig_label_x_tag, fig_label_axis_y=fig_label_y_tag,
                                 fig_dpi=150)

                # group info end
                log_stream.info(' ------> Group "' + group_tag + '" ... DONE')

            # point info end
            log_stream.info(' -----> Point "' + point_tag + '" ... DONE')

        # method end info
        log_stream.info(' ----> View time-series object(s) ... DONE.')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to organize data
    def organize_data(self):

        # method start info
        log_stream.info(' ----> Organize time-series object(s) ... ')

        # get time(s)
        time_run = self.time_run
        # get static object
        point_registry = self.point_obj

        # get path(s)
        file_path_src_data_raw = self.dset_obj_src_data[self.path_name_tag]
        file_path_src_met_raw = self.dset_obj_src_met[self.path_name_tag]
        file_path_anc_raw = self.dset_obj_anc[self.path_name_tag]
        file_path_dst_raw = self.dset_obj_dst_file[self.path_name_tag]

        # get format(s)
        file_format_src_data = self.dset_obj_src_data[self.format_tag]
        file_format_src_met = self.dset_obj_src_met[self.format_tag]

        # get flag(s)
        reset_src, reset_dst = self.reset_src, self.reset_dst

        # iterate over point(s)
        point_data_collections = {}
        for point_fields in point_registry.to_dict(orient="records"):

            # get point information
            point_tag, point_name, point_code = point_fields['tag'], point_fields['name'], point_fields['code']
            # point_geo_x, point_geo_y = point_fields['longitude'], point_fields['latitude']

            # point info start
            log_stream.info(' -----> Point "' + point_tag + '" ... ')

            # define source, ancillary and destination filename(s)
            file_path_src_data_point = self.define_file_name(
                file_path_src_data_raw,
                time_step=time_run, time_start=None, time_end=None,
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            file_path_src_met_point = self.define_file_name(
                file_path_src_met_raw,
                time_step=time_run, time_start=None, time_end=None,
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            file_path_anc_point = self.define_file_name(
                file_path_anc_raw,
                time_step=time_run, time_start=None, time_end=None,
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            file_path_dst_point = self.define_file_name(
                file_path_dst_raw,
                time_step=time_run, time_start=None, time_end=None,
                point_name=point_tag, point_code=point_code, point_tag=point_tag)

            # reset source file
            if reset_src:
                if os.path.exists(file_path_anc_point):
                    os.remove(file_path_anc_point)
                if os.path.exists(file_path_dst_point):
                    os.remove(file_path_dst_point)

            # check ancillary file availability
            if not os.path.exists(file_path_anc_point):

                # get data start
                log_stream.info(' ------> Get datasets ... ')

                # check source reference file availability
                if os.path.exists(file_path_src_data_point) and (os.path.exists(file_path_src_met_point)):

                    # check source data file format
                    if file_format_src_data == 'json':

                        # method to get data
                        point_data_raw = read_file_json(file_path_src_data_point)
                        # method to convert data
                        point_data_converted = convert_data_to_vars(
                            point_data_raw, obj_fields=self.dset_obj_src_data['fields'])

                    elif file_format_src_data == 'csv':

                        # method to get data
                        point_data_raw = read_file_csv(file_path_src_data_point, dframe_index='time')
                        # method to convert data
                        point_data_converted = convert_data_to_vars(
                            point_data_raw, obj_fields=self.dset_obj_src_data['fields'])

                    else:
                        log_stream.error(' ===> Source data type "' + file_format_src_data + '" is not supported.')
                        raise NotImplemented('Case not implemented yet')

                    # check source metrics file format
                    if file_format_src_met == 'csv':

                        # method to get metrics
                        point_met_raw = read_file_csv(file_path_src_met_point, dframe_index='metrics')
                        # method to convert metrics
                        point_met_converted = convert_data_to_vars(
                            point_met_raw, obj_fields=self.dset_obj_src_met['fields'])

                    else:
                        log_stream.error(' ===> Source metrics type "' + file_format_src_met + '" is not supported.')
                        raise NotImplemented('Case not implemented yet')

                else:
                    log_stream.warning(' ===> Datasets file "' + file_path_src_data_point + '" was not available.')
                    point_data_converted, point_met_converted = None, None

                # get data end
                log_stream.info(' ------> Get datasets ... DONE')

                # dump data start
                log_stream.info(' ------> Dump datasets ... ')
                # check data availability
                if (point_data_converted is not None) and (point_met_converted is not None):

                    # method to dump data
                    folder_name_anc_point, file_name_anc_point = os.path.split(file_path_anc_point)
                    make_folder(folder_name_anc_point)
                    # collect data and metrics
                    point_data_obj = {'metrics': point_met_converted, 'data': point_data_converted}
                    write_obj(file_path_anc_point, point_data_obj)

                    # dump data end
                    log_stream.info(' ------> Dump datasets ... DONE')

                else:
                    # dump data end
                    point_data_obj = None
                    log_stream.info(' ------> Dump datasets ... SKIPPED. Datasets is not available')

                # point info end
                log_stream.info(' -----> Point "' + point_tag + '" ... DONE')

            else:

                # point info end
                point_data_obj = read_obj(file_path_anc_point)
                log_stream.info(' -----> Point "' + point_tag + '" ... SKIPPED. Data previously saved')

            point_data_collections[point_tag] = point_data_obj

        # method end info
        log_stream.info(' ----> Organize time-series object(s) ... DONE')

        return point_data_collections

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to check key in a dictionary obj
def check_key_of_obj(key_data, obj_data, value_data_default=None):
    if key_data in list(obj_data.keys()):
        return obj_data[key_data]
    else:
        return value_data_default
# -------------------------------------------------------------------------------------
