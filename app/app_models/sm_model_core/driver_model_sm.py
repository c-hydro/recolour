"""
Class Features

Name:          driver_model_sm
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260615'
Version:       '2.0.0
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd

from lib_data_io_csv import read_datasets_csv, read_metrics_csv, write_datasets_csv, write_metrics_csv

from lib_utils_io import fill_string_with_time, fill_string_with_info
from lib_utils_generic import make_folder

from lib_model_utils import (filter_model_data, organize_model_data, organize_model_parameters,
                             organize_model_results, organize_model_metrics, filter_model_results,
                             plot_model_results_ts, plot_model_results_maps)

from lib_model_core import SMestim_IE_03 as fx_sm_model

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver model
class DriverModel:

    # -------------------------------------------------------------------------------------
    # initialize class
    def __init__(self, time_reference, time_run,
                 alg_data_static, alg_data_dynamic,
                 alg_model, alg_info, alg_template, alg_flags):

        # set time reference
        self.time_reference = time_reference
        self.time_run = time_run

        # set data static object(s)
        self.data_registry = alg_data_static['registry']
        self.data_vars = alg_data_dynamic

        # set algorithm information
        self.alg_flags = alg_flags
        self.alg_info = alg_info
        self.alg_model_data = alg_data_dynamic['destination']
        self.alg_model_results = alg_model['results']
        self.alg_model_metrics = alg_model['metrics']
        self.alg_model_figure = alg_model['figure']
        self.alg_template_time = alg_template['time']
        self.alg_template_datasets = alg_template['datasets']

        # reset flags
        self.reset_model_results = self.alg_flags['reset_model_results']
        self.reset_model_metrics = self.alg_flags['reset_model_metrics']
        self.reset_model_figure = self.alg_flags['reset_model_figure']

        # registry and datasets tag(s)
        self.file_name_tag, self.folder_name_tag = 'file_name', 'folder_name'
        self.time_tag, self.filters_tag = 'time', 'filters'
        self.fields_tag, self.format_tag = 'fields', 'format'

        # model data object(s)
        self.folder_name_data = self.alg_model_data['folder_name']
        self.file_name_data = self.alg_model_data['file_name']
        self.format_data = self.alg_model_data[self.format_tag]
        self.time_data = self.alg_model_data[self.time_tag]
        self.fields_data = self.alg_model_data[self.fields_tag]
        self.file_path_data = os.path.join(self.folder_name_data, self.file_name_data)

        # model results object(s)
        self.folder_name_results = self.alg_model_results['folder_name']
        self.file_name_results = self.alg_model_results['file_name']
        self.format_results = self.alg_model_results[self.format_tag]
        self.time_results = self.alg_model_results[self.time_tag]
        self.fields_results = self.alg_model_results[self.fields_tag]
        self.file_path_results = os.path.join(self.folder_name_results, self.file_name_results)
        # model metrics object(s)
        self.folder_name_metrics = self.alg_model_metrics['folder_name']
        self.file_name_metrics = self.alg_model_metrics['file_name']
        self.format_metrics = self.alg_model_metrics[self.format_tag]
        self.time_metrics = self.alg_model_metrics[self.time_tag]
        self.fields_metrics = self.alg_model_metrics[self.fields_tag]
        self.file_path_metrics = os.path.join(self.folder_name_metrics, self.file_name_metrics)

        # model figure object(s)
        self.folder_name_figure = self.alg_model_figure['folder_name']
        self.file_name_figure = self.alg_model_figure['file_name']
        self.time_figure = self.alg_model_figure[self.time_tag]
        self.fields_figure = self.alg_model_figure[self.fields_tag]
        self.format_figure = self.alg_model_figure[self.format_tag]
        self.file_path_figure = os.path.join(self.folder_name_figure, self.file_name_figure)

        self.mode_figure = self.alg_model_figure.get('mode', 'time_series')
        self.settings_figure = self.alg_model_figure.get('settings', {})

        # model figure options (time_series or maps)
        self.time_series_figure = self.settings_figure["time_series"]
        self.time_series_dpi_figure = self.time_series_figure["dpi"]
        self.time_series_spacing_x_figure = self.time_series_figure["spacing_x"]

        self.maps_figure = self.settings_figure["maps"]
        self.maps_dpi_figure = self.maps_figure["dpi"]
        self.maps_time_select_figure = self.maps_figure["time_select"]

        # add keys to the figure fields
        self.fields_figure["time_reference"] = time_reference

        self.show_figure = False

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get datasets object
    def get_obj_datasets(self, file_name, file_format='csv',
                         file_fields=None, time_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' ------> Read model datasets ... ')

        # check file existence
        if not os.path.exists(file_name):
            log_stream.error(' ===> File does not exist')
            raise IOError('File parameters must be available')

        # check file format
        if file_format == 'csv':

            # time fields
            if time_fields is None:
                time_fields = {}

            # get datasets in ascii format
            fields_obj = read_datasets_csv(
                file_name,
                time_reference=self.time_reference, time_format='%Y-%m-%d %H:%M',
                file_fields=file_fields, registry_fields=registry_fields,
                file_sep=',', file_decimal='.', **time_fields)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Read model datasets ... DONE')

        return fields_obj

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get datasets object
    def get_obj_metrics(self, file_name, file_format='csv', file_type='datasets',
                         file_fields=None, time_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' ------> Read model metrics ... ')

        # check file existence
        if not os.path.exists(file_name):
            log_stream.error(' ===> File does not exist')
            raise IOError('File parameters must be available')

        # check file format
        if file_format == 'csv':

            # get metrics in ascii format
            fields_obj = read_metrics_csv(
                file_name,
                file_fields=file_fields, registry_fields=registry_fields,
                file_sep=',', file_decimal='.')

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Read model metrics ... DONE')

        return fields_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump datasets object
    def dump_obj_datasets(self, file_name, file_dframe, file_format='csv',
                          file_fields=None, time_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' ------> Dump model datasets "' + file_name + '" ... ')

        # check file format
        if file_format == 'csv':

            # dump combined dframe
            folder_name, _ = os.path.split(file_name)
            make_folder(folder_name)

            # write datasets in csv format
            write_datasets_csv(
                file_name, file_dframe, file_fields=file_fields, time_fields=time_fields,
                dframe_index_label='time', dframe_index_format='%Y-%m-%d %H:%M',
                dframe_sep=';', dframe_decimal='.', dframe_float_format='%.3f',
                dframe_index=True, dframe_header=True)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Dump model datasets "' + file_name + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump metrics object
    def dump_obj_metrics(self, file_name, file_dframe, file_format='csv', file_fields=None):

        # info start method
        log_stream.info(' ------> Dump model metrics "' + file_name + '" ... ')

        # check file format
        if file_format == 'csv':

            # dump combined dframe
            folder_name, _ = os.path.split(file_name)
            make_folder(folder_name)

            # write datasets in csv format
            write_metrics_csv(file_name, file_dframe, file_fields=file_fields)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Dump model metrics "' + file_name + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to plot datasets object
    def plot_obj_datasets(self, file_name, dframe_results, dframe_metrics, file_format='png'):

        # info start method
        log_stream.info(' -----> Plot model results "' + file_name + '" ... ')

        # check file format time_series
        if self.mode_figure == "time_series" and file_format == "png":

            # create folder (if needed)
            folder_name, _ = os.path.split(file_name)
            os.makedirs(folder_name, exist_ok=True)

            # plot results in time-series format
            plot_model_results_ts(
                file_name, dframe_results, dframe_metrics,
                fig_spacing_x=self.time_series_spacing_x_figure, fig_fields=self.fields_figure,
                fig_dpi=self.time_series_dpi_figure, fig_show=self.show_figure)

        # check file format maps
        elif self.mode_figure == "maps" and file_format == "png":

            # create folder (if needed)
            folder_name, _ = os.path.split(file_name)
            os.makedirs(folder_name, exist_ok=True)

            # plot results in maps format
            plot_model_results_maps(
                file_name, dframe_results,
                fig_fields=self.fields_figure,
                fig_dpi=self.maps_dpi_figure, fig_show=self.show_figure)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' -----> Plot model results "' + file_name + '" ... DONE')

    # ------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to define file string
    def __define_file_string(self, file_string_tmpl, extended_info=None):

        if extended_info is not None:
            alg_info = {**self.alg_info, **extended_info}
        else:
            alg_info = self.alg_info

        file_string_def = fill_string_with_time(file_string_tmpl, self.time_reference, self.alg_template_time)
        file_string_def = fill_string_with_info(file_string_def, alg_info, self.alg_template_datasets)
        return file_string_def
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to execution model
    def exec(self):

        # method start info
        log_stream.info(' ----> Execution model ... ')

        # get time reference
        time_step_reference = self.time_reference

        # get data object(s)
        data_registry = self.data_registry
        data_vars = self.data_vars

        # get path(s
        file_path_data_tmpl = self.file_path_data
        file_path_results_tmpl = self.file_path_results
        file_path_metrics_tmpl = self.file_path_metrics
        file_path_figure_tmpl = self.file_path_figure

        # get flag(s)
        reset_model_results = self.reset_model_results
        reset_model_metrics = self.reset_model_metrics

        # iterate over geo point(s)
        results_collections, metrics_collections = {}, {}
        for fields_registry in data_registry.to_dict(orient="records"):

            # debug (jesi == 2 in this case
            # fields_registry = data_registry.iloc[2].to_dict()

            # get point information
            point_name, point_tag = fields_registry['name'], fields_registry['tag']

            # info point start
            log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' + point_tag + '" ... ')

            # method to fill the filename(s)
            file_path_data_point = self.__define_file_string(
                file_path_data_tmpl, extended_info={'point_name': point_tag})
            file_path_results_point = self.__define_file_string(
                file_path_results_tmpl, extended_info={'point_name': point_tag})
            file_path_metrics_point = self.__define_file_string(
                file_path_metrics_tmpl, extended_info={'point_name': point_tag})
            file_path_figure_point = self.__define_file_string(
                file_path_figure_tmpl, extended_info={'point_name': point_tag})

            # reset ancillary file if required
            if reset_model_results or reset_model_metrics:
                if os.path.exists(file_path_results_point):
                    os.remove(file_path_results_point)
                if os.path.exists(file_path_metrics_point):
                    os.remove(file_path_metrics_point)
                if os.path.exists(file_path_figure_point):
                    os.remove(file_path_figure_point)

            # check results file availability
            if not os.path.exists(file_path_results_point):

                # check data file availability
                if os.path.exists(file_path_data_point):

                    # get dataframe obj
                    dframe_data = self.get_obj_datasets(
                        file_path_data_point, file_format='csv',
                        time_fields=None,
                        file_fields=None, registry_fields=data_registry)

                    # filter model data
                    dframe_data = filter_model_data(dframe_data, dframe_fields=self.fields_data)
                    # organize model data
                    values_data, values_time = organize_model_data(dframe_data)
                    # organize model parameters
                    values_params = organize_model_parameters(fields_registry)

                    # apply sm model
                    (values_theta, values_ns, values_ns_ln_q, values_ns_rad_q,
                     values_kge, values_rmse, values_rq) = fx_sm_model(values_time, values_data, values_params)

                    # organize result object
                    dframe_result = organize_model_results(
                        dframe_data, values_theta, values_time, dframe_fields=self.fields_results)

                    # dump result object
                    self.dump_obj_datasets(
                        file_path_results_point, dframe_result, file_format=self.format_results,
                        file_fields=self.fields_results, time_fields=self.time_results, registry_fields=fields_registry)
                    # store file results
                    results_collections[point_tag] = file_path_results_point

                    # dump metrics object
                    dframe_metrics = organize_model_metrics(
                        data_metrics={
                            'ns': values_ns, 'ns_ln_q': values_ns_ln_q, 'ns_rad_q': values_ns_rad_q,
                            'kge': values_kge, 'rmse': values_rmse, 'rq': values_rq},
                        data_time={'time': time_step_reference},
                        data_registry=fields_registry,
                        data_fields=self.fields_metrics)

                    # dump metrics object
                    self.dump_obj_metrics(file_path_metrics_point, dframe_metrics, file_format=self.format_metrics)
                    # store file results
                    metrics_collections[point_tag] = file_path_metrics_point

                    # info point end (done)
                    log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                                    point_tag + '" ... DONE')

                else:
                    # info point end (failed)
                    log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                                    point_tag + '" ... FAILED. Datasets are not available')

            else:

                # info point end (skipped)
                log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                                point_tag + '" ... SKIPPED. Results and metrics previously saved.')

                # store file results
                results_collections[point_tag] = file_path_results_point
                # store file results
                metrics_collections[point_tag] = file_path_metrics_point

        # method start info
        log_stream.info(' ----> Execution model ... DONE')

        return results_collections, metrics_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to view results
    def view(self, alg_model_results, alg_model_metrics):

        # method start info
        log_stream.info(' ----> View model ... ')

        # get view mode
        mode_figure = self.alg_model_figure.get('mode', 'time_series')

        # select method to view results
        if mode_figure == 'time_series':

            # create time-series
            self._view_time_series(
                alg_model_results=alg_model_results,
                alg_model_metrics=alg_model_metrics
            )

        elif mode_figure == 'maps':

            # create maps
            self._view_maps(
                alg_model_results=alg_model_results,
                alg_model_metrics=alg_model_metrics
            )

        else:

            log_stream.error(
                ' ===> View model ... SKIPPED. Figure mode "' + str(mode_figure) + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # method end info
        log_stream.info(' ----> View model ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to view time-series results
    def _view_time_series(self, alg_model_results, alg_model_metrics):

        log_stream.info(' ----> View model time-series ... ')

        data_registry = self.data_registry

        file_path_results_tmpl = self.file_path_results
        file_path_metrics_tmpl = self.file_path_metrics
        file_path_figure_tmpl = self.file_path_figure

        reset_model_figure = self.reset_model_figure

        for fields_registry in data_registry.to_dict(orient="records"):

            point_name, point_tag = fields_registry['name'], fields_registry['tag']

            log_stream.info(
                ' -----> Point -- (1) Name: "' + point_name +
                '" :: (2) Tag: "' + point_tag + '" ... '
            )

            file_path_results_point = self.__define_file_string(
                file_path_results_tmpl,
                extended_info={'point_name': point_tag}
            )

            file_path_metrics_point = self.__define_file_string(
                file_path_metrics_tmpl,
                extended_info={'point_name': point_tag}
            )

            file_path_figure_point = self.__define_file_string(
                file_path_figure_tmpl,
                extended_info={'point_name': point_tag}
            )

            if reset_model_figure:
                if os.path.exists(file_path_figure_point):
                    os.remove(file_path_figure_point)

            if os.path.exists(file_path_results_point) and os.path.exists(file_path_metrics_point):

                dframe_results = self.get_obj_datasets(
                    file_path_results_point,
                    file_format='csv',
                    time_fields=None,
                    file_fields=None,
                    registry_fields=fields_registry
                )

                dframe_metrics = self.get_obj_metrics(
                    file_path_metrics_point,
                    file_format='csv',
                    time_fields=None,
                    file_fields=None,
                    registry_fields=fields_registry
                )

                self.plot_obj_datasets(
                    file_path_figure_point,
                    dframe_results,
                    dframe_metrics
                )

                log_stream.info(' -----> Point "' + point_tag + '" ... DONE')

            else:

                log_stream.info(
                    ' -----> Point "' + point_tag +
                    '" ... SKIPPED. Datasets not available'
                )

        log_stream.info(' ----> View model time-series ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to view map results
    def _view_maps(self, alg_model_results, alg_model_metrics):

        # method info start
        log_stream.info(' ----> View model maps ... ')

        data_registry = self.data_registry

        file_path_results_tmpl = self.file_path_results
        file_path_figure_tmpl = self.file_path_figure

        reset_model_figure = self.reset_model_figure

        time_select = self.maps_time_select_figure

        # merge dataset info start
        log_stream.info(' -----> Merge datasets ... ')
        map_collections = []
        for fields_registry in data_registry.to_dict(orient="records"):

            point_name, point_tag = fields_registry['name'], fields_registry['tag']

            log_stream.info(
                ' ------> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' + point_tag + '" ... ')

            file_path_results_point = self.__define_file_string(
                file_path_results_tmpl,
                extended_info={'point_name': point_tag}
            )

            # check file results availability
            if os.path.exists(file_path_results_point):

                # read file results
                dframe_results = self.get_obj_datasets(
                    file_path_results_point,
                    file_format='csv', time_fields=None, file_fields=None, registry_fields=fields_registry)

                # check dframe availability
                if dframe_results is not None and not dframe_results.empty:

                    # filter model results
                    dframe_results = filter_model_results(dframe_results)

                    # compute last strict
                    last_strict = dframe_results.dropna(how='any').index.max()
                    # compute last mod
                    last_mod = (dframe_results.dropna(subset=['rain', 'air_temperature', 'theta_simulated']).index.max())
                    # compute last mod
                    last_obs = (dframe_results.dropna(subset=['rain', 'air_temperature', 'theta_observed']).index.max())

                    # select timestamp
                    if time_select == 'last_strict':
                        time_last = last_strict
                    elif time_select == 'last_mod':
                        time_last = last_mod
                    elif time_select == 'last_obs':
                        time_last = last_obs
                    else:
                        log_stream.error(
                            f' ===> Time selection for mapping results "{time_select}" is not supported'
                        )
                        raise NotImplementedError('Case not implemented yet')

                    # get dataframe row
                    dframe_last = dframe_results.loc[[time_last]]

                    # organize datasets for maps
                    if not dframe_last.empty:

                        for field_key, field_value in fields_registry.items():
                            dframe_last[field_key] = field_value

                        map_collections.append(dframe_last)

                        log_stream.info(
                            ' ------> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                            point_tag + '" ... DONE')

                    else:

                        log_stream.info(
                            ' ------> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                            point_tag + '" ... SKIPPED. Last data not available')

                else:

                    log_stream.info(
                        ' ------> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                        point_tag + '" ... SKIPPED. Results dataframe is empty')


            else:

                log_stream.info(
                    ' ------> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                    point_tag + '" ... SKIPPED. Results file not available')

        # merge dataset info end
        log_stream.info(' -----> Merge datasets ... DONE')

        # plot dataset info start
        log_stream.info(' -----> Plot datasets ... ')

        # check map collections
        if map_collections:

            # concat results
            dframe_map = pd.concat(map_collections, axis=0)

            # define file data
            file_path_figure_map = self.__define_file_string(
                file_path_figure_tmpl,
                extended_info={'point_name': 'map'}
            )

            # delete file data
            if reset_model_figure:
                if os.path.exists(file_path_figure_map):
                    os.remove(file_path_figure_map)

            # plot file data
            self.plot_obj_datasets(
                file_path_figure_map,
                dframe_map, None
            )

            # plot dataset info end
            log_stream.info(' -----> Plot datasets ... DONE')

        else:

            # plot dataset info end
            log_stream.info(' -----> Plot datasets ... SKIPPED. Data not available')

        # method info end
        log_stream.info(' ----> View model maps ... DONE')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
