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
from pathlib import Path

from lib_data_io_csv import read_datasets_csv, read_metrics_csv, write_datasets_csv, write_metrics_csv
from lib_data_io_netcdf import (write_datasets_nc, read_datasets_nc, select_datasets_by_point, check_datasets_nc)

from lib_data_io_netcdf import write_datasets_nc as write_results_nc
from lib_data_io_netcdf import write_datasets_nc as write_ancillary_nc

from lib_utils_io import fill_string_with_time, fill_string_with_info
from lib_utils_generic import make_folder

from lib_model_utils import (filter_model_data, organize_model_data, organize_model_parameters,
                             organize_model_results, organize_model_metrics, filter_model_results,
                             summarize_model_results, organize_model_auxiliary,
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
        self.alg_model_auxiliary = alg_model.get("auxiliary", alg_model.get("metrics"))
        self.alg_model_figure = alg_model['figure']
        self.alg_template_time = alg_template['time']
        self.alg_template_datasets = alg_template['datasets']

        # reset flags
        self.reset_model_results = self.alg_flags['reset_model_results']
        self.reset_model_auxiliary = (
            self.alg_flags)["reset_model_auxiliary"] = self.alg_flags.get(
            "reset_model_auxiliary", self.alg_flags.get("reset_model_metrics", False))
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

        self.no_data_results = self.alg_model_results.get('no_data', -9999)
        self.fill_data_step_results_air_t = self.alg_model_results.get('fill_data_step_air_temperature', 2)
        self.fill_data_step_results_sm = self.alg_model_results.get('fill_data_step_soil_moisture', 2)

        # model metrics object(s)
        self.folder_name_auxiliary = self.alg_model_auxiliary['folder_name']
        self.file_name_auxiliary = self.alg_model_auxiliary['file_name']
        self.format_auxiliary = self.alg_model_auxiliary[self.format_tag]
        self.time_auxiliary = self.alg_model_auxiliary[self.time_tag]
        self.fields_auxiliary = self.alg_model_auxiliary[self.fields_tag]
        self.file_path_auxiliary = os.path.join(self.folder_name_auxiliary, self.file_name_auxiliary)

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

        self.collections_datasets_obj = None
        self.collections_results_obj = None

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get datasets object
    def get_obj_datasets(self, file_name, file_format='csv',
                         point_name='NA', point_tag='NA', point_longitude=None, point_latitude=None,
                         file_fields=None, time_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' ------> Read model datasets ... ')

        # check file existence
        if not os.path.exists(file_name):
            log_stream.error(' ===> File does not exist')
            raise IOError('File datasets must be available')

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

        # check file format
        elif file_format == 'netcdf':

            #`check if collections was previously read
            if self.collections_datasets_obj is None:
                self.collections_datasets_obj = read_datasets_nc(file_name=file_name)

            fields_obj = select_datasets_by_point(self.collections_datasets_obj, point_name=point_name)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Read model datasets ... DONE')

        return fields_obj

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get results object
    def get_obj_results(self, file_name, file_format='csv',
                         point_name='NA', point_tag='NA', point_longitude=None, point_latitude=None,
                         file_fields=None, time_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' ------> Read model results ... ')

        # check file existence
        if not os.path.exists(file_name):
            log_stream.error(' ===> File does not exist')
            raise IOError('File results must be available')

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

        # check file format
        elif file_format == 'netcdf':

            #`check if collections was previously read
            if self.collections_results_obj is None:
                self.collections_results_obj = read_datasets_nc(file_name=file_name)

            fields_obj = select_datasets_by_point(self.collections_results_obj, point_name=point_name)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Read model results ... DONE')

        return fields_obj

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to get datasets object
    def get_obj_metrics(self, file_name, file_format='csv', file_point='NA',
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

        elif file_format == 'netcdf':

            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Read model metrics ... DONE')

        return fields_obj
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to summarize obj datasets
    def summarize_obj_datasets(self, dframe_results):

        # Collect model-result data checks
        log_stream.info(" ------> Summarize model results ... ")

        # summarize results
        summary_results = summarize_model_results(dframe_results)

        # Collect model-result data checks
        log_stream.info(" ------> Summarize model results ... DONE")

        return summary_results

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump datasets object
    def dump_obj_results(self, file_name, file_dframe,
                          point_tag='NA', point_name='NA', point_longitude=-9999, point_latitude=-9999,
                          file_format='csv',
                          file_fields=None, time_fields=None, registry_fields=None):

        # info start method
        log_stream.info(' ------> Dump model results "' + file_name + '" ... ')

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

        elif file_format == 'netcdf':

            # dump combined dframe
            folder_name, _ = os.path.split(file_name)
            make_folder(folder_name)

            # write datasets in netcdf format
            write_results_nc(
                file_name=file_name,
                file_dframe=file_dframe, file_fields=file_fields,
                point_id=point_tag, point_name=point_name,
                longitude=point_longitude, latitude=point_latitude,
            )

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Dump model results "' + file_name + '" ... DONE')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to dump auxiliary object
    def dump_obj_auxiliary(self, file_name, file_dframe,
                           point_tag='NA', point_name='NA', point_longitude=-9999, point_latitude=-9999,
                           file_format='csv', file_fields=None):

        # info start method
        log_stream.info(' ------> Dump model auxiliary "' + file_name + '" ... ')

        # check file format
        if file_format == 'csv':

            # dump combined dframe
            folder_name, _ = os.path.split(file_name)
            make_folder(folder_name)

            # write datasets in csv format
            write_metrics_csv(file_name, file_dframe, file_fields=file_fields)

        elif file_format == 'netcdf':

            # dump combined dframe
            folder_name, _ = os.path.split(file_name)
            make_folder(folder_name)

            # write datasets in netcdf format
            write_ancillary_nc(
                file_name=file_name,
                file_dframe=file_dframe,
                point_id=point_tag, point_name=point_name,
                longitude=point_longitude, latitude=point_latitude,
                file_fields=file_fields)

        else:
            # exit with error if file format is not supported
            log_stream.error(' ===> File format "' + file_format + '" is not supported')
            raise NotImplemented('Case not implemented yet')

        # info end method
        log_stream.info(' ------> Dump model auxiliary "' + file_name + '" ... DONE')

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
        file_path_auxiliary_tmpl = self.file_path_auxiliary
        file_path_figure_tmpl = self.file_path_figure
        # get format
        format_data = self.format_data

        # get flag(s)
        reset_model_results = self.reset_model_results
        reset_model_auxiliary = self.reset_model_auxiliary

        # get first and last point to check dictionary collections
        point_tag_min, point_tag_max = data_registry["tag"].values[0], data_registry["tag"].values[-1]

        # iterate over geo point(s)
        reset_active_by_format = True
        results_collections, auxiliary_collections = {}, {}
        for fields_registry in data_registry.to_dict(orient="records"):

            # debug (jesi == 2 in this case
            # fields_registry = data_registry.iloc[2].to_dict()

            # get point information
            point_name, point_tag = fields_registry['name'], fields_registry['tag']
            point_longitude, point_latitude = fields_registry['longitude'], fields_registry['latitude']

            # method to fill the filename(s)
            if format_data == 'csv':

                # define tags
                point_tag_data = point_tag_results = point_tag_auxiliary = point_tag_figure = point_tag
                # define filenames
                file_path_data_point = self.__define_file_string(
                    file_path_data_tmpl, extended_info={'point_name': point_tag_data})
                file_path_results_point = self.__define_file_string(
                    file_path_results_tmpl, extended_info={'point_name': point_tag_results})
                file_path_auxiliary_point = self.__define_file_string(
                    file_path_auxiliary_tmpl, extended_info={'point_name': point_tag_auxiliary})
                file_path_figure_point = self.__define_file_string(
                    file_path_figure_tmpl, extended_info={'point_name': point_tag_figure})

                # reset ancillary file if required
                if reset_active_by_format:
                    if reset_model_results or reset_model_auxiliary:
                        if os.path.exists(file_path_results_point):
                            os.remove(file_path_results_point)
                        if os.path.exists(file_path_auxiliary_point):
                            os.remove(file_path_auxiliary_point)
                        if os.path.exists(file_path_figure_point):
                            os.remove(file_path_figure_point)
                    reset_active_by_format = True

            elif format_data == 'netcdf':

                # define tags
                point_tag_data = point_tag_results = point_tag_auxiliary = point_tag_figure = 'collections'
                # define filenames
                file_path_data_point = self.__define_file_string(
                    file_path_data_tmpl, extended_info={'point_name': point_tag_data})
                file_path_results_point = self.__define_file_string(
                    file_path_results_tmpl, extended_info={'point_name': point_tag_results})
                file_path_auxiliary_point = self.__define_file_string(
                    file_path_auxiliary_tmpl, extended_info={'point_name': point_tag_auxiliary})
                file_path_figure_point = self.__define_file_string(
                    file_path_figure_tmpl, extended_info={'point_name': point_tag_figure})

                # reset ancillary file if required
                if reset_active_by_format:
                    if reset_model_results or reset_model_auxiliary:
                        if os.path.exists(file_path_results_point):
                            os.remove(file_path_results_point)
                        if os.path.exists(file_path_auxiliary_point):
                            os.remove(file_path_auxiliary_point)
                        if os.path.exists(file_path_figure_point):
                            os.remove(file_path_figure_point)
                        reset_active_by_format = False

            else:
                log_stream.error(' ===> Format data format "' + format_data + '" is not supported')
                raise NotImplemented('Case not implemented yet')

            # define data searching
            if format_data == "csv":
                # One file for each point
                point_exists = os.path.exists(file_path_results_point)
                min_exists, max_exists = False, False

                # Get destination file extension
                file_ext = Path(file_path_results_point).suffix.lower()

                if file_ext != ".csv":
                    raise ValueError(
                        f"CSV results type requires a '.csv' destination file, found '{file_ext}'."
                    )

            elif format_data== "netcdf":
                # One shared NetCDF file containing all point variables
                point_exists, min_exists, max_exists = check_datasets_nc(
                    file_path_results_point, point_tag_results,
                    first_point_expected=point_tag_min, last_point_expected=point_tag_max)

                # Get destination file extension
                file_ext = Path(file_path_results_point).suffix.lower()

                if file_ext != ".nc":
                    raise ValueError(
                        f"NetCDF results type requires a '.nc' destination file, found '{file_ext}'."
                    )

            else:
                raise ValueError(
                    f"Destination file mode '{format_data}' is not supported. "
                    "Supported modes are 'unique' and 'collections'."
                )

            # check if min and max points are available (check for collections)
            if min_exists and max_exists:

                # info data start
                log_stream.info(f' -----> Point -- Collections from {point_tag_min} to {point_tag_max} ... ')
                results_collections[point_tag_results] = file_path_results_point
                log_stream.info(f' -----> Point -- Collections from {point_tag_min} to {point_tag_max} ... DONE')
                break

            else:

                # info data start
                log_stream.info(
                    ' -----> Point -- (1) Name: "' + point_name +
                    '" :: (2) Tag: "' + point_tag_data + '" ... '
                )

                # check results point availability
                if not point_exists:

                    # check data file availability
                    if os.path.exists(file_path_data_point):

                        # get dataframe obj
                        dframe_data = self.get_obj_datasets(
                            file_path_data_point, file_format=format_data,
                            time_fields=None,
                            point_name=point_name, point_tag=point_tag_data,
                            point_longitude=point_longitude, point_latitude=point_latitude,
                            file_fields=None, registry_fields=data_registry)

                        # filter model data
                        dframe_data = filter_model_data(
                            dframe_data, dframe_fields=self.fields_data,
                            interp_limit_sm=self.fill_data_step_results_sm,
                            interp_limit_airt=self.fill_data_step_results_air_t)

                        # organize model data
                        values_data, values_time = organize_model_data(dframe_data)
                        # organize model parameters
                        values_params = organize_model_parameters(fields_registry)

                        # apply sm model
                        (values_theta, values_ns, values_ns_ln_q, values_ns_rad_q,
                         values_kge, values_rmse, values_rq) = fx_sm_model(values_time, values_data, values_params)

                        # organize result object
                        model_result = organize_model_results(
                            dframe_data, values_theta, values_time, dframe_fields=self.fields_results)

                        # summarize result obj
                        model_summary = self.summarize_obj_datasets(model_result)
                        # organize model metrics
                        model_metrics = organize_model_metrics(
                            values_ns,
                            values_ns_ln_q, values_ns_rad_q,
                            values_kge, values_rmse, values_rq)

                        # dump result object
                        self.dump_obj_results(
                            file_path_results_point, model_result,
                            point_tag=point_tag, point_name=point_name,
                            point_longitude=point_longitude, point_latitude=point_latitude,
                            file_format=self.format_results,
                            file_fields=self.fields_results, time_fields=self.time_results, registry_fields=fields_registry)

                        # store dframe in unique collections
                        if point_tag_results not in results_collections.keys():
                            results_collections[point_tag_results] = file_path_results_point

                        # organize auxiliary object
                        dframe_auxiliary = organize_model_auxiliary(
                            data_metrics={**model_metrics, **model_summary},
                            data_time={'time': time_step_reference},
                            data_registry=fields_registry,
                            data_fields=self.fields_auxiliary)

                        # dump auxiliary object
                        self.dump_obj_auxiliary(file_path_auxiliary_point, dframe_auxiliary,
                                                point_tag=point_tag, point_name=point_name,
                                                point_longitude=point_longitude, point_latitude=point_latitude,
                                                file_format=self.format_auxiliary)

                        # store model auxiliary
                        if point_tag_auxiliary not in auxiliary_collections.keys():
                            auxiliary_collections[point_tag_results] = file_path_auxiliary_point

                        # info point end (done)
                        log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                                        point_tag_data + '" ... DONE')

                    else:
                        # info point end (failed)
                        log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                                        point_tag_data + '" ... FAILED. Datasets are not available')

                else:

                    # info point end (skipped)
                    log_stream.info(' -----> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' +
                                    point_tag_data + '" ... SKIPPED. Results and auxiliary previously saved.')

                    # store file results
                    results_collections[point_tag] = file_path_results_point
                    # store file results
                    auxiliary_collections[point_tag] = file_path_auxiliary_point

        # method start info
        log_stream.info(' ----> Execution model ... DONE')

        return results_collections, auxiliary_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to view results
    def view(self, alg_model_results, alg_model_auxiliary):

        # method start info
        log_stream.info(' ----> View model ... ')

        # get view mode
        mode_figure = self.alg_model_figure.get('mode', 'time_series')

        # select method to view results
        if mode_figure == 'time_series':

            # create time-series
            self._view_time_series(
                alg_model_results=alg_model_results,
                alg_model_metrics=alg_model_auxiliary
            )

        elif mode_figure == 'maps':

            # create maps
            self._view_maps(
                alg_model_results=alg_model_results,
                alg_model_metrics=alg_model_auxiliary
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
        file_path_auxiliary_tmpl = self.file_path_auxiliary
        file_path_figure_tmpl = self.file_path_figure

        format_results = self.format_results
        format_auxiliary = self.format_auxiliary
        reset_model_figure = self.reset_model_figure

        for fields_registry in data_registry.to_dict(orient="records"):

            # get point information
            point_name, point_tag = fields_registry['name'], fields_registry['tag']
            point_longitude, point_latitude = fields_registry['longitude'], fields_registry['latitude']

            log_stream.info(
                ' -----> Point -- (1) Name: "' + point_name +
                '" :: (2) Tag: "' + point_tag + '" ... '
            )

            file_path_results_point = self.__define_file_string(
                file_path_results_tmpl,
                extended_info={'point_name': point_tag}
            )

            file_path_auxiliary_point = self.__define_file_string(
                file_path_auxiliary_tmpl,
                extended_info={'point_name': point_tag}
            )

            file_path_figure_point = self.__define_file_string(
                file_path_figure_tmpl,
                extended_info={'point_name': point_tag}
            )

            if reset_model_figure:
                if os.path.exists(file_path_figure_point):
                    os.remove(file_path_figure_point)

            if os.path.exists(file_path_results_point) and os.path.exists(file_path_auxiliary_point):

                dframe_results = self.get_obj_results(
                    file_path_results_point,
                    file_format=format_results,
                    time_fields=None,
                    file_fields=None,
                    registry_fields=fields_registry,
                    point_name=point_name, point_tag=point_tag,
                    point_longitude=point_longitude, point_latitude=point_latitude,
                )

                dframe_auxiliary = self.get_obj_metrics(
                    file_path_auxiliary_point,
                    file_format=format_auxiliary,
                    time_fields=None,
                    file_fields=None,
                    registry_fields=fields_registry
                )

                self.plot_obj_datasets(
                    file_path_figure_point,
                    dframe_results,
                    dframe_auxiliary
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

        format_results = self.format_results
        reset_model_figure = self.reset_model_figure

        time_select = self.maps_time_select_figure

        # merge dataset info start
        log_stream.info(' -----> Merge model results ... ')
        map_collections = []
        for fields_registry in data_registry.to_dict(orient="records"):

            # get point information
            point_name, point_tag = fields_registry['name'], fields_registry['tag']
            point_longitude, point_latitude = fields_registry['longitude'], fields_registry['latitude']

            log_stream.info(
                ' ------> Point -- (1) Name: "' + point_name + '" :: (2) Tag: "' + point_tag + '" ... ')

            # method to fill the filename(s)
            if format_results == 'csv':
                point_tag_results = point_tag
                file_path_results_point = self.__define_file_string(
                    file_path_results_tmpl,
                    extended_info={'point_name': point_tag_results}
                )

            elif format_results == 'netcdf':

                point_tag_results = 'collections'
                file_path_results_point = self.__define_file_string(
                    file_path_results_tmpl,
                    extended_info={'point_name': point_tag_results}
                )

            else:
                log_stream.error(' ===> Format results "' + format_results + '" is not supported')
                raise NotImplemented('Case not implemented yet')

            # check file results availability
            if os.path.exists(file_path_results_point):

                # read file results
                dframe_results = self.get_obj_results(
                    file_path_results_point,
                    file_format=format_results,
                    time_fields=None, file_fields=None, registry_fields=fields_registry,
                    point_name=point_name, point_tag=point_tag,
                    point_longitude=point_longitude, point_latitude=point_latitude,
                )

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
        log_stream.info(' -----> Merge model results ... DONE')

        # plot dataset info start
        log_stream.info(' -----> Plot model results ... ')

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
            log_stream.info(' -----> Plot model results ... DONE')

        else:

            # plot dataset info end
            log_stream.info(' -----> Plot model results ... SKIPPED. Data not available')

        # method info end
        log_stream.info(' ----> View model maps ... DONE')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
