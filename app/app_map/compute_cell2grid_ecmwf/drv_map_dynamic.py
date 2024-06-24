"""
Class Features

Name:          drv_map_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230822'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_info_args import logger_name
from lib_utils_io import fill_path_with_tags
from lib_utils_generic import make_folder, reset_folder, pad_list
from lib_utils_datasets import organize_datasets_cell, organize_datasets_points
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_tiff import organize_file_tiff, write_file_tiff
from lib_data_io_nc import read_file_collection, write_file_collection, organize_file_nc, write_file_nc

from lib_data_analysis import (add_data, get_data, filter_data,
                               resample_data_args, resample_data_fx, interpolate_data, mask_data)

# set logger
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver map
class DrvMap:

    # method to initialize class
    def __init__(self, alg_obj_time, alg_obj_static, alg_obj_settings,
                 tag_section_flags='flags', tag_section_template='template',
                 tag_section_info='info', tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_log='log'):

        self.alg_obj_time = alg_obj_time
        self.alg_obj_static = alg_obj_static

        self.alg_flags = alg_obj_settings[tag_section_flags]
        self.alg_template_time = alg_obj_settings[tag_section_template]['time']
        self.alg_template_dset = alg_obj_settings[tag_section_template]['datasets']
        self.alg_info = alg_obj_settings[tag_section_info]
        self.alg_parameters = alg_obj_settings[tag_section_params]

        self.alg_datasets_src = alg_obj_settings[tag_section_datasets]['dynamic']['source']
        self.alg_datasets_anc_points = alg_obj_settings[tag_section_datasets]['dynamic']['ancillary']['points']
        self.alg_datasets_anc_grid = alg_obj_settings[tag_section_datasets]['dynamic']['ancillary']['grid']
        self.alg_datasets_dst = alg_obj_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_obj_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_file_vars_in = 'variables_in'
        self.tag_file_vars_out = 'variables_out'
        self.tag_file_format = 'format'

        self.alg_datasets_name = self.alg_info['datasets']
        self.alg_datasets_product = self.alg_info['product']

        self.reset_datasets_anc_points = self.alg_flags['reset_ancillary_datasets_points']
        self.reset_datasets_anc_grid = self.alg_flags['reset_ancillary_datasets_grid']
        self.reset_datasets_dst = self.alg_flags['reset_destination_datasets']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_src = self.alg_datasets_src[self.tag_folder_name]
        self.file_name_src = self.alg_datasets_src[self.tag_file_name]
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)
        self.variables_in_src = self.alg_datasets_src[self.tag_file_vars_in]
        self.variables_out_src = self.alg_datasets_src[self.tag_file_vars_out]

        self.folder_name_anc_points = self.alg_datasets_anc_points[self.tag_folder_name]
        self.file_name_anc_points = self.alg_datasets_anc_points[self.tag_file_name]
        self.file_path_anc_points = os.path.join(self.folder_name_anc_points, self.file_name_anc_points)

        self.folder_name_anc_grid = self.alg_datasets_anc_grid[self.tag_folder_name]
        self.file_name_anc_grid = self.alg_datasets_anc_grid[self.tag_file_name]
        self.file_path_anc_grid = os.path.join(self.folder_name_anc_grid, self.file_name_anc_grid)

        self.folder_name_dst = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_datasets_dst[self.tag_file_name]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)
        self.variables_in_dst = self.alg_datasets_dst[self.tag_file_vars_in]
        self.variables_out_dst = self.alg_datasets_dst[self.tag_file_vars_out]
        self.file_format_dst = self.alg_datasets_dst[self.tag_file_format]

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.grid_cells_product = self.alg_obj_static['grid_cells_product']
        self.grid_gpis_product = self.alg_obj_static['grid_gpis_product']
        self.grid_reference_product = self.alg_obj_static['grid_reference_product']
        self.grid_reference_domain = self.alg_obj_static['grid_reference_domain']

        self.vars_name = self.alg_parameters['variables']['name']
        self.vars_min_value = self.alg_parameters['variables']['min_value']
        self.vars_max_value = self.alg_parameters['variables']['max_value']
        self.vars_scale_factor = self.alg_parameters['variables']['scale_factor']
        self.vars_no_data = self.alg_parameters['variables']['no_data']

        self.resampling_max_distance = self.alg_parameters['resample']['max_distance']
        self.resampling_grid_resolution = self.alg_parameters['resample']['grid_resolution']
        self.resampling_min_neighbours = self.alg_parameters['resample']['min_neighbours']
        self.resampling_neighbours = self.alg_parameters['resample']['neighbours']
        self.resampling_variables = self.alg_parameters['resample']['apply_to']

        self.filtering_active = self.alg_parameters['filter']['active']
        self.filtering_method = self.alg_parameters['filter']['method']
        self.filtering_mode = self.alg_parameters['filter']['mode']
        self.filtering_args_width = self.alg_parameters['filter']['args']['width']
        self.filtering_args_stdev = self.alg_parameters['filter']['args']['st_dev']

        self.interpolating_active = self.alg_parameters['interpolate']['active']
        self.interpolating_method = self.alg_parameters['interpolate']['method']
        self.interpolating_max_distance = self.alg_parameters['interpolate']['max_distance']
        self.interpolating_neighbours = self.alg_parameters['interpolate']['neighbours']
        self.interpolating_variables = self.alg_parameters['resample']['apply_to']

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
        alg_logger.info(' ----> Organize dynamic datasets ... ')

        # get file time
        file_time = self.alg_obj_time

        # get grid reference and cells
        grid_cells_product = self.grid_cells_product
        grid_reference_product = self.grid_reference_product
        # get file variables list
        variables_in_src, variables_out_src = self.variables_in_src, self.variables_out_src

        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset
        # get file product
        alg_datasets_name = self.alg_datasets_name

        # get variables info
        vars_name = self.vars_name
        vars_min_value = self.vars_min_value
        vars_max_value = self.vars_max_value
        vars_scale_factor = self.vars_scale_factor
        vars_no_data = self.vars_no_data
        # pad list (using names as list length reference
        vars_min_value, vars_max_value, vars_scale_factor, vars_no_data = pad_list(
            vars_name, [vars_min_value, vars_max_value, vars_scale_factor, vars_no_data])

        # get file path
        file_path_src_raw = self.file_path_src
        file_path_anc_points_raw, file_path_anc_grid_raw = self.file_path_anc_points, self.file_path_anc_grid

        # iterate over grid cells
        file_path_anc_points_list = []
        for grid_cells_name in grid_cells_product:

            # info start get datasets cell
            alg_logger.info(' -----> Get cell "' + str(grid_cells_name) + '" ... ')

            # define source file name
            file_path_src_def = fill_path_with_tags(
                file_path_src_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(grid_cells_name)}, tmpl_tags_dset=alg_template_dset)

            # define ancillary points file name
            file_path_anc_points_def = fill_path_with_tags(
                file_path_anc_points_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(grid_cells_name)}, tmpl_tags_dset=alg_template_dset)

            # clean ancillary points files
            if self.reset_datasets_anc_points:
                if os.path.exists(file_path_anc_points_def):
                    os.remove(file_path_anc_points_def)

            # check file point ancillary availability
            if not os.path.exists(file_path_anc_points_def):

                # organize cell collections and registry
                collections_dframe, registry_dframe = organize_datasets_cell(
                    file_name=file_path_src_def,
                    file_time=file_time,
                    time_frequency_left='D', time_frequency_right='D', time_window_left=5, time_window_right=5,
                    cell_name=grid_cells_name,
                    list_variable_data_in=variables_in_src['datasets'],
                    list_variable_data_out=variables_out_src['datasets'],
                    list_variable_registry_in=variables_in_src['registry'],
                    list_variable_registry_out=variables_out_src['registry'],
                )

                # organize datasets points
                if collections_dframe is not None:

                    # organize cell obj
                    points_dframe = organize_datasets_points(
                        file_time, collections_dframe, registry_dframe)

                    # save cell obj in pickle format
                    # folder_name_anc_points, file_name_anc_points = os.path.split(file_path_anc_points)
                    # make_folder(folder_name_anc_points)
                    # write_file_obj(file_path_anc_points, cell_obj_dict)

                    # save cell obj in netcdf format
                    folder_name_anc_points, file_name_anc_points = os.path.split(file_path_anc_points_def)
                    make_folder(folder_name_anc_points)
                    write_file_collection(file_path_anc_points_def, points_dframe, file_tag_location='location')

                    # info end get datasets cell
                    alg_logger.info(' -----> Get cell "' + str(grid_cells_name) +
                                    '" ... DONE')

                else:
                    # info end get datasets cell
                    alg_logger.info(' -----> Get cell "' + str(grid_cells_name) +
                                    '" ... SKIPPED. Datasets is not available')

            else:

                # info end get datasets cell
                alg_logger.info(' -----> Get cell "' + str(grid_cells_name) +
                                '" ... SKIPPED. Data previously saved')

            # store cell filename(s)
            file_path_anc_points_list.append(file_path_anc_points_def)

        # info start merge cell group
        alg_logger.info(' -----> Merge cell ... ')

        # iterate over cell file(s)
        collections_dframe = None
        for file_path_anc_points_step in file_path_anc_points_list:

            # check file ancillary points availability
            if os.path.exists(file_path_anc_points_step):
                point_dframe_step = read_file_collection(
                    file_path_anc_points_step,
                    file_obj_grid=grid_reference_product,
                    var_dest_list=variables_out_src['datasets'],
                    var_registry_list=variables_out_src['registry'],
                    var_extra_list=['time_start', 'time_end', 'time_reference', 'fx', 'cell'],
                    variable_type=alg_datasets_name,
                    variable_min_value=vars_min_value, variable_max_value=vars_max_value,
                    variable_no_data=vars_no_data, variable_scale_factor=vars_scale_factor,
                    variable_land_filter=True, variable_committed_filter=False)
            else:
                alg_logger.warning(' ===> File ' + file_path_anc_points_step + ' not found')
                point_dframe_step = None

            # merge cell dataset(s)
            if point_dframe_step is not None:
                if collections_dframe is None:
                    collections_dframe = point_dframe_step
                else:
                    collections_dframe = collections_dframe.append(point_dframe_step)
        # info end merge cell group
        alg_logger.info(' -----> Merge cell ... DONE')

        # info end method
        alg_logger.info(' ----> Organize dynamic datasets ... DONE')

        return collections_dframe

    # method to analyze data
    def analyze_data(self, collections_dframe=None):

        # info start method
        alg_logger.info(' ----> Analyze dynamic datasets ... ')

        # get variables info
        vars_name = self.vars_name
        vars_min_value = self.vars_min_value
        vars_max_value = self.vars_max_value
        vars_scale_factor = self.vars_scale_factor
        # pad list (using names as list length reference
        vars_min_value, vars_max_value, vars_scale_factor = pad_list(
            vars_name, [vars_min_value, vars_max_value, vars_scale_factor])

        # get grid reference and cells
        grid_reference_domain = self.grid_reference_domain
        # get file time
        file_time = self.alg_obj_time
        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset

        # get file path
        file_path_anc_grid_raw, file_path_dst_raw = self.file_path_anc_grid, self.file_path_dst

        # define ancillary points file name
        file_path_anc_grid_def = fill_path_with_tags(
            file_path_anc_grid_raw,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)
        # define destination file name
        file_path_dst_def = fill_path_with_tags(
            file_path_dst_raw,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)

        # clean ancillary and destination datasets (if ancillary flag is activated
        if self.reset_datasets_anc_points or self.reset_datasets_anc_grid:
            if os.path.exists(file_path_anc_grid_def):
                os.remove(file_path_anc_grid_def)
            if os.path.exists(file_path_dst_def):
                os.remove(file_path_dst_def)

        # check file ancillary points availability
        if collections_dframe is not None:

            # check file ancillary grid availability
            if not os.path.exists(file_path_anc_grid_def):

                # info start add data
                alg_logger.info(' -----> Add data ... ')
                # method to add data
                collections_dframe = add_data(collections_dframe,
                                              time_reference='time_reference', time_data='time')
                # info end add data
                alg_logger.info(' -----> Add data ... DONE')

                # info start compute data
                alg_logger.info(' -----> Compute data ... ')

                # iterate over variable(s)
                variable_collections = {}
                for var_name, var_min_value, var_max_value, var_scale_factor in zip(
                        vars_name, vars_min_value, vars_max_value, vars_scale_factor):

                    # info start variable
                    alg_logger.info(' ------> Variable "' + var_name + '" ... ')

                    # debug
                    # var_name = 'sm_values'

                    # check limits and scale factor
                    if var_scale_factor is not None:
                        if var_min_value is not None:
                            var_min_value = var_min_value * var_scale_factor
                        if var_max_value is not None:
                            var_max_value = var_max_value * var_scale_factor

                    # check variable availability
                    if var_name in collections_dframe.columns:

                        # method to get data
                        dframe_variable_in = get_data(collections_dframe, variable_name=var_name)
                        # method to filter data
                        dframe_variable_filtered = filter_data(
                            dframe_variable_in,
                            variable_name=var_name, min_value=var_min_value, max_value=var_max_value)

                        # method to create resampling argument(s)
                        resample_kwargs = resample_data_args(
                            grid_reference_domain,
                            resampling_max_distance=self.resampling_max_distance,
                            resampling_grid_resolution=self.resampling_grid_resolution,
                            resampling_min_neighbours=self.resampling_min_neighbours,
                            resampling_neighbours=self.resampling_neighbours,
                            resampling_kernel_active=self.filtering_active,
                            resampling_kernel_method=self.filtering_method,
                            resampling_kernel_mode=self.filtering_mode,
                            resampling_kernel_width=self.filtering_args_width,
                            resampling_kernel_stddev=self.filtering_args_stdev)

                        # method to apply resampling data 1D to data 2D
                        da_variable_resampled = resample_data_fx(
                            dframe_variable_filtered, **resample_kwargs,
                            var_name_data=var_name,
                            var_name_geo_x='longitude', var_name_geo_y='latitude',
                            coord_name_x='longitude', coord_name_y='latitude',
                            dim_name_x='longitude', dim_name_y='latitude'
                        )

                        # method to apply interpolating data 2d over reference (in case of dims mismatch)
                        da_variable_interpolated = interpolate_data(
                            da_variable_resampled, grid_reference_domain, var_name_data=var_name,
                            interpolating_active=self.interpolating_active,
                            interpolating_method=self.interpolating_method,
                            interpolating_max_distance=self.interpolating_max_distance,
                            interpolating_neighbours=self.interpolating_neighbours
                        )

                        # method to mask data
                        da_variable_masked = mask_data(
                            da_variable_interpolated, grid_reference_domain, var_name_data=var_name,
                            mask_domain=True)

                        # save variable in a common collection
                        variable_collections[var_name] = da_variable_masked

                        ''' debug
                        plt.figure()
                        plt.imshow(da_variable_masked.values)
                        plt.colorbar()
                        plt.clim(0, 100)
                        plt.show()
                        '''

                        # info end variable
                        alg_logger.info(' ------> Variable "' + var_name + '" ... DONE')

                    else:
                        # info end variable
                        alg_logger.info(' ------> Variable "' + var_name + '" ... SKIPPED. Data not available')

                # info end compute data
                alg_logger.info(' -----> Compute data ... DONE')

                # save variable collection in pickle format
                folder_name_anc_grid, file_name_anc_grid = os.path.split(file_path_anc_grid_def)
                make_folder(folder_name_anc_grid)
                write_file_obj(file_path_anc_grid_def, variable_collections)

                # info end method
                alg_logger.info(' ----> Analyze dynamic datasets ... DONE')

            else:
                # info end method
                alg_logger.info(' ----> Analyze dynamic datasets ... SKIPPED. Data grid previously saved')

        else:
            # info end method
            alg_logger.info(' ----> Analyze dynamic datasets ... FAILED. Data points not available')

    # method to save data
    def dump_data(self):

        # info start method
        alg_logger.info(' ----> Dump dynamic datasets ... ')

        # get file time
        file_time = self.alg_obj_time

        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset
        # get file variables
        variables_in_dst, variables_out_dst = self.variables_in_dst, self.variables_out_dst

        # get file path
        file_path_anc_grid_raw, file_path_dst_raw = self.file_path_anc_grid, self.file_path_dst

        # define ancillary points file name
        file_path_anc_grid_def = fill_path_with_tags(
            file_path_anc_grid_raw,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)
        # define destination file name
        file_path_dst_def = fill_path_with_tags(
            file_path_dst_raw,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)

        # clean destination datasets
        if self.reset_datasets_anc_grid or self.reset_datasets_dst:
            if os.path.exists(file_path_dst_def):
                os.remove(file_path_dst_def)

        # check file ancillary grid availability
        if os.path.exists(file_path_anc_grid_def):
            # check file destination availability
            if not os.path.exists(file_path_dst_def):

                # method to get data obj
                variable_collection = read_file_obj(file_path_anc_grid_def)

                # check file output format
                if self.file_format_dst == 'netcdf':

                    # method to organize dataset
                    variable_dset = organize_file_nc(variable_collection)

                    # method to write dataset
                    folder_name_dst, file_name_dst = os.path.split(file_path_dst_def)
                    make_folder(folder_name_dst)
                    write_file_nc(file_path_dst_def, variable_dset)

                elif self.file_format_dst == 'tiff':

                    # method to organize dataset
                    variable_dset, attrs_dset = organize_file_tiff(
                        variable_collection, obj_time=file_time,
                        obj_variable_in=variables_in_dst['datasets'], obj_variable_out=variables_out_dst['datasets'])

                    # method to write dataset
                    folder_name_dst, file_name_dst = os.path.split(file_path_dst_def)
                    make_folder(folder_name_dst)
                    write_file_tiff(file_path_dst_def, variable_dset, **attrs_dset)

                else:
                    # exit format not supported
                    alg_logger.error(' ===> File format "' + self.file_format_dst + '" not supported')
                    raise NotImplementedError('Case not implemented yet')

                # info end method
                alg_logger.info(' ----> Dump dynamic datasets ... DONE')

            else:
                # info end method
                alg_logger.info(' ----> Dump dynamic datasets ... SKIPPED. Data previously saved')
        else:
            # info end method
            alg_logger.info(' ----> Dump dynamic datasets ... FAILED. Data grid not available')

# -------------------------------------------------------------------------------------
