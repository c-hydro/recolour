"""
Class Features

Name:          drv_data_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230822'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np

from copy import deepcopy

from lib_info_args import logger_name
from lib_utils_io import fill_path_with_tags

from lib_data_io_nc import organize_file_nc, write_file_nc
from lib_data_io_tiff import organize_file_tiff, write_file_tiff

from lib_utils_generic import reset_folder, pad_list, combine_vars2dict
from lib_cell_datasets import get_datasets_cell, apply_scaling_datasets_cell

from lib_data_analysis import select_data, resample_data, adapt_data, mask_data, filter_data
from lib_data_io_pickle import write_file_obj, read_file_obj

# set logger
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_obj_time, alg_obj_static, alg_obj_settings):

        self.alg_obj_time = alg_obj_time
        self.alg_obj_static = alg_obj_static

        self.alg_flags = alg_obj_settings['flags']
        self.alg_template_time = alg_obj_settings['template']['time']
        self.alg_template_dset = alg_obj_settings['template']['datasets']
        self.alg_info = alg_obj_settings['info']
        self.alg_cells = alg_obj_settings['cells']
        self.alg_parameters = alg_obj_settings['parameters']
        self.alg_log = alg_obj_settings['log']

        self.alg_datasets_src_cell_ref_nrt = alg_obj_settings['datasets']['dynamic'][
            'source']['reference']['cell_datasets_data_nrt']
        self.alg_datasets_src_cell_ref_dr = alg_obj_settings['datasets']['dynamic'][
            'source']['reference']['cell_datasets_data_dr']
        self.alg_datasets_src_cell_other_dr = alg_obj_settings['datasets']['dynamic'][
            'source']['other']['cell_datasets_data_dr']
        self.alg_datasets_anc_cell = alg_obj_settings['datasets']['dynamic']['ancillary']['cell_datasets']
        self.alg_datasets_anc_grid = alg_obj_settings['datasets']['dynamic']['ancillary']['grid_datasets']
        self.alg_datasets_dst_grid = alg_obj_settings['datasets']['dynamic']['destination']['grid_datasets']

        self.alg_datasets_reference = self.alg_info['reference']
        self.alg_datasets_other = self.alg_info['other']

        self.tag_active = 'active'
        self.tag_folder_name, self.tag_file_name = 'folder_name', 'file_name'
        self.tag_file_vars_in, self.tag_file_vars_out = 'variables_in', 'variables_out'
        self.tag_file_format = 'format'

        self.reset_anc_cell = self.alg_flags['reset_ancillary_cell']
        self.reset_anc_grid = self.alg_flags['reset_ancillary_grid']
        self.reset_dst_grid = self.alg_flags['reset_destination_grid']
        self.reset_logs = self.alg_flags['reset_logs']

        # organize source datasets
        self.folder_name_src_cell_ref_nrt = self.alg_datasets_src_cell_ref_nrt[self.tag_folder_name]
        self.file_name_src_cell_ref_nrt = self.alg_datasets_src_cell_ref_nrt[self.tag_file_name]
        self.file_path_src_cell_ref_nrt = os.path.join(self.folder_name_src_cell_ref_nrt, self.file_name_src_cell_ref_nrt)
        self.variables_in_src_cell_ref_nrt = self.alg_datasets_src_cell_ref_nrt[self.tag_file_vars_in]
        self.variables_out_src_cell_ref_nrt = self.alg_datasets_src_cell_ref_nrt[self.tag_file_vars_out]

        self.folder_name_src_cell_ref_dr = self.alg_datasets_src_cell_ref_dr[self.tag_folder_name]
        self.file_name_src_cell_ref_dr = self.alg_datasets_src_cell_ref_dr[self.tag_file_name]
        self.file_path_src_cell_ref_dr = os.path.join(self.folder_name_src_cell_ref_dr, self.file_name_src_cell_ref_dr)
        self.variables_in_src_cell_ref_dr = self.alg_datasets_src_cell_ref_dr[self.tag_file_vars_in]
        self.variables_out_src_cell_ref_dr = self.alg_datasets_src_cell_ref_dr[self.tag_file_vars_out]

        self.folder_name_src_cell_other_dr = self.alg_datasets_src_cell_other_dr[self.tag_folder_name]
        self.file_name_src_cell_other_dr = self.alg_datasets_src_cell_other_dr[self.tag_file_name]
        self.file_path_src_cell_other_dr = os.path.join(self.folder_name_src_cell_other_dr, self.file_name_src_cell_other_dr)
        self.variables_in_src_cell_other_dr = self.alg_datasets_src_cell_other_dr[self.tag_file_vars_in]
        self.variables_out_src_cell_other_dr = self.alg_datasets_src_cell_other_dr[self.tag_file_vars_out]
        self.variables_in_dst_grid = self.alg_datasets_dst_grid[self.tag_file_vars_out]
        self.variables_out_dst_grid = self.alg_datasets_dst_grid[self.tag_file_vars_out]

        # organize ancillary datasets
        self.folder_name_anc_cell = self.alg_datasets_anc_cell[self.tag_folder_name]
        self.file_name_anc_cell = self.alg_datasets_anc_cell[self.tag_file_name]
        self.file_path_anc_cell = os.path.join(self.folder_name_anc_cell, self.file_name_anc_cell)

        self.folder_name_anc_grid = self.alg_datasets_anc_grid[self.tag_folder_name]
        self.file_name_anc_grid = self.alg_datasets_anc_grid[self.tag_file_name]
        self.file_path_anc_grid = os.path.join(self.folder_name_anc_grid, self.file_name_anc_grid)

        # organize destination datasets
        self.folder_name_dst_grid = self.alg_datasets_dst_grid[self.tag_folder_name]
        self.file_name_dst_grid = self.alg_datasets_dst_grid[self.tag_file_name]
        self.file_path_dst_grid = os.path.join(self.folder_name_dst_grid, self.file_name_dst_grid)
        self.variables_in_dst_grid = self.alg_datasets_dst_grid[self.tag_file_vars_in]
        self.variables_out_dst_grid = self.alg_datasets_dst_grid[self.tag_file_vars_out]
        self.format_dst_grid = self.alg_datasets_dst_grid[self.tag_file_format]

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.grid_src_ref_nrt = self.alg_obj_static['grid_src_ref_nrt']
        self.grid_src_ref_dr = self.alg_obj_static['grid_src_ref_dr']
        self.grid_src_other_dr = self.alg_obj_static['grid_src_other_dr']
        self.gpis_src_ref_nrt = self.alg_obj_static['gpis_src_ref_nrt']
        self.cell_src_ref_nrt = self.alg_obj_static['cell_src_ref_nrt']
        self.gpis_src_ref_dr = self.alg_obj_static['gpis_src_ref_dr']
        self.cell_src_ref_dr = self.alg_obj_static['cell_src_ref_dr']
        self.gpis_src_other_dr = self.alg_obj_static['gpis_src_other_dr']
        self.cell_src_other_dr = self.alg_obj_static['cell_src_other_dr']

        grid_dst = self.alg_obj_static['grid_dst']
        self.grid_geo_values_dst = grid_dst.values
        self.grid_geo_x_dst, self.grid_geo_y_dst = np.meshgrid(grid_dst['longitude'], grid_dst['latitude'])
        self.grid_geo_transform, self.grid_geo_epsg = grid_dst.attrs['transform'], grid_dst.attrs['crs']

        self.list_grid_src_common = list(
            set(self.cell_src_ref_nrt).intersection(self.cell_src_ref_dr).intersection(self.cell_src_other_dr))

        self.vars_name = self.alg_parameters['variables']['name']
        self.vars_min_value = self.alg_parameters['variables']['min_value']
        self.vars_max_value = self.alg_parameters['variables']['max_value']
        self.vars_scale_factor = self.alg_parameters['variables']['scale_factor']
        self.vars_no_data = self.alg_parameters['variables']['no_data']

        self.alg_fx_scale = self.alg_parameters['fx']['scale']
        self.alg_fx_select = self.alg_parameters['fx']['select']
        self.alg_fx_filter = self.alg_parameters['fx']['filter']
        self.alg_fx_resample = self.alg_parameters['fx']['resample']
        self.alg_fx_mask = self.alg_parameters['fx']['mask']

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

    # method to dump data
    def dump_data(self):

        # info start method
        alg_logger.info(' ----> Dump dynamic datasets ... ')

        # get file time
        file_time = self.alg_obj_time

        # get destination grid
        grid_geo_values_dst = self.grid_geo_values_dst
        grid_geo_x_dst, grid_geo_y_dst = self.grid_geo_x_dst, self.grid_geo_y_dst
        grid_geo_transform_dst, grid_geo_proj_dst = self.grid_geo_transform, self.grid_geo_epsg

        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset

        # get variable(s) source and destination
        variables_in_dst_grid, variables_out_dst_grid = self.variables_in_dst_grid, self.variables_out_dst_grid

        # get file path
        file_path_tmpl_anc_grid, file_path_tmpl_dst_grid = self.file_path_anc_grid, self.file_path_dst_grid

        # define ancillary grid file name
        file_path_step_anc_grid = fill_path_with_tags(
            file_path_tmpl_anc_grid,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)
        # define destination grid file name
        file_path_step_dst_grid = fill_path_with_tags(
            file_path_tmpl_dst_grid,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)

        # clean ancillary and destination file name
        if self.reset_dst_grid:
            if os.path.exists(file_path_step_dst_grid):
                os.remove(file_path_step_dst_grid)

        # check destination grid file availability
        if not os.path.exists(file_path_step_dst_grid):

            # check ancillary grid file availability
            if os.path.exists(file_path_step_anc_grid):

                # load ancillary obj
                variable_collections = read_file_obj(file_path_step_anc_grid)

                # check destination format
                if self.format_dst_grid == 'netcdf':

                    # method to organize netcdf dataset
                    variable_dset = organize_file_nc(
                        obj_data=variable_collections, obj_time=file_time,
                        obj_variable_in=variables_in_dst_grid, obj_variable_out=variables_out_dst_grid,
                        obj_geo_x=grid_geo_x_dst, obj_geo_y=grid_geo_y_dst)

                    # method to write netcdf dataset
                    folder_name_dst_step, file_name_dst_step = os.path.split(file_path_step_dst_grid)
                    os.makedirs(folder_name_dst_step, exist_ok=True)
                    write_file_nc(file_path_step_dst_grid, variable_dset)

                elif (self.format_dst_grid == 'tiff') or (self.format_dst_grid == 'tif'):

                    # method to organize tiff dataset
                    variable_data, variable_attrs = organize_file_tiff(
                        variable_collections, obj_time=file_time,
                        obj_variable_in=variables_in_dst_grid, obj_variable_out=variables_out_dst_grid,
                        var_attr_description='description', var_attr_time='time',
                        var_name_geo_x='longitude', var_name_geo_y='latitude',
                        obj_transform=grid_geo_transform_dst, obj_proj=grid_geo_proj_dst)

                    # method to write tiff dataset
                    folder_name_dst_step, file_name_dst_step = os.path.split(file_path_step_dst_grid)
                    os.makedirs(folder_name_dst_step, exist_ok=True)
                    write_file_tiff(
                        file_name=file_path_step_dst_grid, file_data=variable_data,
                        **variable_attrs)

                else:
                    alg_logger.error(' ===> Destination grid format "' + self.format_dst_grid + '" is not supported')
                    raise NotImplemented('Case not implemented yet')

                # info end method
                alg_logger.info(' ----> Dump dynamic datasets ... DONE')

            else:
                # info end method
                alg_logger.info(' ----> Dump dynamic datasets ... SKIPPED. Ancillary grid data not available')

        else:
            # info end method
            alg_logger.info(' ----> Dump dynamic datasets ... SKIPPED. Destination grid data previously saved')

    # method to analyze data
    def analyze_data(self):

        # info start method
        alg_logger.info(' ----> Analyze dynamic datasets ... ')

        # get file time
        file_time = self.alg_obj_time

        # get grid cells reference
        list_grid_src_common = self.list_grid_src_common
        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset

        # get destination grid
        grid_geo_values_dst = self.grid_geo_values_dst
        grid_geo_x_dst, grid_geo_y_dst = self.grid_geo_x_dst, self.grid_geo_y_dst

        # set max distance
        settings_select_data = self.alg_fx_select
        settings_resample_data = self.alg_fx_resample
        settings_filter_data = self.alg_fx_filter
        settings_mask_data = self.alg_fx_mask

        # get file path
        file_path_tmpl_anc_cell, file_path_tmpl_anc_grid = self.file_path_anc_cell, self.file_path_anc_grid
        file_path_tmpl_dst_grid = self.file_path_dst_grid

        # define ancillary grid file name
        file_path_step_anc_grid = fill_path_with_tags(
            file_path_tmpl_anc_grid,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)
        # define destination grid file name
        file_path_step_dst_grid = fill_path_with_tags(
            file_path_tmpl_dst_grid,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)

        # clean ancillary and destination file name
        if self.reset_anc_grid:
            if os.path.exists(file_path_step_anc_grid):
                os.remove(file_path_step_anc_grid)
            if os.path.exists(file_path_step_dst_grid):
                os.remove(file_path_step_dst_grid)
        if self.reset_dst_grid:
            if os.path.exists(file_path_step_dst_grid):
                os.remove(file_path_step_dst_grid)

        # check destination grid file availability
        file_path_list_anc_grid = []
        if not os.path.exists(file_path_step_dst_grid):

            # check ancillary grid file availability
            if not os.path.exists(file_path_step_anc_grid):

                # info start merge datasets
                alg_logger.info(' -----> (1) Merge data ... ')

                # iterate over grid cells
                obj_collections_merge, obj_registry_merge = None, None
                for cell_name in list_grid_src_common:

                    # define ancillary cell file name
                    file_path_step_anc_cell = fill_path_with_tags(
                        file_path_tmpl_anc_cell,
                        tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                        tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

                    # check ancillary file cell availability
                    if os.path.exists(file_path_step_anc_cell):

                        # load ancillary obj
                        obj_cell = read_file_obj(file_path_step_anc_cell)
                        # get collections and registry objects
                        obj_collections_cell, obj_registry_cell = obj_cell['collections'], obj_cell['registry']

                        # check collections and registry objects availability
                        if obj_collections_cell is not None and obj_registry_cell is not None:

                            # merge collection cell dataset(s)
                            if obj_collections_merge is None:
                                obj_collections_merge = deepcopy(obj_collections_cell)
                            else:
                                obj_collections_merge = obj_collections_merge.append(obj_collections_cell)

                            # merge registry cell dataset(s)
                            if obj_registry_merge is None:
                                obj_registry_merge = deepcopy(obj_registry_cell)
                            else:
                                obj_registry_merge = obj_registry_merge.append(obj_registry_cell)

                        else:
                            # message warning for NoneType datasets
                            alg_logger.warning(' ===> Cell datasets "' + str(cell_name) + '" is defined by NoneType')

                    else:
                        # message warning for cell datasets ancillary file not found
                        alg_logger.warning(' ===> Cell datasets "' + str(cell_name) + '" filename not found')

                # info end merge datasets
                alg_logger.info(' -----> (1) Merge data ... DONE')

                # info start select datasets
                alg_logger.info(' -----> (2) Select data ... ')
                obj_collections_filter = select_data(
                    obj_collections_merge, **settings_select_data)
                # info end select datasets
                alg_logger.info(' -----> (2) Select data ... DONE')

                # info start adapt datasets
                alg_logger.info(' -----> (3) Adapt data ... ')
                obj_collections_adapt, obj_geo_x_adapt, obj_geo_y_adapt, obj_time_adapt = adapt_data(
                    obj_collections_filter,
                    var_name_geo_x='lon_scaled', var_name_geo_y='lat_scaled', var_name_time='time_scaled')
                # info end adapt datasets
                alg_logger.info(' -----> (3) Adapt data ... DONE')

                # info start resample cell datasets
                alg_logger.info(' -----> (4) Resample data ... ')
                # resample datasets
                obj_collections_resample = resample_data(
                    obj_collections_adapt, obj_geo_x_adapt, obj_geo_y_adapt,
                    grid_geo_x_dst, grid_geo_y_dst, **settings_resample_data)
                # info end resample cell datasets
                alg_logger.info(' -----> (4) Resample data ... DONE')

                # info start filter cell datasets
                alg_logger.info(' -----> (5) Filter data ... ')
                obj_collections_filter = filter_data(obj_collections_resample, **settings_filter_data)
                # info end filter cell datasets
                alg_logger.info(' -----> (5) Filter data ... DONE')

                # info start mask datasets
                alg_logger.info(' -----> (6) Mask data ... ')
                # mask datasets
                obj_collections_mask = mask_data(obj_collections_filter,
                                                 var_values_ref=grid_geo_values_dst, **settings_mask_data)
                # info end mask datasets
                alg_logger.info(' -----> (6) Mask data ... DONE')

                # info start save datasets
                alg_logger.info(' -----> (7) Save data ... ')
                # save data in pickle format
                folder_name_anc, file_name_anc = os.path.split(file_path_step_anc_grid)
                os.makedirs(folder_name_anc, exist_ok=True)
                write_file_obj(file_path_step_anc_grid, obj_collections_mask)
                # info end save datasets
                alg_logger.info(' -----> (7) Save data ... DONE')

                # info end method
                file_path_list_anc_grid.append(file_path_step_anc_grid)
                alg_logger.info(' ----> Analyze dynamic datasets ... DONE')

            else:
                # info end method
                file_path_list_anc_grid.append(file_path_step_anc_grid)
                alg_logger.info(' ----> Analyze dynamic datasets ... SKIPPED. Ancillary grid data previously saved')

        else:
            # info end method
            file_path_list_anc_grid = []
            alg_logger.info(' ----> Analyze dynamic datasets ... SKIPPED. Destination grid data previously saved')

        return file_path_list_anc_grid

    # method to compute data
    def compute_data(self):

        # info start method
        alg_logger.info(' ----> Compute dynamic datasets ... ')

        # get file time
        file_time = self.alg_obj_time

        # get grid cells, obj and reference
        list_grid_src_common = self.list_grid_src_common
        grid_src_ref_nrt, grid_src_ref_dr = self.grid_src_ref_nrt, self.grid_src_ref_dr
        grid_src_other_dr = self.grid_src_other_dr

        # get variable(s) source and destination
        variables_in_src_cell_ref_nrt = self.variables_in_src_cell_ref_nrt
        variables_out_src_cell_ref_nrt = self.variables_out_src_cell_ref_nrt
        variables_in_src_cell_ref_dr = self.variables_in_src_cell_ref_dr
        variables_out_src_cell_ref_dr = self.variables_out_src_cell_ref_dr
        variables_in_src_cell_other_dr = self.variables_in_src_cell_other_dr
        variables_out_src_cell_other_dr = self.variables_out_src_cell_other_dr

        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset

        # get variables info
        vars_name = self.vars_name
        vars_min_value = self.vars_min_value
        vars_max_value = self.vars_max_value
        vars_scale_factor = self.vars_scale_factor
        vars_no_data = self.vars_no_data
        # pad list (using names as list length reference
        vars_min_value, vars_max_value, vars_scale_factor, vars_no_data = pad_list(
            vars_name, [vars_min_value, vars_max_value, vars_scale_factor, vars_no_data])
        # combine variables in dictionary
        vars_parameters = combine_vars2dict(
            var_pivot='name',
            var_dict={'name': vars_name, 'min_value': vars_min_value,
                      'max_value': vars_max_value, 'scale_factor': vars_scale_factor, 'no_data': vars_no_data})

        # set fx settings
        fx_scale_method = self.alg_fx_scale['methods']
        fx_scale_vars_src = self.alg_fx_scale['variables']['source']
        fx_scale_vars_dst = self.alg_fx_scale['variables']['destination']
        fx_scale_max_dist = 25000

        # get file path
        file_path_tmpl_src_cell_ref_nrt = self.file_path_src_cell_ref_nrt
        file_path_tmpl_src_cell_ref_dr = self.file_path_src_cell_ref_dr
        file_path_tmpl_src_cell_other_dr = self.file_path_src_cell_other_dr
        file_path_tmpl_anc_cell, file_path_tmpl_anc_grid = self.file_path_anc_cell, self.file_path_anc_grid
        file_path_tmpl_dst_grid = self.file_path_dst_grid

        # define ancillary grid file name
        file_path_step_anc_grid = fill_path_with_tags(
            file_path_tmpl_anc_grid,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)
        # define destination grid file name
        file_path_step_dst_grid = fill_path_with_tags(
            file_path_tmpl_dst_grid,
            tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
            tmpl_dset_obj=None, tmpl_tags_dset=alg_template_dset)

        # clean ancillary and destination file name
        if self.reset_anc_cell:
            if os.path.exists(file_path_step_anc_grid):
                os.remove(file_path_step_anc_grid)
            if os.path.exists(file_path_step_dst_grid):
                os.remove(file_path_step_dst_grid)
        if self.reset_dst_grid:
            if os.path.exists(file_path_step_dst_grid):
                os.remove(file_path_step_dst_grid)

        # check file destination availability
        if not os.path.exists(file_path_step_dst_grid):

            # iterate over grid cells
            file_path_list_anc_cell = []
            for cell_name in list_grid_src_common:

                # define source cell reference nrt file name
                file_path_step_src_cell_ref_nrt = fill_path_with_tags(
                    file_path_tmpl_src_cell_ref_nrt,
                    tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

                # define source cell reference dr file name
                file_path_step_src_cell_ref_dr = fill_path_with_tags(
                    file_path_tmpl_src_cell_ref_dr,
                    tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

                # define source cell other dr file name
                file_path_step_src_cell_other_dr = fill_path_with_tags(
                    file_path_tmpl_src_cell_other_dr,
                    tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

                # define ancillary file name
                file_path_step_anc_cell = fill_path_with_tags(
                    file_path_tmpl_anc_cell,
                    tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                    tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

                # clean ancillary file name (and destination file name if needed)
                if self.reset_anc_cell:
                    if os.path.exists(file_path_step_anc_cell):
                        os.remove(file_path_step_anc_cell)

                # info start method
                alg_logger.info(' -----> Cell "' + str(cell_name) + '" ... ')

                # check file cell ancillary availability
                if not os.path.exists(file_path_step_anc_cell):

                    # info start get source reference nrt datasets
                    alg_logger.info(' ------> Get cell source reference nrt datasets ... ')

                    # check file source reference nrt datasets availability
                    if os.path.exists(file_path_step_src_cell_ref_nrt):

                        # get collections and registry objects
                        collections_dframe_cell_src_ref_nrt, registry_dframe_cell_src_ref_nrt = (
                            get_datasets_cell(
                                file_name=file_path_step_src_cell_ref_nrt,
                                cell_name=cell_name,
                                list_variable_data_in=variables_in_src_cell_ref_nrt['datasets'],
                                list_variable_data_out=variables_out_src_cell_ref_nrt['datasets'],
                                list_variable_registry_in=variables_in_src_cell_ref_nrt['registry'],
                                list_variable_registry_out=variables_out_src_cell_ref_nrt['registry'],
                                index_name_data_in='time', index_name_data_out='time_ref_nrt',
                                index_name_registry_out='location_id_ref_nrt',
                                index_name_cell_out='cell_ref_nrt'
                        ))

                        # info end get source reference nrt datasets
                        alg_logger.info(' ------> Get cell source reference nrt datasets ... DONE')

                    else:
                        # info end get source reference nrt datasets
                        collections_dframe_cell_src_ref_nrt, registry_dframe_cell_src_ref_nrt = None, None
                        alg_logger.warning(' ===> File not found "' + str(file_path_step_src_cell_ref_nrt) + '"')
                        alg_logger.info(' ------> Get cell source reference nrt datasets ... FAILED.')

                    # info start get source reference dr datasets
                    alg_logger.info(' ------> Get cell source reference dr datasets ... ')

                    # check file source reference dr datasets availability
                    if os.path.exists(file_path_step_src_cell_ref_dr):

                        # get collections and registry objects
                        collections_dframe_cell_src_ref_dr, registry_dframe_cell_src_ref_dr = (
                            get_datasets_cell(
                                file_name=file_path_step_src_cell_ref_dr,
                                cell_name=cell_name,
                                list_variable_data_in=variables_in_src_cell_ref_dr['datasets'],
                                list_variable_data_out=variables_out_src_cell_ref_dr['datasets'],
                                list_variable_registry_in=variables_in_src_cell_ref_dr['registry'],
                                list_variable_registry_out=variables_out_src_cell_ref_dr['registry'],
                                index_name_data_in='time', index_name_data_out='time_ref_dr',
                                index_name_registry_out='location_id_ref_dr',
                                index_name_cell_out='cell_ref_dr'
                        ))

                        # info end get source reference dr datasets
                        alg_logger.info(' ------> Get cell source reference dr datasets ... DONE')

                    else:
                        # info end get source reference dr datasets
                        alg_logger.info(' ------> Get cell source reference dr datasets ... FAILED')
                        alg_logger.error(' ===> File not found "' + str(file_path_step_src_cell_ref_dr) + '"')
                        raise IOError('File not found "' + str(file_path_step_src_cell_ref_dr) + '"')

                    # info start get source other dr datasets
                    alg_logger.info(' ------> Get cell source other dr datasets ... ')

                    # check file source other dr datasets availability
                    if os.path.exists(file_path_step_src_cell_other_dr):

                        # get collections and registry objects
                        collections_dframe_cell_src_other_dr, registry_dframe_cell_src_other_dr = (
                            get_datasets_cell(
                                file_name=file_path_step_src_cell_other_dr,
                                cell_name=cell_name,
                                list_variable_data_in=variables_in_src_cell_other_dr['datasets'],
                                list_variable_data_out=variables_out_src_cell_other_dr['datasets'],
                                list_variable_registry_in=variables_in_src_cell_other_dr['registry'],
                                list_variable_registry_out=variables_out_src_cell_other_dr['registry'],
                                index_name_data_in='time', index_name_data_out='time_other_dr',
                                index_name_registry_out='location_id_other_dr',
                                index_name_cell_out='cell_other_dr'
                        ))

                        # info end get source other dr datasets
                        alg_logger.info(' ------> Get cell source other dr datasets ... DONE')

                    else:
                        # info end get source other dr datasets
                        alg_logger.info(' ------> Get cell source other dr datasets ... FAILED')
                        alg_logger.error(' ===> File not found "' + str(file_path_step_src_cell_other_dr) + '"')
                        raise IOError('File not found "' + str(file_path_step_src_cell_other_dr) + '"')

                    # info start scale cell datasets
                    alg_logger.info(' ------> Scale cell datasets ... ')

                    # check collections and registry objects availability
                    if collections_dframe_cell_src_ref_nrt is not None:

                        # apply metrics to datasets
                        (collections_dframe_out, registry_dframe_out) = apply_scaling_datasets_cell(
                            file_time,
                            collections_dframe_cell_src_ref_nrt, registry_dframe_cell_src_ref_nrt,
                            collections_dframe_cell_src_ref_dr, registry_dframe_cell_src_ref_dr,
                            collections_dframe_cell_src_other_dr, registry_dframe_cell_src_other_dr,
                            var_params=vars_parameters,
                            grid_obj_ref_nrt=grid_src_ref_nrt,
                            grid_obj_ref_dr=grid_src_ref_dr,
                            grid_obj_other_dr=grid_src_other_dr,
                            fx_var_name_in=fx_scale_vars_src, fx_var_name_out=fx_scale_vars_dst,
                            fx_var_methods=fx_scale_method, max_dist=fx_scale_max_dist,
                        )

                        # info end apply metrics to datasets cell
                        alg_logger.info(' ------> Scale cell datasets ... DONE')

                    else:
                        # info end apply metrics to datasets cell
                        collections_dframe_out, registry_dframe_out = None, None
                        alg_logger.warning(' ===> Data not available for scaling')
                        alg_logger.info(' ------> Scale cell datasets ... FAILED')

                    # info start dump datasets cell
                    alg_logger.info(' ------> Save cell datasets ... ')

                    # check collections and registry objects availability
                    if collections_dframe_out is not None and registry_dframe_out is not None:

                        # save cell obj in pickle format
                        folder_name_anc, file_name_anc = os.path.split(file_path_step_anc_cell)
                        os.makedirs(folder_name_anc, exist_ok=True)

                        # organize obj cell
                        obj_cell = {'collections': collections_dframe_out, 'registry': registry_dframe_out}
                        # dump obj cell
                        write_file_obj(file_path_step_anc_cell, obj_cell)

                        # store file path destination cell
                        file_path_list_anc_cell.append(file_path_step_anc_cell)

                        # info end dump datasets cell
                        alg_logger.info(' ------> Save cell datasets ... DONE')

                    else:

                        # info end apply metrics to datasets cell
                        alg_logger.warning(' ===> Datasets are defined by NoneType')
                        alg_logger.info(' ------> Save cell datasets ... SKIPPED')

                    # info end method
                    alg_logger.info(' -----> Cell "' + str(cell_name) + '" ... DONE')

                else:
                    # message data ancillary already available
                    alg_logger.info(' -----> Cell "' + str(cell_name) +
                                    '"  ... SKIPPED. Ancillary cell data previously saved')
                    # store file path destination cell
                    file_path_list_anc_cell.append(file_path_step_anc_cell)

            # info end method
            alg_logger.info(' ----> Compute dynamic datasets ... DONE')

        else:
            # info end method
            file_path_list_anc_cell = []
            alg_logger.info(' ----> Compute dynamic datasets ... SKIPPED. Destination grid data previously saved')

        return file_path_list_anc_cell
# ----------------------------------------------------------------------------------------------------------------------
