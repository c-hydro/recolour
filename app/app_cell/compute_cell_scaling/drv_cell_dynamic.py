"""
Class Features

Name:          drv_map_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230822'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_info_args import logger_name
from lib_utils_io import fill_path_with_tags
from lib_utils_generic import reset_folder, pad_list, remove_duplicates_list
from lib_cell_datasets import get_datasets_cell, apply_scaling_datasets_cell, dump_datasets_cell, dump_datasets_grid

# set logger
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class driver cell
class DrvCell:

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

        self.alg_datasets_src_cell = alg_obj_settings['datasets']['dynamic']['source']['cell_datasets_data']
        self.alg_datasets_ref_cell = alg_obj_settings['datasets']['dynamic']['source']['cell_datasets_ref']
        self.alg_datasets_dst_cell = alg_obj_settings['datasets']['dynamic']['destination']['cell_datasets']

        self.alg_datasets_name = self.alg_info['datasets']
        self.alg_datasets_product = self.alg_info['product']

        self.tag_active = 'active'
        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_file_vars_in = 'variables_in'
        self.tag_file_vars_out = 'variables_out'
        self.tag_file_format = 'format'

        self.reset_dst_cell = self.alg_flags['reset_destination_cell']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_src_cell = self.alg_datasets_src_cell[self.tag_folder_name]
        self.file_name_src_cell = self.alg_datasets_src_cell[self.tag_file_name]
        self.file_path_src_cell = os.path.join(self.folder_name_src_cell, self.file_name_src_cell)
        self.variables_in_src_cell = self.alg_datasets_src_cell[self.tag_file_vars_in]
        self.variables_out_src_cell = self.alg_datasets_src_cell[self.tag_file_vars_out]

        self.folder_name_ref_cell = self.alg_datasets_ref_cell[self.tag_folder_name]
        self.file_name_ref_cell = self.alg_datasets_ref_cell[self.tag_file_name]
        self.file_path_ref_cell = os.path.join(self.folder_name_ref_cell, self.file_name_ref_cell)
        self.variables_in_ref_cell = self.alg_datasets_ref_cell[self.tag_file_vars_in]
        self.variables_out_ref_cell = self.alg_datasets_ref_cell[self.tag_file_vars_out]

        self.folder_name_dst_cell = self.alg_datasets_dst_cell[self.tag_folder_name]
        self.file_name_dst_cell = self.alg_datasets_dst_cell[self.tag_file_name]
        self.file_path_dst_cell = os.path.join(self.folder_name_dst_cell, self.file_name_dst_cell)
        self.variables_in_dst_cell = self.alg_datasets_dst_cell[self.tag_file_vars_in]
        self.variables_out_dst_cell = self.alg_datasets_dst_cell[self.tag_file_vars_out]

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.grid_data = self.alg_obj_static['grid_cell_data']
        self.grid_ref = self.alg_obj_static['grid_cell_ref']
        self.list_gpis_data = self.alg_obj_static['list_gpis_data']
        self.list_grid_data = self.alg_obj_static['list_cell_data']
        self.list_gpis_ref = self.alg_obj_static['list_gpis_ref']
        self.list_grid_ref = self.alg_obj_static['list_cell_ref']

        self.list_grid_common = list(set(self.list_grid_data).intersection(self.list_grid_ref))

        self.vars_name = self.alg_parameters['variables']['name']
        self.vars_min_value = self.alg_parameters['variables']['min_value']
        self.vars_max_value = self.alg_parameters['variables']['max_value']
        self.vars_scale_factor = self.alg_parameters['variables']['scale_factor']
        self.vars_no_data = self.alg_parameters['variables']['no_data']

        self.fx_active = self.alg_parameters['fx']['active']
        self.fx_method = self.alg_parameters['fx']['methods']
        self.fx_vars_src = self.alg_parameters['fx']['variables']['source']
        self.fx_vars_dst = self.alg_parameters['fx']['variables']['destination']

        # grid file
        self.file_name_grid = 'grid.nc'
        self.folder_name_grid = self.folder_name_dst_cell
        self.file_path_grid = os.path.join(self.folder_name_grid, self.file_name_grid)

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

    # method to compute data
    def compute_data(self):

        # info start method
        alg_logger.info(' ----> Compute dynamic datasets ... ')

        # get file time
        file_time = self.alg_obj_time

        # get grid cells, obj and reference
        list_grid_data, list_grid_ref, list_grid_common = self.list_grid_data, self.list_grid_ref, self.list_grid_common
        grid_data, grid_ref = self.grid_data, self.grid_ref

        # get variable(s) source and destination
        variables_in_src_cell, variables_out_src_cell = self.variables_in_src_cell, self.variables_out_src_cell
        variables_in_ref_cell, variables_out_ref_cell = self.variables_in_ref_cell, self.variables_out_ref_cell
        variables_in_dst_cell, variables_out_dst_cell = self.variables_in_dst_cell, self.variables_out_dst_cell
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
        # remove duplicates
        var_min_value = remove_duplicates_list(vars_min_value)[0]
        var_max_value = remove_duplicates_list(vars_max_value)[0]
        var_no_data_value = remove_duplicates_list(vars_no_data)[0]
        max_dist = 25000

        # get file path
        file_path_src_cell_raw, file_path_ref_cell_raw = self.file_path_src_cell, self.file_path_ref_cell
        file_path_dst_cell_raw = self.file_path_dst_cell
        file_path_grid_raw = self.file_path_grid

        # iterate over grid cells
        file_path_dst_cell_list = []
        for cell_name in list_grid_common:

            # define source cell file name
            file_path_src_cell_def = fill_path_with_tags(
                file_path_src_cell_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

            # define reference cell file name
            file_path_ref_cell_def = fill_path_with_tags(
                file_path_ref_cell_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

            # define destination cell file name
            file_path_dst_cell_def = fill_path_with_tags(
                file_path_dst_cell_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

            # define destination grid file name
            file_path_grid_def = fill_path_with_tags(
                file_path_grid_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(cell_name)}, tmpl_tags_dset=alg_template_dset)

            # clean destination file name
            if self.reset_dst_cell:
                if os.path.exists(file_path_dst_cell_def):
                    os.remove(file_path_dst_cell_def)

            # info start method
            alg_logger.info(' -----> Cell "' + str(cell_name) + '" ... ')

            # check fx active or not
            if self.fx_active:

                # check file cell destination availability
                if not os.path.exists(file_path_dst_cell_def):

                    # check file cell source availability
                    if os.path.exists(file_path_src_cell_def):

                        # info start get datasets cell
                        alg_logger.info(' ------> Get cell source datasets ... ')

                        # check file source cell datasets availability
                        if os.path.exists(file_path_src_cell_def):

                            # get collections and registry objects
                            collections_dframe_cell_src, registry_dframe_cell_src = get_datasets_cell(
                                file_name=file_path_src_cell_def,
                                cell_name=cell_name,
                                list_variable_data_in=variables_in_src_cell['datasets'],
                                list_variable_data_out=variables_out_src_cell['datasets'],
                                list_variable_registry_in=variables_in_src_cell['registry'],
                                list_variable_registry_out=variables_out_src_cell['registry'],
                                index_name_data_in='time', index_name_data_out='time_data',
                                index_name_registry_out='location_id_data', index_name_cell_out='cell_data'
                            )

                            # info end get datasets cell
                            alg_logger.info(' ------> Get cell source datasets ... DONE')

                        else:
                            # info end get datasets cell
                            collections_dframe_cell_src, registry_dframe_cell_src = None, None
                            alg_logger.warning(' ===> File not found "' + str(file_path_src_cell_def) + '"')
                            alg_logger.info(' ------> Get cell source datasets ... FAILED.')

                        # info start get datasets cell
                        alg_logger.info(' ------> Get cell reference datasets ... ')

                        # check file reference cell datasets availability
                        if os.path.exists(file_path_ref_cell_def):

                            # get collections and registry objects
                            collections_dframe_cell_ref, registry_dframe_cell_ref = get_datasets_cell(
                                file_name=file_path_ref_cell_def,
                                cell_name=cell_name,
                                list_variable_data_in=variables_in_ref_cell['datasets'],
                                list_variable_data_out=variables_out_ref_cell['datasets'],
                                list_variable_registry_in=variables_in_ref_cell['registry'],
                                list_variable_registry_out=variables_out_ref_cell['registry'],
                                index_name_data_in='time', index_name_data_out='time_ref',
                                index_name_registry_out='location_id_ref', index_name_cell_out='cell_ref'
                            )

                            # info end get datasets cell
                            alg_logger.info(' ------> Get cell reference datasets ... DONE')

                        else:
                            # info end get datasets cell
                            collections_dframe_cell_ref, registry_dframe_cell_ref = None, None
                            alg_logger.warning(' ===> File not found "' + str(file_path_ref_cell_def) + '"')
                            alg_logger.info(' ------> Get cell reference datasets ... FAILED.')

                        # info start apply scaling to datasets cell
                        alg_logger.info(' ------> Apply scaling to cell datasets ... ')

                        # check collections and registry objects availability
                        if collections_dframe_cell_src is not None and collections_dframe_cell_ref is not None:

                            # apply metrics to datasets
                            (collections_dframe_out, registry_dframe_out,
                             collections_attrs_out, registry_attrs_out) = apply_scaling_datasets_cell(
                                file_time,
                                collections_dframe_cell_src, registry_dframe_cell_src,
                                collections_dframe_cell_ref, registry_dframe_cell_ref,
                                grid_obj_data=grid_data, grid_obj_ref=grid_ref,
                                fx_var_name_in=self.fx_vars_src,
                                fx_var_name_out=self.fx_vars_dst,
                                fx_var_methods=self.fx_method,
                                max_dist=max_dist,
                                var_min=var_min_value, var_max=var_max_value, var_no_data=var_no_data_value,
                                index_name_time_data='time_data', index_name_time_ref='time_ref',
                                index_name_location_data='location_id_data', index_name_location_ref='location_id_ref',
                                index_name_geo_x_data='lon_data', index_name_geo_y_data='lat_data',
                                index_name_cell_data='cell_data', index_name_row_size_data='row_size_data',
                                index_name_geo_x_ref='lon_ref', index_name_geo_y_ref='lat_ref',
                                index_name_cell_ref='cell_ref', index_name_row_size_ref='row_size_ref'
                            )

                            # info end apply metrics to datasets cell
                            alg_logger.info(' ------> Apply scaling to cell datasets ... DONE')

                        else:
                            # info end apply metrics to datasets cell
                            collections_dframe_out, registry_dframe_out = None, None
                            collections_attrs_out, registry_attrs_out = None, None
                            alg_logger.warning(' ===> Data not available for scaling')
                            alg_logger.info(' ------> Apply scaling to cell datasets ... FAILED')

                        # info start dump datasets cell
                        alg_logger.info(' ------> Dump cell datasets ... ')

                        # check collections and registry objects availability
                        if collections_dframe_out is not None and registry_dframe_out is not None:

                            # dump collections and registry objects
                            dump_datasets_cell(file_path_dst_cell_def,
                                               collections_dframe_out, registry_dframe_out,
                                               collections_attrs_out, registry_attrs_out,
                                               list_variable_data_in=variables_in_dst_cell['datasets'],
                                               list_variable_data_out=variables_out_dst_cell['datasets'],
                                               list_variable_registry_in=variables_in_dst_cell['registry'],
                                               list_variable_registry_out=variables_out_dst_cell['registry'])

                            # info end dump datasets cell
                            alg_logger.info(' ------> Dump cell datasets ... DONE')

                            # info start dump grid cell
                            alg_logger.info(' ------> Dump cell grid ... ')

                            # dump grid object
                            dump_datasets_grid(file_path_grid_def, grid_data)

                            # info end dump grid cell
                            alg_logger.info(' ------> Dump cell grid ... DONE')

                            # store file path destination cell
                            file_path_dst_cell_list.append(file_path_dst_cell_def)

                        else:

                            # info end apply metrics to datasets cell
                            alg_logger.warning(' ===> Data not available for dumping')
                            alg_logger.info(' ------> Dump cell datasets ... FAILED')

                        # info end method
                        alg_logger.info(' -----> Cell "' + str(cell_name) + '" ... DONE')

                    else:
                        # message data source not available
                        alg_logger.warning(' ===> Source data "' + str(file_path_src_cell_def) + '" not available')
                        alg_logger.info(' -----> Cell "' + str(cell_name) + '"  ... SKIPPED. Source data not available')
                else:
                    # message data destination already available
                    alg_logger.info(' -----> Cell "' + str(cell_name) + '"  ... SKIPPED. Data previously saved')
            else:
                # message fx not active
                alg_logger.info(' -----> Cell "' + str(cell_name) + '"  ... SKIPPED. Compute data not active')

        # info end method
        alg_logger.info(' ----> Compute dynamic datasets ... DONE')

        return file_path_dst_cell_list
# ----------------------------------------------------------------------------------------------------------------------
