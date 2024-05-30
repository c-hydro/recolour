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
from lib_utils_generic import reset_folder, pad_list
from lib_utils_datasets import get_datasets_cell, apply_rzsm_datasets_cell, dump_datasets_cell, dump_datasets_grid

# set logger
alg_logger = logging.getLogger(logger_name)

# debug
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class driver cell
class DrvCell:

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
        self.alg_datasets_dst = alg_obj_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_obj_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_file_vars_in = 'variables_in'
        self.tag_file_vars_out = 'variables_out'
        self.tag_file_format = 'format'

        self.alg_datasets_name = self.alg_info['datasets']
        self.alg_datasets_product = self.alg_info['product']

        self.reset_dst = self.alg_flags['reset_destination']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_src = self.alg_datasets_src[self.tag_folder_name]
        self.file_name_src = self.alg_datasets_src[self.tag_file_name]
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)
        self.variables_in_src = self.alg_datasets_src[self.tag_file_vars_in]
        self.variables_out_src = self.alg_datasets_src[self.tag_file_vars_out]

        self.folder_name_dst = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_datasets_dst[self.tag_file_name]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)
        self.variables_in_dst = self.alg_datasets_dst[self.tag_file_vars_in]
        self.variables_out_dst = self.alg_datasets_dst[self.tag_file_vars_out]

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.grid_cells = self.alg_obj_static['grid_cells']
        self.grid_gpis = self.alg_obj_static['grid_gpis']
        self.grid_reference = self.alg_obj_static['grid_reference']
        self.grid_obj = self.alg_obj_static['grid_obj']

        self.vars_name = self.alg_parameters['variables']['name']
        self.vars_mode = self.alg_parameters['variables']['mode']
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
        self.folder_name_grid = self.folder_name_dst
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
        grid_cells_product = self.grid_cells
        grid_obj_product = self.grid_obj
        grid_reference_product = self.grid_reference

        # get file variables
        variables_in_src, variables_out_src = self.variables_in_src, self.variables_out_src
        variables_in_data_dst, variables_out_data_dst = self.variables_in_dst, self.variables_out_dst
        # get file template
        alg_template_time, alg_template_dset = self.alg_template_time, self.alg_template_dset

        # get variables info
        vars_name = self.vars_name
        vars_mode = self.vars_mode
        vars_min_value = self.vars_min_value
        vars_max_value = self.vars_max_value
        vars_scale_factor = self.vars_scale_factor
        vars_no_data = self.vars_no_data
        # pad list (using names as list length reference
        vars_mode, vars_min_value, vars_max_value, vars_scale_factor, vars_no_data = pad_list(
            vars_name, [vars_mode, vars_min_value, vars_max_value, vars_scale_factor, vars_no_data])

        # get file path
        file_path_src_raw, file_path_dst_raw = self.file_path_src, self.file_path_dst
        file_path_grid_raw = self.file_path_grid

        # iterate over grid cells
        file_path_dst_list = []
        for grid_cells_name in grid_cells_product:

            # define source file name
            file_path_src_def = fill_path_with_tags(
                file_path_src_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(grid_cells_name)}, tmpl_tags_dset=alg_template_dset)

            # define destination file name
            file_path_dst_def = fill_path_with_tags(
                file_path_dst_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(grid_cells_name)}, tmpl_tags_dset=alg_template_dset)

            # define grid file name
            file_path_grid_def = fill_path_with_tags(
                file_path_grid_raw,
                tmpl_time_obj=file_time, tmpl_tags_time=alg_template_time,
                tmpl_dset_obj={'cell_n': str(grid_cells_name)}, tmpl_tags_dset=alg_template_dset)

            # clean destination file name
            if self.reset_dst:
                if os.path.exists(file_path_dst_def):
                    os.remove(file_path_dst_def)

            # info start method
            alg_logger.info(' -----> Cell "' + str(grid_cells_name) + '" ... ')

            # check fx active or not
            if self.fx_active:

                # check file cell source availability
                if os.path.exists(file_path_src_def):

                    # check file point ancillary availability
                    if not os.path.exists(file_path_dst_def):

                        # info start get datasets cell
                        alg_logger.info(' ------> Get datasets ... ')

                        # get collections and registry objects
                        collections_obj_in, registry_obj_in = get_datasets_cell(
                            file_name=file_path_src_def,
                            cell_name=grid_cells_name,
                            list_variable_data_in=variables_in_src['datasets'],
                            list_variable_data_out=variables_out_src['datasets'],
                            list_variable_registry_in=variables_in_src['registry'],
                            list_variable_registry_out=variables_out_src['registry'],
                            index_name_data='time', index_name_cell='cell', index_name_registry='location_id'
                        )

                        # info end get datasets cell
                        alg_logger.info(' ------> Get datasets ... DONE')

                        # info start apply fx dataset cell
                        alg_logger.info(' ------> Filter datasets ... ')

                        # filter collections and registry objects
                        (collections_obj_out, registry_obj_out,
                         collections_attrs_out, registry_attrs_out) = apply_rzsm_datasets_cell(
                            file_time, collections_obj_in, registry_obj_in,
                            fx_var_name_in=self.fx_vars_src,
                            fx_var_name_out=self.fx_vars_dst,
                            fx_var_methods=self.fx_method,
                            data_var_name=vars_name, data_var_mode=vars_mode,
                            data_var_min_value=vars_min_value, data_var_max_value=vars_max_value,
                            data_var_scale_factor=vars_scale_factor, data_var_undef=vars_no_data
                        )

                        # info start apply fx to dataset cell
                        alg_logger.info(' ------> Filter datasets ... DONE')

                        # info start dump datasets cell
                        alg_logger.info(' ------> Dump datasets ... ')

                        # dump collections and registry objects
                        dump_datasets_cell(file_path_dst_def,
                                           collections_obj_out, registry_obj_out,
                                           collections_attrs_out, registry_attrs_out,
                                           list_variable_data_in=variables_in_data_dst['datasets'],
                                           list_variable_data_out=variables_out_data_dst['datasets'],
                                           list_variable_registry_in=variables_in_data_dst['registry'],
                                           list_variable_registry_out=variables_out_data_dst['registry'])

                        # info end dump datasets cell
                        alg_logger.info(' ------> Dump datasets ... DONE')

                        # info start dump grid cell
                        alg_logger.info(' ------> Dump grid ... ')

                        # dump grid object
                        dump_datasets_grid(file_path_grid_def, grid_obj_product, grid_reference_product)

                        # info end dump grid cell
                        alg_logger.info(' ------> Dump grid ... DONE')

                        # info start method
                        alg_logger.info(' -----> Cell "' + str(grid_cells_name) + '" ... DONE')
                    else:

                        # info end get datasets cell
                        alg_logger.info(' -----> Cell "' + str(grid_cells_name) +
                                        '"  ... SKIPPED. Data previously saved')

                    # store cell filename(s) -- if compute data not active return destination path file
                    file_path_dst_list.append(file_path_dst_def)

                else:
                    # message data source not available
                    alg_logger.warning(' ===> Source data "' + str(file_path_src_def) + '" not available')
                    alg_logger.info(' -----> Cell "' + str(grid_cells_name) +
                                    '"  ... SKIPPED. Source data not available')

            else:
                # message fx not active
                alg_logger.info(' -----> Cell "' + str(grid_cells_name) + '"  ... SKIPPED. Compute data not active')

        # info end method
        alg_logger.info(' ----> Compute dynamic datasets ... DONE')

        return file_path_dst_list
# ----------------------------------------------------------------------------------------------------------------------
