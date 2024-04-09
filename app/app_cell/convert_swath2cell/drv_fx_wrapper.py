"""
Class Features

Name:          drv_fx_wrapper
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230804'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import shutil

from lib_utils_generic import reset_folder
from lib_utils_time import parse_time_string, fill_time_string

from drv_fx_runner import DrvFxRunner
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm readers list
product_list = ['h16', 'h101', 'h102', 'h103', 'h104', 'h105']
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver fx
class DrvFxWrapper:

    # method to initialize class
    def __init__(self, alg_settings,
                 alg_time_reference=None, alg_time_start=None, alg_time_end=None,
                 tag_section_flags='flags',
                 tag_section_product='product', tag_section_parameters='parameters',
                 tag_section_template='template',
                 tag_section_datasets='datasets',
                 tag_section_time='time', tag_section_log='log'):

        self.alg_time_reference = alg_time_reference
        self.alg_time_start = alg_time_start
        self.alg_time_end = alg_time_end

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_product = alg_settings[tag_section_product]
        self.alg_parameters = alg_settings[tag_section_parameters]
        self.alg_template = alg_settings[tag_section_template]
        self.alg_time = alg_settings[tag_section_time]
        self.alg_datasets_grid = alg_settings[tag_section_datasets]['static']
        self.alg_datasets_swath = alg_settings[tag_section_datasets]['dynamic']['swath']
        self.alg_datasets_cell = alg_settings[tag_section_datasets]['dynamic']['ancillary']['cell']
        self.alg_datasets_chunk = alg_settings[tag_section_datasets]['dynamic']['ancillary']['chunk']
        self.alg_datasets_workspace = alg_settings[tag_section_datasets]['dynamic']['ancillary']['workspace']
        self.alg_datasets_ts = alg_settings[tag_section_datasets]['dynamic']['ts']
        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_flag_active = 'active'
        self.tag_flag_initialize = 'initialize'

        self.reset_ancillary_cell = self.alg_flags['reset_ancillary_cell']
        self.reset_ancillary_chunk = self.alg_flags['reset_ancillary_chunk']
        self.reset_ancillary_workspace = self.alg_flags['reset_ancillary_workspace']
        self.reset_datasets_ts = self.alg_flags['reset_datasets_ts']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_grid = self.alg_datasets_grid[self.tag_folder_name]
        self.file_name_grid = self.alg_datasets_grid[self.tag_file_name]
        self.file_path_grid = os.path.join(self.folder_name_grid, self.file_name_grid)

        self.folder_name_swath = self.alg_datasets_swath[self.tag_folder_name]
        self.file_name_swath = self.alg_datasets_swath[self.tag_file_name]
        self.file_path_swath = os.path.join(self.folder_name_swath, self.file_name_swath)

        self.folder_name_ts = self.alg_datasets_ts[self.tag_folder_name]
        self.file_name_ts = self.alg_datasets_ts[self.tag_file_name]
        self.file_path_ts = os.path.join(self.folder_name_ts, self.file_name_ts)

        self.folder_name_cell = self.alg_datasets_cell[self.tag_folder_name]
        self.file_name_cell = self.alg_datasets_cell[self.tag_file_name]
        self.flag_active_cell = self.alg_datasets_cell[self.tag_flag_active]
        self.flag_initialize_cell = self.alg_datasets_cell[self.tag_flag_initialize]
        self.file_path_cell = os.path.join(self.folder_name_cell, self.file_name_cell)

        self.folder_name_chunk = self.alg_datasets_chunk[self.tag_folder_name]
        self.file_name_chunk = self.alg_datasets_chunk[self.tag_file_name]
        self.flag_active_chunk = self.alg_datasets_chunk[self.tag_flag_active]
        self.flag_initialize_chunk = self.alg_datasets_chunk[self.tag_flag_initialize]
        self.file_path_chunk = os.path.join(self.folder_name_chunk, self.file_name_chunk)

        self.folder_name_workspace = self.alg_datasets_workspace[self.tag_folder_name]
        self.file_name_workspace = self.alg_datasets_workspace[self.tag_file_name]
        self.flag_active_workspace = self.alg_datasets_workspace[self.tag_flag_active]
        self.flag_initialize_workspace = self.alg_datasets_workspace[self.tag_flag_initialize]
        self.file_path_workspace = os.path.join(self.folder_name_workspace, self.file_name_workspace)

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root=self.tag_folder_name, dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.drv_fx_obj = None

    # method to check datasets (if clean or not)
    def _clean_ancillary_folders(self, dset_obj,
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

    # method to organize fx args
    def organize_fx_args(self):

        # info method start
        logging.info(' ---> Organize fx args ... ')

        alg_product_name = None
        if 'name' in list(self.alg_product.keys()):
            alg_product_name = self.alg_product['name']
        alg_product_bbox = None
        if 'bbox' in list(self.alg_product.keys()):
            alg_product_bbox = self.alg_product['bbox']
        alg_write_mode = 'w'
        if 'writing_mode' in list(self.alg_parameters.keys()):
            alg_write_mode = self.alg_parameters['writing_mode']
        alg_write_n_resampled = 1000
        if 'write_n_resampled' in list(self.alg_parameters.keys()):
            alg_write_n_resampled = self.alg_parameters['write_n_resampled']
        alg_spatial_resolution = 25000
        if 'spatial_resolution' in list(self.alg_parameters.keys()):
            alg_spatial_resolution = self.alg_parameters['spatial_resolution']
        alg_weight_function = 'hamming'
        if 'weight_function' in list(self.alg_parameters.keys()):
            alg_weight_function = self.alg_parameters['weight_function']

        # get time information
        tmp_time_reference = self.alg_time_reference
        alg_time_reference = parse_time_string(tmp_time_reference)
        tmp_time_start = self.alg_time_start
        alg_time_start = parse_time_string(tmp_time_start)
        tmp_time_end = self.alg_time_end
        alg_time_end = parse_time_string(tmp_time_end)

        if (alg_product_name is None) or (alg_product_name not in product_list):
            logging.error(' ===> Product name "' + str(alg_product_name) + '" is not supported')
            raise RuntimeError('Product name must be defined by one of "' + str(product_list) + '" names')

        alg_sub_path_tmpl_swath = None
        if 'sub_path_swath' in list(self.alg_template.keys()):
            alg_sub_path_tmpl_swath = self.alg_template['sub_path_swath']

        alg_datetime_tmpl_cell = None
        if 'datetime_cell' in list(self.alg_template.keys()):
            alg_datetime_tmpl_cell = self.alg_template['datetime_cell']
        alg_sub_path_tmpl_cell = None
        if 'sub_path_cell' in list(self.alg_template.keys()):
            alg_sub_path_tmpl_cell = self.alg_template['sub_path_cell']

        alg_datetime_tmpl_chunk_start = None
        if 'datetime_chunk_start' in list(self.alg_template.keys()):
            alg_datetime_tmpl_chunk_start = self.alg_template['datetime_chunk_start']
        alg_datetime_tmpl_chunk_end = None
        if 'datetime_chunk_end' in list(self.alg_template.keys()):
            alg_datetime_tmpl_chunk_end = self.alg_template['datetime_chunk_end']
        alg_sub_path_tmpl_chunk = None
        if 'sub_path_chunk' in list(self.alg_template.keys()):
            alg_sub_path_tmpl_chunk = self.alg_template['sub_path_chunk']

        alg_datetime_tmpl_workspace_start = None
        if 'datetime_workspace_start' in list(self.alg_template.keys()):
            alg_datetime_tmpl_workspace_start = self.alg_template['datetime_workspace_start']
        alg_datetime_tmpl_workspace_end = None
        if 'datetime_workspace_end' in list(self.alg_template.keys()):
            alg_datetime_tmpl_workspace_end = self.alg_template['datetime_workspace_end']
        alg_sub_path_tmpl_workspace = None
        if 'sub_path_workspace' in list(self.alg_template.keys()):
            alg_sub_path_tmpl_workspace = self.alg_template['sub_path_workspace']

        alg_datetime_tmpl_ts = None
        if 'datetime_ts' in list(self.alg_template.keys()):
            alg_datetime_tmpl_ts = self.alg_template['datetime_ts']
        alg_sub_path_tmpl_ts = None
        if 'sub_path_ts' in list(self.alg_template.keys()):
            alg_sub_path_tmpl_ts = self.alg_template['sub_path_ts']

        # fill path cell
        file_path_cell = fill_time_string(
            self.file_path_cell, alg_time_reference,
            ['datetime_cell', 'sub_path_cell'],
            [alg_datetime_tmpl_cell, alg_sub_path_tmpl_cell])
        folder_name_cell, file_name_cell = os.path.split(file_path_cell)
        # fill path chunk
        file_path_chunk = fill_time_string(
            self.file_path_chunk, alg_time_reference,
            [],
            [alg_sub_path_tmpl_chunk])
        folder_name_chunk, file_name_chunk = os.path.split(file_path_chunk)
        # fill path workspace
        file_path_workspace = fill_time_string(
            self.file_path_workspace,
            [alg_time_start, alg_time_end, alg_time_reference],
            ['datetime_workspace_start',  'datetime_workspace_end', 'sub_path_workspace'],
            [alg_datetime_tmpl_workspace_start, alg_datetime_tmpl_workspace_end, alg_sub_path_tmpl_workspace])
        folder_name_workspace, file_name_workspace = os.path.split(file_path_workspace)
        # fill path ts
        file_path_ts = fill_time_string(
            self.file_path_ts, alg_time_reference,
            ['datetime_ts', 'sub_path_ts'],
            [alg_datetime_tmpl_ts, alg_sub_path_tmpl_ts])
        folder_name_ts, file_name_ts = os.path.split(file_path_ts)

        fx_kwargs = {
            'product_name': alg_product_name, 'product_bbox': alg_product_bbox,
            'writing_mode': alg_write_mode, 'write_n_resampled': alg_write_n_resampled,
            'time_reference': alg_time_reference, 'time_start': alg_time_start, 'time_end': alg_time_end,
            "folder_name_swath": self.folder_name_swath, "file_name_swath": self.file_name_swath,
            'folder_name_ts': folder_name_ts, 'file_name_ts': file_name_ts,
            'folder_name_grid': self.folder_name_grid, 'file_name_grid': self.file_name_grid,
            'folder_name_workspace': folder_name_workspace, 'file_name_workspace': file_name_workspace,
            'folder_name_chunk': folder_name_chunk, 'file_name_chunk': file_name_chunk,
            'folder_name_cell': folder_name_cell, 'file_name_cell': file_name_cell,
            'spatial_resolution': alg_spatial_resolution, 'weight_function': alg_weight_function,
            'sub_path_swath': alg_sub_path_tmpl_swath,
            'sub_path_cell': alg_sub_path_tmpl_cell, 'sub_path_workspace': alg_sub_path_tmpl_workspace,
            'sub_path_ts': alg_sub_path_tmpl_ts,
            'reset_ancillary_cell': self.reset_ancillary_cell,
            'reset_ancillary_chunk': self.reset_ancillary_chunk,
            'reset_ancillary_workspace': self.reset_ancillary_workspace,
            'reset_datasets_ts': self.reset_datasets_ts}

        # info method end
        logging.info(' ---> Organize fx args ... DONE')

        return fx_kwargs

    # method to organize fx classes
    def organize_fx_classes(self, fx_kwargs):

        # info method start
        logging.info(' ---> Organize fx classes ... ')

        # set driver fx runner
        self.drv_fx_obj = DrvFxRunner(**fx_kwargs)

        # select product classes source and destination
        product_class_obj_src, product_class_obj_dst = self.drv_fx_obj.select_class_io()
        # select fx class
        fx_class_obj = self.drv_fx_obj.select_class_fx()
        # define product driver srt
        product_class_driver_src, product_time_intervals, \
            product_dt_delta, product_dt_buffer = self.drv_fx_obj.organize_class_io_src(product_class_obj_src)
        # define product driver dst
        product_class_driver_dst = self.drv_fx_obj.organize_class_io_dst(product_class_obj_dst)
        # define fx driver
        fx_class_driver = self.drv_fx_obj.organize_class_fx(
            fx_class_obj,
            product_class_driver_src, product_dt_delta, product_dt_buffer,
            product_class_driver_dst)

        # define class time obj
        fx_class_time = {'time_intervals': product_time_intervals,
                         'time_dt_delta': product_dt_delta, 'time_dt_buffer': product_dt_buffer}

        # info method end
        logging.info(' ---> Organize fx classes... DONE')

        return fx_class_driver, fx_class_time

    # method to execute fx
    def execute_fx(self, fx_class_driver, fx_class_time):

        # info method start
        logging.info(' ---> Execute fx ... ')

        # get fx time
        fx_time_intervals = fx_class_time['time_intervals']
        # run fx class
        self.drv_fx_obj.execute_class_fx_resampler(
            fx_class_driver, fx_time_intervals,
            init_file_ws=self.flag_initialize_workspace, use_file_ws=self.flag_active_workspace,
            init_file_chunk=self.flag_initialize_chunk, use_file_chunk=self.flag_active_chunk,
            init_file_cell=self.flag_initialize_cell, use_file_cell=self.flag_active_cell
        )

        # info method end
        logging.info(' ---> Execute fx ... DONE')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------

