"""
Class Features

Name:          drv_fx_runner
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import os
import numpy as np

from datetime import timedelta

from ascat.read_native.cdr import load_grid

from lib_utils_generic import reset_folder

from lib_fx_datasets_source import AscatNrtBufrFileList, AscatEpsBufrFileList
from lib_fx_datasets_destination import GriddedNcIndexedRaggedTs as GriddedNcRaggedTs
from lib_fx_datasets_destination import GriddedNcContiguousRaggedTs as GriddedNcRaggedTs
#from time_series import GriddedNcContiguousRaggedTs as GriddedNcRaggedTs
#from time_series import GriddedNcIndexedRaggedTs as GriddedNcRaggedTs
from lib_fx_datasets_generic import convert_sub_path_str_2_dict, convert_sub_path_2_root_path

from lib_fx_orbit_resampler_ascat import OrbitResamplerAscat

# suppress warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# product list supported
product_list = ['h16', 'h101', 'h102', 'h103', 'h104', 'h105']
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver fx
class DrvFxRunner:

    # method to initialize class
    def __init__(self, product_name,
                 time_start, time_end,
                 sub_path_swath, sub_path_cell, sub_path_workspace, sub_path_ts,
                 folder_name_swath, file_name_swath,
                 folder_name_ts, file_name_ts,
                 folder_name_cell, file_name_cell,
                 folder_name_workspace, file_name_workspace,
                 folder_name_grid, file_name_grid='TUW_WARP5_grid_info_2_3.nc',
                 writing_mode='w', write_n_resampled=1000,
                 spatial_resolution=25000, weight_function='hamming',
                 reset_ancillary_cell=False, reset_ancillary_workspace=False, reset_datasets_ts=False,
                 **kwargs):

        self.product_name = product_name
        self.time_start = time_start
        self.time_end = time_end

        self.sub_path_swath = sub_path_swath
        self.sub_path_cell = sub_path_cell
        self.sub_path_workspace = sub_path_workspace
        self.sub_path_ts = sub_path_ts

        self.folder_name_swath = folder_name_swath
        self.file_name_swath = file_name_swath
        self.folder_name_cell = folder_name_cell
        self.file_name_cell = file_name_cell
        self.folder_name_workspace = folder_name_workspace
        self.file_name_workspace = file_name_workspace
        self.folder_name_ts = folder_name_ts
        self.file_name_ts = file_name_ts

        self.folder_name_grid = folder_name_grid
        self.file_name_grid = file_name_grid

        self.writing_mode = writing_mode
        self.write_n_resampled = write_n_resampled

        self.spatial_resolution = spatial_resolution
        self.weight_function = weight_function

        if os.path.exists(os.path.join(self.folder_name_grid, self.file_name_grid)):
            self.target_grid = load_grid(os.path.join(self.folder_name_grid, self.file_name_grid))
        else:
            logging.error(' ===> File grid "' + os.path.join(self.folder_name_grid, self.file_name_grid) + "' not found")
            raise FileNotFoundError('File grid must be defined to correctly run the algorithm')

        self.reset_ancillary_cell = reset_ancillary_cell
        self.reset_ancillary_workspace = reset_ancillary_workspace
        self.reset_datasets_ts = reset_datasets_ts

    # method to select io class source and destination
    def select_class_io(self):

        # info method start
        logging.info(' ---> Select io class ... ')

        if self.product_name not in product_list:
            logging.error(' ===> Product "' + self.product_name + '" is not supported by the algorithm')
            raise NotImplemented('Case not implemented yet')

        product_class_lut_src = {
            'h16': AscatNrtBufrFileList, 'h103': AscatNrtBufrFileList,
            'h101': AscatNrtBufrFileList, 'h102': AscatNrtBufrFileList,
            'h104': AscatEpsBufrFileList, 'h105': AscatEpsBufrFileList
        }

        product_class_lut_dst = {
            'h16': GriddedNcRaggedTs, 'h103': GriddedNcRaggedTs,
            'h101': GriddedNcRaggedTs, 'h102': GriddedNcRaggedTs,
            'h104': GriddedNcRaggedTs, 'h105': GriddedNcRaggedTs
        }

        product_class_obj_src = product_class_lut_src[self.product_name]
        product_class_obj_dst = product_class_lut_dst[self.product_name]

        # info method end
        logging.info(' ---> Select io class ... DONE')

        return product_class_obj_src, product_class_obj_dst

    # method to select fx class
    def select_class_fx(self):

        # info method start
        logging.info(' ---> Select fx class ... ')

        if self.product_name not in product_list:
            logging.error(' ===> Product "' + self.product_name + '" is not supported by the algorithm')
            raise NotImplemented('Case not implemented yet')

        fx_class_lut = {
            'h16': OrbitResamplerAscat, 'h103': OrbitResamplerAscat,
            'h101': OrbitResamplerAscat, 'h102': OrbitResamplerAscat,
            'h104': OrbitResamplerAscat, 'h105': OrbitResamplerAscat
        }

        fx_class_obj = fx_class_lut[self.product_name]

        # info method end
        logging.info(' ---> Select fx class ... DONE')

        return fx_class_obj

    # method to organize source class
    def organize_class_io_src(self, product_class_obj):

        # info method start
        logging.info(' ---> Organize io class source ... ')

        time_start, time_end = self.time_start,  self.time_end

        if self.product_name == 'h16':

            # product_sub_folder_tmpl = {'years': '{year}', 'months': '{month}', 'days': '{day}', 'hours': '{hour}'}
            product_sub_folder_tmpl = convert_sub_path_str_2_dict(self.sub_path_swath)
            product_sub_folder_root = convert_sub_path_2_root_path(self.folder_name_swath)
            product_filename_tmpl = '{product_id}_{date}*.buf'

            product_id = self.product_name

            product_class_driver = product_class_obj(
                root_path=product_sub_folder_root, product_id=product_id,
                filename_template=product_filename_tmpl,
                subfolder_template=product_sub_folder_tmpl)

            product_time_stamps, product_time_intervals = np.array(
                product_class_driver.tstamps_for_daterange(time_start, time_end))
            product_file_names = product_class_driver.search_period(time_start, time_end)

            product_dt_delta = timedelta(minutes=3)
            product_dt_buffer = timedelta(days=0)

        elif self.product_name == 'h104':

            product_sub_folder_tmpl = None
            product_filename_tmpl = 'ASCA_{product_id}_02_M0{sat_id}_{date}Z_*_*_*_*.nat.gz'

            product_id = 'SMR'
            product_class_driver = product_class_obj(
                root_path=self.folder_name_swath, product_id=product_id, sat_name='C',
                filename_template=product_filename_tmpl,
                subfolder_template=product_sub_folder_tmpl)

            product_time_stamps, product_time_intervals = np.array(
                product_class_driver.tstamps_for_daterange(time_start, time_end))

            product_dt_delta = timedelta(minutes=100)
            product_dt_buffer = timedelta(hours=0)

            # ascat_file_names = ascat_io_swaths.search_period(start_dt, end_dt)
            # file_datasets = ascat_io_swaths.read_period(start_dt, end_dt)

        elif self.product_name == 'h105':

            product_sub_folder_tmpl = None
            product_filename_tmpl = 'ASCA_{product_id}_02_M0{sat_id}_{date}Z_*_*_*_*.nat.gz'

            product_id = 'SMO'
            product_class_driver = product_class_obj(
                root_path=self.folder_name_swath, product_id=product_id, sat_name='C',
                filename_template=product_filename_tmpl,
                subfolder_template=product_sub_folder_tmpl)

            product_time_stamps, product_time_intervals = np.array(
                product_class_driver.tstamps_for_daterange(time_start, time_end))

            product_dt_delta = timedelta(minutes=100)
            product_dt_buffer = timedelta(hours=0)

            # ascat_file_names = ascat_io_swaths.search_period(start_dt, end_dt)
            # file_datasets = ascat_io_swaths.read_period(start_dt, end_dt)

        else:
            logging.error(' ===> Product "' + self.product_name + '" is not supported by the algorithm')
            raise NotImplemented('Case not implemented yet')

        # info method end
        logging.info(' ---> Organize io class source ... DONE')

        return product_class_driver, product_time_intervals, product_dt_delta, product_dt_buffer

    # method to organize destination class
    def organize_class_io_dst(self, product_class_obj):

        # info method start
        logging.info(' ---> Organize io class destination ... ')

        folder_name_ts = self.folder_name_ts
        file_name_ts = self.file_name_ts

        product_class_driver = product_class_obj(
            folder_name_ts, file_name_ts=file_name_ts, grid=self.target_grid, mode=self.writing_mode,
            ioclass_kws={'time_units': "days since 1858-11-17 00:00:00"})

        # info method end
        logging.info(' ---> Organize io class destination ... DONE')

        return product_class_driver

    # method to organize fx class
    def organize_class_fx(self,
                          fx_class_obj,
                          product_class_driver_src, product_dt_delta, product_dt_buffer,
                          product_class_driver_dst):

        # info method start
        logging.info(' ---> Organize fx class ... ')

        fx_class_driver = fx_class_obj(
            product_class_driver_src, product_class_driver_dst,
            spatial_res=self.spatial_resolution,
            dt=15,
            wfunc=self.weight_function,
            write_orbit_buffer=True, to_xarray=True,
            dt_delta=product_dt_delta, dt_buffer=product_dt_buffer)

        # info method end
        logging.info(' ---> Organize fx class ... DONE')

        return fx_class_driver

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to execute fx class resampler
    def execute_class_fx_resampler(self, fx_class_driver, fx_time_intervals,
                                   init_file_ws=False, use_file_ws=False,
                                   init_file_cell=False, use_file_cell=False):

        # info method start
        logging.info(' ---> Execute fx resampler ... ')

        # organize workspace folders and filename(s)
        reset_folder(self.folder_name_workspace, self.reset_ancillary_workspace)
        file_path_workspace = os.path.join(self.folder_name_workspace, self.file_name_workspace)
        os.makedirs(self.folder_name_workspace, exist_ok=True)
        # organize cell folders and filename(s)
        reset_folder(self.folder_name_cell, self.reset_ancillary_cell)
        file_path_cell = os.path.join(self.folder_name_cell, self.file_name_cell)
        os.makedirs(self.folder_name_cell, exist_ok=True)
        # organize ts folders and filename(s)
        reset_folder(self.folder_name_ts, self.reset_datasets_ts)
        file_path_ts = os.path.join(self.folder_name_ts, self.file_name_ts)
        os.makedirs(self.folder_name_ts, exist_ok=True)

        # run resample method
        fx_class_driver.resample(
            fx_time_intervals,
            write_n_resampled=self.write_n_resampled,
            init_file_ws=init_file_ws, use_file_ws=use_file_ws, name_file_ws=file_path_workspace,
            init_file_cell=init_file_cell, use_file_cell=use_file_cell, name_file_cell=file_path_cell,
            name_file_ts=file_path_ts)

        # info method end
        logging.info(' ---> Execute fx resampler ... DONE')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
