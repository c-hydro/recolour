"""
Class Features

Name:          drv_process_main
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import shutil
import netCDF4
import multiprocessing as mp

from datetime import datetime
from itertools import repeat

from cpl_process_metrics import CplMetrics
from cpl_process_datasets import CplDatasets
from cpl_process_time import CplTime
from cpl_process_analysis import CplAnalysis
from cpl_process_mode import CplMode
from cpl_process_logs import CplLog

from lib_utils_log import open_file_log, dump_message_to_file, close_file_log

from pygeogrids.grids import CellGrid
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver process
class DrvProcess:

    # method to initialize class
    def __init__(self, alg_cells, alg_settings,
                 tag_section_mode='mode', tag_section_flags='flags',
                 tag_section_domain='domain',
                 tag_section_params='parameters', tag_section_datasets='datasets',
                 tag_section_time='time', tag_section_log='log'):

        self.alg_cells = alg_cells

        self.alg_mode = alg_settings[tag_section_mode]
        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_params = alg_settings[tag_section_params]
        self.alg_domain = alg_settings[tag_section_domain]
        self.alg_time = alg_settings[tag_section_time]
        self.alg_datasets_src = alg_settings[tag_section_datasets]['source']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['destination']

        if 'swi' in alg_settings[tag_section_datasets]['source']['k2']:
            self.swi_option = alg_settings[tag_section_datasets]['source']['k2']['swi']
        else:
            logging.warning(' ===> SWI option not activated for k2 datasets in the configuration file')
            self.swi_option = None

        self.alg_log = alg_settings[tag_section_log]

        self.coupler_logs = None
        self.coupler_time = None
        self.coupler_datasets = None
        self.coupler_metrics = None
        self.coupler_analysis = None
        self.coupler_mode = None

        self.dset_process_groups = None
        self.dset_process_execution = None
        self.dset_process_cpu = None
        self.dset_process_count = 0

        self.file_log_generic = None
        self.file_log_groups = None
        self.file_log_process = None

        self.dset_process_status = False

        self.dset_path_dst = None

        self.reset_datasets_src = self.alg_flags['reset_datasets_src']
        self.reset_datasets_dst = self.alg_flags['reset_datasets_dst']
        self.reset_logs = self.alg_flags['reset_logs']

        self._clean_ancillary_folders(self.alg_datasets_src,
                                      dset_key_root='tmp', dset_key_sub='path_tmp',
                                      dset_clean=self.reset_datasets_src)
        self._clean_ancillary_folders(self.alg_datasets_dst,
                                      dset_key_root='path_analysis', dset_key_sub=None,
                                      dset_clean=self.reset_datasets_dst)
        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

    # method to check datasets (if clean or not)
    def _clean_ancillary_folders(self, dset_obj, dset_key_root=None, dset_key_sub=None,
                        dset_clean=True):

        if (dset_key_root is not None) and (dset_key_sub is not None):
            for dset_key, dset_fields in dset_obj.items():
                dset_path = dset_fields[dset_key_root][dset_key_sub]
                self._remove_datasets(dset_path, dset_clean=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is not None):
            dset_path = dset_obj[dset_key_sub]
            self._remove_datasets(dset_path, dset_clean=dset_clean)
        elif (dset_key_root is not None) and (dset_key_sub is None):
            dset_path = dset_obj[dset_key_root]
            self._remove_datasets(dset_path, dset_clean=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is None):
            pass
        else:
            pass

    # method to remove datasets
    @staticmethod
    def _remove_datasets(dset_path, dset_clean=True):
        if os.path.exists(dset_path):
            if dset_clean:
                shutil.rmtree(dset_path)
        os.makedirs(dset_path, exist_ok=True)

    # method to setup process
    def setup_process(self):

        # info setup process start
        logging.info(' ---> Setup process ... ')

        # method to set and organize logs information
        self.coupler_logs = CplLog(self.alg_log)
        self.file_log_generic, self.file_log_groups, self.file_log_process = self.coupler_logs.setup_logs()

        # method to set and organize time information
        self.coupler_time = CplTime(self.alg_time, time_window_hours=self.alg_params['temporal_window'])
        time_period, time_window_fraction = self.coupler_time.setup_time()

        # method(s) to set and organize datasets information
        self.coupler_datasets = CplDatasets(self.alg_datasets_src, self.alg_datasets_dst)
        dset_interfaces, dset_modes = self.coupler_datasets.setup_datasets_src()
        self.dset_path_dst = self.coupler_datasets.setup_datasets_dst()

        if len(dset_interfaces) == 3:
            metrics_type = 'extended'
        elif len(dset_interfaces) == 2:
            metrics_type = 'basic'
        else:
            raise NotImplemented('Case not implemented yet.')

        # method to set and organize metrics information
        self.coupler_metrics = CplMetrics(dset_interfaces,
                                          metrics_seasonal=self.alg_params['seasonal_metrics'],
                                          metrics_anomaly=self.alg_params['seasonal_anomaly'],
                                          metrics_type=metrics_type,
                                          swi_option=self.swi_option)
        dset_metrics = self.coupler_metrics.setup_metrics()

        # method to set and organize analysis information
        self.coupler_analysis = CplAnalysis(dset_interfaces, dset_metrics, time_period, time_window_fraction)

        # method to set and organize mode information
        self.coupler_mode = CplMode(self.alg_mode)
        self.dset_process_groups, self.dset_process_execution, self.dset_process_cpu = self.coupler_mode.setup_mode(
            cells=self.alg_cells, file_log_groups=self.file_log_groups)

        self.dset_process_status = True

        # info setup process end
        logging.info(' ---> Setup process ... DONE')

    # method to execute process
    def execute_process(self):

        # info execute process start
        logging.info(' ---> Execute process in "' + self.dset_process_execution + '" ... ')

        process_status = self.dset_process_status
        process_cells = self.alg_cells
        process_groups = self.dset_process_groups

        # wrap process for sequential mode
        if self.dset_process_execution == 'sequential':

            # iterate over cell(s)
            for cell_group, cell_n in enumerate(process_cells):

                # run process in sequential mode
                process_results = self.exec_process_wrap([cell_group, cell_n])

        # wrap process for parallel mode
        elif self.dset_process_execution == 'parallel':

            # open parallel process(es) using pool
            mp_pool = mp.Pool(processes=self.dset_process_cpu)

            # iterate over cell group(s)
            for cell_group, cell_list in enumerate(process_groups):

                # update args for parallel mode
                process_args = zip(repeat(cell_group), cell_list)
                # run process in parallel mode
                process_results = mp_pool.map_async(self.exec_process_wrap, process_args)

            # close and join parallel process(es)
            mp_pool.close()
            mp_pool.join()

        else:
            logging.error(' ===> Process mode "' + self.dset_process_execution + '" is not expected by the algorithm')
            raise NotImplemented('Case not implemented yet')

        # info wrap process end
        logging.info(' ---> Execute process in "' + self.dset_process_execution + '" ... DONE')

    # method to wrap process execution
    def exec_process_wrap(self, dset_settings):

        # get process time start
        process_time_start = datetime.now()
        # parse the fx arguments
        dset_group, dset_cell = dset_settings[0], dset_settings[1]

        # create log file to handle process
        file_log_handle, file_log_group, file_log_cell = open_file_log(
            dset_group, dset_cell, execution_mode=self.dset_process_execution,
            file_path=self.file_log_process, file_update=True)

        # info start process execution
        logging.info(' ----> Execute process ... ')
        logging.info(' ----> Start time "' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '" ')
        logging.info('  (1) Process Group: "' + file_log_group + '"')
        logging.info('  (2) Process Cell(s): "' + file_log_cell + '"')

        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='Execute process ... ')
        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='Start time "' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '"')
        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='  (1) Process Group: "' + file_log_group + '"')
        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='  (2) Process Cell(s): "' + file_log_cell + '"')

        # set and define process obj and jobs
        dset_process_obj, dset_process_jobs = self.coupler_analysis.setup_analysis(dset_cell)
        dset_job_n = dset_process_jobs.__len__()

        # iterate over jobs
        for dset_job_idx, dset_job_step in enumerate(dset_process_jobs):

            # compute job idx
            file_job_idx = str(dset_job_idx + 1).zfill(len(str(dset_job_n)))
            # compute job percentage
            dset_job_percentage = ((dset_job_idx + 1) / dset_job_n) * 100
            file_job_percentage = "{:.2f}".format(dset_job_percentage)

            # info start process job
            logging.info(
                ' -----> Process Info :: IDX: ' + file_job_idx + ' of ' + str(dset_job_n) + ' :: JOB: "' +
                str(dset_job_step) + '" :: STEP: ' + file_job_percentage + ' [%] ... ')

            dump_message_to_file(
                file_log_handle, line_arrow='---->', line_break=True,
                line_msg='Process Info :: IDX: ' + file_job_idx + ' of ' + str(dset_job_n) + ' :: JOB: ' +
                str(dset_job_step) + '" :: STEP: ' + file_job_percentage + ' [%] ... ')

            # debug
            # dset_job_step = dset_process_jobs[43]
            # dset_job_step = dset_process_jobs[10]

            # ***********************************************************
            # execute process core                                      #
            # ***********************************************************
            results = self.exec_process_core(dset_job_step, dset_process_obj, self.dset_path_dst)
            # ***********************************************************
            #                                                           #
            # ***********************************************************

            # info start process job
            logging.info(
                ' -----> Process Info :: IDX: ' + file_job_idx + ' of ' + str(dset_job_n) + ' :: JOB: "' +
                str(dset_job_step) + '" :: STEP: ' + file_job_percentage + ' [%] ... DONE')

            dump_message_to_file(
                file_log_handle, line_arrow='---->', line_break=True,
                line_msg='Process Info :: IDX: ' + file_job_idx + ' of ' + str(dset_job_n) + ' :: JOB: ' +
                str(dset_job_step) + ' :: STEP: ' + file_job_percentage + ' [%] ... DONE')

        # get process time end, time_diff and time_elapsed
        process_time_end = datetime.now()
        process_time_diff = process_time_end - process_time_start
        process_time_elapsed = divmod(process_time_diff.total_seconds(), 3600)[0]

        # info end process execution
        logging.info(' ----> End time "' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '" \n')
        logging.info(' ----> Elapsed time "' + str(process_time_elapsed) + '" [hours] \n')
        logging.info(' ----> Execute process ... DONE')

        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='End time "' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '"')
        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='Elapsed time "' + str(process_time_elapsed) + '" [hours]')
        dump_message_to_file(file_log_handle, line_arrow='---->', line_break=True,
                             line_msg='Execute process ... DONE')

        # method to close log file
        close_file_log(file_log_handle)

    # method to execute process core
    @staticmethod
    def exec_process_core(job, process, save_path=None, write_to_cells=True, handle_errors='ignore'):

        # run computation(s) over datasets
        if not write_to_cells:
            # return results
            return process.calc(*job, handle_errors=handle_errors)
        else:
            # dump results in separate files per cell
            result = process.calc(*job, handle_errors=handle_errors)
            if save_path is None:
                return result

            if bool(result) is True:

                if process.data_manager.reference_grid is not None:
                    if type(process.data_manager.reference_grid) is CellGrid:
                        grid = process.data_manager.reference_grid
                        try:
                            cell = grid.gpi2cell(job[0])
                        except:
                            cell = grid.gpi2cell(job[0][0])

                        gpis, lon, lat = grid.grid_points_for_cell(cell)
                        for key in result.keys():
                            fname = '{:04d}'.format(cell) + '.nc'
                            filename = os.path.join(save_path, fname)

                            if not os.path.exists(filename):
                                ncfile = netCDF4.Dataset(filename, 'w')

                                global_attr = {}
                                s = "%Y-%m-%d %H:%M:%S"
                                global_attr['date_created'] = datetime.now().strftime(s)
                                ncfile.setncatts(global_attr)

                                ncfile.createDimension('dim', None)
                            else:
                                ncfile = netCDF4.Dataset(filename, 'a')

                            index = len(ncfile.dimensions['dim'])
                            for field in result[key]:

                                if field in ncfile.variables.keys():
                                    var = ncfile.variables[field]
                                else:
                                    var_type = result[key][field].dtype
                                    kwargs = {'fill_value': -9999.0}
                                        # [martina]: was -99999 before, which made small sense to me
                                    # if dtype is a object the assumption is that the data is a
                                    # string
                                    if var_type == object:
                                        var_type = str
                                        kwargs = {}
                                    var = ncfile.createVariable(field, var_type,
                                                                'dim', **kwargs)
                                var[index:] = result[key][field]

                            ncfile.close()

                            return result
                    else:
                        return result
                else:
                    return result
            else:
                return result

# -------------------------------------------------------------------------------------
