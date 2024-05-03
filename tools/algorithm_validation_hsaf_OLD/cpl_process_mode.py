"""
Class Features

Name:          cpl_process_mode
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import multiprocessing as mp

from lib_utils_generic import slice_list
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class coupler process mode
class CplMode:

    # method to initialize class
    def __init__(self, mode_collections,
                 tag_mode_type='mp_mode', tag_mode_cpu='mp_cpu_n',
                 tag_mode_parallel='parallel', tag_mode_seq='sequential'):

        self.mode_type = mode_collections[tag_mode_type]

        if tag_mode_cpu in list(mode_collections.keys()):
            self.mode_cpu = mode_collections[tag_mode_cpu]
        else:
            self.mode_cpu = mp.cpu_count() - 1

        self.tag_mode_parallel = tag_mode_parallel
        self.tag_mode_seq = tag_mode_seq

    # method to setup process mode
    def setup_mode(self, cells=1394, file_log_groups='file_log_group.txt'):

        # format cells obj
        if not isinstance(cells, list):
            cells = [cells]

        # method to select mode execution
        mode_execution = self.__select_mode_execution()
        # method to select mode cpu(s)
        mode_cpu = self.__select_mode_cpu(mode_execution)

        # get group cells
        mode_groups = self.get_group_cells(cells, mode_cpu)
        # write groups cell to file
        self.write_group_cells(mode_groups, execution_mode=mode_execution, file_groups=file_log_groups)

        return mode_groups, mode_execution, mode_cpu

    # method to select mode cpu(s)
    def __select_mode_cpu(self, mode_execution):

        if mode_execution == self.tag_mode_seq:
            cpu_n = 1
        elif mode_execution == self.tag_mode_parallel:
            cpu_n = self.mode_cpu
        else:
            logging.error()
            raise RuntimeError()

        if cpu_n < 1:
            cpu_n = 1

        return cpu_n

    # method to select mode execution
    def __select_mode_execution(self):
        mode_type = self.mode_type
        if mode_type:
            mode_execution = self.tag_mode_parallel
        else:
            mode_execution = self.tag_mode_seq
        return mode_execution

    # method to save group cells
    @staticmethod
    def write_group_cells(dset_groups, execution_mode='NA', file_groups='file_groups.txt',
                          file_clean=True,
                          file_header=' ### process groups ### \n',
                          file_line_default=' group: {group_idx} -- cells: {group_data} \n'):

        # define file group
        file_groups = file_groups.format(execution_mode=execution_mode)

        # clean old file group
        if file_clean:
            if os.path.exists(file_groups):
                os.remove(file_groups)

        # open and write new file group
        with open(file_groups, 'w') as file_handle:
            file_handle.write(file_header)
            for group_idx, group_data in enumerate(dset_groups):
                file_line_step = file_line_default.format(group_idx=group_idx, group_data=group_data)
                file_handle.write(file_line_step)

    # method to divide cells in groups
    @staticmethod
    def get_group_cells(obj_cells, mp_cpu=1):

        if not isinstance(obj_cells, list):
            obj_cells = [obj_cells]

        n_cells = obj_cells.__len__()
        n_groups = n_cells / mp_cpu

        if n_groups < 1:
            n_groups = 1

        n_elem = int(np.floor(n_cells / n_groups))

        # debug
        # obj_cells, n_elem = [1394, 1395, 1396], 2

        obj_groups = slice_list(obj_cells, n_elem)

        return obj_groups

# -------------------------------------------------------------------------------------
