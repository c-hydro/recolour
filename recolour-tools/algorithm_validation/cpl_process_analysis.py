"""
Class Features

Name:          cpl_process_analysis
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np

from pytesmo.validation_framework.validation import Validation

from lib_utils_generic import get_dataset_modes, get_dataset_names
from cpl_process_metrics import BasicMetrics, ExtendedMetrics
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class coupler process analysis
class CplAnalysis:

    # method to initialize class
    def __init__(self, dset_interfaces, dset_metrics_obj, time_period, time_window_fraction=0.5):

        self.dset_interfaces = dset_interfaces
        self.dset_metrics_obj = dset_metrics_obj
        self.time_period = time_period
        self.time_window_fraction = time_window_fraction

        # TODO: check if we want to introduce scaling
        self.scaling = None

        # set datasets reference
        self.dset_name_reference = get_dataset_modes(self.dset_interfaces, dset_mode='reference')
        # set datasets names
        self.dset_name_list = get_dataset_names(self.dset_name_reference, self.dset_interfaces)
        self.dset_n = self.dset_name_list.__len__()

    # method to setup process analysis
    def setup_analysis(self, process_analysis_cell=1394):

        # create process metrics obj
        process_metrics_obj = self.__create_process_metrics()

        # ********************************************************* #
        # set process analysis obj                                  #
        # ********************************************************* #
        process_analysis_obj = Validation(
            datasets=self.dset_interfaces,
            spatial_ref=self.dset_name_reference,
            metrics_calculators=process_metrics_obj,
            temporal_window=self.time_window_fraction,
            scaling=self.scaling,
            period=self.time_period,
        )
        # ********************************************************* #

        # set process analysis jobs
        process_analysis_jobs = self.create_process_jobs(process_analysis_obj, process_analysis_cell)

        return process_analysis_obj, process_analysis_jobs

    # method to create process jobs obj
    def create_process_jobs(self, process_obj, process_cell):

        # get process analysis jobs
        # process_jobs = process_obj.get_processing_jobs()

        # get reference reader
        process_dset_reader_reference = process_obj.data_manager.datasets[self.dset_name_reference]['class']

        # get reference datasets gpi, lons and lats
        cell_gpi, cell_lons, cell_lats = process_dset_reader_reference.grid.grid_points_for_cell(process_cell)

        # define reference datasets points
        cell_points = np.zeros((len(cell_gpi), 3))
        cell_points[:, 0], cell_points[:, 1], cell_points[:, 2] = cell_gpi, cell_lons, cell_lats

        # define jobs
        process_jobs = []
        for point in cell_points:
            process_jobs.append((int(point[0]), point[1], point[2]))

        return process_jobs

    # method to create process metrics obj
    def __create_process_metrics(self):

        if self.dset_n == 3:
            process_metrics_obj = {(3, 3): self.dset_metrics_obj.calc_metrics}
        elif self.dset_n == 2:
            process_metrics_obj = {(2, 2): self.dset_metrics_obj.calc_metrics}
        else:
            logging.error(' ===> Datasets elements expected can be 2 or 3. '
                          'Actually tha algorithm found "' + str(self.dset_n) + '" datasets')
            raise RuntimeError('Check your settings to set correctly the datasets')

        return process_metrics_obj


# -------------------------------------------------------------------------------------
