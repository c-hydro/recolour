"""
Class Features

Name:          cpl_process_time
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
from datetime import datetime
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class coupler process time
class CplTime:

    # method to initialize class
    def __init__(self, time_collections, time_format='%Y-%m-%d', time_window_hours=12):

        self.time_collections = time_collections
        self.time_format = time_format
        self.time_window_hours = time_window_hours

    # method to setup process time
    def setup_time(self):

        time_start = self.__check_key(self.time_collections, 'time_start')
        time_end = self.__check_key(self.time_collections, 'time_end')
        time_format = self.__check_key(self.time_collections, 'time_format')

        time_period = [datetime.strptime(time_start, time_format), datetime.strptime(time_end, time_format)]
        time_window_fraction = float(self.time_window_hours) / 24.

        return time_period, time_window_fraction

    # method to check key availability
    @staticmethod
    def __check_key(dset_obj, dset_key='time_start'):
        if dset_key in list(dset_obj.keys()):
            dset_value = dset_obj[dset_key]
        else:
            logging.error(' ===> Dataset key "' + dset_key + '" is not available in the time obj')
            raise RuntimeError('The key is needed by the algorithm. Please check your settings file')
        return dset_value
# -------------------------------------------------------------------------------------
