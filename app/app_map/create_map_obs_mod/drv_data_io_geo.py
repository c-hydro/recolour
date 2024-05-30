"""
Class Features

Name:          drv_data_io_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20190729'
Version:       '1.0.0'
"""

#######################################################################################
# Library
import logging

from os import makedirs, remove
from os.path import join, split, exists

from lib_info_args import logger_name
from lib_data_io_geo import select_grid

# Log
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Class to compute geographical data
class DataGeo:

    # -------------------------------------------------------------------------------------
    # Method class initialization
    def __init__(self, file_grid_global, file_grid_domain, file_domain, file_updating=False):

        # -------------------------------------------------------------------------------------
        # Get path(s) and filename(s)
        filepath_grid_global, filename_grid_global = split(file_grid_global)
        filepath_grid_domain, filename_grid_domain = split(file_grid_domain)
        filepath_domain, filename_domain = split(file_domain)

        # Store information in global workspace
        self.filepath_grid_global = filepath_grid_global
        self.filename_grid_global = filename_grid_global
        self.filepath_grid_domain = filepath_grid_domain
        self.filename_grid_domain = filename_grid_domain
        self.filepath_domain = filepath_domain
        self.filename_domain = filename_domain

        if not exists(join(self.filepath_grid_global, self.filename_grid_global)):
            log_stream.error(' =====> ERROR: file ' + join(self.filepath_grid_global, self.filename_grid_global) +
                             ' does not exist!', exc_info=True)
            raise RuntimeError('File grid global does not exist')
        if not exists(join(self.filepath_domain, self.filename_domain)):
            log_stream.error(' =====> ERROR: file ' + join(self.filepath_domain, self.filename_domain) +
                             ' does not exist!', exc_info=True)
            raise RuntimeError('File grid domain does not exist')

        if not exists(self.filepath_grid_domain):
            makedirs(self.filepath_grid_domain)

        if file_updating:
            if exists(join(self.filepath_grid_domain, self.filename_grid_domain)):
                remove(join(self.filepath_grid_domain, self.filename_grid_domain))
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get geographical data
    def getDataGeo(self):

        log_stream.info(' ---> Reference grid ... ')
        if not exists(join(self.filepath_grid_domain, self.filename_grid_domain)):
            select_grid(join(self.filepath_grid_global, self.filename_grid_global),
                        join(self.filepath_grid_domain, self.filename_grid_domain),
                        join(self.filepath_domain, self.filename_domain))
            log_stream.info(' ---> Reference grid ... DONE')
        else:
            log_stream.info(' ---> Reference grid ... PREVIOUSLY CREATED')

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
