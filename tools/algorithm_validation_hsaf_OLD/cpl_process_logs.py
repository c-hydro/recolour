"""
Class Features

Name:          cpl_process_logs
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import tempfile

from lib_utils_generic import make_folder
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class coupler process logs
class CplLog:

    # method to initialize class
    def __init__(self, log_collections,
                 tag_path_log='path_log', tag_file_log='file_log'):

        self.log_collections = log_collections

        self.path_log = None
        if tag_path_log in list(log_collections.keys()):
            self.path_log = log_collections[tag_path_log]
        if (self.path_log == "") or (self.path_log is None):
            self.path_log = tempfile.mkdtemp()
            logging.warning(' ===> Log folder is null or defined by NoneType. Create temporary folder.'
                            ' Log will be saved in "' + self.path_log + '"')
        self.file_log_generic = 'file_log_generic.txt'
        if tag_file_log in list(log_collections.keys()):
            self.file_log_generic = log_collections[tag_file_log]
        if (self.file_log_generic == "") or (self.file_log_generic is None):
            temp_obj = tempfile.NamedTemporaryFile(suffix='.txt')
            _, self.file_log_generic = os.path.split(temp_obj.name)
            logging.warning(' ===> Log file is null or defined by NoneType. Create temporary file.'
                            ' Log file will be named as "' + self.file_log_generic + '"')

        self.file_log_groups = 'file_log_groups_{execution_mode}.txt'
        self.file_log_process = 'file_log_process_{execution_mode}_group_{group_n}_cell_{cell_n}.txt'

    # method to setup process logs
    def setup_logs(self):

        # define file log(s)
        file_path_generic = os.path.join(self.path_log, self.file_log_generic)
        file_path_groups = os.path.join(self.path_log, self.file_log_groups)
        file_path_process = os.path.join(self.path_log, self.file_log_process)

        # make analysis folder
        make_folder(self.path_log)

        return file_path_generic, file_path_groups, file_path_process

# -------------------------------------------------------------------------------------
