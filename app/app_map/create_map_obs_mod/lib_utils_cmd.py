"""
Library Features:

Name:          lib_utils_cmd
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240209'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import subprocess

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Method to execute process
def exec_process(command_line=None, command_path=None):

    try:

        # Execute command-line
        os.chdir(command_path)
        process_handle = subprocess.Popen(
            command_line, shell=True,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Read standard output
        while True:
            std_out = process_handle.stdout.readline()
            if isinstance(std_out, bytes):
                std_out = std_out.decode('UTF-8')

            if std_out == '' and process_handle.poll() is not None:

                if process_handle.poll() == 0:
                    break
                else:
                    log_stream.error(' ===> Process failed!')
                    raise RuntimeError('Errors occurred in running process!')

            if std_out:
                std_out = str(std_out.strip())

        # Collect stdout and stderr and exitcode
        std_out, std_error = process_handle.communicate()
        std_exit = process_handle.poll()

        if std_out == b'' or std_out == '':
            std_out = None
        if std_error == b'' or std_error == '':
            std_error = None

        # Return variable(s)
        return std_out, std_error, std_exit

    except subprocess.CalledProcessError:
        # Exit code for process error
        log_stream.error(' ===> Process error in calling command line!')
        raise IOError('Errors occurred in running process!')

    except OSError:
        # Exit code for os error
        log_stream.error(' ===> Process error in running command line!')
        raise OSError('Errors occurred in running process!')

# ----------------------------------------------------------------------------------------------------------------------
