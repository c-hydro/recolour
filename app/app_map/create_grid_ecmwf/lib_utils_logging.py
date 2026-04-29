"""
Library Features:

Name:          lib_utils_logging
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260429'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# method to set logging information
def set_logging(
        logger_name='algorithm_logger',
        logger_folder=None,
        logger_file='log.txt',
        logger_format=None,
        logger_level=logging.INFO,
        logger_mode='w'):

    import os
    import sys
    import logging
    from copy import deepcopy

    if logger_format is None:
        logger_format = deepcopy(logger_format_def)

    if logger_file is None:
        logger_file = 'log.txt'

    # define logger path
    if logger_folder is not None:
        logger_path = os.path.join(logger_folder, logger_file)
    else:
        logger_path = logger_file

    # if only filename is given, save in script folder
    if os.path.dirname(logger_path) == '':
        logger_folder_name = os.path.dirname(os.path.abspath(sys.argv[0]))
        logger_path = os.path.join(logger_folder_name, logger_path)
    else:
        logger_path = os.path.abspath(logger_path)

    # create log folder
    os.makedirs(os.path.dirname(logger_path), exist_ok=True)

    # reset logging configuration
    logging.basicConfig(level=logger_level, force=True)

    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.setLevel(logger_level)

    # create named logger
    logger = logging.getLogger(logger_name)
    logger.handlers.clear()
    logger.setLevel(logger_level)
    logger.propagate = False

    formatter = logging.Formatter(logger_format)

    # file handler
    file_handler = logging.FileHandler(logger_path, mode=logger_mode)
    file_handler.setLevel(logger_level)
    file_handler.setFormatter(formatter)

    # console handler
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logger_level)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    # suppress noisy third-party debug logs
    logging.getLogger("findlibs").setLevel(logging.WARNING)
    logging.getLogger("gribapi").setLevel(logging.WARNING)
    logging.getLogger("eccodes").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    # monkey patch logging module
    logging.info = logger.info
    logging.warning = logger.warning
    logging.error = logger.error
    logging.debug = logger.debug
    logging.critical = logger.critical
# -------------------------------------------------------------------------------------
