# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import json

from lib_utils_info import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to setup configuration
def setup_config(file_path):
    if os.path.exists(file_path):
        with open(file_path) as file_handle:
            file_cfg = json.load(file_handle)
    else:
        alg_logger.error(' ===> Error in reading settings file "' + file_path + '"')
        raise IOError('File not found')
    return file_cfg
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to setup logging
def setup_logging(log_path):

    # manage log path
    log_folder, log_file = os.path.split(log_path)
    os.makedirs(log_folder, exist_ok=True)

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()  # Avoid duplicates if called multiple times

    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console (prompt) handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

# ----------------------------------------------------------------------------------------------------------------------
