#!/usr/bin/python3
"""
RECOLOUR APP - TIME-SERIES GRID TIFF 2 CSV  - REprocess paCkage for sOiL mOistUre pRoducts
__date__ = '20231016'
__version__ = '1.0.0'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org',
        'Francesco Silvestro (francesco.silvestro@cimafoundation.org)'

__library__ = 'SM'

General command line:
python3 app_sm_ts_analysis.py -settings_file configuration.json

Version(s):
20231016 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
from argparse import ArgumentParser

import pandas as pd

import os

# default logger information
logger_name = 'app_join_ts_logger'
logger_file = 'app_join_ts.txt'
logger_handle = 'file'  # 'file' or 'stream'
logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                '%(message)-80s %(filename)s:[%(lineno)-6s - %(funcName)-20s()] '
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for joining time-series k1, k2, k3'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2023-10-27'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # get algorithm settings
    file_path_k1, file_path_k2, file_path_k3, file_path_out = get_args()

    # set logging
    set_logging(logger_name=logger_name, logger_format=logger_format,
                logger_folder=os.getcwd(),
                logger_file=logger_file)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ')

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # get ts 1
    logging.info(' ----> Get time-series k1 ... ')
    if file_path_k1 is not None:
        if os.path.exists(file_path_k1):
            file_data_k1 = read_file_csv(file_path_k1)
        else:
            logging.error(' ===> File "' + file_path_k1 + '" does not exists')
            raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    else:
        logging.error(' ===> File k1 is not defined by the user')
        raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    logging.info(' ----> Get time-series k1 ... DONE')

    # get ts 2
    logging.info(' ----> Get time-series k2 ... ')
    if file_path_k2 is not None:
        if os.path.exists(file_path_k2):
            file_data_k2 = read_file_csv(file_path_k2)
        else:
            logging.error(' ===> File "' + file_path_k2 + '" does not exists')
            raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    else:
        logging.error(' ===> File k2 is not defined by the user')
        raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    logging.info(' ----> Get time-series k2 ... DONE')

    # get ts 3
    logging.info(' ----> Get time-series k3 ... ')
    if file_path_k3 is not None:
        if os.path.exists(file_path_k3):
            file_data_k3 = read_file_csv(file_path_k3)
        else:
            logging.warning(' ===> File "' + file_path_k2 + '" does not exists')
            file_data_k3 = None
    else:
        logging.warning(' ===> File k3 is not defined by the user')
        file_data_k3 = None
        logging.info(' ----> Get time-series k3 ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # join ts
    logging.info(' ----> Join time-series k1-k2-k3 ... ')
    if (file_data_k1 is not None) and (file_data_k2 is not None) and (file_data_k3 is not None):
        file_data_out = pd.concat([file_data_k1, file_data_k2, file_data_k3])
    elif (file_data_k1 is not None) and (file_data_k2 is not None) and (file_data_k3 is None):
        file_data_out = pd.concat([file_data_k1, file_data_k2])
    else:
        logging.error(' ===> Datasets formats are not supported in this mode')
        raise NotImplemented('Case not implemented yet')

    # remove dups keeping the first occurrence
    file_data_out = file_data_out[~file_data_out.index.duplicated(keep='first')]

    # adjust the dataframe
    file_data_out = file_data_out.reset_index()
    file_data_out = file_data_out.set_index('time')
    file_data_out.index.name = 'time'

    logging.info(' ----> Join time-series k1-k2-k3 ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # write ts
    logging.info(' ----> Write time-series k1-k2-k3 ... ')
    folder_name_out, file_name_out = os.path.split(file_path_out)
    os.makedirs(folder_name_out, exist_ok=True)

    write_file_csv(file_path_out, file_data_out,
                   dframe_index=True, dframe_float_format='%.3f')
    logging.info(' ----> Write time-series k1-k2-k3 ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file csv
def read_file_csv(file_name, dframe_index='time', dframe_date_format='%Y-%m-%d %H:%M',
                  dframe_sep=',', dframe_decimal='.', dframe_float_precision='legacy'):

    file_dframe = pd.read_csv(file_name, decimal=dframe_decimal, sep=dframe_sep,
                              float_precision=dframe_float_precision)

    if dframe_index not in list(file_dframe.columns):
        logging.error(' ===> Index column "' + dframe_index +
                         '"  must be available in the source dataframe. Check the source file')
        raise RuntimeError('Including the index column in the source file for skipping this error.')

    file_dframe[dframe_index] = pd.DatetimeIndex(file_dframe[dframe_index].values).strftime(dframe_date_format)
    file_dframe[dframe_index] = pd.DatetimeIndex(file_dframe[dframe_index])

    file_dframe = file_dframe.reset_index()
    if 'index' in list(file_dframe.columns):
        file_dframe = file_dframe.drop(['index'], axis=1)
    file_dframe = file_dframe.set_index(dframe_index)
    file_dframe.sort_index(inplace=True)

    return file_dframe

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write file csv
def write_file_csv(file_name, file_dframe,
                   dframe_sep=',', dframe_decimal='.', dframe_float_format='%.2f',
                   dframe_index=False, dframe_header=True,
                   dframe_index_label='time'):

    file_dframe.to_csv(
        file_name, mode='w',
        index=dframe_index, sep=dframe_sep, decimal=dframe_decimal,
        index_label=dframe_index_label,
        header=dframe_header, float_format=dframe_float_format,  quotechar='"')

# ----------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_name='algorithm_logger', logger_folder=None,
                logger_file='log.txt', logger_format=None):

    if logger_format is None:
        logger_format = deepcopy(logger_format_def)
    if logger_file is None:
        logger_file = deepcopy(logger_file_def)

    if logger_folder is not None:
        logger_path = os.path.join(logger_folder, logger_file)
    else:
        logger_path = deepcopy(logger_file)

    logger_loc = os.path.split(logger_path)
    if logger_loc[0] == '' or logger_loc[0] == "":
        logger_folder_name, logger_file_name = os.path.dirname(os.path.abspath(sys.argv[0])), logger_loc[1]
    else:
        logger_folder_name, logger_file_name = logger_loc[0], logger_loc[1];

    os.makedirs(logger_folder_name, exist_ok=True)

    # define logger path
    logger_path = os.path.join(logger_folder_name, logger_file_name)

    # Remove old logging file
    if os.path.exists(logger_path):
        os.remove(logger_path)

    # Open logger
    logging.getLogger(logger_name)
    logging.root.setLevel(logging.DEBUG)

    # Open logging basic configuration
    logging.basicConfig(level=logging.DEBUG, format=logger_format, filename=logger_file, filemode='w')

    # Set logger handle
    logger_handle_1 = logging.FileHandler(logger_path, 'w')
    logger_handle_2 = logging.StreamHandler()
    # Set logger level
    logger_handle_1.setLevel(logging.DEBUG)
    logger_handle_2.setLevel(logging.DEBUG)
    # Set logger formatter
    logger_formatter = logging.Formatter(logger_format)
    logger_handle_1.setFormatter(logger_formatter)
    logger_handle_2.setFormatter(logger_formatter)

    # Add handle to logging
    logging.getLogger('').addHandler(logger_handle_1)
    logging.getLogger('').addHandler(logger_handle_2)

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-k1', action="store", dest="k1")
    parser_handle.add_argument('-k2', action="store", dest="k2")
    parser_handle.add_argument('-k3', action="store", dest="k3")
    parser_handle.add_argument('-out', action="store", dest="out")
    parser_values = parser_handle.parse_args()

    file_path_k1, file_path_k2, file_path_k3, file_path_out = None, None, None, 'out.csv'
    if parser_values.k1:
        file_path_k1 = parser_values.k1
    if parser_values.k2:
        file_path_k2 = parser_values.k2
    if parser_values.k3:
        file_path_k3 = parser_values.k3
    if parser_values.out:
        file_path_out = parser_values.out

    return file_path_k1, file_path_k2, file_path_k3, file_path_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Call script from external library
if __name__ == '__main__':

    main()
# ----------------------------------------------------------------------------------------------------------------------
