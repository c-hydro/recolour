#!/usr/bin/python3
"""
RECOLOUR APP - XML Renderer  - REprocess paCkage for sOiL mOistUre pRoducts
__date__ = '20231110'
__version__ = '1.0.0'
__author__ = 'Fabio Delogu (fabio.delogu@cimafoundation.org'

__library__ = 'recolour'

General command line:
python3 app_xml_renderer.py -product file_name.tiff -xml_in file_name_in.xml -xml_out file_name_out.xml

Version(s):
20231110 (1.0.0) --> Beta release
"""

# ----------------------------------------------------------------------------------------------------------------------
# Complete library
import logging
import os
import time
import json
from argparse import ArgumentParser
from copy import deepcopy

from lib_utils_io_tiff import read_file_tiff
from lib_utils_io_xml import read_file_xml, update_file_xml, write_file_xml

from lib_utils_envs import organize_envs_obj
from lib_utils_generic import (get_geographical_bounds, extract_time_from_filename, extract_parts_from_filename,
                               get_time_xml, get_time_product, convert_time_stamp_to_string, select_filename)

# default logger information
logger_name = 'app_xml_renderer_logger'
logger_file = 'app_xml_renderer.txt'
logger_handle = 'file'  # 'file' or 'stream'
logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                '%(message)-80s %(filename)s:[%(lineno)-6s - %(funcName)-20s()] '
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for rendering xml'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2023-11-10'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # get algorithm settings
    file_path_product, file_path_xml_renderer, file_path_xml_in, file_path_xml_out = get_args()
    # read file renderer
    xml_renderer = read_file_settings(file_path_xml_renderer)

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
    # get product
    logging.info(' ----> Get product ... ')
    if file_path_product is not None:
        if os.path.exists(file_path_product):
            file_data, file_attrs, file_lons, file_lats, common_attrs = read_file_tiff(file_path_product)
        else:
            logging.error(' ===> File "' + file_path_product + '" does not exists')
            raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    else:
        logging.error(' ===> File product is not defined by the user')
        raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    logging.info(' ----> Get product ... DONE')

    # get xml
    logging.info(' ----> Get xml ... ')
    if file_path_xml_in is not None:
        if os.path.exists(file_path_xml_in):
            xml_data_root, xml_tree = read_file_xml(file_path_xml_in)
        else:
            logging.error(' ===> File "' + file_path_xml_in + '" does not exists')
            raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    else:
        logging.error(' ===> File xml is not defined by the user')
        raise FileNotFoundError('File is mandatory to correctly run the algorithm')
    logging.info(' ----> Get xml ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # collect information start
    logging.info(' ----> Collect information ... ')

    # get xml template, keys and values
    xml_flags = xml_renderer['flag']
    xml_pattern, xml_template = xml_renderer['pattern'], xml_renderer['template']
    xml_keys, xml_values = xml_renderer['xml_keys'], xml_renderer['xml_values']

    # collect information end
    logging.info(' ----> Collect information ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # organize information start
    logging.info(' ----> Organize information ... ')

    # get geographical information
    lon_west, lon_east, lat_south, lat_north = get_geographical_bounds(file_lons, file_lats)

    # get time product begin and end
    time_stamp_product_begin, time_stamp_product_end = extract_time_from_filename(
        file_path_product,
        time_format=xml_template['time_format_file_product'],
        time_pattern=xml_pattern['time_format_file_product'],
        time_sep=xml_pattern['time_separator_file_product'])

    time_product_begin_file = convert_time_stamp_to_string(
        time_stamp_product_begin, time_format=xml_template['time_format_period_begin_file'])
    time_product_end_file = convert_time_stamp_to_string(
        time_stamp_product_end, time_format=xml_template['time_format_period_end_file'])

    time_product_begin_tag = convert_time_stamp_to_string(
        time_stamp_product_begin, time_format=xml_template['time_format_period_begin_tag'])
    time_product_end_tag = convert_time_stamp_to_string(
        time_stamp_product_end, time_format=xml_template['time_format_period_end_tag'])

    # organize file name tags
    file_name_tags = {'time_format_period_begin_file': time_product_begin_file,
                      'time_format_period_end_file': time_product_end_file}

    # select product filename
    folder_name_product, file_name_product_args = extract_parts_from_filename(file_path_product)
    file_name_product_default = xml_template['name_format_product']
    file_name_product_filled = file_name_product_default.format(**file_name_tags)

    file_name_product = select_filename(
        file_name_product_args, file_name_product_filled,
        flag_expected=xml_flags['use_expected_file_name_product'], flag_mandatory=False)
    file_path_product = os.path.join(folder_name_product, file_name_product)

    # select xml filename
    folder_name_xml, file_name_xml_out_args = extract_parts_from_filename(file_path_xml_out)
    file_name_xml_default = xml_template['name_format_xml']
    file_name_xml_out_filled = file_name_xml_default.format(**file_name_tags)

    file_name_xml_out = select_filename(
        file_name_xml_out_args, file_name_xml_out_filled,
        flag_expected=xml_flags['use_expected_file_name_xml'], flag_mandatory=False)
    file_path_xml_out = os.path.join(folder_name_xml, file_name_xml_out)

    # get time creation xml
    time_creation_xml = get_time_xml(time_format=xml_template['time_format_creation_xml'])
    # get time creation and publication product
    time_creation_product = get_time_product(
        file_path_product, time_format=xml_template['time_format_creation_product'])
    time_publication_product = get_time_product(
        file_path_product, time_format=xml_template['time_format_publication_product'])

    # organize envs obj
    xml_envs = organize_envs_obj(
        filename_product=file_name_product, filename_xml=file_name_xml_out,
        creation_time_product=time_creation_product, publication_time_product=time_publication_product,
        creation_time_xml=time_creation_xml,
        lon_east=lon_east, lon_west=lon_west, lat_south=lat_south, lat_north=lat_north,
        begin_time_product=time_product_begin_tag, end_time_product=time_product_end_tag)

    # organize information end
    logging.info(' ----> Organize information ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # update xml
    logging.info(' ----> Update xml ... ')

    # update xml keys and values
    xml_data = update_file_xml(xml_data_root,
                               xml_keys=xml_keys, xml_values=xml_values, xml_envs=xml_envs,
                               xml_field_error_update=xml_flags['xml_field_not_update_error'],
                               xml_field_not_declared_text=xml_flags['xml_field_not_declared_text'],
                               xml_field_not_declared_attributes=xml_flags['xml_field_not_declared_attributes'])

    logging.info(' ----> Update xml ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # write xml
    logging.info(' ----> Write xml ... ')

    folder_name_out, file_name_out = os.path.split(file_path_xml_out)
    os.makedirs(folder_name_out, exist_ok=True)
    write_file_xml(file_path_xml_out, xml_data)

    logging.info(' ----> Write xml ... DONE')
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

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file settings
def read_file_settings(file_name):
    if os.path.exists(file_name):
        with open(file_name) as file_handle:
            data_settings = json.load(file_handle)
    else:
        logging.error(' ===> Error in reading settings file "' + file_name + '"')
        raise IOError('File not found')
    return data_settings
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set logging information
def set_logging(logger_name='algorithm_logger', logger_folder=None,
                logger_file='log.txt', logger_format=None):

    if logger_format is None:
        logger_format = deepcopy(logger_format)
    if logger_file is None:
        logger_file = deepcopy(logger_file)

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

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():

    # parser argument(s)
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-product', action="store", dest="product")
    parser_handle.add_argument('-xml_renderer', action="store", dest="xml_renderer")
    parser_handle.add_argument('-xml_in', action="store", dest="xml_in")
    parser_handle.add_argument('-xml_out', action="store", dest="xml_out")
    parser_values = parser_handle.parse_args()

    file_path_product, file_path_xml_renderer = None, 'app_xml_renderer.json'
    file_path_xml_in, file_path_xml_out = None, 'out.xml'
    if parser_values.product:
        file_path_product = parser_values.product
    if parser_values.xml_renderer:
        file_path_xml_renderer = parser_values.xml_renderer
    if parser_values.xml_in:
        file_path_xml_in = parser_values.xml_in
    if parser_values.xml_out:
        file_path_xml_out = parser_values.xml_out

    return file_path_product, file_path_xml_renderer, file_path_xml_in, file_path_xml_out

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Call script from external library
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------------------------------
