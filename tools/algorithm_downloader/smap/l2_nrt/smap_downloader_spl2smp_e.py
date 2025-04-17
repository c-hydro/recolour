#!/usr/bin/python3
"""
RECOLOUR Downloading Tool - SATELLITE SMAP (modified and extended NSIDC Data Download Script)

__date__ = '20250410'
__version__ = '1.3.0'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org)',
        'Andrea Libertino (andrea.libertino@cimafoundation.org)',
        'Martina Natali (martina01.natali@edu.unife.it)'
__library__ = 'RECOLOUR'

If you wish, you may store your Earthdata username/password in a .netrc
file in your $HOME directory and the script will automatically attempt to
read this file. The .netrc file should have the following format:
   machine urs.earthdata.nasa.gov login myusername password mypassword
where 'myusername' and 'mypassword' are your Earthdata credentials.

General command line:
python smap_downloader_spl2smp_e.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20250410 (1.3.0) --> Add new version of the downloader
20230929 (1.2.0) --> Download L2 product SPL2SMP_E (low timeliness, half-orbit product)
20230910 (1.1.1) --> Maintain original grid, subset on geo file
20230713 (1.1.0) --> Download L3 product SPL3SMP_E (high timeliness, global product)
20200511 (1.0.2) --> Fix bugs 
20200510 (1.0.1) --> Add multiprocessing mode and cleaning procedure(s) for ancillary and source file(s)
20200504 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import time
import os

from lib_info_args import logger_name, logger_format, time_format_algorithm
from lib_info_settings import get_args, read_file_json, parse_file_json

from lib_utils_time import set_data_time, set_run_time

from lib_utils_log import set_logging
from lib_utils_geo import set_geo_info

from lib_utils_wget import create_file_url, search_file_url, create_file_wget, download_file_url
from lib_utils_wget import get_credentials

from lib_processing_data import collect_data_info, set_data_destination, set_data_ancillary, set_data_source
from lib_processing_fx import reproj_data

# logger stream
logger_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'RECOLOUR DOWNLOADING TOOL - SATELLITE SMAP'
alg_version = '1.3.0'
alg_release = '2025-04-10'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'
# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# script main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # read algorithm settings
    data_settings = read_file_json(alg_settings)
    # parse algorithm settings
    (log_settings, lock_settings, time_settings, product_settings,
     alg_ancillary_settings, alg_template_settings, alg_flags_settings,
     data_static_settings, data_dynamic_settings) = parse_file_json(data_settings)

    # Set algorithm logging
    os.makedirs(log_settings['folder'], exist_ok=True)
    set_logging(logger_file=os.path.join(log_settings['folder'], log_settings['filename']),
                logger_format=logger_format)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Info algorithm
    logger_stream.info(' ============================================================================ ')
    logger_stream.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logger_stream.info(' ==> START ... ')
    logger_stream.info(' ')

    # Time algorithm information
    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get user and password
    user, password = get_credentials()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get the time information
    time_run, time_range = set_run_time(alg_time, time_settings)
    # get the geo information
    geo_grid, geo_proj_wkt, geo_geotrans, grid_sm_2d = set_geo_info(
        file_name=os.path.join(data_static_settings['geo_file']['folder'], data_static_settings['geo_file']['filename']))
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # iterate over time steps
    for time_step in time_range[1:]:

        # info time start
        logger_stream.info(' ----> TIME RUN: "' + str(time_run) + '" TIME STEP: "' + str(time_step) + '" ... ')

        # get data time information
        # time_range, time_start, time_end = set_data_time(time_step, data_dynamic_settings['time'])

        # Collect product information
        (short_name_list, version_list, template_root_list,\
            template_vars_data_list, template_group_data_list,\
            bounding_box,
            remote_url, remote_folder) = collect_data_info(
            info_time=time_step,
            info_product=product_settings,
            info_bbox=data_static_settings['bounding_box'],
            info_url=data_dynamic_settings['url'])

        # define remote url
        remote_url = create_file_url(remote_url=remote_url, remote_folder=remote_folder)
        # search remote file
        remote_file_found = search_file_url(
            remote_url,
            file_lock=os.path.join(lock_settings['folder'], lock_settings['filename']),
            remove_lock=True)

        # iterate over product ids
        for product_id, (short_name_step, version_step) in enumerate(zip(short_name_list, version_list)):

            # info product start
            logger_stream.info(' -----> PRODUCT "' + str(short_name_step) + ' ' + str(version_step) + '" ... ')

            # iterate over file urls
            for remote_file_url in remote_file_found:

                # info url start
                logger_stream.info(' ------> FILE URL "' + str(remote_file_url) + '" ... ')

                # prepare source object
                data_source_obj = set_data_source(
                    product_id, remote_file_url,
                    file_obj=data_dynamic_settings['source'],
                    variable_obj=template_vars_data_list, root_obj=template_root_list,
                    ancillary_obj=alg_ancillary_settings, template_obj=alg_template_settings,
                    flag_remove_source=alg_flags_settings['remove_dynamic_data_source'])

                # prepare ancillary object
                data_ancillary_obj = set_data_ancillary(
                    product_id,
                    source_obj=data_source_obj,
                    file_obj=data_dynamic_settings['ancillary'],
                    variable_obj=template_vars_data_list,
                    ancillary_obj=alg_ancillary_settings, template_obj=alg_template_settings,
                    flag_remove_ancillary_global=alg_flags_settings['remove_dynamic_data_ancillary_global'],
                    flag_remove_ancillary_domain=alg_flags_settings['remove_dynamic_data_ancillary_domain'])

                # prepare destination object
                data_destination_obj = set_data_destination(
                    product_id,
                    source_obj=data_source_obj,
                    file_obj=data_dynamic_settings['destination'],
                    variable_obj=template_vars_data_list, group_obj=template_group_data_list,
                    ancillary_obj=alg_ancillary_settings, template_obj=alg_template_settings,
                    flag_remove_destination=alg_flags_settings['remove_dynamic_data_destination'])

                # define wget command
                file_name_url, folder_name_save = data_source_obj['file_url'], data_source_obj['folder_name']
                wget_cmd = create_file_wget(user=user, password=password, save_path=folder_name_save)

                # execute wget command
                file_path_found = download_file_url(wget_cmd, file_name_url, path_download=folder_name_save)
                folder_path_found, file_name_source = os.path.split(file_path_found)
                # info file found
                logger_stream.info(' ==== FOLDER SOURCE: "' + str(folder_path_found) + '" ')
                logger_stream.info(' ==== FILE SOURCE: "' + str(file_path_found) + '" ')

                # execute data translate and reproject
                reproj_data(
                    data_source_obj, data_ancillary_obj, data_destination_obj,
                    geo_grid, geo_proj_wkt, geo_geotrans, grid_sm_2d,
                    flag_remove_data_out_of_domain=alg_flags_settings['remove_data_out_of_domain'])

                # info url end
                logger_stream.info(' ------> FILE URL "' + str(remote_file_url) + '" ... DONE')

            # info product end
            logger_stream.info(' -----> PRODUCT "' + str(short_name_step) + ' ' + str(version_step) + '" ... DONE')

        # info time end
        logger_stream.info(' ----> TIME RUN: "' + str(time_run) + '" TIME STEP: "' + str(time_step) + '" ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logger_stream.info(' ')
    logger_stream.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logger_stream.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logger_stream.info(' ==> ... END')
    logger_stream.info(' ==> Bye, Bye')
    logger_stream.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Call script from external library
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------------------------------
