"""
RECOLOUR APP - PROCESSING TOOL - SOIL MOISTURE PRODUCT ASCAT OBS DR [ASCAT - H16, H103]

__date__ = '20240202'
__version__ = '5.1.0'
__author__ = 'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'recolour'

General command line:
python3 app_sm_obs_dr_ascat.py -settingfile configuration_product.json

Version:
20240201 (5.1.0) --> Refactor to recolour package
20190805 (5.0.0) --> Beta release using pytesmo library and time-series data format
20161007 (4.0.1) --> Fix bug(s) and update code(s)
20160606 (4.0.0) --> First release 4.0.0
20140708 (3.0.8) --> Last release based on previously H07 and H14 script(s) versions 1-2-3
"""

# ----------------------------------------------------------------------------------------------------------------------
# Library
import time
import argparse

from lib_info_logging import setLoggingFile

from lib_info_args import logger_formatter, logger_handle, logger_name

from drv_configuration_time import DataTime
from drv_data_io_geo import DataGeo
from drv_configuration_algorithm import DataAlgorithm

from drv_data_io_obs import DataProductTime, DataProductWrapper
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # Version and algorithm information
    project_name = 'Recolour'
    alg_version = '5.1.0'
    alg_type = 'Processing Tool'
    alg_name = 'ASCAT OBS - Data Record'
    # Time algorithm information
    alg_start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get script argument(s)
    [script_name, file_settings, time_arg] = get_args()

    # Set algorithm configuration
    driver_algorithm_settings = DataAlgorithm(file_settings)
    [data_algorithm_settings, data_algorithm_path,
     data_algorithm_flags, data_algorithm_colormap] = driver_algorithm_settings.getDataSettings()

    # Set logging file
    log_stream = setLoggingFile(data_algorithm_path['log'], logger_name, logger_handle, logger_formatter)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Start Program
    log_stream.info('[' + project_name + ' ' + alg_type + ' - ' + alg_name + ' (Version ' + alg_version + ')]')
    log_stream.info('[' + project_name + '] Start Program ... ')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get data time
    log_stream.info(' --> Set algorithm time ... ')
    driver_algorithm_time = DataTime(
        time_arg,
        time_settings=data_algorithm_settings['time']['time_now'],
        time_period_past=int(data_algorithm_settings['time']['time_period']),
        time_frequency=data_algorithm_settings['time']['time_frequency'],
        time_rounding=data_algorithm_settings['time']['time_rounding'])
    data_algorithm_time = driver_algorithm_time.getDataTime()
    log_stream.info(' --> Set algorithm time ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get data geo
    log_stream.info(' --> Set grid domain ... ')
    data_algorithm_geo = DataGeo(
        data_algorithm_path['file_grid_ref_global'], data_algorithm_path['file_grid_ref_domain'],
        data_algorithm_path['file_domain'], file_updating=data_algorithm_flags['cleaning_static_data'])
    data_algorithm_geo.getDataGeo()
    log_stream.info(' --> Set grid domain ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Iterate over time steps
    for time_step in data_algorithm_time['time_steps']:

        # -------------------------------------------------------------------------------------
        # Compute dynamic time
        log_stream.info(' --> Compute analysis time ... ')
        driver_dynamic_time = DataProductTime(
            time_step=time_step,
            time_run=data_algorithm_time['time_run'],
            time_settings=data_algorithm_settings['data']['dynamic']['time'])
        data_dynamic_time = driver_dynamic_time.computeDataTime()
        log_stream.info(' --> Compute analysis time ... DONE')
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Compute dynamic product
        log_stream.info(' --> Compute data record product ... ')
        driver_dynamic_product = DataProductWrapper(
            var_in=data_algorithm_settings['variables']['input']['sm_data'],
            var_out=data_algorithm_settings['variables']['outcome']['sm_data'],
            settings=data_algorithm_settings,
            flags=data_algorithm_flags,
            time=data_dynamic_time,
            sm_domain=data_algorithm_path['file_domain'],
            sm_grid_ref_global=data_algorithm_path['file_grid_ref_global'],
            sm_grid_ref_domain = data_algorithm_path['file_grid_ref_domain'],
            rzsm_grid=data_algorithm_path['file_grid_mod_rzsm'],
            sm_data_in_ref=[
                data_algorithm_path['sm_data_in_ref1'],
                data_algorithm_path['sm_data_in_ref2']
            ],
            sm_data_raw_ref=[
                data_algorithm_path['sm_data_raw_ref1'],
                data_algorithm_path['sm_data_raw_ref2']
            ],
            sm_data_out_ref=[
                data_algorithm_path['sm_data_out_ref1'],
                data_algorithm_path['sm_data_out_ref2']
            ],
            sm_data_out=data_algorithm_path['sm_data_out'],
            sm_data_analysis=data_algorithm_path['sm_data_analysis'],
            sm_data_points=data_algorithm_path['sm_data_points'],
            sm_data_maps=data_algorithm_path['sm_data_maps'],
            sm_data_dr_global=data_algorithm_path['sm_data_dr_global'],
            sm_data_dr_domain=data_algorithm_path['sm_data_dr_domain'],
            rzsm_data_dr=data_algorithm_path['rzsm_data_dr'],
            sm_updating_out_idx=data_algorithm_flags['cleaning_dynamic_out_idx'],
            sm_updating_out_contiguous=data_algorithm_flags['cleaning_dynamic_out_contiguous'],
            sm_updating_analysis=data_algorithm_flags['cleaning_dynamic_analysis'],
            sm_updating_points=data_algorithm_flags['cleaning_dynamic_points'],
            sm_updating_maps=data_algorithm_flags['cleaning_dynamic_maps'],
        )

        # Method to resample data product (from swath(s) to indexed time-series)
        driver_dynamic_product.resampleDataProduct(dr=True)
        # Method to create data product
        driver_dynamic_product.execDataRecord()
        log_stream.info(' --> Compute data record product ... DONE')
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Clean ancillary archive file(s)
        driver_dynamic_product.cleanDataProduct(
            sm_cleaning_out=data_algorithm_flags['cleaning_dynamic_tmp'],
            sm_cleaning_analysis=data_algorithm_flags['cleaning_dynamic_tmp'],
            sm_cleaning_points=data_algorithm_flags['cleaning_dynamic_tmp'])
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Note about script parameter(s)
    log_stream.info('NOTE - Algorithm parameter(s)')
    log_stream.info('Script: ' + str(script_name))
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # End Program
    alg_time_elapsed = round(time.time() - alg_start_time, 1)

    log_stream.info('[' + project_name + ' ' + alg_type + ' - ' + alg_name + ' (Version ' + alg_version + ')]')
    log_stream.info('End Program - Time elapsed: ' + str(alg_time_elapsed) + ' seconds')
    # ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get script argument(s)
def get_args():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('-settings_file', action="store", dest="settings_file")
    args_parser.add_argument('-time', action="store", dest="time_arg")
    args_values = args_parser.parse_args()

    script_name = args_parser.prog

    settings_file = 'configuration.json'
    if 'settings_file' in args_values:
        settings_file = args_values.settings_file
    time_arg = None
    if 'time_arg' in args_values:
        time_arg = args_values.time_arg

    return script_name, settings_file, time_arg
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
