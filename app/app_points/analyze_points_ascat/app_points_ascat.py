# app_points_ascat.py
"""
RECOLOUR APPS - COMPUTE METRICS - REprocess paCkage for sOiL mOistUre pRoducts

Application to read SSM points, plot them, and generate GeoTIFF.
Time is rounded to the nearest 00:00 or 12:00 on the date.

__date__ = '20250612'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

Usage:
    python app_points_ascat.py -settings_file config.json -time "YYYY-MM-DD HH:MM"

Version(s):
20250612 (1.0.0) --> Beta release
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import argparse
import os
import time
from datetime import datetime, timedelta

from lib_info_args import logger_name
from lib_info_logging import set_logging

from lib_utils_time import get_time, round_time
from lib_utils_io import read_json, read_data_grid, read_data_points, write_data_grid
from lib_utils_generic import build_filepath
from lib_utils_analysis import resample_points_to_grid
from lib_utils_plot import plot_data_points, plot_data_grid

# set logger
alg_logger = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for analyzing SSM points ascat and generating points and GeoTIFF maps'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2025-06-12'
# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm main
def main() -> None:

    # ------------------------------------------------------------------------------------------------------------------
    # method to get script argument(s)
    algorithm_args = get_args()
    # method to read configuration file
    algorithm_cfg = read_json(algorithm_args.settings_file)

    # set logging
    set_logging(logger_name=logger_name,
                logger_folder=algorithm_cfg['log']['folder_name'],
                logger_file=algorithm_cfg['log']['file_name'])
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (start)
    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    # time algorithm
    start_time = time.time()

    # info configuration start
    alg_logger.info(' ===> Get configurations ... ')

    # configuration settings parameters
    params_cfg = algorithm_cfg.get('parameters', None)
    if params_cfg is None:
        alg_logger.error(' ===> Parameters configuration is missing in the settings file')
        raise ValueError('Check parameters configuration in the settings file')

    # configuration settings geo
    geo_cfg = algorithm_cfg.get('geo', None)
    if geo_cfg is None:
        alg_logger.error(' ===> Geo configuration is missing in the settings file')
        raise ValueError('Check geo configuration in the settings file')

    # configuration settings data
    data_cfg_root = algorithm_cfg.get('data', None)
    if data_cfg_root is None:
        alg_logger.error(' ===> Datasets root configuration is missing in the settings file')
        raise ValueError('Check datasets root configuration in the settings file')
    data_cfg_src = data_cfg_root.get('source', None)
    if data_cfg_src is None:
        alg_logger.error(' ===> Datasets source configuration is missing in the settings file')
        raise ValueError('Check datasets source configuration in the settings file')
    # configuration settings destination
    data_cfg_dst = data_cfg_root.get('destination', None)
    if data_cfg_dst is None:
        alg_logger.error(' ===> Datasets destination configuration is missing in the settings file')
        raise ValueError('Check datasets destination configuration in the settings file')

    # configuration settings plot
    plot_cfg_root = algorithm_cfg.get('plot', None)
    if plot_cfg_root is None:
        alg_logger.error(' ===> Plot configuration root is missing in the settings file')
        raise ValueError('Check plot configuration root in the settings file')
    plot_cfg_point = plot_cfg_root.get('point', None)
    if plot_cfg_point is None:
        alg_logger.error(' ===> Plot configuration point is missing in the settings file')
        raise ValueError('Check plot configuration point in the settings file')
    plot_cfg_map = plot_cfg_root.get('map', None)
    if plot_cfg_map is None:
        alg_logger.error(' ===> Plot configuration map is missing in the settings file')
        raise ValueError('Check plot configuration map in the settings file')

    # info configuration end
    alg_logger.info(' ===> Get configurations ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info time start
    alg_logger.info(' ===> Get time ... ')

    # time settings
    dt_run = get_time(algorithm_args.time, format='%Y-%m-%d %H:%M')
    dt_rounded = round_time(dt_run)
    alg_logger.info(f" ===> Original time: {dt_run}, rounded to {dt_rounded}")

    # info time end
    alg_logger.info(' ===> Get time ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info data source start
    alg_logger.info(' ===> Get geo ... ')

    # define geo
    folder_name_geo, file_name_geo = geo_cfg['folder_name'], geo_cfg['file_name']
    path_name_geo = build_filepath(folder_name_geo, file_name_geo, dt=dt_rounded)
    # get geo darray
    geo_da = read_data_grid(path_name_geo)
    geo_attrs = geo_da.attrs

    # info data source start
    alg_logger.info(' ===> Get geo  ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info data source start
    alg_logger.info(' ===> Get data ... ')

    # define csv source
    folder_name_csv, file_name_csv = data_cfg_src['folder_name'], data_cfg_src['file_name']
    path_name_csv = build_filepath(folder_name_csv, file_name_csv, dt=dt_rounded)
    # get csv points
    points_lons, points_lats, points_data, points_dframe = read_data_points(path_name_csv)

    # info data source start
    alg_logger.info(' ===> Get data  ... DONE')

    # info data resampling start
    alg_logger.info(' ===> Resample data ... ')

    # apply resampling points to grid
    grid_da = resample_points_to_grid(
        points_dframe, geo_da['longitude'], geo_da['latitude'],
        var_name_data='sm', var_name_geo_x='lon', var_name_geo_y='lat',
        resampling_max_distance= params_cfg.get('resampling_max_distance', 18000.0),
        resampling_neighbours=params_cfg.get('resampling_neighbours', 8),
        resampling_method=params_cfg.get('resampling_method', 'nearest'),
        filtering_active=params_cfg.get('filtering_active', True),
    )
    # extract grid data and coordinates
    grid_data, grid_lons, grid_lats = grid_da.values, grid_da['longitude'].values, grid_da['latitude'].values

    # info data resampling end
    alg_logger.info(' ===> Resample data ... DONE')

    # info data destination start
    alg_logger.info(' ===> Dump data ... ')

    # define tif destination
    folder_name_tif, file_name_tif = data_cfg_dst['folder_name'], data_cfg_dst['file_name']
    path_name_tif = build_filepath(folder_name_tif, file_name_tif, dt=dt_rounded)
    os.makedirs(os.path.dirname(path_name_tif), exist_ok=True)

    # write grid data to GeoTIFF
    grid_height, grid_width = grid_data.shape
    grid_proj = geo_attrs.get('proj', 'EPSG:4326')
    grid_transform = geo_attrs.get('transform', None)

    write_data_grid(file_name=path_name_tif, file_data=grid_data,
                    file_wide=grid_width, file_high=grid_height,
                    file_proj=grid_proj, file_transform=grid_transform)

    # info data destination end
    alg_logger.info(' ===> Dump data ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info plot points start
    alg_logger.info(' ===> Plot points ... ')

    # define destination points
    folder_name_points, file_name_points = plot_cfg_point['folder_name'], plot_cfg_point['file_name']
    path_name_points = build_filepath(folder_name_points, file_name_points, dt=dt_rounded)
    os.makedirs(os.path.dirname(path_name_points), exist_ok=True)

    # plot raw points
    plot_data_points(
        points_lons, points_lats, points_data,
        output_path=path_name_points,
        vmin=plot_cfg_point.get('vmin', 0.0),
        vmax=plot_cfg_point.get('vmax', 100.0),
        cmap=plot_cfg_point.get('cmap', 'viridis'),
        marker_size=plot_cfg_point.get('marker_size', 10)
    )

    # info plot points end
    alg_logger.info(' ===> Plot points ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info plot map start
    alg_logger.info(' ===> Plot map ... ')

    # define destination map
    folder_name_map, file_name_map = plot_cfg_map['folder_name'], plot_cfg_map['file_name']
    path_name_map = build_filepath(folder_name_map, file_name_map, dt=dt_rounded)
    os.makedirs(os.path.dirname(path_name_map), exist_ok=True)

    # plot grid data
    plot_data_grid(
        grid_lons, grid_lats, grid_data,
        output_path=path_name_map,
        cmap=plot_cfg_map.get('cmap', 'viridis'),
        vmin=plot_cfg_map.get('vmin', 0.0),
        vmax=plot_cfg_map.get('vmax', 100.0)
    )

    # info plot map end
    alg_logger.info(' ===> Plot map ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(' ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    alg_logger.info(' ==> ... END')
    alg_logger.info(' ==> Bye, Bye')
    alg_logger.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to get script argument(s)
def get_args():
    parser = argparse.ArgumentParser(description='Analyze SSM points ASCAT')
    parser.add_argument('-settings_file', required=True,
                        help='Path to JSON config file')
    parser.add_argument('-time', required=True,
                        help='Datetime string e.g. "2025-06-10 12:43"')
    args = parser.parse_args()

    return args
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# call main function
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------------------------------
