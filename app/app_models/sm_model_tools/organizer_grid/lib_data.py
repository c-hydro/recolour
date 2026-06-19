"""
Library Features:

Name:          lib_data
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260618'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations

import logging
import os
import numpy as np

from lib_utils_geo import read_file_grid
from lib_utils_io import read_file_registry, collect_data, write_file_geotiff
from lib_utils_analysis import interpolate_points2grid
from config_info import LOGGER_NAME, VALUE_NODATA_DEFAULT

# logger stream
logger_stream = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to execute data process
def process(settings, time_reference):

    # ------------------------------------------------------------------------------------------------------------------
    # start message - script
    logger_stream.info(" ----> Process execution ... ")

    # get registry info
    folder_name_reg = settings['data_static']["registry"]["folder"]
    file_name_reg = settings['data_static']["registry"]["filename"]
    path_name_reg = os.path.join(folder_name_reg, file_name_reg)

    params_reg = {
        "delimiter": settings['data_static']["registry"].get("delimiter", ","),
        "tag_col": settings['data_static']["registry"].get("tag_col", "tag"),
        "lon_col": settings['data_static']["registry"].get("lon_col", "longitude"),
        "lat_col": settings['data_static']["registry"].get("lat_col", "latitude"),
        "valid_col": settings['data_static']["registry"].get("valid_col", "valid"),
        "use_only_valid": bool(settings['data_static']["registry"].get("use_only_valid", True)),
    }

    # get grid info
    folder_name_grid = settings['data_static']["grid"]["folder"]
    file_name_grid = settings['data_static']["grid"]["filename"]
    path_name_grid = os.path.join(folder_name_grid, file_name_grid)

    params_grid = {
        "band": settings['data_static']["grid"].get("band", 1),
        "no_data": settings['data_static']["grid"].get("no_data", -9999)
    }

    # get source info
    folder_name_src = settings['data_dynamic']["source"]["folder"]
    file_name_src = settings['data_dynamic']["source"]["filename"]
    path_name_src = os.path.join(folder_name_src, file_name_src)

    params_src = {
        "time_format": settings['data_dynamic']["source"].get("time_format", "%Y-%m-%d %H:%M"),
        "step_hours": float( settings['data_dynamic']["source"].get("step_hours", 1)),
        "max_previous_steps": int(settings['data_dynamic']["source"].get("max_previous_steps", 10)),
        "nodata_values": settings['data_dynamic']["source"].get("nodata_values",[VALUE_NODATA_DEFAULT, -9998]),
        "valid_min": settings['data_dynamic']["source"].get("valid_min", 0),
        "valid_max": settings['data_dynamic']["source"].get("valid_max", 1),
        "delimiter": settings['data_dynamic']["source"].get("delimiter", ","),
        "time_col": settings['data_dynamic']["source"].get("time_col", "time"),
        "value_col": settings['data_dynamic']["source"].get("value_col", "theta_simulated"),
        "registry_lon_col": settings['data_dynamic']["source"].get("registry_lon_col", "longitude"),
        "registry_lat_col": settings['data_dynamic']["source"].get("registry_lat_col", "latitude"),
    }

    # get destination info
    folder_name_dst = settings['data_dynamic']["destination"]["folder"]
    file_name_dst = settings['data_dynamic']["destination"]["filename"]
    path_name_dst = os.path.join(folder_name_dst, file_name_dst)

    # get interpolation info
    params_interp = {
        "no_data": settings['interpolation'].get("no_data", VALUE_NODATA_DEFAULT),
        "method": settings['interpolation'].get("method", "idw"),
        "roi_m": settings['interpolation'].get("roi_m", 10000),
        "min_points": settings['interpolation'].get("min_points", 1),
        "max_points": settings['interpolation'].get("max_points", 8),
        "idw_power": settings['interpolation'].get("idw_power", 8),
    }
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # read file registry
    logger_stream.info(f' ----> Read registry {path_name_reg} ... ')
    registry_obj = read_file_registry(path_name_reg, params_reg)
    logger_stream.info(f' ----> Read registry {path_name_reg} ... DONE')

    # read file grid
    logger_stream.info(f' ----> Read grid {path_name_grid} ... ')
    grid_mask, grid_lon, grid_lat, grid_profile = read_file_grid(path_name_grid, params_grid)
    logger_stream.info(f' ----> Read grid {path_name_grid} ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # collect data points
    logger_stream.info(f' ----> Collect data points ... ')
    points_data, points_times, points_lags = collect_data(path_name_src, registry_obj, time_reference, params_src)

    # check valid points (for interpolation)
    if len(points_data) == 0:
        logger_stream.warning(f' ====> No valid point values available for interpolation')
        logger_stream.info(f' ----> Collect data points ... FAILED')
        return None
    else:
        logger_stream.info(f' ----> Valid points: {len(points_data)}')
        logger_stream.info(f' ----> Collect data points ... DONE')

    # interpolate data points
    logger_stream.info(f' ----> Interpolate data points ... ')
    grid_data = interpolate_points2grid(points_data, grid_lon, grid_lat, grid_mask, params_interp)
    logger_stream.info(f' ----> Interpolate data points ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message - script
    logger_stream.info(" ----> Process execution ... DONE")

    return grid_data, grid_mask, grid_profile
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to save data
def save(grid_data, grid_mask, grid_profile, settings, time_reference):

    # ------------------------------------------------------------------------------------------------------------------
    # start message - script
    logger_stream.info(" ----> Dumping execution ... ")
    logger_stream.info(f" ----> Reference time: {time_reference}")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get destination info
    folder_name_dst = settings['data_dynamic']["destination"]["folder"]
    file_name_dst = settings['data_dynamic']["destination"]["filename"]
    path_name_dst = os.path.join(folder_name_dst, file_name_dst)

    name_dst = settings['data_dynamic']["destination"].get("name", "soil_moisture"),
    compress_dst = settings['data_dynamic']["destination"].get("compress", "deflate")
    no_data_dst = settings['data_dynamic']["destination"].get("no_data", np.nan)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # check no data
    if no_data_dst is None:
        no_data_dst = np.nan

    # define path destination
    path_name_dst = path_name_dst.format(
        time_now=time_reference,
        var_name=name_dst,
        yyyy=time_reference.strftime("%Y"),
        mm=time_reference.strftime("%m"),
        dd=time_reference.strftime("%d"),
        yyyymmddhh=time_reference.strftime("%Y%m%d%H"),
    )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # dump data
    write_file_geotiff(path_name_dst, grid_data, grid_profile, no_data_dst, compress=compress_dst)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message - script
    logger_stream.info(" ----> Dumping execution ... DONE")
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------