#!/usr/bin/env python3
"""
RECOLOUR APPS - SOIL_GRIDS DOWNLOADER

General command line:
python soil_grids_downloader_layers.py -settings_file configuration.json
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import sys
import argparse
import json
import time
import logging

from pathlib import Path
from datetime import datetime

import rasterio
from rasterio.warp import transform_bounds
from owslib.wcs import WebCoverageService
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for downloading soil_grids files'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2026-05-21'

# get logger
alg_logger = logging.getLogger('app_downloader')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
## basic utils
# helper to configure logging
def get_logger(settings):

    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    log_folder = log_settings.get("folder", "./log")
    log_filename = log_settings.get("filename", "downloader.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = datetime.now().strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    alg_logger.handlers = []
    alg_logger.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    fh = logging.FileHandler(log_path)
    fh.setLevel(level)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)

    alg_logger.addHandler(fh)
    alg_logger.addHandler(ch)
# helper to get args
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        type=str,
        required=True
    )
    return parser.parse_args()

# helper to read file json
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f' ===> File "{file_name}" not found. Exit')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# i/o utils
# method to read geotiff
def read_bbox_from_geotiff(file_path):
    with rasterio.open(file_path) as src:
        bbox = transform_bounds(
            src.crs,
            "EPSG:4326",
            *src.bounds,
            densify_pts=21
        )

    return {"xmin": bbox[0], "ymin": bbox[1], "xmax": bbox[2], "ymax": bbox[3],}

# helper to apply bbox
def apply_bbox_buffer(bbox, buffer_deg):
    return {
        "xmin": bbox["xmin"] - buffer_deg,
        "ymin": bbox["ymin"] - buffer_deg,
        "xmax": bbox["xmax"] + buffer_deg,
        "ymax": bbox["ymax"] + buffer_deg,
    }

# helper to convert bbox
def bbox_to_tuple(bbox):
    return (
        bbox["xmin"],
        bbox["ymin"],
        bbox["xmax"],
        bbox["ymax"],
    )

# helper to check geotiff
def is_geotiff(content):
    return content[:4] in [b"II*\x00", b"MM\x00*"]

# helper to get bbox
def get_domain_bbox(folder_name: str, file_name: str, bbox_buffer_deg: float = 0.0) -> tuple[dict, dict]:

    path_file = os.path.join(folder_name, file_name)

    bbox_raw = read_bbox_from_geotiff(path_file)
    bbox_def = apply_bbox_buffer(bbox_raw, bbox_buffer_deg)

    return bbox_raw, bbox_def
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to download data
def download_soilgrids(config, bbox):

    # get info
    folder_name = config["soilgrids"]["folder_name"]
    os.makedirs(folder_name, exist_ok=True)
    file_template = config["soilgrids"].get("file_name", "{variable}_{depth}_{statistic}.tif")
    overwrite = int(config["soilgrids"].get("overwrite", True))

    # variables info
    template_url = config["soilgrids"]["service_url"]
    variables = config["soilgrids"]["variables"]
    depths = config["soilgrids"]["depths"]
    statistic = config["soilgrids"]["statistic"]
    width = int(config["soilgrids"].get("width", 1024))
    height = int(config["soilgrids"].get("height", 1024))
    wcs_version = config["soilgrids"].get("wcs_version", "1.0.0")
    output_format = config["soilgrids"].get("format", "GEOTIFF_INT16")

    bbox_wcs = bbox_to_tuple(bbox)

    # iterate over variables
    for variable in variables:

        # info variable start
        alg_logger.info(f' ----> Variable {variable} ... ')

        # service url
        service_url = template_url.format(variable=variable)
        # service wcs
        service_wcs = WebCoverageService(service_url,version=wcs_version)

        # info connection
        alg_logger.info(" ::: Connecting to %s", service_url)

        # iterate over depth(s)
        for depth in depths:

            # info depth start
            alg_logger.info(f' -----> Depth {depth} ... ')

            # define covarage and output name
            coverage_id = f"{variable}_{depth}_{statistic}"

            # define path name
            file_name = file_template.format(variable=variable, depth=depth, statistic=statistic)
            path_name = os.path.join(folder_name, file_name)

            # check file availability
            if os.path.exists(path_name):
                if overwrite:
                    alg_logger.info(f" ::: File: {path_name} REMOVED")
                    os.remove(path_name)
                else:
                    alg_logger.info(f" ::: File: {path_name} FOUND")
                    continue

            if coverage_id not in service_wcs.contents:
                alg_logger.warning("Coverage not found: %s", coverage_id)
                continue

            # info download
            alg_logger.info(" ::: Download %s", coverage_id)

            # get service response
            service_response = service_wcs.getCoverage(
                identifier=coverage_id, bbox=bbox_wcs, crs="EPSG:4326", format=output_format,
                width=width,height=height
            )
            # get service content
            service_content = service_response.read()

            # check geotiff format
            if not is_geotiff(service_content):

                service_preview = service_content[:1000].decode("utf-8", errors="ignore")

                alg_logger.error(" ::: Server response is not a GeoTIFF")
                alg_logger.error(" ::: Response preview:\n%s", service_preview)
                raise RuntimeError(f" ::: Downloaded file is not a GeoTIFF: {coverage_id}")

            # write file
            with open(path_name, "wb") as file_handle:
                file_handle.write(service_content)

            # info save
            alg_logger.info(" ::: Save: %s", path_name)

            # info depth end
            alg_logger.info(f' -----> Depth {depth} ... DONE')

        # info variable end
        alg_logger.info(f' ----> Variable {variable} ... DONE')

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get algorithm args
    algorithm_args = get_args()
    # read settings
    algorithm_settings = read_file_json(algorithm_args.settings_file)
    # configure logger
    get_logger(algorithm_settings)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info start
    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    alg_logger.info(' ---> Settings file: ' + algorithm_args.settings_file)

    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get domain reference
    alg_logger.info(' ---> Get domain bbox ... ')
    bbox_raw, bbox_download = get_domain_bbox(
        folder_name= algorithm_settings["domain"]["folder_name"],
        file_name=algorithm_settings["domain"]["file_name"],
        bbox_buffer_deg=algorithm_settings["domain"]["bbox_buffer_deg"]
    )

    # info bbox
    alg_logger.info(" ::: BBox - reference domain: %s", bbox_raw)
    alg_logger.info(" ::: BBox - download datasets: %s", bbox_download)

    alg_logger.info(' ---> Get domain bbox ... DONE')

    # method to download data
    alg_logger.info(' ---> Download domain datasets ... ')
    download_soilgrids(algorithm_settings, bbox_download)
    alg_logger.info(' ---> Download domain datasets ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info end
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
# entrypoint
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
