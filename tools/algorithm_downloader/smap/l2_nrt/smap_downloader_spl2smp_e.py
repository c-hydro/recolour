#!/usr/bin/python3
"""
HyDE Downloading Tool - SATELLITE SMAP (modified and extended NSIDC Data Download Script)

__date__ = '20200511'
__version__ = '1.0.2'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org)',
        'Andrea Libertino (andrea.libertino@cimafoundation.org)',
        'Martina Natali (martina01.natali@edu.unife.it)'
__library__ = 'HyDE'

If you wish, you may store your Earthdata username/password in a .netrc
file in your $HOME directory and the script will automatically attempt to
read this file. The .netrc file should have the following format:
   machine urs.earthdata.nasa.gov login myusername password mypassword
where 'myusername' and 'mypassword' are your Earthdata credentials.

General command line:
python smap_downloader_spl3smp_e.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20230929 (1.2.0) --> Download L2 product SPL2SMP_E (low timeliness, half-orbit product)
20230910 (1.1.1) --> Maintain original grid, subset on geo file
20230713 (1.1.0) --> Download L3 product SPL3SMP_E (high timeliness, global product)
20200511 (1.0.2) --> Fix bugs 
20200510 (1.0.1) --> Add multiprocessing mode and cleaning procedure(s) for ancillary and source file(s)
20200504 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
import base64
import itertools
import json
import netrc
import re
import os
import ssl
import sys
import math
import h5py

import pandas as pd
import numpy as np
import xarray as xr
import pyresample
import pyproj
import rasterio
from rasterio.crs import CRS

import cartopy
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import AxesGrid
from pygeogrids.grids import BasicGrid, genreg_grid

# import matplotlib.pylab as plt
import matplotlib.pyplot as plt

from multiprocessing import Pool, cpu_count
from contextlib import contextmanager
from argparse import ArgumentParser
from copy import deepcopy
from datetime import datetime
from osgeo import gdal, gdalconst
from osgeo.gdal import osr

from urllib.parse import urlparse
from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
from urllib.error import HTTPError, URLError

from lib_utils_time import set_data_time, set_run_time, set_run_time_modified
from lib_utils_cmr import get_credentials, build_version_query_params, build_cmr_query_url, cmr_request,\
    cmr_download_mp, cmr_download_seq, cmr_filter_urls, cmr_search
from lib_utils_log import set_logging
from lib_utils_io import clean_data_tmp, set_data_outcome, set_data_ancillary, set_data_source, collect_product_info
from lib_data_io import read_file_geo, read_grid_data
from lib_data_smap import read_smap_l2, process_smap_l2
from lib_data_io_tiff import read_file_tiff, write_file_tiff
from lib_utils_generic import fill_tags2string, create_darray_2d, make_folder, clean_data_tmp, set_data_outcome,\
    set_data_ancillary, set_data_source, collect_product_info
from lib_info_settings import read_file_json, get_args
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'HYDE DOWNLOADING TOOL - SATELLITE SMAP'
alg_version = '1.2.0'
alg_release = '2023-09-29'

# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to keep the working directory
@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to process cmr
def process_cmr(filename_obj_source, fileroot_obj_source,
                filename_obj_ancillary_global, filename_obj_ancillary_domain,
                filename_obj_outcome,
                grid_obj, geo_proj, geo_geotrans, grid_sm_2d,
                template_vars_data_list):

    filename_n = filename_obj_source.__len__()
    template_vars_data_list = [x.replace("__","_") if "__" in x else x for x in template_vars_data_list]
    template_vars_data_obj = template_vars_data_list * filename_n

    for (fn_source_time, fp_source), fr_source_list, fp_anc_global_list, fp_anc_domain_list, \
        fp_outcome_list, var_list in zip(
            filename_obj_source.items(), fileroot_obj_source.values(),
            filename_obj_ancillary_global.values(), filename_obj_ancillary_domain.values(),
            filename_obj_outcome.values(), template_vars_data_obj):

        if isinstance(fp_source, list):
            fp_source = fp_source[0]
        ff_source, fn_source = os.path.split(fp_source[0])
        logging.info(' ----> Process file: ' + fn_source + ' ... ')

        for fr_source_step, fp_anc_global_step, fp_anc_domain_step, fp_outcome_step, var_step in zip(
                fr_source_list, fp_anc_global_list, fp_anc_domain_list, fp_outcome_list, var_list):

            # Info
            logging.info(' -----> Process variable: ' + var_step + ' ... ')

            # Create folder(s) for ancillary and outcome file(s)
            ff_anc_global_step, fn_anc_global_step = os.path.split(fp_anc_global_step)
            make_folder(ff_anc_global_step)
            fp_anc_global_step_1 = fp_anc_global_step.replace('_global', '')
            ff_anc_domain_step, fn_anc_domain_step = os.path.split(fp_anc_domain_step)
            make_folder(ff_anc_domain_step)
            ff_outcome_step, fn_outcome_step = os.path.split(fp_outcome_step)
            make_folder(ff_outcome_step)

            # Open file and read it to build the dataset before translation
            logging.info(' ------> Translating and reproiecting over domain ... ')
            one_down, variable = os.path.split(fr_source_step)
            two_down, subgroup = os.path.split(one_down)

            # Reproject from global tiff file in epsg:6933 modified to domain tiff file in epsg:4326
            # bbox information are used to subset the destination image
            # while original resolution is preserved
            smap_data, smap_meta, smap_geo_x, smap_geo_y, smap_col, smap_row =\
                read_smap_l2(fp_source, subgroup, variable)

            smap_domain_2d = []
            for data, meta in zip(smap_data, smap_meta):
                data_2d = process_smap_l2(data, smap_geo_x, smap_geo_y, smap_col, smap_row, grid_obj, grid_sm_2d, meta)
                data_2d = np.nan_to_num(data_2d, nan=-9999.)
                smap_domain_2d.append(data_2d)

            smap_domain_high, smap_domain_wide = smap_domain_2d[0].shape

            if np.all(np.where(smap_domain_2d[0]==-9999., True, False)):
                logging.warning(' =======> No valid soil_moisture data over domain. The product will be discarded ... SKIPPING')
                os.remove(fp_source)
                continue

            logging.info(' ------> Translating and reproiecting over domain ... DONE')

            # Write domain tiff file in epsg:4326
            logging.info(' ------> Saving ' + fp_outcome_step + ' ... ')
            write_file_tiff(fp_outcome_step, smap_domain_2d, smap_meta, smap_domain_wide,
                            smap_domain_high, geo_geotrans, geo_proj, flip=True)
            logging.info(' ------> Saving ' + fp_outcome_step + ' ... DONE')

            # Info
            logging.info(' -----> Process variable: ' + var_step + ' ... DONE')

        logging.info(' ----> Process file: ' + fn_source + ' ... DONE')
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    # data_settings corresponds to the dictionary which is stored in the json
    # file, so refer to that for properly setting these data

    # Set algorithm logging
    make_folder(data_settings['data']['log']['folder'])
    set_logging(logger_file=os.path.join(data_settings['data']['log']['folder'],
                                         data_settings['data']['log']['filename']))
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
    # Get algorithm time range information
    # Time run is needed only for logging information
    time_run, time_run_range = set_run_time(alg_time, data_settings['time'])

    # Get algorithm geographical information
    # include only if polygon is provided as file, which must be a tif
    geo_settings = data_settings['data']['static']['geo_file']
    grid_name = os.path.join(geo_settings['folder'], geo_settings['filename'])

    # -----------------------------------------------------------------------
    # Static grid data reading
    grid_data, grid_attrs = read_grid_data(grid_name)
    geo_proj_obj, geo_geotrans_obj = grid_attrs['crs'], grid_attrs['transform']
    geo_geotrans = (geo_geotrans_obj.xoff, geo_geotrans_obj.a, 0, geo_geotrans_obj.yoff, 0, geo_geotrans_obj.e)
    geo_proj_wkt = geo_proj_obj.wkt

    grid_lons_1d, grid_lats_1d = grid_data['longitude'].values, np.flipud(grid_data['latitude'].values)

    # define grid mesh
    grid_lons_2d, grid_lats_2d = np.meshgrid(grid_lons_1d, grid_lats_1d)
    # define grid data 2d
    grid_sm_2d = np.zeros(shape=(grid_lons_2d.shape[0], grid_lats_2d.shape[1]))
    grid_sm_2d[:, :] = np.nan

    # grid obj using pyresample definition
    grid_obj = pyresample.geometry.GridDefinition(lats=grid_lats_2d, lons=grid_lons_2d)

    # -----------------------------------------------------------------------

    # Starting info
    logging.info(' --> TIME RUN: ' + str(time_run))

    # Iterate over time steps
    for time_run_step in time_run_range[1:]:

        # Starting info
        logging.info(' ---> TIME STEP: ' + str(time_run_step) + ' ... ')

        # Get data time range
        time_range, time_start, time_end = set_data_time(time_run_step,
                                                         data_settings['data'][
                                                             'dynamic'][
                                                             'time'])

        # Collect product information
        short_name_list, version_list, template_root_list,\
            template_vars_data_list, template_group_data_list,\
            bounding_box, cmr_url_list, urs_url_list, cmr_page_size_list,\
            url_filename_search_list, cmr_file_url_list, polygon_list,\
            filename_filter_list = collect_product_info(
            info_product=data_settings['product'],
            info_bbox=data_settings['data']['static']['bounding_box'],
            info_url=data_settings['data']['dynamic']['url'])

        # Iterate over product(s) type
        for file_id, (short_name_step, version_step, cmr_url_step, urs_url_step, cmr_page_size_step,
                      url_filename_search_step, cmr_file_url_step,
        polygon_step, filename_filter_step) in enumerate(
            zip(short_name_list, version_list, cmr_url_list, urs_url_list,
                cmr_page_size_list, url_filename_search_list,
                cmr_file_url_list, polygon_list, filename_filter_list)):
            # Retrieve url(s)
            if not url_filename_search_step:
                url_filename_search_step = cmr_search(
                    short_name_step, version_step, time_start, time_end,
                    bounding_box=bounding_box, polygon=polygon_step,
                    filename_filter=filename_filter_step,
                    cmr_page_size=cmr_page_size_step, cmr_url=cmr_url_step,
                    cmr_file_url=cmr_file_url_step)

            # Prepare folder(s) and filename(s)
            time_stamp, filename_list_url, filename_obj_source, fileroot_obj_source = set_data_source(
                file_id, url_filename_search_step,
                file_obj=data_settings['data']['dynamic']['source'],
                variable_obj=template_vars_data_list,
                root_obj=template_root_list,
                ancillary_obj=data_settings['algorithm']['ancillary'],
                template_obj=data_settings['algorithm']['template'],
                flag_cleaning_source=data_settings['algorithm']['flags']['cleaning_dynamic_data_source'])

            # Prepare ancillary filename(s)
            filename_obj_ancillary_global, filename_obj_ancillary_domain = set_data_ancillary(
                time_stamp, file_id,
                file_obj=data_settings['data']['dynamic']['ancillary'],
                variable_obj=template_vars_data_list,
                ancillary_obj=data_settings['algorithm']['ancillary'],
                template_obj=data_settings['algorithm']['template'],
                flag_cleaning_ancillary_global=data_settings['algorithm']['flags']['cleaning_dynamic_data_ancillary_global'],
                flag_cleaning_ancillary_domain=data_settings['algorithm']['flags']['cleaning_dynamic_data_ancillary_domain'])

            # Prepare outcome filename(s)
            filename_obj_outcome = set_data_outcome(
                time_stamp, file_id,
                file_obj=data_settings['data']['dynamic']['outcome'],
                variable_obj=template_vars_data_list,
                group_obj=template_group_data_list,
                ancillary_obj=data_settings['algorithm']['ancillary'],
                template_obj=data_settings['algorithm']['template'],
                flag_cleaning_outcome=data_settings['algorithm']['flags']['cleaning_dynamic_data_ancillary_global'])

            # Download file(s)
            if data_settings['algorithm']['flags']['downloading_mp']:
                cmr_download_mp(filename_list_url, filename_obj_source,
                                process_n=data_settings['algorithm']['ancillary']['process_mp'])
            else:
                cmr_download_seq(filename_list_url, filename_obj_source)

            # Process file(s)
            process_cmr(filename_obj_source, fileroot_obj_source,
                        filename_obj_ancillary_global,
                        filename_obj_ancillary_domain, filename_obj_outcome,
                        grid_obj, geo_proj_wkt, geo_geotrans, grid_sm_2d,
                        template_vars_data_list)

            # Clean files source and ancillary (if needed)
            clean_data_tmp(filename_obj_source, filename_obj_ancillary_global,
                           filename_obj_ancillary_domain, flag_cleaning_source=
                           data_settings['algorithm']['flags'][
                               'cleaning_dynamic_data_source'],
                           flag_cleaning_tmp=data_settings['algorithm']['flags']['cleaning_dynamic_data_tmp'])

        # Ending info
        logging.info(' ---> TIME STEP: ' + str(time_run_step) + ' ... DONE')


    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(
        ' ============================================================================ ')  # -------------------------------------------------------------------------------------


    # -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------
