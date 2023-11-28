"""
Library Features:

Name:          lib_data_io_h5
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20231110'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
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

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to read smap l2 data
def read_smap_l2(file_name, subgroup, variable):

    with h5py.File(file_name, mode='r') as f:

        data_list = []
        meta_list = []

        # Get original data in valid interval
        # is layer
        name = f'/{subgroup}/{variable}'
        data = f[name][:]
        units = f[name].attrs['units']
        units = units.decode('ascii', 'replace')
        long_name = f[name].attrs['long_name']
        long_name = long_name.decode('ascii', 'replace')
        _FillValue = f[name].attrs['_FillValue']
        valid_max = f[name].attrs['valid_max']
        valid_min = f[name].attrs['valid_min']
        invalid = np.logical_or(data > valid_max, data < valid_min)
        invalid = np.logical_or(invalid, data == _FillValue)
        data_raw = data.copy()
        data_raw[invalid] = np.nan
        data_meta = {'description': 'soil_moisture', 'units': '[cm**3cm**-3]'}

        # Get quality mask
        # quality mask itself is layer
        retrieval_qual_flag = f[f'/{subgroup}/retrieval_qual_flag'][:]
        qual_meta = {'description': 'retrieval_qual_flag', 'units': '[-]'}

        low_quality = np.logical_and(retrieval_qual_flag != 0, retrieval_qual_flag != 8)
        data_quality = data_raw.copy()
        data_quality[low_quality] = np.nan
        data_qual_meta = {'description': 'soil_moisture_qual', 'units': '[cm**3cm**-3]'}

        # Get bulk density to calculate porosity and get sm as sat degree
        # is layer
        bulk_density = f[f'/{subgroup}/bulk_density'][:]
        bulk_density[bulk_density==-9999.] = np.nan
        bulk_meta = {'description': 'bulk_density', 'units': 'gcm**-3'}

        porosity = 1. - bulk_density/2.65
        # data_quality = data_quality/porosity

        # Deepen layer with exponential filter for SWI
        # WIP WIP WIP WIP

        # Get vegetation data
        # is layer
        surf_temp = f[f'/{subgroup}/surface_temperature'][:]
        surf_temp_meta = {'description': 'surface_temperature', 'units':'K'}

        # Get vegetation data
        # is layer
        veg_water = f[f'/{subgroup}/vegetation_water_content'][:]
        veg_meta = {'description': 'vegetation_water_content', 'units':'kgm**-2'}

        # Get vegetation data
        # is layer
        veg_opacity = f[f'/{subgroup}/vegetation_opacity'][:]
        veg_opacity_meta = {'description': 'vegetation_opacity', 'units':'[-]'}


        # Build the layers. PAY ATTENTION TO THE ORDER, MUST BE EQUAL TO L3:
        # to ensure consistency, check the json file for l3
        data_list.append(data_raw);             meta_list.append(data_meta) # soil_moisture
        data_list.append(bulk_density);         meta_list.append(bulk_meta) # bulk_density
        data_list.append(retrieval_qual_flag.astype(np.float32));  meta_list.append(qual_meta) # retrieval_qual_flag
        data_list.append(surf_temp);            meta_list.append(surf_temp_meta)  # surface_temperature
        data_list.append(veg_water);            meta_list.append(veg_meta)  # vegetation_water_content
        data_list.append(veg_opacity);          meta_list.append(veg_opacity_meta)  # vegetation_opacity
        data_list.append(data_quality);         meta_list.append(data_qual_meta) # soil_moisture_qual

        # Get the geolocation data
        latitude = f[f'/{subgroup}/latitude'][:]
        longitude = f[f'/{subgroup}/longitude'][:]

        col = f[f'/{subgroup}/EASE_column_index'][:]
        row = f[f'/{subgroup}/EASE_row_index'][:]

        return data_list, meta_list, longitude, latitude, col, row

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to resample smap l2 data
def process_smap_l2(file_data, file_geo_x, file_geo_y, file_col, file_row, grid_obj, grid_sm_2d, file_metadata):

    file_reference = {
        'epsg': 6933, 'x_min': -17367530.45, 'y_max': 7314540.83, 'res': 9008.05, 'n_cols': 3856,
        'n_rows': 1624}
    file_proj = pyproj.Proj(file_reference['epsg'])

    file_x = file_reference['x_min'] + file_col * file_reference['res'] + file_reference['res'] / 2
    file_y = file_reference['y_max'] - file_row * file_reference['res'] - file_reference['res'] / 2

    # get data, lons and lats
    lons_raw, lats_raw = file_proj(file_x, file_y, inverse=True)
    idx_finite = np.argwhere(np.isfinite(file_data))[:, 0]

    # select data finite
    values_finite = file_data[idx_finite]
    lons_finite = lons_raw[idx_finite]
    lats_finite = lats_raw[idx_finite]

    file_obj = pyresample.geometry.SwathDefinition(lons=lons_finite, lats=lats_finite)

    # if there are no finite values, use raw data (will be empy band)
    if len(values_finite)<1:
        return np.full_like(grid_sm_2d, -9999.)

    # join point index to grid index
    values_finite_resampled = \
        pyresample.kd_tree.resample_nearest(
        source_geo_def=file_obj, data=values_finite, target_geo_def=grid_obj, radius_of_influence=9000, fill_value=np.nan)


    return values_finite_resampled


# -------------------------------------------------------------------------------------
