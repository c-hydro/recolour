"""
Library Features:

Name:          lib_info_settings
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
# Method to read file json
# json file contains a dict of relevant data for the execution of the script
# such as all the product's characteristics
# Path to the json file must be stored in the alg_settings argument from the
# running script or configuration file

def read_file_json(file_name):
    env_ws = {}
    for env_item, env_value in os.environ.items():
        env_ws[env_item] = env_value

    with open(file_name, "r") as file_handle:
        json_block = []
        for file_row in file_handle:

            for env_key, env_value in env_ws.items():
                env_tag = '$' + env_key
                if env_tag in file_row:
                    env_value = env_value.strip("'\\'")
                    file_row = file_row.replace(env_tag, env_value)
                    file_row = file_row.replace('//', '/')

            # Add the line to our JSON block
            json_block.append(file_row)

            # Check whether we closed our JSON block
            if file_row.startswith('}'):
                # Do something with the JSON dictionary
                json_dict = json.loads(''.join(json_block))
                # Start a new block
                json_block = []

    return json_dict

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get script argument(s)
# Arguments and the execution of the script are defined in a bash file but
# could also be defined in the configuration file (in Pycharm)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-time', action="store", dest="alg_time")
    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    if parser_values.alg_time:
        alg_time = parser_values.alg_time
    else:
        alg_time = None

    return alg_settings, alg_time

# -------------------------------------------------------------------------------------

