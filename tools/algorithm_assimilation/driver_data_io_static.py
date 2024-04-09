"""
Class Features

Name:          driver_data_io_static
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20200515'
Version:       '1.0.0'
"""

######################################################################################
# Library
import logging
import os
import time
import numpy as np
import pandas as pd

from copy import deepcopy

import shapely
from shapely.geometry import shape, Point
from scipy.spatial import cKDTree

from multiprocessing import Pool, cpu_count

from lib_data_io_ascii import read_file_raster
from lib_data_io_nc import read_file_geo
from lib_data_io_generic import write_obj, read_obj
from lib_utils_system import make_folder

from lib_utils_geo import create_map_idx

# Debug
# import matplotlib.pylab as plt
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverStatic
class DriverStatic:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, src_dict, dst_dict=None, parameters_dict=None,
                 alg_ancillary=None, alg_template_tags=None,
                 flag_datasets_terrain='terrain', flag_datasets_cn='cn',
                 flag_datasets_channel_network='channel_network',
                 flag_datasets_sm='soil_moisture', flag_datasets_map='map',
                 flag_cleaning_static=True):

        self.flag_datasets_terrain = flag_datasets_terrain
        self.flag_datasets_cn = flag_datasets_cn
        self.flag_datasets_channel_network = flag_datasets_channel_network
        self.flag_datasets_sm = flag_datasets_sm
        self.flag_datasets_map = flag_datasets_map

        self.alg_ancillary = alg_ancillary

        self.alg_template_tags = alg_template_tags
        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'

        self.folder_name_terrain = src_dict[self.flag_datasets_terrain][self.folder_name_tag]
        self.file_name_terrain = src_dict[self.flag_datasets_terrain][self.file_name_tag]
        self.file_path_terrain = os.path.join(self.folder_name_terrain, self.file_name_terrain)

        self.folder_name_cn = src_dict[self.flag_datasets_cn][self.folder_name_tag]
        self.file_name_cn = src_dict[self.flag_datasets_cn][self.file_name_tag]
        self.file_path_cn = os.path.join(self.folder_name_cn, self.file_name_cn)

        self.folder_name_channel_network = src_dict[self.flag_datasets_channel_network][self.folder_name_tag]
        self.file_name_channel_network = src_dict[self.flag_datasets_channel_network][self.file_name_tag]
        self.file_path_channel_network = os.path.join(self.folder_name_channel_network, self.file_name_channel_network)

        self.folder_name_sm = src_dict[self.flag_datasets_sm][self.folder_name_tag]
        self.file_name_sm = src_dict[self.flag_datasets_sm][self.file_name_tag]
        self.file_path_sm = os.path.join(self.folder_name_sm, self.file_name_sm)

        self.folder_name_dst = dst_dict[self.folder_name_tag]
        self.file_name_dst = dst_dict[self.file_name_tag]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

        self.flag_cleaning_static = flag_cleaning_static

        self.chunks_n = parameters_dict['chunks']
        self.cpu_n = parameters_dict['cpu_datasets']

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize minimum distance
    @staticmethod
    def organize_minimum_distance(pnt_dict, pnt_idx='sm_n'):

        global pnt_collections
        if pnt_collections is None:
            pnt_collections = deepcopy(pnt_dict)
        else:

            pnt_keys_expected = list(pnt_collections.keys())

            for pnt_key, pnt_value in pnt_dict.items():
                if pnt_key in pnt_keys_expected:
                    pnt_value_tmp = pnt_collections[pnt_key]
                    if not isinstance(pnt_value_tmp, list):
                        pnt_value_tmp = [pnt_value_tmp]
                    pnt_value_tmp.append(pnt_value)
                    pnt_collections[pnt_key] = pnt_value_tmp
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to compute minimum distance
    def find_minimum_distance(self, geo_list):

        pnt_n = geo_list[0]
        pnt_x = geo_list[1]
        pnt_y = geo_list[2]
        pnt_id = geo_list[3]

        abs_y = np.abs(terrain_geo_ys - pnt_y)
        abs_x = np.abs(terrain_geo_xs - pnt_x)

        c = np.maximum(abs_x, abs_y)
        idx_xy = np.argmin(c)
        terrain_value = terrain_values.flat[idx_xy]
        terrain_idx = terrain_idxs.flat[idx_xy]
        terrain_geo_x = terrain_geo_xs.flat[idx_xy]
        terrain_geo_y = terrain_geo_ys.flat[idx_xy]

        terrain_idx_x, terrain_idx_y = np.where(c == np.min(c))

        pnt_dict = {'sm_geo_x': pnt_x, 'sm_geo_y': pnt_y, 'sm_id': pnt_id, 'sm_n': pnt_n,
                    'terrain_geo_x': terrain_geo_x, 'terrain_geo_y': terrain_geo_y,
                    'terrain_id': terrain_idx, 'terrain_value': terrain_value,
                    'terrain_idx_x': terrain_idx_x[0], 'terrain_idx_y': terrain_idx_y[0]}

        return pnt_dict

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to map geographical information
    def map_grid(self, obj_terrain, obj_sm):

        global terrain_values
        global terrain_geo_xs
        global terrain_geo_ys
        global terrain_idxs

        global pnt_collections

        logging.info(' ----> Create grid mapping terrain/variable ... ')
        time_start = time.time()

        terrain_values = obj_terrain['values']
        terrain_geo_xs = obj_terrain['longitude']
        terrain_geo_ys = obj_terrain['latitude']

        terrain_idxs = create_map_idx(terrain_values)

        sm_geo_x = obj_sm['longitude']
        sm_geo_y = obj_sm['latitude']
        sm_id = obj_sm['id']

        pnt_collections = None
        fx_pool = Pool(processes=self.cpu_n)
        for pnt_i, (pnt_x, pnt_y, pnt_id) in enumerate(zip(sm_geo_x, sm_geo_y, sm_id)):
            sm_geo = [pnt_i, pnt_x, pnt_y, pnt_id]

            # SERIAL
            # pnt_dict = self.find_minimum_distance(sm_geo)
            # self.organize_minimum_distance(pnt_dict)

            # ASYNC MULTIPROCESSING
            fx_pool.apply_async(self.find_minimum_distance, args=(sm_geo, ), callback=self.organize_minimum_distance)

        fx_pool.close()
        fx_pool.join()

        if pnt_collections is not None:
            pnt_dframe = pd.DataFrame(pnt_collections, index=[pnt_collections['sm_n']])
        else:
            pnt_dframe = None

        time_end = time.time()
        time_elapsed = time_end - time_start

        logging.info(' ----> Create grid mapping terrain/variable ... DONE [ELAPSED: ' +
                     str(np.floor(time_elapsed)) + ' seconds]')

        return pnt_dframe

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize geographical data
    def organize_static(self):

        # Info start
        logging.info(' ---> Organize static datasets ... ')

        flag_cleaning_static = self.flag_cleaning_static
        file_path_dst = self.file_path_dst

        if flag_cleaning_static:
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)

        if not os.path.exists(file_path_dst):

            # Read data terrain
            obj_terrain = read_file_raster(self.file_path_terrain)
            # Read data cn
            obj_cn = read_file_raster(self.file_path_cn)
            # Read data channel network
            obj_channel_network = read_file_raster(self.file_path_channel_network)

            # Read data soil moisture
            obj_sm = read_file_geo(self.file_path_sm)

            # Map terrain and soil moisture geographical information
            obj_map = self.map_grid(obj_terrain, obj_sm)

            # Data collection object
            static_data_collections = {
                self.flag_datasets_terrain: obj_terrain,
                self.flag_datasets_cn: obj_cn,
                self.flag_datasets_channel_network: obj_channel_network,
                self.flag_datasets_sm: obj_sm,
                self.flag_datasets_map: obj_map}

            folder_name, file_name = os.path.split(file_path_dst)
            make_folder(folder_name)

            write_obj(file_path_dst, static_data_collections)

            # Info end
            logging.info(' ---> Organize static datasets ... DONE')

        else:
            # Info end
            static_data_collections = read_obj(file_path_dst)
            logging.info(' ---> Organize static datasets ... LOADED. Datasets previously computed')

        return static_data_collections
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
