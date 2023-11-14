"""
Library Features:

Name:          lib_data_io_json
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to join data point
def join_data_point(point_time, point_data, point_collection, name_index='time'):

    # join point collections
    if point_collection is None:
        # create dataframe collection
        point_collection = pd.DataFrame(data=point_data, index=[point_time])
        point_collection.index.name = name_index
    else:
        # create tmp collection
        tmp_collection = pd.DataFrame(data=point_data, index=[point_time])
        tmp_collection.index.name = name_index
        # append new line to dataframe collection
        point_collection = pd.concat([point_collection, tmp_collection])

    return point_collection
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to extract grid values to point values
def extract_data_grid2point(
        grid_obj_da, grid_obj_geo, point_obj_geo, method_spatial_operation='average',
        method_spatial_mask=True, value_spatial_mask=0):

    grid_data = grid_obj_da.values.ravel()
    grid_geo = grid_obj_geo.values.ravel()

    point_obj_dict = {}
    for point_id, point_row in point_obj_geo.iterrows():
        point_name, point_tag, point_idx_1d_list = point_row['name'], point_row['tag'], point_row['point_idx_1d']

        log_stream.info(' -------> Point "' + point_tag + '" ...')

        if not isinstance(point_idx_1d_list, list):
            point_idx_1d_list = [point_idx_1d_list]

        point_data_list_unmasked, point_data_list_masked = [], []
        for point_idx_1d_step in point_idx_1d_list:

            point_data, point_geo = grid_data[point_idx_1d_step], grid_geo[point_idx_1d_step]

            if method_spatial_mask:
                if point_geo != value_spatial_mask:
                    point_data_list_unmasked.append(point_data)
                else:
                    point_data_list_masked.append(point_data)

        if method_spatial_operation == 'average':
            point_data_array = np.array(point_data_list_unmasked, dtype=float)
            point_data_num = float(np.mean(point_data_array))
        else:
            log_stream.error(' ===> Spatial data operation "' + method_spatial_operation + '" is not supported')
            raise NotImplementedError('Case not implemented yet. Only "average" method is available.')

        point_obj_dict[point_tag] = point_data_num

        log_stream.info(' -------> Point "' + point_tag + '" ... DONE')

    return point_obj_dict
# ----------------------------------------------------------------------------------------------------------------------
