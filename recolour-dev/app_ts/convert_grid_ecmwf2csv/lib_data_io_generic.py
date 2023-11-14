"""
Library Features:

Name:          lib_data_io_generic
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
# method to join data point
def join_data_point(point_time, point_data, point_collection, name_index='time'):

    # join point collections
    if point_collection is None:
        # iterate over variable(s)
        point_collection = {}
        for var_name, var_value in point_data.items():
            # create dataframe collection
            var_df = pd.DataFrame(data=var_value, index=[point_time])
            var_df.index.name = name_index
            # store dataframe in a common obj
            point_collection[var_name] = var_df
    else:

        # iterate over variable(s)
        for var_name, var_value in point_data.items():
            # get tmp collection
            collection_df = point_collection[var_name]
            # create dataframe collection
            tmp_df = pd.DataFrame(data=var_value, index=[point_time])
            tmp_df.index.name = name_index

            # append new line to dataframe collection
            collection_df = pd.concat([collection_df, tmp_df])
            # store in the collection
            point_collection[var_name] = collection_df

    return point_collection
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extract grid values to point values
def extract_data_grid2point(
        grid_obj_data, grid_obj_geo, point_obj_geo,
        method_spatial_operation='average',
        method_spatial_mask=True, value_spatial_mask=0, value_decimal_round=4):

    # iterate over variable(s)
    point_obj_collections = {}
    for grid_name, grid_da in grid_obj_data.items():

        # info variable start
        log_stream.info(' -------> Variable "' + grid_name + '" ...')

        # get dataset and geo
        grid_data = grid_da.values.ravel()
        grid_geo = grid_obj_geo.values.ravel()

        # iterate over point(s)
        point_obj_dict = {}
        for point_id, point_row in point_obj_geo.iterrows():
            point_name, point_tag, point_idx_1d_list = point_row['name'], point_row['tag'], point_row['point_idx_1d']

            # info point start
            log_stream.info(' --------> Point "' + point_tag + '" ...')

            if not isinstance(point_idx_1d_list, list):
                point_idx_1d_list = [point_idx_1d_list]

            # iterate over indexes
            point_data_list_unmasked, point_data_list_masked = [], []
            for point_idx_1d_step in point_idx_1d_list:

                point_data, point_geo = grid_data[point_idx_1d_step], grid_geo[point_idx_1d_step]

                if method_spatial_mask:
                    if point_geo != value_spatial_mask:
                        point_data_list_unmasked.append(point_data)
                    else:
                        point_data_list_masked.append(point_data)

            # select methods to get point(s) values
            if method_spatial_operation == 'average':
                point_data_array = np.array(point_data_list_unmasked, dtype=float)
                point_data_num = float(np.mean(point_data_array))
            else:
                log_stream.error(' ===> Spatial data operation "' + method_spatial_operation + '" is not supported')
                raise NotImplementedError('Case not implemented yet. Only "average" method is available.')

            # round point data
            point_data_num = round(point_data_num, value_decimal_round)

            # store point
            point_obj_dict[point_tag] = point_data_num

            # info point end
            log_stream.info(' --------> Point "' + point_tag + '" ... DONE')

        # store all variable(s)
        point_obj_collections[grid_name] = point_obj_dict

        # info variable end
        log_stream.info(' -------> Variable "' + grid_name + '" ... DONE')

    return point_obj_collections
# ----------------------------------------------------------------------------------------------------------------------
