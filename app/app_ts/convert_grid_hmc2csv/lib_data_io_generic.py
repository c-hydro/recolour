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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to combine data over the expected time range
def combine_data_point_by_time(point_data_collections_src, point_geo,
                               time_start_expected, time_end_expected,
                               time_frequency_expected='D', time_reverse=True):
    # create time range expected
    time_range_expected = pd.date_range(start=time_start_expected, end=time_end_expected, freq=time_frequency_expected)

    # iterate over variable(s)
    point_data_collections_dst = {}
    if point_data_collections_src is not None:
        for point_var, point_data_raw in point_data_collections_src.items():

            # create null data
            null_data = np.zeros(shape=(time_range_expected.shape[0], point_geo.shape[0]))
            null_data[:, :] = np.nan
            null_dict = {}
            for point_id, point_label in enumerate(point_data_raw.columns):
                null_dict[point_label] = null_data[:, point_id]
            point_data_expected = pd.DataFrame(data=null_dict, index=time_range_expected)

            # update expected data with raw data
            point_data_expected.update(point_data_raw)
            # time reverse flag
            if time_reverse:
                point_data_expected = point_data_expected.sort_index(ascending=False)

            # store data in a common obj
            point_data_collections_dst[point_var] = point_data_expected

    else:
        log_stream.warning(' ===> No data available to combine over the expected time range')
        point_data_collections_dst = None

    return point_data_collections_dst

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to join data point
def join_data_point(point_time, point_data, point_collection, name_index='time'):

    # join point collections
    if point_collection is None:
        # iterate over variable(s)
        if point_data is not None:
            point_collection = {}
            for var_name, var_value in point_data.items():
                # create dataframe collection
                var_df = pd.DataFrame(data=var_value, index=[point_time])
                var_df.index.name = name_index
                # store dataframe in a common obj
                point_collection[var_name] = var_df
        else:
            log_stream.warning(' ===> No data available to join')
    else:

        # iterate over variable(s)
        if point_data is not None:
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
        else:
            log_stream.warning(' ===> No data available to join')

    return point_collection
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to extract grid values to point values
def extract_data_grid2point(
        grid_obj_data, grid_obj_geo, point_obj_geo,
        method_spatial_operation='average',
        method_spatial_mask=True, value_no_data=-9999.0,
        value_spatial_mask=0, value_decimal_round=4):

    # iterate over variable(s)
    point_obj_collections = {}
    for grid_name, grid_da in grid_obj_data.items():

        # info variable start
        log_stream.info(' -------> Variable "' + grid_name + '" ...')

        # get dataset and geo
        grid_data = grid_da.values.ravel()
        grid_geo = grid_obj_geo.values.ravel()

        ''' debug
        grid_values = grid_da.values[:, :, 0]
        grid_values[grid_values < 0] = np.nan
        plt.figure()
        plt.imshow(grid_values)
        plt.colorbar()
        plt.show()
        '''

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

                if point_data == value_no_data:
                    point_data = np.nan

                if method_spatial_mask:
                    if point_geo != value_spatial_mask:
                        point_data_list_unmasked.append(point_data)
                    else:
                        point_data_list_masked.append(point_data)

            if np.isnan(point_data_list_unmasked).any():
                log_stream.warning(' ===> Some values are defined by NaN')
            elif np.isnan(point_data_list_unmasked).all():
                log_stream.warning(' ===> All values are defined by NaN')

            # select methods to get point(s) values
            if method_spatial_operation == 'average':
                point_data_array = np.array(point_data_list_unmasked, dtype=float)
                point_data_num = float(np.nanmean(point_data_array))
            elif method_spatial_operation == 'maximum':
                point_data_array = np.array(point_data_list_unmasked, dtype=float)
                point_data_num = float(np.nanmax(point_data_array))
            else:
                log_stream.error(' ===> Spatial data operation "' + method_spatial_operation + '" is not supported')
                raise NotImplementedError('Case not implemented yet. Methods "average" or "maximum" are available.')

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
