"""
Library Features:

Name:          lib_utils_envs
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize envs obj
def organize_envs_obj(filename_product=None, filename_xml=None,
                      creation_time_product=None, publication_time_product=None, creation_time_xml=None,
                      lon_east=None, lon_west=None, lat_south=None, lat_north=None,
                      begin_time_product=None, end_time_product=None):
    envs_obj = {
        'ENV_FILENAME_PRODUCT': filename_product,
        'ENV_CREATION_TIME_XML': creation_time_xml,
        'ENV_CREATION_TIME_PRODUCT': creation_time_product,
        'ENV_PUBLICATION_TIME_PRODUCT': publication_time_product,
        'ENV_FILENAME_XML': filename_xml,
        'ENV_LON_EAST': lon_east,
        'ENV_LON_WEST': lon_west,
        'ENV_LAT_SOUTH': lat_south,
        'ENV_LAT_NORTH': lat_north,
        'ENV_BEGIN_TIME_PRODUCT': begin_time_product,
        'ENV_END_TIME_PRODUCT': end_time_product
    }

    return envs_obj
# ----------------------------------------------------------------------------------------------------------------------
