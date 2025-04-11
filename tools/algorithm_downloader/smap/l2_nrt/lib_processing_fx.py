"""
Library Features:

Name:          lib_utils_processing
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20250410'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np

from lib_data_smap import read_smap_l2, process_smap_l2
from lib_data_io_tiff import read_file_tiff, write_file_tiff
from lib_info_args import logger_name

# logger stream
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to reproject data
def reproj_data(source_obj, ancillary_obj, destination_obj,
                geo_grid, geo_proj, geo_geotrans, grid_sm_2d,
                flag_remove_data_out_of_domain: bool = False):

    # get variable information
    var_name_list = source_obj['variable_name']
    var_path_list = source_obj['variable_path']

    # get file path(s)
    file_path_src = source_obj['file_path']
    file_path_anc_domain_list = ancillary_obj['file_path_domain']
    file_path_anc_global_list = ancillary_obj['file_path_global']
    file_path_dst_list = destination_obj['file_path']

    ff_source, fn_source = os.path.split(file_path_src)

    # info file start
    logger_stream.info(' ------> Analyze file "' + fn_source + '" ... ')

    for file_path_anc_global_step, file_path_anc_domain_step, file_path_dst_step, var_name_step, var_path_step in zip(
            file_path_anc_global_list, file_path_anc_domain_list, file_path_dst_list,
            var_name_list, var_path_list):

        # info variable start
        logger_stream.info(' -------> Variable "' + var_name_step + '" ... ')

        # get variable format
        one_down, variable = os.path.split(var_path_step)
        two_down, subgroup = os.path.split(one_down)

        # create folder(s) for ancillary and outcome file(s)
        ff_anc_global_step, fn_anc_global_step = os.path.split(file_path_anc_global_step)
        os.makedirs(ff_anc_global_step, exist_ok=True)

        ff_anc_domain_step, fn_anc_domain_step = os.path.split(file_path_anc_domain_step)
        os.makedirs(ff_anc_domain_step, exist_ok=True)

        ff_destination_step, fn_destination_step = os.path.split(file_path_dst_step)
        os.makedirs(ff_destination_step, exist_ok=True)

        # info start processing
        logger_stream.info(' --------> Translate and reproject over domain ... ')

        # Reproject from global tiff file in epsg:6933 modified to domain tiff file in epsg:4326
        # bbox information are used to subset the destination image
        # while original resolution is preserved
        smap_data, smap_meta, smap_geo_x, smap_geo_y, smap_col, smap_row =\
            read_smap_l2(file_path_src, subgroup, variable)

        smap_domain_2d = []
        for data, meta in zip(smap_data, smap_meta):
            data_2d = process_smap_l2(data, smap_geo_x, smap_geo_y, smap_col, smap_row, geo_grid, grid_sm_2d, meta)
            data_2d = np.nan_to_num(data_2d, nan=-9999.)
            smap_domain_2d.append(data_2d)
        # get file dimensions
        smap_domain_high, smap_domain_wide = smap_domain_2d[0].shape

        # check data over domain
        if np.all(np.where(smap_domain_2d[0] == -9999., True, False)):

            # info no data available
            logger_stream.warning(' ===> No valid soil_moisture data over domain.')

            if flag_remove_data_out_of_domain:
                logger_stream.info(' -------> Cleaning source file ... ')
                os.remove(file_path_src)
                logger_stream.info(' -------> Cleaning source file ... DONE')

            # info end translate and reproject (skip)
            logger_stream.info(' --------> Translate and reproject over domain ... SKIPPED. '
                               'Datasets are not available over the domain ')
        else:

            # info end translate and reproject (done)
            logger_stream.info(' --------> Translate and reproject over domain ... DONE')

            # info start save
            logger_stream.info(' --------> Save "' + file_path_dst_step + '" ... ')
            write_file_tiff(file_path_dst_step, smap_domain_2d, smap_meta, smap_domain_wide,
                            smap_domain_high, geo_geotrans, geo_proj, flip=True)
            # info end save
            logger_stream.info(' --------> Save "' + file_path_dst_step + '" ... DONE')

        # info variable end
        logger_stream.info(' -------> Variable "' + var_name_step + '" ... DONE')

    # info file end
    logger_stream.info(' ------> Analyze file "' + fn_source + '" ... DONE')
# ----------------------------------------------------------------------------------------------------------------------
