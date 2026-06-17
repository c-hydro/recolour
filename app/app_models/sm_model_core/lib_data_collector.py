"""
Library Features:

Name:          lib_data_collector
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd

from copy import deepcopy

from lib_utils_generic import fill_tags2string
from lib_utils_obj import map_vars_dframe, sanitize_string, fill_tags_time
from lib_utils_time import replace_time_part
from lib_utils_io import fill_string_with_time, fill_string_with_info
from lib_info_args import logger_name, time_format_algorithm

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# method to define file string
def __define_file_string(file_string_tmpl, extended_info=None):
    if extended_info is not None:
        alg_info = {**self.alg_info, **extended_info}
    else:
        alg_info = self.alg_info

    file_string_def = fill_string_with_time(file_string_tmpl, self.time_reference, self.alg_template_time)
    file_string_def = fill_string_with_info(file_string_def, alg_info, self.alg_template_datasets)
    return file_string_def


# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to collect data one file one poiont
def collect_data_one_file_one_point():

    # iterate over geo point(s)
    obj_collections = {}
    for fields_data in data_registry.to_dict(orient="records"):

        # get point information
        point_name, point_tag = fields_data['name'], fields_data['tag']

        # info point start
        log_stream.info(' -----> Point -- (1) Name: "' + point_tag + '" :: (2) Tag: "' + point_tag + '" ... ')

        # method to fill the filename(s)
        file_path_src_rain_point = self.__define_file_string(
            file_path_src_rain_tmpl, extended_info={'point_name': point_tag})
        file_path_src_airt_point = self.__define_file_string(
            file_path_src_airt_tmpl, extended_info={'point_name': point_tag})
        file_path_src_sm_point = self.__define_file_string(
            file_path_src_sm_tmpl, extended_info={'point_name': point_tag})

        file_path_dst_point = self.__define_file_string(
            file_path_dst_tmpl, extended_info={'point_name': point_tag})

        # reset ancillary file if required
        if reset_data_dynamic:
            if os.path.exists(file_path_dst_point):
                os.remove(file_path_dst_point)

        # check ancillary file availability
        if not os.path.exists(file_path_dst_point):

            # get rain dataframe
            dframe_rain = self.get_obj_datasets(
                file_path_src_rain_point,
                file_format=self.format_rain, file_delimiter=self.delimiter_rain, file_mandatory=True,
                time_fields=self.time_rain, file_fields=self.fields_rain, registry_fields=data_registry)

            # get air temperature dataframe
            dframe_airt = self.get_obj_datasets(
                file_path_src_airt_point,
                file_format=self.format_airt, file_delimiter=self.delimiter_airt, file_mandatory=True,
                time_fields=self.time_airt, file_fields=self.fields_airt, registry_fields=data_registry)

            # get soil moisture dataframe
            dframe_sm = self.get_obj_datasets(
                file_path_src_sm_point,
                file_format=self.format_sm, file_delimiter=self.delimiter_sm, file_mandatory=False,
                time_fields=self.time_sm, file_fields=self.fields_sm, registry_fields=data_registry)

            # create combined dataframe
            dframe_combined = combine_data_point_by_time(
                dframe_k1=dframe_rain, dframe_k2=dframe_airt, dframe_k3=dframe_sm)

            # check combined dataframe
            if dframe_combined is not None:

                # dump combined dataframe
                self.dump_obj_datasets(
                    file_path_dst_point, dframe_combined, file_format=self.format_dst,
                    file_fields=self.fields_dst, time_fields=self.time_dst, registry_fields=data_registry)

                # store file path
                obj_collections[point_tag] = file_path_dst_point

                # info point end
                log_stream.info(' -----> Point -- (1) Name: "' + point_tag + '" :: (2) Tag: "' + point_tag +
                                '" ... DONE')
            else:
                log_stream.info(' -----> Point -- (1) Name: "' + point_tag + '" :: (2) Tag: "' + point_tag +
                                '" ... SKIPPED. Datasets not available')

        else:

            # store file path
            obj_collections[point_tag] = file_path_dst_point

            # info point end
            log_stream.info(' -----> Point -- (1) Name: "' + point_tag + '" :: (2) Tag: "' + point_tag +
                            '" ... SKIPPED. Datasets previously saved')

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to wrap parameters in ascii format
def read_parameters_ascii(file_name, file_column='param_{:}'):

    with open(file_name, 'r') as file:
        file_content = file.readlines()
    file_content = [x.strip() for x in file_content]
    file_content = [float(x) for x in file_content if x]

    file_object = {}
    for i, x in enumerate(file_content):
        file_object[file_column.format(i + 1)] = x

    return file_object
# ----------------------------------------------------------------------------------------------------------------------
