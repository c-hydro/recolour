"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20231110'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
import os
import pandas as pd

from lib_utils_io import fill_tags2string
from lib_info_args import logger_name

# logger stream
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to collect data info
def collect_data_info(info_time=None,
        info_product=None, info_bbox=None, info_url=None):

    pr_short_name = info_product['short_name']
    pr_version = info_product['version']
    pr_template_root = info_product['template_root']
    pr_template_vars_data = info_product['template_vars_data']
    pr_template_group_data = info_product['template_group_data']

    pr_bbox_lon_right = info_bbox['lon_right']
    pr_bbox_lon_left = info_bbox['lon_left']
    pr_bbox_lat_top = info_bbox['lat_top']
    pr_bbox_lat_bottom = info_bbox['lat_bottom']

    # wget variable(s)
    pr_remote_url = info_url['remote_url']
    pr_remote_folder = info_url['remote_folder']

    if info_time is not None:
        if isinstance(info_time, str):
            info_time = pd.Timestamp(info_time)
        info_time = info_time.strftime("%Y.%m.%d")

    pr_remote_folder = pr_remote_folder.format(remote_sub_path_time=info_time)

    # bounding_box [min_lon min_lat max_lon max_lat]
    pr_bbox_tmp = [str(pr_bbox_lon_right), str(pr_bbox_lat_bottom), str(pr_bbox_lon_left), str(pr_bbox_lat_top)]
    pr_bbox_ref = ','.join(pr_bbox_tmp)

    return (pr_short_name, pr_version,
            pr_template_root, pr_template_vars_data, pr_template_group_data,
            pr_bbox_ref, pr_remote_url, pr_remote_folder)

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create data destination object
def set_data_destination(file_id, source_obj, variable_obj=None, group_obj=None, file_obj=None,
                         ancillary_obj=None, template_obj=None, flag_remove_destination=False):

    time_stamp = source_obj['time_stamp']

    folder_name_raw = file_obj['folder'][file_id]
    file_name_raw = file_obj['filename'][file_id]

    group_list = group_obj[file_id]
    variable_list = variable_obj[file_id][:]
    domain = ancillary_obj['domain']

    time_def = time_stamp.to_pydatetime()

    for var in variable_list:
        i = variable_list.index(var)

        if "AM" in var: timing = "0600"
        elif "PM" in var: timing = "1800"
        else: timing = ""

        if "/" in var: var = var.split("/")[-1]
        if "_pm" in var: var = var.replace("_pm", "")
        if timing: var = "_".join([timing, var])
        var = var.replace("__", "_") if "__" in var else var
        variable_list[i] = var

    var_list = []
    for variable, group in zip(variable_list, group_list):

        template_values = {"domain": domain,
                           "var_name": variable,
                           "group_name": group,
                           "destination_sub_path_time": time_def,
                           "destination_datetime": time_def}

        folder_name_def = fill_tags2string(folder_name_raw, template_obj, template_values)
        file_name_def = fill_tags2string(file_name_raw, template_obj, template_values)
        path_name_def = os.path.join(folder_name_def, file_name_def)

        if flag_remove_destination:
            if os.path.exists(path_name_def):
                os.remove(path_name_def)

        var_list.append(path_name_def)

    destination_info = {
        "domain": domain,
        "timing": timing,
        "time_datetime": time_def,
        "time_stamp": time_stamp,
        "variable_name": variable_list,
        "group_name": group_list,
        "file_path": var_list,
    }

    return destination_info
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create data ancillary object
def set_data_ancillary(file_id, source_obj, variable_obj=None, file_obj=None,
                       ancillary_obj=None, template_obj=None,
                       flag_remove_ancillary_global=False, flag_remove_ancillary_domain=False):

    time_stamp = source_obj['time_stamp']

    folder_global_raw = file_obj['global']['folder'][file_id]
    filename_global_raw = file_obj['global']['filename'][file_id]
    folder_domain_raw = file_obj['domain']['folder'][file_id]
    filename_domain_raw = file_obj['domain']['filename'][file_id]

    variable_list = variable_obj[file_id][:]
    domain = ancillary_obj['domain']

    time_def = time_stamp.to_pydatetime()

    # l3 product compatibility
    for var in variable_list:
        i = variable_list.index(var)

        if "AM" in var: timing = "0600"
        elif "PM" in var: timing = "1800"
        else: timing=""

        if "/" in var: var = var.split("/")[-1]
        if "_pm" in var: var = var.replace("_pm", "")
        if timing is not "": var = "_".join([timing, var])
        var = var.replace("__","_") if "__" in var else var
        variable_list[i] = var

    var_global_list, var_domain_list = [], []
    for variable in variable_list:

        template_values = {"domain": domain,
                           "var_name": variable,
                           "ancillary_sub_path_time": time_def,
                           "ancillary_datetime": time_def}

        folder_global_def = fill_tags2string(folder_global_raw, template_obj, template_values)
        filename_global_def = fill_tags2string(filename_global_raw, template_obj, template_values)
        path_global_def = os.path.join(folder_global_def, filename_global_def)

        folder_domain_def = fill_tags2string(folder_domain_raw, template_obj, template_values)
        filename_domain_def = fill_tags2string(filename_domain_raw, template_obj, template_values)
        path_domain_def = os.path.join(folder_domain_def, filename_domain_def)

        if flag_remove_ancillary_global:
            if os.path.exists(path_global_def):
                os.remove(path_global_def)
        if flag_remove_ancillary_domain:
            if os.path.exists(path_domain_def):
                os.remove(path_domain_def)

        var_global_list.append(path_global_def)
        var_domain_list.append(path_domain_def)

    ancillary_info = {
        "domain": domain,
        "timing": timing,
        "time_datetime": time_def,
        "time_stamp": time_stamp,
        "variable_name": variable_list,
        "file_path_domain": var_domain_list,
        "file_path_global": var_global_list,
    }

    return ancillary_info
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create data source object
def set_data_source(file_id, file_url,
                    file_obj=None, variable_obj=None, root_obj=None, ancillary_obj=None, template_obj=None,
                    file_suffix='.h5', flag_remove_source=False):

    folder_name_raw, file_name_raw= file_obj['folder'][file_id], file_obj['filename'][file_id]
    domain = ancillary_obj['domain']

    variable_name_list = variable_obj[file_id][:]
    file_variable_raw = root_obj[file_id]

    source_info = {}
    if file_url.endswith(file_suffix):

        # manage file information (from name)
        match_time = re.search(r'\d{4}\d{2}\d{2}T\d{2}\d{2}\d{2}', file_url)
        orbit     = str(re.search(r'_E_\d{5}', file_url).group())
        dir       = str(re.search(r'_A_|_D_', file_url).group())
        releaseid = str(re.search(r'_R\d{5}_', file_url).group())
        version   = str(re.search(r'\d{3}.h5', file_url).group())

        # manage time information
        time_str = match_time.group()
        time_stamp = pd.Timestamp(time_str)

        time_def = time_stamp.to_pydatetime()
        template_values = {
            "domain": domain,
            "source_sub_path_time": time_def,
            "source_datetime": time_def,
            "orbit": orbit,
            "direction": dir,
            "release_id": releaseid,
            "version": version,
        }

        # create file information
        folder_name_def = fill_tags2string(folder_name_raw, template_obj, template_values)
        file_name_def = fill_tags2string(file_name_raw, template_obj, template_values)
        path_name_def = os.path.join(folder_name_def, file_name_def)

        if flag_remove_source:
            if os.path.exists(path_name_def):
                try:
                    os.remove(path_name_def)
                except BaseException as e:
                    logger_stream.warning(' ===> Try to delete directly the folder "' + str(e) + '"')
                    os.rmdir(path_name_def)

        # create variable info
        file_variable_list = []
        for variable in variable_name_list:

            template_values = {
                "domain": domain,
                "file_name": path_name_def,
                "var_name": variable,
                "source_sub_path_time": time_def,
                "source_datetime": time_def}

            file_variable_def = fill_tags2string(file_variable_raw, template_obj, template_values)
            file_variable_list.append(file_variable_def)

        source_info = {
            "domain": domain,
            "time_datetime": time_def,
            "time_stamp": time_stamp,
            "time_string": time_str,
            "orbit": orbit,
            "direction": dir,
            "release_id": releaseid,
            "version": version,
            "file_name": file_name_def,
            "file_path": path_name_def,
            "file_url": file_url,
            "folder_name": folder_name_def,
            "variable_name": variable_name_list,
            "variable_path": file_variable_list,
        }

    return source_info

# ----------------------------------------------------------------------------------------------------------------------
