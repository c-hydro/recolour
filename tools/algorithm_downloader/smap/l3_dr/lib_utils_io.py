"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20231110'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import re
import os
import pandas as pd

from lib_utils_generic import fill_tags2string
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to clean tmp data such as ancillary or global (if domain is set)
def clean_data_tmp(filename_obj_source, filename_obj_ancillary_global, filename_obj_ancillary_domain,
                   flag_cleaning_source=False, flag_cleaning_tmp=False):
    if flag_cleaning_source:
        for data_key, data_value in filename_obj_source.items():
            if not isinstance(data_value, list):
                data_value = list(data_value)
            for data_step in data_value:
                if os.path.exists(data_step):
                    os.remove(data_step)
    if flag_cleaning_tmp:
        for data_key, data_value in filename_obj_ancillary_global.items():
            for data_step in data_value:
                if os.path.exists(data_step):
                    os.remove(data_step)
        for data_key, data_value in filename_obj_ancillary_domain.items():
            for data_step in data_value:
                if os.path.exists(data_step):
                    os.remove(data_step)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create data outcome list
def set_data_outcome(time_stamp_list, file_id, variable_obj=None, group_obj=None, file_obj=None,
                       ancillary_obj=None, template_obj=None, flag_cleaning_outcome=False):

    folder_raw = file_obj['folder'][file_id]
    filename_raw = file_obj['filename'][file_id]

    group_list = group_obj[file_id]
    variable_list = variable_obj[file_id][:]
    domain = ancillary_obj['domain']

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

    filename_list = {}
    for time_stamp in time_stamp_list:

        var_list = []
        for variable, group in zip(variable_list, group_list):

            time_step = time_stamp.to_pydatetime()
            template_values = {"domain": domain,
                               "var_name": variable,
                               "group_name": group,
                               "outcome_sub_path_time": time_step,
                               "outcome_datetime": time_step}

            folder_step = fill_tags2string(folder_raw, template_obj, template_values)
            filename_step = fill_tags2string(filename_raw, template_obj, template_values)
            path_step = os.path.join(folder_step, filename_step)

            if flag_cleaning_outcome:
                if os.path.exists(path_step):
                    os.remove(path_step)

            var_list.append(path_step)

        filename_list[time_stamp] = var_list

    return filename_list
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create data ancillary list
def set_data_ancillary(time_stamp_list, file_id, variable_obj=None, file_obj=None,
                       ancillary_obj=None, template_obj=None,
                       flag_cleaning_ancillary_global=False, flag_cleaning_ancillary_domain=False):

    folder_global_raw = file_obj['global']['folder'][file_id]
    filename_global_raw = file_obj['global']['filename'][file_id]
    folder_domain_raw = file_obj['domain']['folder'][file_id]
    filename_domain_raw = file_obj['domain']['filename'][file_id]

    variable_list = variable_obj[file_id][:]
    domain = ancillary_obj['domain']

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

    filename_global_list = {}
    filename_domain_list = {}
    for time_stamp in time_stamp_list:

        var_global_list = []
        var_domain_list = []
        for variable in variable_list:

            time_step = time_stamp.to_pydatetime()
            template_values = {"domain": domain,
                               "var_name": variable,
                               "ancillary_sub_path_time": time_step,
                               "ancillary_datetime": time_step}

            folder_global_step = fill_tags2string(folder_global_raw, template_obj, template_values)
            filename_global_step = fill_tags2string(filename_global_raw, template_obj, template_values)
            path_global_step = os.path.join(folder_global_step, filename_global_step)

            folder_domain_step = fill_tags2string(folder_domain_raw, template_obj, template_values)
            filename_domain_step = fill_tags2string(filename_domain_raw, template_obj, template_values)
            path_domain_step = os.path.join(folder_domain_step, filename_domain_step)

            if flag_cleaning_ancillary_global:
                if os.path.exists(path_global_step):
                    os.remove(path_global_step)
            if flag_cleaning_ancillary_domain:
                if os.path.exists(path_domain_step):
                    os.remove(path_domain_step)

            var_global_list.append(path_global_step)
            var_domain_list.append(path_domain_step)

        filename_global_list[time_stamp] = var_global_list
        filename_domain_list[time_stamp] = var_domain_list

    return filename_global_list, filename_domain_list
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create data source list
def set_data_source(file_id, filename_url,
                    file_obj=None, variable_obj=None, root_obj=None, ancillary_obj=None, template_obj=None,
                    filename_suffix='.h5', flag_cleaning_source=False):

    if not isinstance(filename_url, list):
        filename_url = [filename_url]

    folder_raw = file_obj['folder'][file_id]
    filename_raw = file_obj['filename'][file_id]
    domain = ancillary_obj['domain']

    variable_list = variable_obj[file_id][:]
    fileroot_raw = root_obj[file_id]

    time_stamp_list = []
    filename_list_url = []
    filename_obj_source = {}
    fileroot_obj_source = {}
    for filename_url_step in filename_url:

        if filename_url_step.endswith(filename_suffix):

            # match_time = re.search(r'\d{4}\d{2}\d{2}', filename_url_step) # TODO: check time format of source file
            # original:
            # match_time = re.search(r'\d{4}\d{2}\d{2}\w\d{2}\d{2}\d{2}', filename_url_step)
            match_time = re.search(r'\d{4}\d{2}\d{2}T\d{2}\d{2}\d{2}', filename_url_step)

            # -------------------------------
            orbit     = str(re.search(r'_E_\d{5}', filename_url_step).group())
            dir       = str(re.search(r'_A_|_D_', filename_url_step).group())
            releaseid = '18290' # str(re.search(r'\d{5}', filename_url_step).group())
            version   = str(re.search(r'\d{3}.h5', filename_url_step).group())

            # -------------------------------

            time_str = match_time.group()
            time_stamp = pd.Timestamp(time_str)

            time_step = time_stamp.to_pydatetime()
            template_values = {"domain": domain,
                               "source_sub_path_time": time_step,
                               "source_datetime": time_step,
                               "orbit": orbit,
                               "direction": dir,
                               "release_id": releaseid,
                               "version": version,
                               }

            folder_step = fill_tags2string(folder_raw, template_obj, template_values)
            filename_step = fill_tags2string(filename_raw, template_obj, template_values)

            path_step = os.path.join(folder_step, filename_step)

            if flag_cleaning_source:
                if os.path.exists(path_step):
                    os.remove(path_step)

            time_stamp_list.append(time_stamp)
            filename_list_url.append(filename_url_step)

            if time_step not in list(filename_obj_source.keys()):
                filename_obj_source[time_step] = [path_step]
            else:
                logging.error(' ===> Time is always set in source obj')
                raise NotImplementedError('Merge filename(s) is not implemented yet')

            var_fileroot_list = []
            for variable in variable_list:
                template_values = {"domain": domain,
                                   "file_name": path_step,
                                   "var_name": variable,
                                   "source_sub_path_time": time_step,
                                   "source_datetime": time_step}

                fileroot_step = fill_tags2string(fileroot_raw, template_obj, template_values)
                var_fileroot_list.append(fileroot_step)
            fileroot_obj_source[time_stamp] = var_fileroot_list

    return time_stamp_list, filename_list_url, filename_obj_source, fileroot_obj_source

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to collect product info
def collect_product_info(info_product=None, info_bbox=None, info_url=None):

    pr_short_name = info_product['short_name']
    pr_version = info_product['version']
    pr_template_root = info_product['template_root']
    pr_template_vars_data = info_product['template_vars_data']
    pr_template_group_data = info_product['template_group_data']

    pr_bbox_lon_right = info_bbox['lon_right']
    pr_bbox_lon_left = info_bbox['lon_left']
    pr_bbox_lat_top = info_bbox['lat_top']
    pr_bbox_lat_bottom = info_bbox['lat_bottom']

    pr_cmr_url = info_url['cmr_url']
    pr_urs_url = info_url['urs_url']
    pr_cmr_page_size = info_url['cmr_page_size']
    pr_url_list = info_url['url_list']
    pr_cmr_file_url = info_url['cmr_file_url']
    pr_polygon = info_url['polygon']
    pr_filename_filter = info_url['filename_filter']

    # bounding_box [min_lon min_lat max_lon max_lat]
    pr_bbox_tmp = [str(pr_bbox_lon_right), str(pr_bbox_lat_bottom), str(pr_bbox_lon_left), str(pr_bbox_lat_top)]
    pr_bbox_ref = ','.join(pr_bbox_tmp)

    return pr_short_name, pr_version, pr_template_root, pr_template_vars_data, pr_template_group_data, \
        pr_bbox_ref, \
        pr_cmr_url, pr_urs_url, pr_cmr_page_size, pr_url_list, pr_cmr_file_url, pr_polygon, pr_filename_filter

# -------------------------------------------------------------------------------------
