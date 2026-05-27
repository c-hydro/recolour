"""
Library Features:

Name:           lib_utils_basic
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import re
import numpy as np

from glob import glob

from lib_utils_time import iter_time_steps, resolve_time_tags
from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to make folder
def make_folder(path_name):
    os.makedirs(path_name, exist_ok=True)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to format cell tags
def format_cell_tags(cell_value, cell_digits=4):
    cell_int = int(cell_value)
    return {
        "cell": str(cell_int),
        "cell_n": f"{cell_int:0{int(cell_digits)}d}"
    }
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to resolve generic tags
def resolve_generic_tags(template_string, tags):
    template_resolved = resolve_time_tags(template_string, tags)

    for tag_key, tag_value in tags.items():
        if isinstance(tag_value, (str, int, float)):
            template_resolved = template_resolved.replace("{" + str(tag_key) + "}", str(tag_value))

    return template_resolved
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to discover porosity files
def discover_porosity_files(porosity_settings, cell_list, cell_digits=4, strict=True):

    porosity_folder = porosity_settings["folder"]
    porosity_filename_template = porosity_settings.get("filename", "soilgrids_cell_{cell:04d}.nc")

    selected_files = []
    missing_files = []
    seen_files = set()

    unique_cells = sorted(np.unique(cell_list).astype(np.int32))

    for cell_id in unique_cells:

        cell_tags = format_cell_tags(cell_id, cell_digits=cell_digits)

        try:
            porosity_filename = porosity_filename_template.format(
                cell=int(cell_id),
                cell_n=cell_tags["cell_n"]
            )
        except KeyError:
            porosity_filename = resolve_generic_tags(
                porosity_filename_template,
                cell_tags
            )

        file_path = os.path.join(porosity_folder, porosity_filename)

        if os.path.exists(file_path):
            if file_path not in seen_files:
                seen_files.add(file_path)
                selected_files.append(file_path)
        else:
            missing_files.append(file_path)

    if strict and missing_files:
        msg = "\n".join(missing_files)
        raise RuntimeError(
            "Porosity conversion is enabled, but some required porosity files are missing:\n"
            f"{msg}"
        )

    return sorted(selected_files)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to discover source files
def discover_source_files(source_settings, time_settings,
                          time_start, time_end, reference_time,
                          cell_list=None, cell_digits=4):

    source_folder_template = source_settings["folder"]
    source_filename_template = source_settings.get("filename", "h122_cell_{cell_n}.nc")

    selected_files = []
    seen_files = set()

    for time_step in iter_time_steps(time_settings, time_start, time_end):
        common_tags = {
            "time_step": time_step,
            "time_start": time_start,
            "time_end": time_end,
            "time_run": reference_time
        }

        src_folder = resolve_time_tags(source_folder_template, common_tags)

        logger.info(f" ----> Search source folder: {src_folder}")

        if not os.path.exists(src_folder):
            continue

        if cell_list is not None:
            for cell_id in cell_list:
                cell_tags = format_cell_tags(cell_id, cell_digits=cell_digits)

                source_filename = resolve_generic_tags(
                    source_filename_template,
                    {
                        **common_tags,
                        **cell_tags
                    }
                )

                file_path = os.path.join(src_folder, source_filename)

                if os.path.exists(file_path) and file_path not in seen_files:
                    seen_files.add(file_path)
                    selected_files.append(file_path)
        else:
            search_pattern = source_filename_template
            search_pattern = search_pattern.replace("{cell}", "*")
            search_pattern = search_pattern.replace("{cell_n}", "*")
            search_pattern = resolve_time_tags(search_pattern, common_tags)

            file_list = sorted(
                f for f in glob(os.path.join(src_folder, search_pattern))
                if re.search(r"\b\d{4}\b", os.path.splitext(os.path.basename(f))[0])
            )

            for file_path in file_list:
                if file_path not in seen_files:
                    seen_files.add(file_path)
                    selected_files.append(file_path)

    return sorted(selected_files)
# ----------------------------------------------------------------------------------------------------------------------
