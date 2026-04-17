#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SSM H122 TRANSFER - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20260415'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_transfer_ssm_h122.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20260415 (1.0.0) --> First development
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import time
import os
import re
import sys
import json
import shutil
import argparse

from glob import glob
from datetime import datetime, timedelta

import numpy as np
import xarray as xr
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for transfering ssm h122 files over a configured geographic domain'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2026-04-15'

# algorithm globals
logger = logging.getLogger("app_transfer")
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# basic utils
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f'File "{file_name}" not found. Exit')

def make_folder(path_name):
    os.makedirs(path_name, exist_ok=True)

def resolve_time_template(path_template, time_obj):
    return time_obj.strftime(path_template)

def normalize_time_to_hour(time_obj):
    return time_obj.replace(minute=0, second=0, microsecond=0)

def parse_reference_time(time_string=None):
    if time_string is None:
        return normalize_time_to_hour(datetime.now())

    try:
        time_obj = datetime.strptime(time_string, "%Y-%m-%d %H:%M")
    except ValueError as exc:
        raise ValueError(
            'Argument -time must have format "YYYY-MM-DD HH:MM"'
        ) from exc

    return normalize_time_to_hour(time_obj)

def is_template_hourly(path_template):
    return any(token in path_template for token in ["%H", "%M", "%S"])

def is_template_daily(path_template):
    return any(token in path_template for token in ["%Y", "%m", "%d", "%j"])

def is_recent_by_mtime(file_path, time_limit):
    try:
        return datetime.fromtimestamp(os.path.getmtime(file_path)) >= time_limit
    except OSError:
        return False
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to parse the file name and extract the time information
def parse_file_time(file_name, filename_regex, filename_time_format=None):

    match = filename_regex.match(file_name)
    if not match:
        return None

    groups = match.groupdict()

    if "stamp" in groups and groups["stamp"] is not None:
        if filename_time_format is None:
            raise ValueError(
                "filename_time_format is required when filename_regex uses group 'stamp'"
            )
        return datetime.strptime(groups["stamp"], filename_time_format)

    if "date" in groups and "hour" in groups:
        return datetime.strptime(groups["date"] + groups["hour"], "%Y%m%d%H")

    raise ValueError(
        "filename_regex must define either named group 'stamp' "
        "or named groups 'date' and 'hour'"
    )
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# method to find variables
def find_var_name(dataset, candidates):
    for name in candidates:
        if name in dataset.variables:
            return name
        if name in dataset.coords:
            return name
    return None

# method to normalize longitudes to the range [-180, 180]
def normalize_longitudes(lon_values):
    lon = np.asarray(lon_values, dtype=float)
    return ((lon + 180.0) % 360.0) - 180.0

# method to compute a boolean mask of points that are inside the specified lat/lon domain
def compute_domain_mask(lat_values, lon_values, lat_min, lat_max, lon_min, lon_max):
    lat_ok = (lat_values >= lat_min) & (lat_values <= lat_max)

    if lon_min <= lon_max:
        lon_ok = (lon_values >= lon_min) & (lon_values <= lon_max)
    else:
        lon_ok = (lon_values >= lon_min) | (lon_values <= lon_max)

    return lat_ok & lon_ok

# method to check if the file contains data points that are over the specified geographic domain
def check_file_over_domain(file_path, settings):

    domain_settings = settings["domain"]
    netcdf_settings = settings["netcdf"]

    # get geographic domain settings
    lat_min = domain_settings.get("lat_min")
    lat_max = domain_settings.get("lat_max")
    lon_min = domain_settings.get("lon_min")
    lon_max = domain_settings.get("lon_max")
    require_mode = domain_settings.get("require", "any").lower()

    # check if are defined or not
    if None in (lat_min, lat_max, lon_min, lon_max):
        info = {
            "n_points": None,
            "n_inside": None,
            "lat_range": (None, None),
            "lon_range": (None, None),
        }
        return True, info

    # convert AFTER check
    lat_min = float(lat_min)
    lat_max = float(lat_max)
    lon_min = float(lon_min)
    lon_max = float(lon_max)

    dataset = None
    try:
        dataset = xr.open_dataset(file_path)

        lat_name = find_var_name(dataset, netcdf_settings["lat_var_candidates"])
        lon_name = find_var_name(dataset, netcdf_settings["lon_var_candidates"])

        if lat_name is None or lon_name is None:
            raise RuntimeError(
                f'Latitude/longitude variables not found in "{os.path.basename(file_path)}"'
            )

        lat_values = np.asarray(dataset[lat_name].values, dtype=float).ravel()
        lon_values = normalize_longitudes(dataset[lon_name].values).ravel()

        valid_mask = np.isfinite(lat_values) & np.isfinite(lon_values)
        lat_values = lat_values[valid_mask]
        lon_values = lon_values[valid_mask]

        if lat_values.size == 0:
            info = {
                "n_points": 0,
                "n_inside": 0,
                "lat_range": (None, None),
                "lon_range": (None, None),
            }
            return False, info

        inside_mask = compute_domain_mask(
            lat_values, lon_values, lat_min, lat_max, lon_min, lon_max
        )

        n_points = int(lat_values.size)
        n_inside = int(np.count_nonzero(inside_mask))

        if require_mode == "all":
            over_domain = n_inside == n_points
        elif require_mode == "any":
            over_domain = n_inside > 0
        else:
            raise ValueError('domain.require must be "any" or "all"')

        info = {
            "n_points": n_points,
            "n_inside": n_inside,
            "lat_range": (float(np.nanmin(lat_values)), float(np.nanmax(lat_values))),
            "lon_range": (float(np.nanmin(lon_values)), float(np.nanmax(lon_values))),
        }

        return over_domain, info

    finally:
        if dataset is not None:
            dataset.close()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to iterate over source patterns based on the settings and reference time, yielding unique search patterns for files to process
def iter_source_patterns(settings, reference_time):
    src_dir = settings["src_dir"]
    look_back_days = int(settings.get("look_back_days", 3))
    look_back_hours = int(settings.get("look_back_hours", look_back_days * 24))

    seen = set()

    # hourly template
    if is_template_hourly(src_dir):
        for hour_step in range(look_back_hours + 1):
            time_step = reference_time - timedelta(hours=hour_step)
            src_path = resolve_time_template(src_dir, time_step)
            search_pattern = os.path.join(src_path, "*H122_C_LIIB_*")
            if search_pattern not in seen:
                seen.add(search_pattern)
                yield search_pattern

    # daily template
    elif is_template_daily(src_dir):
        for day_step in range(look_back_days + 1):
            time_step = reference_time - timedelta(days=day_step)
            src_path = resolve_time_template(src_dir, time_step)
            search_pattern = os.path.join(src_path, "*H122_C_LIIB_*")
            if search_pattern not in seen:
                seen.add(search_pattern)
                yield search_pattern

    # flat folder
    else:
        for day_step in range(look_back_days + 1):
            time_step = reference_time - timedelta(days=day_step)
            date_step = time_step.strftime("%Y%m%d")
            search_pattern = os.path.join(src_dir, f"*H122_C_LIIB_{date_step}*")
            if search_pattern not in seen:
                seen.add(search_pattern)
                yield search_pattern

# method to build the destination path for a file based on the destination template and the file time
def build_destination_path(dst_dir_template, file_time, file_name):
    dst_dir = resolve_time_template(dst_dir_template, file_time)
    return os.path.join(dst_dir, file_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to copy / move the file, creating destination folder if needed
def transfer_file(src_path, dst_path, action="copy", preserve_metadata=True):
    make_folder(os.path.dirname(dst_path))

    if action == "move":
        shutil.move(src_path, dst_path)
    elif action == "copy":
        if preserve_metadata:
            shutil.copy2(src_path, dst_path)
        else:
            shutil.copy(src_path, dst_path)
    else:
        raise ValueError('copy.action must be "copy" or "move"')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to transfer files (copy or move)
def transfer(settings, reference_time):

    stats = {
        "processed_ok": 0,
        "skipped_existing": 0,
        "skipped_old": 0,
        "skipped_no_match": 0,
        "skipped_outside": 0,
        "errors": 0,
    }

    transferred_files = []
    skipped_files = []
    error_files = []

    visited_files = set()

    for source_pattern in iter_source_patterns(settings, reference_time):
        file_list = sorted(glob(source_pattern))

        for src_path in file_list:
            if src_path in visited_files:
                continue
            visited_files.add(src_path)

            if not os.path.isfile(src_path):
                continue

            if not is_recent_by_mtime(src_path, reference_time - timedelta(days=int(settings.get("look_back_days", 3)))):
                stats["skipped_old"] += 1
                skipped_files.append({"file": src_path, "reason": "old"})
                continue

            file_name = os.path.basename(src_path)

            try:
                file_time = parse_file_time(
                    file_name=file_name,
                    filename_regex=re.compile(settings["filename_regex"]),
                    filename_time_format=settings.get("filename_time_format")
                )
            except Exception as exc:
                stats["skipped_no_match"] += 1
                skipped_files.append({"file": src_path, "reason": f"bad_filename_datetime: {exc}"})
                continue

            if file_time is None:
                stats["skipped_no_match"] += 1
                skipped_files.append({"file": src_path, "reason": "regex_no_match"})
                continue

            dst_path = build_destination_path(settings["dst_dir"], file_time, file_name)

            try:
                over_domain, info = check_file_over_domain(src_path, settings)
            except Exception as exc:
                stats["errors"] += 1
                error_files.append({"file": src_path, "error": str(exc)})
                continue

            if not over_domain:
                stats["skipped_outside"] += 1
                skipped_files.append({"file": src_path, "reason": "outside_domain"})
                continue

            if os.path.exists(dst_path) and not settings.get("copy", {}).get("overwrite", False):
                stats["skipped_existing"] += 1
                skipped_files.append({"file": src_path, "reason": "already_exists", "destination": dst_path})
                continue

            try:
                if not settings.get("copy", {}).get("dry_run", False):
                    transfer_file(
                        src_path,
                        dst_path,
                        action=settings.get("copy", {}).get("action", "copy"),
                        preserve_metadata=settings.get("copy", {}).get("preserve_metadata", True)
                    )

                stats["processed_ok"] += 1
                transferred_files.append({"source": src_path, "destination": dst_path})

            except Exception as exc:
                stats["errors"] += 1
                error_files.append({"file": src_path, "error": str(exc)})

    return stats, transferred_files, skipped_files, error_files
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to collect report
def collect_report(
        settings, reference_time, stats, start_time,
        transferred_files=None, skipped_files=None, error_files=None):

    if transferred_files is None:
        transferred_files = []
    if skipped_files is None:
        skipped_files = []
    if error_files is None:
        error_files = []

    include_file_lists = settings.get("report", {}).get("include_file_lists", True)
    if not include_file_lists:
        transferred_files = []
        skipped_files = []
        error_files = []

    elapsed = round(time.time() - start_time, 1)

    report = {
        "algorithm": {
            "name": alg_name,
            "version": alg_version,
            "release": alg_release
        },
        "run": {
            "reference_time": reference_time.strftime("%Y-%m-%d %H:%M:%S"),
            "elapsed_seconds": elapsed
        },
        "paths": {
            "src_dir": settings.get("src_dir"),
            "dst_dir": settings.get("dst_dir")
        },
        "domain": settings.get("domain", {}),
        "copy": settings.get("copy", {}),
        "stats": stats,
        "files": {
            "transferred": transferred_files,
            "skipped": skipped_files,
            "errors": error_files
        }
    }

    return report

# method to write report
def write_report(report, file_path):
    folder_name = os.path.dirname(file_path)
    if folder_name:
        make_folder(folder_name)

    with open(file_path, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=4)

# method to save report
def save_report(settings, report):
    report_settings = settings.get("report", {})
    if not report_settings.get("enabled", False):
        return

    report_folder = report_settings.get("folder", "./report")
    report_filename_template = report_settings.get(
        "filename", "sync_h122_report_%Y%m%d%H%M.json"
    )

    # resolve timestamp directly from template
    execution_time = datetime.now()
    report_filename = execution_time.strftime(report_filename_template)
    report_path = os.path.join(report_folder, report_filename)

    write_report(report, report_path)

    logger.info(' ----> Report saved: ' + report_path)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to get logger
def get_logger(settings):
    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    log_folder = log_settings.get("folder", "./log")
    log_filename = log_settings.get("filename", "transfer.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = datetime.now().strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    # reset handlers (important if script reused)
    logger.handlers = []
    logger.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    # file handler
    fh = logging.FileHandler(log_path)
    fh.setLevel(level)
    fh.setFormatter(formatter)

    # console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to parse the args (cli)
def get_args():
    parser = argparse.ArgumentParser(
        description="Sync H122 files over a configured geographic domain"
    )

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        required=True,
        help="Path to JSON settings file",
    )

    parser.add_argument(
        "-time",
        dest="time_run",
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM". Minutes are rounded down to 00.',
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# main of the script
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()
    # get settings
    settings = read_file_json(args.settings_file)
    # initialize logger
    get_logger(settings)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (start)
    logger.info(' ============================================================================ ')
    logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logger.info(' ==> START ... ')
    logger.info(' ')

    # info settings
    logger.info(f" ---> Settings file:        {args.settings_file}")
    logger.info(f" ---> Source location:      {settings.get('src_dir')}")
    logger.info(f" ---> Destination location: {settings.get('dst_dir')}")

    # time algorithm
    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get reference time - start
    logger.info(f" ---> Set time ... ")
    try:

        # get reference time
        reference_time = parse_reference_time(args.time_run)

        logger.info(f" ----> Get reference time: {reference_time.strftime('%Y-%m-%d %H:%M:%S')}")
        if args.time_run:
            logger.info(f" ----> Set by user: {args.time_run} (rounded to hour)")
        else:
            logger.info(" ----> Set by system (rounded to hour)")

        # get reference time - end (done)
        logger.info(f" ---> Set time ... DONE")

    except Exception as exc:

        # get reference time - end (failed)
        logger.error(f" ===> ERROR: parsing time: {exc}")
        logger.info(f" ---> Set time ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # execute transfer - start
    logger.info(' ---> Execute transfer ... ')
    try:

        # run workflow to sync datasets
        stats, transferred_files, skipped_files, error_files = transfer(
            settings, reference_time
        )

        # execute transfer - end (done)
        logger.info(' ---> Execute transfer ... DONE')

    except Exception as exc:

        # execute transfer - end (failed)
        logger.error(f" ===> ERROR in executing transfer algorithm: {exc}")
        logger.info(' ---> Execute transfer ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # save report - start
    logger.info(" ---> Save report ...")
    try:
        # collect report
        report = collect_report(
            settings=settings,
            reference_time=reference_time,
            stats=stats,
            start_time=start_time,
            transferred_files=transferred_files,
            skipped_files=skipped_files,
            error_files=error_files
        )

        # save report
        save_report(settings, report)

        # save report - end (done)
        logger.info(" ---> Save report ... DONE")

    except Exception as exc:

        # save report - end (failed)
        logger.error(f" ===> ERROR saving report: {exc}")
        logger.info(" ---> Save report ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    logger.info(' ')
    logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    logger.info(' ==> ... END')
    logger.info(' ==> Bye, Bye')
    logger.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper for running as script
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
