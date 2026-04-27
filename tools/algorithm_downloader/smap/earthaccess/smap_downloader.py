#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SMAP SPL2SMP_E DOWNLOADER

General command line:
python smap_downloader_spl2smp_e.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import warnings
import re
import sys
import time
import json
import atexit
import signal
import logging
import argparse

from pathlib import Path
from datetime import datetime, timedelta

import earthaccess

# suppress earthaccess warnings
warnings.filterwarnings("ignore",category=FutureWarning,module="earthaccess.store")
warnings.filterwarnings("ignore",message=r".*DataGranule\.size.*",category=FutureWarning)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for downloading SMAP SPL2SMP_E files'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2026-04-24'

# algorithm globals
alg_logger = logging.getLogger('app_downloader')
acquired_lock = None
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# basic utils

# method to read file json
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f' ===> File "{file_name}" not found. Exit')


# method to get args
def get_args():

    parser = argparse.ArgumentParser(
        description='Download SMAP SPL2SMP_E granules using earthaccess'
    )

    parser.add_argument(
        '-settings_file',
        dest='settings_file',
        required=True,
        help='Path to JSON settings file'
    )

    parser.add_argument(
        '-f', '-force',
        dest='force_lock',
        action='store_true',
        help='Force removal of stale lock files'
    )

    parser.add_argument(
        '-time',
        dest='time_run',
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM"'
    )

    return parser.parse_args()


# method to get logger
def get_logger(settings):
    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    log_folder = log_settings.get("folder", "./log")
    log_filename = log_settings.get("filename", "smap_downloader.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = datetime.now().strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    # reset handlers
    alg_logger.handlers = []
    alg_logger.setLevel(level)

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

    alg_logger.addHandler(fh)
    alg_logger.addHandler(ch)


# make folder
def make_folder(folder_name):

    if folder_name is not None and folder_name != '':
        os.makedirs(folder_name, exist_ok=True)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# time utils

# method to parse time string
def parse_time_string(time_string, time_format='%Y-%m-%d %H:%M'):

    if time_string is None:
        return None

    try:
        time_obj = datetime.strptime(time_string, time_format)
    except ValueError as exc:
        raise ValueError(
            f'time "{time_string}" must have format "{time_format}"'
        ) from exc

    return time_obj


# method to parse date string
def parse_date_string(date_string):

    if date_string is None:
        return None

    try:
        date_obj = datetime.strptime(date_string, '%Y-%m-%d')
    except ValueError as exc:
        raise ValueError(
            f'date "{date_string}" must have format "YYYY-MM-DD"'
        ) from exc

    return date_obj


# method to round time to day
def round_time_to_day(time_obj):

    time_obj = time_obj.replace(hour=0, minute=0, second=0, microsecond=0)

    return time_obj


# get reference time
def get_reference_time(args, settings):

    downloader_settings = settings.get('downloader', {})
    round_to_day = downloader_settings.get('reference_time_round_to_day', True)

    if args.time_run is not None:
        reference_time = parse_time_string(args.time_run)
    else:
        reference_time = datetime.now()

    if round_to_day:
        reference_time = round_time_to_day(reference_time)

    return reference_time


# get time window
def get_time_window(settings):

    time_settings = settings.get('time', {})

    time_start_raw = time_settings.get('time_start', None)
    time_end_raw = time_settings.get('time_end', None)

    time_start = parse_date_string(time_start_raw) if time_start_raw is not None else None
    time_end = parse_date_string(time_end_raw) if time_end_raw is not None else None

    return time_start, time_end


# method to iterate over dates
def daterange(start_date, end_date):

    current_date = start_date

    while current_date <= end_date:
        yield current_date
        current_date += timedelta(days=1)

# method to extract datetime from filename
def extract_datetime_from_filename(filename):

    match = re.search(r'(\d{8}T\d{6})', filename)
    if match is None:
        return None

    return datetime.strptime(match.group(1), "%Y%m%dT%H%M%S")

# method to fill time template
def fill_time_template(path_template, time_obj):

    if "{time_step" not in path_template:
        return path_template

    pattern = r"\{time_step:\s*(.*?)\}"

    def replace(match):
        time_format = match.group(1)
        return time_obj.strftime(time_format)

    return re.sub(pattern, replace, path_template)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# run utils

# method to select run mode
def select_run_mode(args, settings):

    time_start, time_end = get_time_window(settings)
    reference_time = get_reference_time(args, settings)

    downloader_settings = settings.get('downloader', {})
    n_days = int(downloader_settings.get('n_days', 3))

    if (time_start is not None) and (time_end is not None):
        mode = 'date_range'
    else:
        mode = 'rolling_window'
        time_end = reference_time
        time_start = reference_time - timedelta(days=n_days)

    return mode, reference_time, time_start, time_end, n_days


# method to acquire lock
def acquire_lock(settings, force_lock=False):

    global acquired_lock

    lock_settings = settings.get('lock', {})
    lock_dir = lock_settings['folder']
    lock_basename = lock_settings['basename']
    max_slots = int(lock_settings.get('max_slots', 1))

    make_folder(lock_dir)

    for lock_id in range(1, max_slots + 1):

        lock_file = os.path.join(lock_dir, f'{lock_basename}{lock_id}')

        if os.path.exists(lock_file) and force_lock:
            alg_logger.warning(f' ===> Forcing removal of stale lock: {lock_file}')
            try:
                os.remove(lock_file)
            except OSError as exc:
                alg_logger.error(f' ===> Error removing stale lock "{lock_file}": {exc}')
                continue

        if os.path.exists(lock_file):
            continue

        try:
            file_descriptor = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            with os.fdopen(file_descriptor, 'w') as file_handle:
                file_handle.write(str(os.getpid()) + '\n')

            acquired_lock = lock_file
            alg_logger.info(f' ----> Acquired slot {lock_id} using {lock_file}')
            return True

        except FileExistsError:
            continue
        except OSError as exc:
            alg_logger.error(f' ===> Error creating lock "{lock_file}": {exc}')

    alg_logger.error(f' ===> All {max_slots} slots are currently in use. Exiting.')
    return False


# method to release lock
def release_lock():

    global acquired_lock

    if acquired_lock is not None:
        try:
            if os.path.exists(acquired_lock):
                os.remove(acquired_lock)
                alg_logger.info(f' ----> Released slot ({acquired_lock})')
        except OSError as exc:
            alg_logger.error(f' ===> Error releasing lock "{acquired_lock}": {exc}')
        finally:
            acquired_lock = None


# method to handle signal
def handle_signal(signum, frame):
    alg_logger.warning(f' ===> Caught signal {signum}. Releasing lock and exiting.')
    release_lock()
    sys.exit(1)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# smap utils

# method to check earthdata session
def check_earthdata_session(settings):

    auth_settings = settings.get('earthdata', {})
    login_strategy = auth_settings.get('strategy', None)

    alg_logger.info(' ')
    alg_logger.info(' ::: CHECK SESSION: EARTHDATA AUTHENTICATION')
    alg_logger.info(f' ::: Strategy: {login_strategy if login_strategy is not None else "default"}')
    alg_logger.info(' ')

    try:
        if login_strategy is None:
            auth = earthaccess.login()
        else:
            auth = earthaccess.login(strategy=login_strategy)
    except Exception as exc:
        alg_logger.error(f' ===> ERROR: Earthdata authentication failed: {exc}')
        return False

    if auth is None:
        alg_logger.error(' ===> ERROR: Earthdata authentication failed.')
        return False

    alg_logger.info(' ----> OK: Earthdata authentication completed.')
    alg_logger.info(' ----> CHECK SESSION PASSED')

    return True


# method to extract filename
def extract_filename(granule):

    try:
        links = granule.data_links(access="external")
        if links:
            return Path(links[0]).name
    except Exception:
        pass

    try:
        links = granule.data_links()
        if links:
            return Path(links[0]).name
    except Exception:
        pass

    try:
        return granule["meta"]["native-id"]
    except Exception:
        return None


# method to search smap granules
def search_smap_granules(settings, day):

    product_settings = settings.get('product', {})
    spatial_settings = settings.get('spatial', {})

    short_name = product_settings.get('short_name', 'SPL2SMP_E')
    version = product_settings.get('version', '006')
    bounding_box = tuple(spatial_settings.get('bounding_box', [6.0, 36.0, 19.0, 47.5]))

    next_day = day + timedelta(days=1)

    results = earthaccess.search_data(
        short_name=short_name,
        version=version,
        temporal=(
            day.strftime('%Y-%m-%d 00:00:00'),
            next_day.strftime('%Y-%m-%d 00:00:00'),
        ),
        bounding_box=bounding_box,
    )

    return results


# method to define data folder
def define_data_folder(settings, time_obj):

    downloader_settings = settings.get('downloader', {})
    data_template = downloader_settings['data_folder']

    # apply time formatting
    data_folder = fill_time_template(data_template, time_obj)

    output_folder = Path(data_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    return output_folder


# method to filter existing granules
def filter_existing_granules(granules, output_folder):

    selected_granules = []
    skipped_granules = []

    for granule in granules:

        filename = extract_filename(granule)

        if filename is None:
            selected_granules.append(granule)
            continue

        target_file = output_folder / filename

        if target_file.exists():
            skipped_granules.append(filename)
        else:
            selected_granules.append(granule)

    return selected_granules, skipped_granules

# method to download smap granules
def download_smap_granules(settings, granules):

    downloaded_files = []
    skipped_files = []

    for granule in granules:

        filename = extract_filename(granule)

        if filename is None:
            alg_logger.warning(' ===> Filename not found, skipping granule')
            continue

        file_time = extract_datetime_from_filename(filename)

        if file_time is None:
            alg_logger.warning(f' ===> Time not found in filename: {filename}')
            continue

        # round file time to hour
        file_time = file_time.replace(minute=0, second=0, microsecond=0)

        # define output folder using file time
        output_folder = define_data_folder(settings, file_time)

        target_file = output_folder / filename

        if target_file.exists():
            alg_logger.info(f' ----> Skip existing: {filename}')
            skipped_files.append(filename)
            continue

        alg_logger.info(f' ----> Downloading: {filename}')
        alg_logger.info(f' ----> Output folder: {output_folder}')

        file_list = earthaccess.download(
            [granule],
            str(output_folder)
        )

        for file_path in file_list:
            alg_logger.info(f' ----> Downloaded: {Path(file_path).name}')
            downloaded_files.append(file_path)

    return downloaded_files, skipped_files

# method to download date range
def download_mode_date_range(settings, time_start, time_end):

    product_settings = settings.get('product', {})
    short_name = product_settings.get('short_name', 'SPL2SMP_E')
    version = product_settings.get('version', '006')

    total_found = 0
    total_downloaded = 0
    total_skipped = 0

    alg_logger.info(' ----> Running date-range mode')
    alg_logger.info(f' ::: Product: {short_name}.{version}')
    alg_logger.info(f' ::: Time start: {time_start.strftime("%Y-%m-%d")}')
    alg_logger.info(f' ::: Time end  : {time_end.strftime("%Y-%m-%d")}')

    for day in daterange(time_start, time_end):

        day_str = day.strftime('%Y-%m-%d')

        alg_logger.info(' ')
        alg_logger.info(f' ::: DAY: {day_str}')
        alg_logger.info(f' ::: Searching {short_name}.{version} granules ...')

        granules = search_smap_granules(settings, day)

        n_found = len(granules)
        total_found += n_found

        alg_logger.info(f' ::: Granules found: {n_found}')

        if not granules:
            continue

        downloaded_files, skipped_files = download_smap_granules(
            settings=settings,
            granules=granules
        )

        total_downloaded += len(downloaded_files)
        total_skipped += len(skipped_files)

    alg_logger.info(' ')
    alg_logger.info(' ::: DOWNLOAD SUMMARY')
    alg_logger.info(f' ::: Total granules found   : {total_found}')
    alg_logger.info(f' ::: Total downloaded       : {total_downloaded}')
    alg_logger.info(f' ::: Total skipped existing : {total_skipped}')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # set global variables
    global alg_name, alg_version, alg_release
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()
    # get settings
    settings = read_file_json(args.settings_file)
    # initialize logger
    get_logger(settings)

    # set algorithm information
    alg_info = settings.get('algorithm', {})
    alg_name = alg_info.get('name', alg_name)
    alg_version = alg_info.get('version', alg_version)
    alg_release = alg_info.get('release', alg_release)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm start
    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    # info settings
    alg_logger.info(' ---> Settings file: ' + args.settings_file)

    # time algorithm
    start_time = time.time()

    # register release handlers
    atexit.register(release_lock)
    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # select download mode - start
    alg_logger.info(' ---> Select download mode ... ')
    try:
        mode, reference_time, time_start, time_end, n_days = select_run_mode(args, settings)
    except Exception as exc:
        alg_logger.error(f' ===> Error in selecting download mode: {exc}')
        alg_logger.info(' ---> Select download mode ... FAILED')
        sys.exit(1)

    alg_logger.info(' ::: Reference time: ' + reference_time.strftime('%Y-%m-%d %H:%M:%S'))
    alg_logger.info(' ::: Download mode: ' + mode)
    alg_logger.info(' ::: Time start: ' + time_start.strftime('%Y-%m-%d'))
    alg_logger.info(' ::: Time end:   ' + time_end.strftime('%Y-%m-%d'))

    if mode == 'rolling_window':
        alg_logger.info(' ::: Days back: ' + str(n_days))

    alg_logger.info(' ---> Select download mode ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # initialize download mode - start
    alg_logger.info(' ---> Initialize download mode ... ')
    alg_logger.info(' ')

    if not acquire_lock(settings=settings, force_lock=args.force_lock):
        sys.exit(1)

    if not check_earthdata_session(settings=settings):
        alg_logger.error(' ===> Aborting before download because session validation failed.')
        alg_logger.info(' ---> Initialize download mode ... FAILED')
        sys.exit(1)
    else:
        alg_logger.info(' ---> Initialize download mode ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # execute download mode - start
    alg_logger.info(' ---> Execute download mode ... ')
    alg_logger.info(' ::: Timestamp Start -- ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    try:

        if time_start > time_end:
            raise RuntimeError('time_start must be less than or equal to time_end')

        download_mode_date_range(
            settings=settings,
            time_start=time_start,
            time_end=time_end
        )

        alg_logger.info(' ::: Timestamp End -- ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        alg_logger.info(' ---> Execute download mode ... DONE')

    except Exception as exc:

        alg_logger.error(f' ===> ERROR in executing download algorithm: {exc}')
        alg_logger.info(' ---> Execute download mode ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm end
    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(' ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    alg_logger.info(' ==> ... END')
    alg_logger.info(' ==> Bye, Bye')
    alg_logger.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# call script from external library
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------------------------------