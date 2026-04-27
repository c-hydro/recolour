#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SSM H122 DOWNLOADER - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20260415'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python ecmwf_downloader_rzsm.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20260415 (1.0.0) --> First development
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import re
import sys
import time
import json
import atexit
import signal
import logging
import argparse
import subprocess

from datetime import datetime
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for downloading ssm h122 files'
alg_type = 'Package'
alg_version = '1.0.0'
alg_release = '2026-04-15'

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
        description='Mirror or date-filter FTP products using lftp'
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
    log_filename = log_settings.get("filename", "downloader.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = datetime.now().strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    # reset handlers (important if script reused)
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

# execute command--
def execute_command(command_string):

    process = subprocess.run(
        command_string,
        shell=True,
        text=True,
        capture_output=True
    )

    return process.returncode, process.stdout, process.stderr
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# time utils
# method to parse time string
def parse_time_string(time_string):

    if time_string is None:
        return None

    try:
        time_obj = datetime.strptime(time_string, '%Y-%m-%d %H:%M')
    except ValueError as exc:
        raise ValueError(
            f'time "{time_string}" must have format "YYYY-MM-DD HH:MM"'
        ) from exc

    return time_obj

# method to round time to hour
def round_time_to_hour(time_obj):

    time_obj = time_obj.replace(minute=0, second=0, microsecond=0)

    return time_obj

# get reference time
def get_reference_time(args, settings):

    mirror_settings = settings.get('mirror', {})
    round_to_hour = mirror_settings.get('reference_time_round_to_hour', True)

    if args.time_run is not None:
        reference_time = parse_time_string(args.time_run)
    else:
        reference_time = datetime.now()

    if round_to_hour:
        reference_time = round_time_to_hour(reference_time)

    return reference_time

# get time window
def get_time_window(settings):

    time_settings = settings.get('time', {})

    time_start_raw = time_settings.get('time_start', None)
    time_end_raw = time_settings.get('time_end', None)

    time_start = parse_time_string(time_start_raw) if time_start_raw is not None else None
    time_end = parse_time_string(time_end_raw) if time_end_raw is not None else None

    return time_start, time_end

# method to parse file time
def parse_file_time(file_name, settings):

    filter_settings = settings.get('filter', {})
    file_name_regex = re.compile(filter_settings['filename_regex'])
    file_time_format = filter_settings['filename_time_format']

    match = file_name_regex.match(file_name)
    if match is None:
        return None

    file_time = datetime.strptime(match.group('stamp'), file_time_format)

    return file_time
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# run utils

# method to select run mode
def select_run_mode(args, settings):

    time_start, time_end = get_time_window(settings)
    reference_time = get_reference_time(args, settings)
    n_days = int(settings.get('mirror', {}).get('n_days', 3))

    if (time_start is not None) and (time_end is not None):
        mode = 'date_filter'
    else:
        mode = 'mirror'

    return mode, reference_time, time_start, time_end, n_days

# method to acquire lock-
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
# ftp utils

# method to check ftp session
def check_ftp_session(settings):

    ftp_settings = settings.get('ftp', {})

    ftp_url = ftp_settings['ftp_url']
    netrc_machine_label = ftp_settings['netrc_machine_label']
    remote_check_path = ftp_settings['check_path']
    timeout = ftp_settings.get('timeout', 20)
    max_retries = ftp_settings.get('max_retries', 2)
    reconnect_interval_base = ftp_settings.get('reconnect_interval_base', 5)
    ssl_verify = ftp_settings.get('ssl_verify_certificate', False)

    ssl_verify_value = 'yes' if ssl_verify else 'no'

    home_path = os.path.expanduser('~')
    netrc_file = os.path.join(home_path, '.netrc')
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    alg_logger.info(' ')
    alg_logger.info(' ::: CHECK SESSION: FTP AUTHENTICATION + REMOTE PATH ACCESS')
    alg_logger.info(f' ::: Host           : {ftp_url}')
    alg_logger.info(f' ::: Netrc machine  : {netrc_machine_label}')
    alg_logger.info(f' ::: Check path     : {remote_check_path}')
    alg_logger.info(f' ::: HOME           : {home_path}')
    alg_logger.info(f' ::: Netrc file     : {netrc_file}')
    alg_logger.info(f' ::: Timestamp      : {timestamp}')
    alg_logger.info(' ')

    if not os.path.isfile(netrc_file):
        alg_logger.error(f' ===> ERROR: .netrc file not found at {netrc_file}')
        return False

    lftp_script = f"""
set net:timeout {timeout}
set net:max-retries {max_retries}
set net:reconnect-interval-base {reconnect_interval_base}
set ssl:verify-certificate {ssl_verify_value}
set cmd:fail-exit yes

pwd
cls /products
cls /products/h122
cd {remote_check_path}
pwd
bye
""".strip()

    command = f"lftp {ftp_url} <<'EOF'\n{lftp_script}\nEOF"

    return_code, stdout, stderr = execute_command(command)

    if stdout:
        for line in stdout.splitlines():
            alg_logger.info(line)

    if stderr:
        for line in stderr.splitlines():
            alg_logger.error(line)

    alg_logger.info(' ')

    if return_code != 0:
        alg_logger.error(' ===> ERROR: FTP session check failed.')
        alg_logger.error(' ===> Possible causes:')
        alg_logger.error(' ===>   - wrong username/password in .netrc')
        alg_logger.error(' ===>   - .netrc not being used by lftp')
        alg_logger.error(' ===>   - remote folder access denied')
        return False

    if remote_check_path in stdout:
        alg_logger.info(' ----> OK: FTP session established and remote path is accessible.')
        alg_logger.info(' ----> CHECK SESSION PASSED')
        return True

    alg_logger.error(' ===> ERROR: Authentication may have worked, but remote path access was not confirmed.')
    return False

# method to list remote files
def list_remote_files(settings, remote_folder):

    ftp_settings = settings.get('ftp', {})

    ftp_url = ftp_settings['ftp_url']
    timeout = ftp_settings.get('timeout', 30)
    max_retries = ftp_settings.get('max_retries', 5)
    reconnect_interval_base = ftp_settings.get('reconnect_interval_base', 5)
    ssl_verify = ftp_settings.get('ssl_verify_certificate', False)

    ssl_verify_value = 'yes' if ssl_verify else 'no'

    lftp_script = f"""
open {ftp_url}
set net:timeout {timeout}
set net:max-retries {max_retries}
set net:reconnect-interval-base {reconnect_interval_base}
set ssl:verify-certificate {ssl_verify_value}
set cmd:fail-exit yes

cls -1 "{remote_folder}"
quit
""".strip()

    command = f"lftp <<'EOF'\n{lftp_script}\nEOF"

    return_code, stdout, stderr = execute_command(command)

    if stderr:
        for line in stderr.splitlines():
            alg_logger.error(line)

    if return_code != 0:
        raise RuntimeError(f' ===> Unable to list remote folder: {remote_folder}')

    file_list = []
    if stdout:
        for line in stdout.splitlines():
            file_name = line.strip()
            if file_name != '':
                file_list.append(file_name)

    return file_list

# method to filter remote files by time window
def filter_remote_files_by_time(settings, file_list, time_start, time_end):

    selected_files = []
    skipped_files = []

    for file_name in file_list:

        file_time = parse_file_time(file_name, settings)
        if file_time is None:
            skipped_files.append(file_name)
            continue

        if time_start <= file_time <= time_end:
            selected_files.append((file_name, file_time))

    return selected_files, skipped_files

# method to download remote files
def download_remote_files(settings, remote_folder, selected_files):

    ftp_settings = settings.get('ftp', {})
    mirror_settings = settings.get('mirror', {})

    ftp_url = ftp_settings['ftp_url']
    timeout = ftp_settings.get('timeout', 30)
    max_retries = ftp_settings.get('max_retries', 5)
    reconnect_interval_base = ftp_settings.get('reconnect_interval_base', 5)
    ssl_verify = ftp_settings.get('ssl_verify_certificate', False)

    local_folder_mirror = mirror_settings['local_folder_mirror']

    ssl_verify_value = 'yes' if ssl_verify else 'no'

    product_name = os.path.basename(remote_folder.rstrip('/'))
    local_dir = os.path.join(local_folder_mirror.rstrip('/'), product_name)
    make_folder(local_dir)

    if len(selected_files) == 0:
        alg_logger.warning(f' ===> No files selected for product: {product_name}')
        return

    lftp_lines = [
        f'open {ftp_url}',
        f'set net:timeout {timeout}',
        f'set net:max-retries {max_retries}',
        f'set net:reconnect-interval-base {reconnect_interval_base}',
        f'set ssl:verify-certificate {ssl_verify_value}',
        ''
    ]

    start_ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    lftp_lines.extend([
        'echo =================================================================',
        f'echo PRODUCT: {product_name}',
        f'echo Remote : {remote_folder}',
        f'echo Local  : {local_dir}',
        f'echo Start  : {start_ts}',
        'echo -----------------------------------------------------------------'
    ])

    for file_name, file_time in selected_files:
        remote_file = f'{remote_folder.rstrip("/")}/{file_name}'
        local_file = os.path.join(local_dir, file_name)

        lftp_lines.append(f'echo GET FILE: {file_name}')
        lftp_lines.append(f'get -c "{remote_file}" -o "{local_file}"')

    lftp_lines.extend([
        'echo -----------------------------------------------------------------',
        f'echo PRODUCT DONE: {product_name}',
        'echo =================================================================',
        '',
        'quit'
    ])

    lftp_script = '\n'.join(lftp_lines)
    command = f"lftp <<'EOF'\n{lftp_script}\nEOF"

    return_code, stdout, stderr = execute_command(command)

    if stdout:
        for line in stdout.splitlines():
            alg_logger.info(line)

    if stderr:
        for line in stderr.splitlines():
            alg_logger.error(line)

    if return_code != 0:
        raise RuntimeError(f' ===> Download failed for remote folder: {remote_folder}')

# method to download mirror mode
def download_mode_mirror(settings, n_days):

    alg_logger.info(' ----> Running mirror mode')

    ftp_settings = settings.get('ftp', {})
    mirror_settings = settings.get('mirror', {})
    product_settings = settings.get('products', {})

    ftp_url = ftp_settings['ftp_url']
    remote_folders = product_settings['remote_folders']
    local_folder_mirror = mirror_settings['local_folder_mirror']

    timeout = ftp_settings.get('timeout', 30)
    max_retries = ftp_settings.get('max_retries', 5)
    reconnect_interval_base = ftp_settings.get('reconnect_interval_base', 5)
    ssl_verify = ftp_settings.get('ssl_verify_certificate', False)

    parallel_transfer_count = mirror_settings.get('parallel_transfer_count', 4)
    use_pget_n = mirror_settings.get('use_pget_n', 4)

    ssl_verify_value = 'yes' if ssl_verify else 'no'

    mirror_flags = [f'--newer-than={n_days}d']

    if mirror_settings.get('only_newer', True):
        mirror_flags.append('--only-newer')
    if mirror_settings.get('continue_download', True):
        mirror_flags.append('--continue')
    if mirror_settings.get('no_empty_dirs', True):
        mirror_flags.append('--no-empty-dirs')
    if mirror_settings.get('verbose', True):
        mirror_flags.append('--verbose')

    mirror_flags_string = ' '.join(mirror_flags)

    lftp_lines = [
        f'open {ftp_url}',
        f'set net:timeout {timeout}',
        f'set net:max-retries {max_retries}',
        f'set net:reconnect-interval-base {reconnect_interval_base}',
        f'set ssl:verify-certificate {ssl_verify_value}',
        f'set mirror:parallel-transfer-count {parallel_transfer_count}',
        f'set mirror:use-pget-n {use_pget_n}',
        ''
    ]

    for remote_folder in remote_folders:

        product_name = os.path.basename(remote_folder.rstrip('/'))
        local_dir = os.path.join(local_folder_mirror.rstrip('/'), product_name)
        make_folder(local_dir)

        start_ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        lftp_lines.extend([
            'echo =================================================================',
            f'echo PRODUCT: {product_name}',
            f'echo Remote : {remote_folder}',
            f'echo Local  : {local_dir}',
            f'echo Start  : {start_ts}',
            'echo -----------------------------------------------------------------',
            f'mirror {mirror_flags_string} "{remote_folder}" "{local_dir}"',
            'echo -----------------------------------------------------------------',
            f'echo PRODUCT DONE: {product_name}',
            'echo =================================================================',
            ''
        ])

    lftp_lines.append('quit')

    lftp_script = '\n'.join(lftp_lines)
    command = f"lftp <<'EOF'\n{lftp_script}\nEOF"

    return_code, stdout, stderr = execute_command(command)

    if stdout:
        for line in stdout.splitlines():
            alg_logger.info(line)

    if stderr:
        for line in stderr.splitlines():
            alg_logger.error(line)

    if return_code != 0:
        raise RuntimeError('Mirror mode failed')

# method to download date filter mode
def download_mode_date_filter(settings, time_start, time_end):

    product_settings = settings.get('products', {})
    remote_folders = product_settings['remote_folders']

    alg_logger.info(' ----> Running date-filter mode')
    alg_logger.info(f' ::: Time start: {time_start.strftime("%Y-%m-%d %H:%M:%S")}')
    alg_logger.info(f' ::: Time end  : {time_end.strftime("%Y-%m-%d %H:%M:%S")}')

    for remote_folder in remote_folders:

        product_name = os.path.basename(remote_folder.rstrip('/'))

        alg_logger.info(' ')
        alg_logger.info(f' ::: PRODUCT: {product_name}')
        alg_logger.info(f' ::: LIST REMOTE FILES: {remote_folder}')

        file_list = list_remote_files(settings, remote_folder)

        alg_logger.info(f' ::: Remote files found: {len(file_list)}')

        selected_files, skipped_files = filter_remote_files_by_time(
            settings=settings,
            file_list=file_list,
            time_start=time_start,
            time_end=time_end
        )

        alg_logger.info(f' ::: Remote files selected: {len(selected_files)}')
        alg_logger.info(f' ::: Remote files skipped (no valid timestamp): {len(skipped_files)}')

        download_remote_files(
            settings=settings,
            remote_folder=remote_folder,
            selected_files=selected_files
        )
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
    # info algorithm (start)
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
    alg_logger.info(f" ---> Select download mode ... ")
    try:
        # select download mode
        mode, reference_time, time_start, time_end, n_days = select_run_mode(args, settings)
    except Exception as exc:
        # select download mode - end (failed)
        alg_logger.error(f' ===> Error in selecting download mode: {exc}')
        alg_logger.info(f" ---> Select download mode ... FAILED")
        sys.exit(1)

    # info run mode
    alg_logger.info(' ::: Reference time: ' + reference_time.strftime('%Y-%m-%d %H:%M:%S'))
    if mode == 'date_filter':
        alg_logger.info(' ::: Download mode: date_filter')
        alg_logger.info(' ::: Time start: ' + time_start.strftime('%Y-%m-%d %H:%M:%S'))
        alg_logger.info(' ::: Time end:   ' + time_end.strftime('%Y-%m-%d %H:%M:%S'))
        alg_logger.info(' ::: Time window taken from JSON settings')
    else:
        alg_logger.info(' ::: Download mode: mirror')
        alg_logger.info(' ::: Days back: ' + str(n_days))
        if args.time_run is not None:
            alg_logger.info(' ::: Reference time taken from command line -time')
        else:
            alg_logger.info(' ::: Reference time taken from system time NOW')

    # select download mode - end (done)
    alg_logger.info(f" ---> Select download mode ... DONE")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # initialize download mode - start
    ftp_settings = settings.get('ftp', {})
    alg_logger.info(' ---> Initialize download mode ... ')
    alg_logger.info(' ')

    # acquire lock
    if not acquire_lock(settings=settings, force_lock=args.force_lock):
        sys.exit(1)

    # check ftp session
    if not check_ftp_session(settings=settings):
        # initialize download mode - end (failed - errors in ftp session validation)
        alg_logger.error(' ===> Aborting before download because session validation failed.')
        alg_logger.info(' ---> Initialize download mode ... FAILED')
        sys.exit(1)
    else:
        # initialize download mode - end (done)
        alg_logger.info(' ---> Initialize download mode ... DONE')

    # execute the download mode - start
    alg_logger.info(' ---> Execute download mode ... ')
    alg_logger.info(' ::: Timestamp Start -- ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    alg_logger.info(' ::: Host -- ' + ftp_settings.get('ftp_url', 'undefined'))
    try:

        # execute mode
        if mode == 'date_filter':

            # info about mode
            alg_logger.info(' ::: Mode -- filter by dates')
            if time_start > time_end:
                raise RuntimeError('time_start must be less than or equal to time_end')
            # execute date filter
            download_mode_date_filter(settings=settings, time_start=time_start, time_end=time_end)

        else:
            # info about mode
            alg_logger.info(' ::: Mode -- mirror')
            # execute mirror
            download_mode_mirror(settings=settings, n_days=n_days)

        # execute the download mode - end (done)
        alg_logger.info(' ::: Timestamp End -- ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        alg_logger.info(' ---> Execute download mode ... DONE')

    except Exception as exc:

        # execute the download mode - end (failed)
        alg_logger.error(f" ===> ERROR in executing download algorithm: {exc}")
        alg_logger.info(' ---> Execute download mode ... FAILED')
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # info algorithm (end)
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
