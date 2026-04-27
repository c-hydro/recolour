#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SSM H122 TRANSFER - EXECUTION WRAPPER

__date__ = '20260420'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_transfer_ssm_h122_wrapper.py
python app_transfer_ssm_h122_wrapper.py -time_now "2026-04-20 00:00" -time_period 2
python app_transfer_ssm_h122_wrapper.py -time_start "2026-04-18 00:00" -time_end "2026-04-20 00:00"

Rules:
- if both -time_start and -time_end are provided, they are used
- otherwise the script uses:
    -time_now (default: system time, rounded to hour)
    -time_period (default: 2 days)
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import warnings
import sys
import time
import argparse
import subprocess
from datetime import datetime, timedelta

# suppress warnings
warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated as an API"
)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Execution wrapper for ssm h122 transfer over a time period'
alg_type = 'Exec'
alg_version = '1.1.0'
alg_release = '2026-04-20'
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# default script and settings file (can be overridden by command line arguments)
default_script_file = "app_transfer_ssm_h122.py"
default_settings_file = "app_transfer_ssm_h122.json"

# default environment settings
virtual_env_folder = '/home/fabio/Documents/Work_Area/Code_Development/Workspace/recolour/conda_p311/bin/'
virtual_env_name = 'recolour_converter_libraries'
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# time utils
# helper to normalize time according to selected round option
def normalize_time(time_obj, round_option="day"):
    round_option = round_option.lower()

    if round_option == "day":
        return time_obj.replace(hour=0, minute=0, second=0, microsecond=0)

    elif round_option in ["hour", "hour_floor"]:
        return time_obj.replace(minute=0, second=0, microsecond=0)

    elif round_option == "hour_ceil":
        time_floor = time_obj.replace(minute=0, second=0, microsecond=0)
        if time_obj == time_floor:
            return time_floor
        return time_floor + timedelta(hours=1)

    else:
        raise ValueError(
            'round_option must be one of: "day", "hour", "hour_floor", "hour_ceil"'
        )


# helper to parse time string
def parse_time(time_string, round_option="day"):
    try:
        time_obj = datetime.strptime(time_string, "%Y-%m-%d %H:%M")
    except ValueError as exc:
        raise ValueError('Time must have format "YYYY-MM-DD HH:MM"') from exc
    return normalize_time(time_obj, round_option=round_option)


# helper to get current time or parse provided time string
def get_time_now(time_string=None, round_option="day"):
    if time_string is None:
        return normalize_time(datetime.now(), round_option=round_option)
    return parse_time(time_string, round_option=round_option)


# helper to compute time range
def get_time_range(time_now=None, time_period=2, time_start=None, time_end=None,
                   round_option="day"):
    """
    Priority:
    1. if time_start and time_end are both set -> use them
    2. otherwise use time_now and time_period

    round_option:
    - "day"        -> round to midnight
    - "hour"       -> round down to hour
    - "hour_floor" -> round down to hour
    - "hour_ceil"  -> round up to next hour
    """
    if (time_start is not None) and (time_end is not None):
        time_from = parse_time(time_start, round_option=round_option)
        time_to = parse_time(time_end, round_option=round_option)

    elif (time_start is None) and (time_end is None):
        time_to = get_time_now(time_now, round_option=round_option)
        time_from = time_to - timedelta(days=int(time_period))

    else:
        raise RuntimeError(
            'time_start and time_end must be both provided, or neither provided'
        )

    if time_from > time_to:
        raise RuntimeError('time_start must be less than or equal to time_end')

    return time_from, time_to

# helper to build time steps
def iter_time_steps(time_start, time_end, step_hours=24, reverse=False):
    if step_hours <= 0:
        raise ValueError("step_hours must be > 0")

    time_steps = []
    time_step = time_start

    while time_step <= time_end:
        time_steps.append(time_step)
        time_step += timedelta(hours=int(step_hours))

    if reverse:
        time_steps = list(reversed(time_steps))

    return time_steps
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# environment utils

# helper to update environment variables for subprocess
def update_env(base_env, script_folder, env_folder=None, env_name=None):
    env = base_env.copy()

    # add script folder to PYTHONPATH
    pypath_old = env.get("PYTHONPATH", "")
    if pypath_old.strip():
        env["PYTHONPATH"] = pypath_old + os.pathsep + script_folder
    else:
        env["PYTHONPATH"] = script_folder

    # add environment bin folder to PATH
    if env_folder is not None:
        path_old = env.get("PATH", "")
        env["PATH"] = env_folder + os.pathsep + path_old if path_old else env_folder

    # optional conda env name metadata
    if env_name is not None:
        env["CONDA_DEFAULT_ENV"] = env_name

    return env

# helper to get python executable from environment
def get_python_executable(script_folder, env_folder=None, use_env_python=True):
    if use_env_python and env_folder is not None:
        python_env = os.path.join(env_folder, "python")
        if os.path.exists(python_env):
            return python_env

    return sys.executable
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# execution utils

# helper to run command
def run_command(script_path, settings_path, time_step, python_exe=None):
    if python_exe is None:
        python_exe = sys.executable

    time_string = time_step.strftime("%Y-%m-%d %H:%M")

    cmd = [
        python_exe,
        script_path,
        "-settings_file", settings_path,
        "-time", time_string
    ]

    print(f' ---> Execute time step: {time_string}')
    print(' ----> Command: ' + ' '.join(cmd))

    result = subprocess.run(cmd)

    if result.returncode == 0:
        print(f' ---> Execute time step: {time_string} ... DONE')
    else:
        print(f' ---> Execute time step: {time_string} ... FAILED')

    return result.returncode

# helper to get script arguments
def get_args():

    # define wrapper folder (based on execution file location)
    wrapper_folder = os.path.dirname(os.path.abspath(__file__))
    # define default script and settings file paths based on wrapper folder
    wrapper_script_file = os.path.join(wrapper_folder, default_script_file)
    wrapper_settings_file = os.path.join(wrapper_folder, default_settings_file)

    # parse arguments
    parser = argparse.ArgumentParser(
        description=f"Wrapper to execute {wrapper_script_file} over a time period"
    )

    parser.add_argument(
        "-script_file",
        dest="script_file",
        default=wrapper_script_file,
        help=f"Path to {wrapper_script_file} (default: same folder as wrapper)",
    )

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        default=wrapper_settings_file,
        help="Path to JSON settings file (default: same folder as wrapper)",
    )

    parser.add_argument(
        "-time_now",
        dest="time_now",
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM" '
             '(default: system time rounded to hour)',
    )

    parser.add_argument(
        "-time_period",
        dest="time_period",
        type=int,
        default=2,
        help="Number of days backward from time_now (default: 2)",
    )

    parser.add_argument(
        "-time_start",
        dest="time_start",
        default=None,
        help='Start time in format "YYYY-MM-DD HH:MM"',
    )

    parser.add_argument(
        "-time_end",
        dest="time_end",
        default=None,
        help='End time in format "YYYY-MM-DD HH:MM"',
    )

    parser.add_argument(
        "-step_hours",
        dest="step_hours",
        type=int,
        default=24,
        help="Step in hours between runs (default: 24)",
    )

    parser.add_argument(
        "-reverse",
        dest="reverse",
        action="store_true",
        help="Run time steps from latest to oldest",
    )

    parser.add_argument(
        "-round_option",
        dest="round_option",
        default="day",
        choices=["day", "hour", "hour_floor", "hour_ceil"],
        help='Time rounding option: "day"=midnight, '
             '"hour" or "hour_floor"=round down to hour, '
             '"hour_ceil"=round up to next hour'
    )

    parser.add_argument(
        "-time_utc",
        dest="time_utc",
        action="store_true",
        help="Use UTC time for default time_now",
    )

    parser.add_argument(
        "-stop_on_error",
        dest="stop_on_error",
        action="store_true",
        help="Stop execution on first failing time step",
    )

    parser.add_argument(
        "-virtual_env_folder",
        dest="virtual_env_folder",
        default=virtual_env_folder,
        help="Path to conda/venv bin folder",
    )

    parser.add_argument(
        "-virtual_env_name",
        dest="virtual_env_name",
        default=virtual_env_name,
        help="Name of conda/venv environment",
    )

    parser.add_argument(
        "-use_env_python",
        dest="use_env_python",
        action="store_true",
        help="Use python executable found in virtual_env_folder",
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()

    # get script file (or exit if not found)
    script_file = os.path.abspath(args.script_file)
    if not os.path.exists(script_file):
        raise FileNotFoundError(f'Script file "{script_file}" not found')

    # get settings file (or exit if not found)
    settings_file = os.path.abspath(args.settings_file)
    if not os.path.exists(settings_file):
        raise FileNotFoundError(f'Settings file "{settings_file}" not found')

    # get script folder
    script_folder = os.path.dirname(script_file)

    # get time parameters and validate
    if args.time_period < 0:
        raise RuntimeError("time_period must be >= 0")
    if args.step_hours <= 0:
        raise RuntimeError("step_hours must be > 0")

    # get virtual environment folder (or None if not provided)
    if args.virtual_env_folder is not None:
        env_folder = os.path.abspath(args.virtual_env_folder)
    else:
        env_folder = None
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # compute time range
    time_start, time_end = get_time_range(
        time_now=args.time_now,
        time_period=args.time_period,
        time_start=args.time_start,
        time_end=args.time_end,
        round_option=args.round_option
    )

    # build steps
    time_steps = iter_time_steps(
        time_start=time_start,
        time_end=time_end,
        step_hours=args.step_hours,
        reverse=args.reverse
    )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # update environment for subprocesses
    run_env = update_env(
        base_env=os.environ,
        script_folder=script_folder,
        env_folder=env_folder,
        env_name=args.virtual_env_name
    )

    # get python executable to use
    python_exe = get_python_executable(
        script_folder=script_folder,
        env_folder=env_folder,
        use_env_python=args.use_env_python
    )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start info
    print(' ============================================================================ ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> START ... ')
    print(' ')
    print(f' ---> Wrapper folder:     {os.path.dirname(os.path.abspath(__file__))}')
    print(f' ---> Application folder: {script_folder}')
    print(f' ---> Script file:        {script_file}')
    print(f' ---> Settings file:      {settings_file}')
    print(f' ---> Time start:         {time_start.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f' ---> Time end:           {time_end.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f' ---> Step hours:         {args.step_hours}')
    print(f' ---> Reverse:            {args.reverse}')
    print(f' ---> Stop on error:      {args.stop_on_error}')
    print(f' ---> Time UTC:           {args.time_utc}')
    print(f' ---> Number steps:       {len(time_steps)}')
    print(' ')
    print(f' ---> Virtual env folder: {env_folder}')
    print(f' ---> Virtual env name:   {args.virtual_env_name}')
    print(f' ---> Use env python:     {args.use_env_python}')
    print(f' ---> Python executable:  {python_exe}')
    print(f' ---> PYTHONPATH:         {run_env.get("PYTHONPATH", "")}')
    print(' ')

    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # iterate over time steps and execute script
    n_done, n_failed = 0, 0
    for time_step in time_steps:

        # execute command for current time step
        return_code = run_command(
            script_path=script_file,
            settings_path=settings_file,
            time_step=time_step,
            python_exe=sys.executable
        )

        # update counters and check stop_on_error flag
        if return_code == 0:
            n_done += 1
        else:
            n_failed += 1
            if args.stop_on_error:
                print(' ===> Stop execution due to stop_on_error flag')
                break
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message
    elapsed = round(time.time() - start_time, 1)

    print(' ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(f' ==> TIME ELAPSED: {elapsed} seconds')
    print(f' ==> RUNS DONE: {n_done}')
    print(f' ==> RUNS FAILED: {n_failed}')
    print(' ==> ... END')
    print(' ==> Bye, Bye')
    print(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# call script
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
