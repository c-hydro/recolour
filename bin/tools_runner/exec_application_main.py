#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - EXECUTION WRAPPER
"""

import os
import sys
import json
import time
import argparse
import warnings
import subprocess
from datetime import datetime, timedelta

warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated as an API"
)

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Execution wrapper for applications and tools'
alg_type = 'Exec'
alg_version = '1.3.0'
alg_release = '2026-04-20'
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# default script and settings file
default_script_file = "../../tools/algorithm_transfer/h122/app_transfer_ssm_h122.py"
default_settings_file = "../../tools/algorithm_transfer/h122/app_transfer_ssm_h122.json"

# default environment settings
virtual_env_folder = '/home/fabio/Documents/Work_Area/Code_Development/Workspace/recolour/conda_p311/bin/'
virtual_env_name = 'recolour_converter_libraries'

# default summary settings
default_input_summary_folder = None
default_input_summary_name = None

default_output_summary_folder = "./summary"
default_output_summary_name = None

default_following_algorithm = "next_algorithm"
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# generic utils
def is_none(value):
    if value is None:
        return True
    if str(value).lower() in ["none", "null", ""]:
        return True
    return False
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# time utils
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


def parse_time(time_string, round_option="day"):

    try:
        time_obj = datetime.strptime(time_string, "%Y-%m-%d %H:%M")
    except ValueError as exc:
        raise ValueError('Time must have format "YYYY-MM-DD HH:MM"') from exc

    return normalize_time(time_obj, round_option=round_option)


def get_time_now(time_string=None, round_option="day", time_utc=False):

    if time_string is None:

        if time_utc:
            time_obj = datetime.utcnow()
        else:
            time_obj = datetime.now()

        return normalize_time(time_obj, round_option=round_option)

    return parse_time(time_string, round_option=round_option)


def get_time_range(time_now=None, time_period=2, time_start=None, time_end=None,
                   round_option="day", time_utc=False):

    if (time_start is not None) and (time_end is not None):

        time_from = parse_time(time_start, round_option=round_option)
        time_to = parse_time(time_end, round_option=round_option)

    elif (time_start is None) and (time_end is None):

        time_to = get_time_now(
            time_string=time_now,
            round_option=round_option,
            time_utc=time_utc
        )

        time_from = time_to - timedelta(days=int(time_period))

    else:
        raise RuntimeError(
            'time_start and time_end must be both provided, or neither provided'
        )

    if time_from > time_to:
        raise RuntimeError('time_start must be less than or equal to time_end')

    return time_from, time_to


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
# summary utils
def build_input_summary_path(summary_folder=None, summary_name=None):

    if is_none(summary_folder) or is_none(summary_name):
        return None

    return os.path.join(summary_folder, summary_name)


def build_output_summary_path(summary_folder=None, summary_name=None,
                              time_end=None, following_algorithm=None):

    if is_none(summary_folder):
        return None

    if not os.path.exists(summary_folder):
        os.makedirs(summary_folder, exist_ok=True)

    if not is_none(summary_name):
        return os.path.join(summary_folder, summary_name)

    if time_end is None:
        time_tag = datetime.now().strftime("%Y%m%d%H%M")
    else:
        time_tag = time_end.strftime("%Y%m%d%H%M")

    if is_none(following_algorithm):
        following_algorithm = default_following_algorithm

    filename = f"{following_algorithm}_{time_tag}.json"

    return os.path.join(summary_folder, filename)


def check_summary_file(summary_file):

    if is_none(summary_file):
        return True

    if not os.path.exists(summary_file):
        print(f' ===> Input summary file not found: {summary_file}')
        return False

    try:
        with open(summary_file, "r") as file_handle:
            summary_data = json.load(file_handle)
    except Exception as exc:
        print(f' ===> Input summary file is not readable: {summary_file}')
        print(f' ===> Error: {exc}')
        return False

    summary_status = summary_data.get("status", None)

    if summary_status == "DONE":
        return True

    print(f' ===> Input summary status is not DONE: {summary_status}')
    return False


def save_summary_file(summary_file, time_now, time_now_utc, round_option,
                      time_start, time_end, time_steps,
                      n_done, n_failed, status="DONE",
                      algorithm="ssm_h122_transfer",
                      following_algorithm=None):

    if is_none(summary_file):
        return

    summary_folder = os.path.dirname(os.path.abspath(summary_file))
    if summary_folder:
        os.makedirs(summary_folder, exist_ok=True)

    summary_data = {
        "project": project_name,
        "algorithm": algorithm,
        "following_algorithm": following_algorithm,
        "status": status,

        "time_now": time_now.strftime("%Y-%m-%d %H:%M") if time_now is not None else None,
        "time_now_utc": time_now_utc.strftime("%Y-%m-%d %H:%M") if time_now_utc is not None else None,
        "round_option": round_option,

        "time_start": time_start.strftime("%Y-%m-%d %H:%M"),
        "time_end": time_end.strftime("%Y-%m-%d %H:%M"),
        "time_steps": [
            time_step.strftime("%Y-%m-%d %H:%M") for time_step in time_steps
        ],

        "runs_expected": len(time_steps),
        "runs_done": n_done,
        "runs_failed": n_failed,
        "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }

    with open(summary_file, "w") as file_handle:
        json.dump(summary_data, file_handle, indent=4)

    print(f' ---> Output summary saved: {summary_file}')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# environment utils
def update_env(base_env, script_folder, env_folder=None, env_name=None):

    env = base_env.copy()

    pypath_old = env.get("PYTHONPATH", "")
    if pypath_old.strip():
        env["PYTHONPATH"] = pypath_old + os.pathsep + script_folder
    else:
        env["PYTHONPATH"] = script_folder

    if not is_none(env_folder):
        path_old = env.get("PATH", "")
        env["PATH"] = env_folder + os.pathsep + path_old if path_old else env_folder

    if not is_none(env_name):
        env["CONDA_DEFAULT_ENV"] = env_name

    return env


def get_python_executable(env_folder=None, use_env_python=True):

    if use_env_python and not is_none(env_folder):

        python_env = os.path.join(env_folder, "python")

        if os.path.exists(python_env):
            return python_env

    return sys.executable
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# execution utils
def run_command(script_path, settings_path, time_step, python_exe=None, env=None):

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

    result = subprocess.run(cmd, env=env)

    if result.returncode == 0:
        print(f' ---> Execute time step: {time_string} ... DONE')
    else:
        print(f' ---> Execute time step: {time_string} ... FAILED')

    return result.returncode
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# argument parser
def get_args():

    wrapper_folder = os.path.dirname(os.path.abspath(__file__))
    wrapper_script_file = os.path.join(wrapper_folder, default_script_file)
    wrapper_settings_file = os.path.join(wrapper_folder, default_settings_file)

    parser = argparse.ArgumentParser(
        description=f"Wrapper to execute {wrapper_script_file} over a time period"
    )

    parser.add_argument("-script_file", dest="script_file", default=wrapper_script_file)
    parser.add_argument("-settings_file", dest="settings_file", default=wrapper_settings_file)

    parser.add_argument("-time_now", dest="time_now", default=None)
    parser.add_argument("-time_period", dest="time_period", type=int, default=2)
    parser.add_argument("-time_start", dest="time_start", default=None)
    parser.add_argument("-time_end", dest="time_end", default=None)
    parser.add_argument("-step_hours", dest="step_hours", type=int, default=24)
    parser.add_argument("-reverse", dest="reverse", action="store_true")

    parser.add_argument(
        "-round_option",
        dest="round_option",
        default="day",
        choices=["day", "hour", "hour_floor", "hour_ceil"]
    )

    parser.add_argument("-time_utc", dest="time_utc", action="store_true")
    parser.add_argument("-stop_on_error", dest="stop_on_error", action="store_true")

    parser.add_argument("-virtual_env_folder", dest="virtual_env_folder", default=virtual_env_folder)
    parser.add_argument("-virtual_env_name", dest="virtual_env_name", default=virtual_env_name)
    parser.add_argument("-use_env_python", dest="use_env_python", action="store_true")

    parser.add_argument(
        "-input_summary_folder",
        dest="input_summary_folder",
        default=default_input_summary_folder
    )

    parser.add_argument(
        "-input_summary_name",
        dest="input_summary_name",
        default=default_input_summary_name
    )

    parser.add_argument(
        "-output_summary_folder",
        dest="output_summary_folder",
        default=default_output_summary_folder
    )

    parser.add_argument(
        "-output_summary_name",
        dest="output_summary_name",
        default=default_output_summary_name
    )

    parser.add_argument(
        "-following_algorithm",
        dest="following_algorithm",
        default=default_following_algorithm
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    args = get_args()

    script_file = os.path.abspath(args.script_file)
    if not os.path.exists(script_file):
        raise FileNotFoundError(f'Script file "{script_file}" not found')

    settings_file = os.path.abspath(args.settings_file)
    if not os.path.exists(settings_file):
        raise FileNotFoundError(f'Settings file "{settings_file}" not found')

    script_folder = os.path.dirname(script_file)

    if args.time_period < 0:
        raise RuntimeError("time_period must be >= 0")

    if args.step_hours <= 0:
        raise RuntimeError("step_hours must be > 0")

    if not is_none(args.virtual_env_folder):
        env_folder = os.path.abspath(args.virtual_env_folder)
    else:
        env_folder = None

    time_now_obj = get_time_now(
        time_string=args.time_now,
        round_option=args.round_option,
        time_utc=args.time_utc
    )

    time_now_utc_obj = get_time_now(
        time_string=args.time_now,
        round_option=args.round_option,
        time_utc=True
    )

    time_start, time_end = get_time_range(
        time_now=time_now_obj.strftime("%Y-%m-%d %H:%M"),
        time_period=args.time_period,
        time_start=args.time_start,
        time_end=args.time_end,
        round_option=args.round_option,
        time_utc=args.time_utc
    )

    time_steps = iter_time_steps(
        time_start=time_start,
        time_end=time_end,
        step_hours=args.step_hours,
        reverse=args.reverse
    )

    input_summary_file = build_input_summary_path(
        summary_folder=args.input_summary_folder,
        summary_name=args.input_summary_name
    )

    output_summary_file = build_output_summary_path(
        summary_folder=args.output_summary_folder,
        summary_name=args.output_summary_name,
        time_end=time_end,
        following_algorithm=args.following_algorithm
    )

    run_env = update_env(
        base_env=os.environ,
        script_folder=script_folder,
        env_folder=env_folder,
        env_name=args.virtual_env_name
    )

    python_exe = get_python_executable(
        env_folder=env_folder,
        use_env_python=args.use_env_python
    )

    print(' ============================================================================ ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> START ... ')
    print(' ')
    print(f' ---> Wrapper folder:          {os.path.dirname(os.path.abspath(__file__))}')
    print(f' ---> Application folder:      {script_folder}')
    print(f' ---> Script file:             {script_file}')
    print(f' ---> Settings file:           {settings_file}')
    print(f' ---> Time now:                {time_now_obj.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f' ---> Time now UTC:            {time_now_utc_obj.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f' ---> Round option:            {args.round_option}')
    print(f' ---> Time start:              {time_start.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f' ---> Time end:                {time_end.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f' ---> Step hours:              {args.step_hours}')
    print(f' ---> Reverse:                 {args.reverse}')
    print(f' ---> Stop on error:           {args.stop_on_error}')
    print(f' ---> Time UTC:                {args.time_utc}')
    print(f' ---> Number steps:            {len(time_steps)}')
    print(' ')
    print(f' ---> Input summary file:      {input_summary_file}')
    print(f' ---> Output summary file:     {output_summary_file}')
    print(f' ---> Following algorithm:     {args.following_algorithm}')
    print(' ')
    print(f' ---> Virtual env folder:      {env_folder}')
    print(f' ---> Virtual env name:        {args.virtual_env_name}')
    print(f' ---> Use env python:          {args.use_env_python}')
    print(f' ---> Python executable:       {python_exe}')
    print(f' ---> PYTHONPATH:              {run_env.get("PYTHONPATH", "")}')
    print(' ')

    start_time = time.time()

    n_done, n_failed = 0, 0
    summary_status = "DONE"

    if not check_summary_file(input_summary_file):

        summary_status = "SKIPPED"

        print(' ===> Upstream summary check failed')
        print(' ===> Time iterations will not start')

    else:

        print(' ---> Upstream summary check ... DONE')
        print(' ')

        for time_step in time_steps:

            return_code = run_command(
                script_path=script_file,
                settings_path=settings_file,
                time_step=time_step,
                python_exe=python_exe,
                env=run_env
            )

            if return_code == 0:
                n_done += 1

            else:
                n_failed += 1
                summary_status = "FAILED"

                if args.stop_on_error:
                    print(' ===> Stop execution due to stop_on_error flag')
                    break

    if n_failed > 0:
        summary_status = "FAILED"

    save_summary_file(
        summary_file=output_summary_file,
        time_now=time_now_obj,
        time_now_utc=time_now_utc_obj,
        round_option=args.round_option,
        time_start=time_start,
        time_end=time_end,
        time_steps=time_steps,
        n_done=n_done,
        n_failed=n_failed,
        status=summary_status,
        algorithm="ssm_h122_transfer",
        following_algorithm=args.following_algorithm
    )

    elapsed = round(time.time() - start_time, 1)

    print(' ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(f' ==> TIME ELAPSED: {elapsed} seconds')
    print(f' ==> SUMMARY STATUS: {summary_status}')
    print(f' ==> RUNS DONE: {n_done}')
    print(f' ==> RUNS FAILED: {n_failed}')
    print(' ==> ... END')
    print(' ==> Bye, Bye')
    print(' ============================================================================ ')


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------