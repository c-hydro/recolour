#!/bin/bash -e

# -----------------------------------------------------------------------------------------
# Script information
script_name='SM MODEL RUNNER - REALTIME'
script_version="1.0.5"
script_date='2026/06/17'

export TZ="Europe/Rome"

# Environment file
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python script and settings
fp_script="/hydro/library/fp_package_connectors/sm_model_core/app_model_sm_main.py"
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/runner/app_model_sm_grid_realtime.json"

# Lock folder/file
folder_lock="/hydro/lock/sm_cnr"
fp_lock="${folder_lock}/run_sm_model_grid_realtime.lock"

# Arguments
time_arg=""
run_hour=""
force_run=false
previous_day=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)
            force_run=true
            shift
            ;;
        -H|--hour)
            if [ -z "${2:-}" ]; then
                echo "ERROR: --hour requires a value between 0 and 23"
                exit 1
            fi
            run_hour="$2"
            shift 2
            ;;
        --previous-day|--previous_day)
            previous_day=true
            shift
            ;;
        -*)
            echo "ERROR: Unknown option '$1'"
            exit 1
            ;;
        *)
            if [ -n "${time_arg}" ]; then
                echo "ERROR: Multiple time arguments provided"
                echo "First : ${time_arg}"
                echo "Second: $1"
                exit 1
            fi
            time_arg="$1"
            shift
            ;;
    esac
done

if [ -n "${run_hour}" ]; then
    if ! [[ "${run_hour}" =~ ^([01]?[0-9]|2[0-3])$ ]]; then
        echo "ERROR: invalid hour '${run_hour}' (allowed values: 00-23)"
        exit 1
    fi
    run_hour=$(printf "%02d" "${run_hour}")
fi

if [ -z "${time_arg}" ]; then
    base_date=$(date "+%Y-%m-%d")
    time_now=$(date "+%Y-%m-%d %H:%M")
else
    base_date=$(date -d "${time_arg}" "+%Y-%m-%d")
    time_now=$(date -d "${time_arg}" "+%Y-%m-%d %H:%M")
fi

if ${previous_day}; then
    base_date=$(date -d "${base_date} -1 day" "+%Y-%m-%d")
fi

if [ -n "${run_hour}" ]; then
    time_now="${base_date} ${run_hour}:00"
else
    if ${previous_day}; then
        time_hm=$(date -d "${time_now}" "+%H:%M")
        time_now="${base_date} ${time_hm}"
    fi
fi

# -----------------------------------------------------------------------------------------

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo " ====> SYSTEM TIME   : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo " ====> REFERENCE TIME: ${time_now}"
echo " ====> HOUR OVERRIDE : ${run_hour:-none}"
echo " ====> PREVIOUS DAY  : ${previous_day}"
echo " ====> FORCE RUN     : ${force_run}"

echo " ====> CREATE LOCK FOLDER ..."
mkdir -p "${folder_lock}"
echo " ====> CREATE LOCK FOLDER ... DONE: ${folder_lock}"

echo " ====> LOAD PYTHON ENVIRONMENT ..."

if [ -f "${fp_env_file}" ]; then
    source "${fp_env_file}"
else
    echo " ====> LOAD PYTHON ENVIRONMENT ... FAILED"
    echo " ====> ENVIRONMENT FILE NOT FOUND: ${fp_env_file}"
    exit 1
fi

echo " ====> LOAD PYTHON ENVIRONMENT ... DONE"

if ${force_run}; then
    echo " ====> FORCE MODE ENABLED"

    if command -v fuser >/dev/null 2>&1; then
        echo " ====> KILL PROCESS USING LOCK FILE, IF ANY: ${fp_lock}"
        fuser -k "${fp_lock}" >/dev/null 2>&1 || true
    fi

    echo " ====> REMOVE LOCK FILE: ${fp_lock}"
    rm -f "${fp_lock}"
fi

echo " ====> RUN MODEL ..."
echo " ====> PYTHON SCRIPT : ${fp_script}"
echo " ====> SETTINGS FILE : ${fp_settings}"
echo " ====> LOCK FILE     : ${fp_lock}"

(
    flock -n 9 || {
        echo " ====> RUN MODEL ... SKIPPED"
        echo " ====> ANOTHER INSTANCE IS ALREADY RUNNING"
        echo " ====> USE -f OR --force TO OVERRIDE"
        exit 1
    }

    cmd=(
        python "${fp_script}"
        -settings_file "${fp_settings}"
        -time "${time_now}"
    )

    echo " ====> COMMAND: ${cmd[*]}"
    "${cmd[@]}"

) 9>"${fp_lock}"

exit_code=$?

if [ ${exit_code} -ne 0 ]; then
    echo " ====> RUN MODEL ... FAILED"
    exit ${exit_code}
fi

echo " ====> RUN MODEL ... DONE"
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
