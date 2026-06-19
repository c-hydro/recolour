#!/bin/bash -e

script_name='SM MODEL POINTS 2 GRID - REALTIME'
script_version="1.0.0"
script_date='2026/06/18'

export TZ="Europe/Rome"

fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

fp_script="/hydro/library/fp_package_connectors/sm_model_tools/organizer_grid/sm_model_points2grid.py"
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/organizer_grid/sm_model_points2grid_realtime.json"

folder_lock="/hydro/lock/sm_cnr"
fp_lock="${folder_lock}/run_points2grid_realtime.lock"

time_arg=""
force_run=false
hour_override=""
previous_day=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)
            force_run=true
            shift
            ;;
        --hour|-hour)
            if [ -z "${2:-}" ]; then
                echo "ERROR: --hour requires a value between 0 and 23"
                exit 1
            fi
            hour_override="$2"
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

if [ -n "${hour_override}" ]; then
    if ! [[ "${hour_override}" =~ ^[0-9]+$ ]] || \
       [ "${hour_override}" -lt 0 ] || \
       [ "${hour_override}" -gt 23 ]; then
        echo "ERROR: Invalid hour '${hour_override}'. Expected 0-23"
        exit 1
    fi
fi

if [ -z "${time_arg}" ]; then
    base_date=$(date "+%Y-%m-%d")
    base_minute=$(date "+%M")
else
    base_date=$(date -d "${time_arg}" "+%Y-%m-%d")
    base_minute=$(date -d "${time_arg}" "+%M")
fi

if ${previous_day}; then
    base_date=$(date -d "${base_date} -1 day" "+%Y-%m-%d")
fi

if [ -n "${hour_override}" ]; then
    time_now="${base_date} $(printf "%02d" "${hour_override}"):00"
else
    if [ -z "${time_arg}" ]; then
        time_now=$(date "+%Y-%m-%d %H:%M")
    else
        time_now=$(date -d "${time_arg}" "+%Y-%m-%d %H:%M")
    fi
fi

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo " ====> SYSTEM TIME   : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo " ====> REFERENCE TIME: ${time_now}"
echo " ====> HOUR OVERRIDE : ${hour_override:-none}"
echo " ====> PREVIOUS DAY  : ${previous_day}"
echo " ====> FORCE RUN     : ${force_run}"

mkdir -p "${folder_lock}"

echo " ====> LOAD PYTHON ENVIRONMENT ..."
if [ -f "${fp_env_file}" ]; then
    source "${fp_env_file}"
else
    echo " ====> ENVIRONMENT FILE NOT FOUND: ${fp_env_file}"
    exit 1
fi
echo " ====> LOAD PYTHON ENVIRONMENT ... DONE"

if ${force_run}; then
    echo " ====> FORCE MODE ENABLED"

    if command -v fuser >/dev/null 2>&1; then
        fuser -k "${fp_lock}" >/dev/null 2>&1 || true
    fi

    rm -f "${fp_lock}"
fi

echo " ====> RUN POINTS2GRID ..."
echo " ====> PYTHON SCRIPT : ${fp_script}"
echo " ====> SETTINGS FILE : ${fp_settings}"
echo " ====> LOCK FILE     : ${fp_lock}"

(
    flock -n 9 || {
        echo " ====> RUN POINTS2GRID ... SKIPPED"
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
    echo " ====> RUN POINTS2GRID ... FAILED"
    exit ${exit_code}
fi

echo " ====> RUN POINTS2GRID ... DONE"
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
