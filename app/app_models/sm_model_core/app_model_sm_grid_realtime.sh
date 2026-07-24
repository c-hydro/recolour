#!/bin/bash -e

# -----------------------------------------------------------------------------------------
# Script information
script_name='SM MODEL RUNNER - REALTIME'
script_version="1.2.0"
script_date='2026/07/21'

export TZ="Europe/Rome"

# Environment file
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python script and settings
fp_script="/hydro/library/fp_package_connectors/sm_model_core/app_model_sm_main.py"
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/runner/app_model_sm_grid_realtime.json"

# Lock folder/file
folder_lock="/hydro/lock/sm_cnr"


# -----------------------------------------------------------------------------------------
# Arguments
#
# Realtime behavior:
#   no --time-start and no --time-end
#
# History behavior:
#   both --time-start and --time-end are provided

time_arg=""
time_start_arg=""
time_end_arg=""

run_hour=""
force_run=false
previous_day=false

execution_mode=""

# -----------------------------------------------------------------------------------------
# Usage

usage() {
    echo "Usage:"
    echo ""
    echo "Realtime behavior:"
    echo "  $0 [\"YYYY-MM-DD HH:MM\"] \\"
    echo "     [--hour \"HH:MM\"] [--previous-day] [-f|--force]"
    echo ""
    echo "History behavior:"
    echo "  $0 --time-start \"YYYY-MM-DD\" \\"
    echo "     --time-end \"YYYY-MM-DD\" \\"
    echo "     [--hour \"HH:MM\"] [-f|--force]"
    echo ""
    echo "Realtime examples:"
    echo "  $0"
    echo "  $0 --hour \"12:00\""
    echo "  $0 --previous-day --hour \"00:00\""
    echo "  $0 \"2026-07-21 15:30\""
    echo "  $0 \"2026-07-21\" --hour \"12:00\""
    echo ""
    echo "History example:"
    echo "  $0 --time-start \"2026-07-01\" \\"
    echo "     --time-end \"2026-07-21\" \\"
    echo "     --hour \"12:00\""
}

# -----------------------------------------------------------------------------------------
# Parse arguments

while [[ $# -gt 0 ]]; do
    case "$1" in
        --time-start|--time_start)
            if [ -z "${2:-}" ]; then
                echo "ERROR: $1 requires a value"
                exit 1
            fi

            time_start_arg="$2"
            shift 2
            ;;

        --time-end|--time_end)
            if [ -z "${2:-}" ]; then
                echo "ERROR: $1 requires a value"
                exit 1
            fi

            time_end_arg="$2"
            shift 2
            ;;

        -H|--hour)
            if [ -z "${2:-}" ]; then
                echo "ERROR: --hour requires a value such as 00:00 or 12:00"
                exit 1
            fi

            run_hour="$2"
            shift 2
            ;;

        --previous-day|--previous_day)
            previous_day=true
            shift
            ;;

        -f|--force)
            force_run=true
            shift
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        -*)
            echo "ERROR: Unknown option '$1'"
            usage
            exit 1
            ;;

        *)
            if [ -n "${time_arg}" ]; then
                echo "ERROR: Multiple positional time arguments provided"
                echo "First : ${time_arg}"
                echo "Second: $1"
                exit 1
            fi

            time_arg="$1"
            shift
            ;;
    esac
done

# -----------------------------------------------------------------------------------------
# Detect execution behavior

if [ -n "${time_start_arg}" ] || [ -n "${time_end_arg}" ]; then
    execution_mode="history"
else
    execution_mode="realtime"
fi

# -----------------------------------------------------------------------------------------
# Validate and normalize hour

if [ -n "${run_hour}" ]; then
    if ! [[ "${run_hour}" =~ ^([01]?[0-9]|2[0-3]):[0-5][0-9]$ ]]; then
        echo "ERROR: Invalid hour '${run_hour}'"
        echo "ERROR: Use HH:MM, for example 00:00, 06:00, 12:00 or 23:00"
        exit 1
    fi

    hour_value="${run_hour%%:*}"
    minute_value="${run_hour##*:}"

    run_hour=$(printf "%02d:%02d" \
        "$((10#${hour_value}))" \
        "$((10#${minute_value}))")
fi

# -----------------------------------------------------------------------------------------
# Configure realtime behavior

if [ "${execution_mode}" = "realtime" ]; then

    if [ -n "${time_arg}" ]; then
        if ! date -d "${time_arg}" >/dev/null 2>&1; then
            echo "ERROR: Invalid reference time '${time_arg}'"
            echo "ERROR: Expected YYYY-MM-DD or YYYY-MM-DD HH:MM"
            exit 1
        fi

        base_time=$(date -d "${time_arg}" "+%Y-%m-%d %H:%M")
    else
        base_time=$(date "+%Y-%m-%d %H:%M")
    fi

    base_date=$(date -d "${base_time}" "+%Y-%m-%d")
    base_hour=$(date -d "${base_time}" "+%H:%M")

    if ${previous_day}; then
        base_date=$(date -d "${base_date} -1 day" "+%Y-%m-%d")
    fi

    if [ -n "${run_hour}" ]; then
        time_now="${base_date} ${run_hour}"
    else
        time_now="${base_date} ${base_hour}"
    fi

# -----------------------------------------------------------------------------------------
# Configure history behavior

else

    if [ -n "${time_arg}" ]; then
        echo "ERROR: Positional time cannot be used with --time-start/--time-end"
        exit 1
    fi

    if ${previous_day}; then
        echo "ERROR: --previous-day cannot be used with history behavior"
        exit 1
    fi

    if [ -z "${time_start_arg}" ]; then
        echo "ERROR: --time-start is mandatory for history behavior"
        usage
        exit 1
    fi

    if [ -z "${time_end_arg}" ]; then
        echo "ERROR: --time-end is mandatory for history behavior"
        usage
        exit 1
    fi

    if ! date -d "${time_start_arg}" >/dev/null 2>&1; then
        echo "ERROR: Invalid time start '${time_start_arg}'"
        exit 1
    fi

    if ! date -d "${time_end_arg}" >/dev/null 2>&1; then
        echo "ERROR: Invalid time end '${time_end_arg}'"
        exit 1
    fi

    if [ -z "${run_hour}" ]; then
        run_hour="00:00"
    fi

    time_start_date=$(date -d "${time_start_arg}" "+%Y-%m-%d")
    time_end_date=$(date -d "${time_end_arg}" "+%Y-%m-%d")

    time_start_epoch=$(date -d "${time_start_date} 00:00" "+%s")
    time_end_epoch=$(date -d "${time_end_date} 00:00" "+%s")

    if [ "${time_start_epoch}" -gt "${time_end_epoch}" ]; then
        echo "ERROR: time-start must be before or equal to time-end"
        echo "ERROR: TIME START: ${time_start_date}"
        echo "ERROR: TIME END  : ${time_end_date}"
        exit 1
    fi

    number_days=$(((time_end_epoch - time_start_epoch) / 86400 + 1))
fi

# -----------------------------------------------------------------------------------------
# Information

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo " ====> SYSTEM TIME   : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo " ====> EXECUTION MODE: ${execution_mode}"
echo " ====> FORCE RUN     : ${force_run}"

if [ "${execution_mode}" = "realtime" ]; then
    echo " ====> REFERENCE TIME: ${time_now}"
    echo " ====> HOUR OVERRIDE : ${run_hour:-none}"
    echo " ====> PREVIOUS DAY  : ${previous_day}"
else
    echo " ====> TIME START    : ${time_start_date} ${run_hour}"
    echo " ====> TIME END      : ${time_end_date} ${run_hour}"
    echo " ====> REFERENCE HOUR: ${run_hour}"
    echo " ====> FREQUENCY     : daily"
    echo " ====> NUMBER OF DAYS: ${number_days}"
fi

# -----------------------------------------------------------------------------------------
# Create lock folder

echo " ====> CREATE LOCK FOLDER ..."
mkdir -p "${folder_lock}"
echo " ====> CREATE LOCK FOLDER ... DONE: ${folder_lock}"

# Time tag (history mode)
if [[ -n "$time_start_arg" && -n "$time_end_arg" ]]; then
    time_tag="_$(date -d "$time_start_arg" +"%Y%m%d%H%M")_$(date -d "$time_end_arg" +"%Y%m%d%H%M")"
else
    time_tag=""
fi

# Lock file
fp_lock="${folder_lock}/run_sm_model_grid_realtime${time_tag}.lock"

# -----------------------------------------------------------------------------------------
# Load Python environment

echo " ====> LOAD PYTHON ENVIRONMENT ..."

if [ -f "${fp_env_file}" ]; then
    source "${fp_env_file}"
else
    echo " ====> LOAD PYTHON ENVIRONMENT ... FAILED"
    echo " ====> ENVIRONMENT FILE NOT FOUND: ${fp_env_file}"
    exit 1
fi

echo " ====> LOAD PYTHON ENVIRONMENT ... DONE"

# -----------------------------------------------------------------------------------------
# Force mode

if ${force_run}; then
    echo " ====> FORCE MODE ENABLED"

    if command -v fuser >/dev/null 2>&1; then
        echo " ====> KILL PROCESS USING LOCK FILE, IF ANY: ${fp_lock}"
        fuser -k "${fp_lock}" >/dev/null 2>&1 || true
    fi

    echo " ====> REMOVE LOCK FILE: ${fp_lock}"
    rm -f "${fp_lock}"
fi

# -----------------------------------------------------------------------------------------
# Run model

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

    # -------------------------------------------------------------------------------------
    # Realtime behavior: run once

    if [ "${execution_mode}" = "realtime" ]; then

        cmd=(
            python "${fp_script}"
            -settings_file "${fp_settings}"
            -time "${time_now}"
        )

        echo " ====> COMMAND: ${cmd[*]}"

        "${cmd[@]}"

    # -------------------------------------------------------------------------------------
    # History behavior: iterate daily

    else

        current_date="${time_start_date}"
        run_index=0
        failed_runs=0

        while true; do

            run_index=$((run_index + 1))
            time_now="${current_date} ${run_hour}"

            echo " -----------------------------------------------------------------------------------"
            echo " ====> RUN ${run_index}/${number_days}"
            echo " ====> REFERENCE TIME: ${time_now}"

            cmd=(
                python "${fp_script}"
                -settings_file "${fp_settings}"
                -time "${time_now}"
            )

            echo " ====> COMMAND: ${cmd[*]}"

            if "${cmd[@]}"; then
                echo " ====> RUN ${run_index}/${number_days} ... DONE"
            else
                run_exit_code=$?

                echo " ====> RUN ${run_index}/${number_days} ... FAILED"
                echo " ====> REFERENCE TIME: ${time_now}"
                echo " ====> EXIT CODE    : ${run_exit_code}"

                failed_runs=$((failed_runs + 1))
            fi

            # Include time_end in the processing
            if [ "${current_date}" = "${time_end_date}" ]; then
                break
            fi

            current_date=$(date -d "${current_date} +1 day" "+%Y-%m-%d")
        done

        echo " -----------------------------------------------------------------------------------"
        echo " ====> TOTAL RUNS     : ${number_days}"
        echo " ====> SUCCESSFUL RUNS: $((number_days - failed_runs))"
        echo " ====> FAILED RUNS    : ${failed_runs}"

        if [ "${failed_runs}" -gt 0 ]; then
            exit 1
        fi
    fi

) 9>"${fp_lock}"

exit_code=$?

# -----------------------------------------------------------------------------------------
# Check result

if [ "${exit_code}" -ne 0 ]; then
    echo " ====> RUN MODEL ... FAILED"
    echo " ====> EXIT CODE: ${exit_code}"
    exit "${exit_code}"
fi

echo " ====> RUN MODEL ... DONE"
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
