#!/bin/bash -e

# -----------------------------------------------------------------------------------------
# Script information
script_name='SM MODEL GRID - ORGANIZER TS - SWI T01 - REALTIME (24h)'
script_version="1.2.0"
script_date='2026/07/20'

# Environment file
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python script and settings
fp_script="/hydro/library/fp_package_connectors/sm_model_tools/organizer_ts/sm_model_points_extractor.py"
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/organizer_ts/sm_model_points_extractor_swi_t01_realtime_24h.json"

# Lock folder/file
folder_lock="/hydro/lock/sm_cnr"
fp_lock="${folder_lock}/run_extract_points_swi_t01_realtime_24h.lock"

# Default daily execution hour
time_hour_default="00:00"

# -----------------------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------------------
print_usage() {

    echo
    echo "Usage:"
    echo
    echo "  Realtime mode:"
    echo
    echo "    $0"
    echo "    $0 \"2026-06-10 12:00\""
    echo "    $0 --time \"2026-06-10 12:00\""
    echo "    $0 --force"
    echo "    $0 --time \"2026-06-10\" --hour 23"
    echo "    $0 --time \"2026-07-20\" --previous_day --hour 23"
    echo "    $0 --time \"2026-06-10 12:00\" --force"
    echo
    echo "  History mode:"
    echo
    echo "    $0 --time-start \"2018-01-01\" --time-end \"2018-01-03\""
    echo
    echo "    $0 \\"
    echo "        --time-start \"2018-01-01\" \\"
    echo "        --time-end \"2018-01-03\" \\"
    echo "        --hour \"00:00\""
    echo
    echo "Options:"
    echo
    echo "  --time TIME"
    echo "      Single execution time."
    echo
    echo "  --time-start TIME"
    echo "      First date of the history period."
    echo
    echo "  --time-end TIME"
    echo "      Last date of the history period."
    echo
    echo "  --hour HH or HH:MM"
    echo "      Hour applied to the selected date."
    echo "      In realtime mode it overrides the hour contained in --time."
    echo "      In history mode it is applied to every date."
    echo "      Default for history mode: ${time_hour_default}"
    echo
    echo "  --previous_day"
    echo "      Realtime mode only. Subtract one day from the selected date."
    echo "      Example: --time \"2026-07-20\" --previous_day --hour 23"
    echo "               produces 2026-07-19 23:00."
    echo
    echo "  -f, --force"
    echo "      Force execution by removing the existing lock."
    echo
    echo "  -h, --help"
    echo "      Show this help."
    echo
}

# -----------------------------------------------------------------------------------------
# Run one time step
# -----------------------------------------------------------------------------------------
run_time_step() {

    local time_step="$1"

    cmd=(
        python -u "${fp_script}"
        -settings_file "${fp_settings}"
        -time "${time_step}"
    )

    echo
    echo " -------------------------------------------------------------------------------"
    echo " ====> REFERENCE TIME: ${time_step}"
    echo " ====> COMMAND       : ${cmd[*]}"
    echo " -------------------------------------------------------------------------------"

    "${cmd[@]}"

    step_exit_code=$?

    if [ ${step_exit_code} -ne 0 ]; then
        echo " ====> TIME STEP FAILED: ${time_step}"
        echo " ====> EXIT CODE      : ${step_exit_code}"
        return ${step_exit_code}
    fi

    echo " ====> TIME STEP COMPLETED: ${time_step}"

    return 0
}

# -----------------------------------------------------------------------------------------
# Parse arguments
#
# Realtime:
#   ./script.sh
#   ./script.sh "2026-06-10 12:00"
#   ./script.sh --time "2026-06-10 12:00"
#
# History:
#   ./script.sh --time-start "2018-01-01" --time-end "2018-01-03"
#   ./script.sh --time-start "2018-01-01" --time-end "2018-01-03" --hour "00:00"
# -----------------------------------------------------------------------------------------
time_now=""
time_start=""
time_end=""
time_hour=""
previous_day=false
force_run=false

while [[ $# -gt 0 ]]; do

    case "$1" in

        --time)
            if [ $# -lt 2 ]; then
                echo " ====> ERROR: missing value for --time"
                exit 1
            fi

            time_now="$2"
            shift 2
            ;;

        --time-start)
            if [ $# -lt 2 ]; then
                echo " ====> ERROR: missing value for --time-start"
                exit 1
            fi

            time_start="$2"
            shift 2
            ;;

        --time-end)
            if [ $# -lt 2 ]; then
                echo " ====> ERROR: missing value for --time-end"
                exit 1
            fi

            time_end="$2"
            shift 2
            ;;

        --hour)
            if [ $# -lt 2 ]; then
                echo " ====> ERROR: missing value for --hour"
                exit 1
            fi

            time_hour="$2"
            shift 2
            ;;

        --previous_day)
            previous_day=true
            shift
            ;;

        -f|--force)
            force_run=true
            shift
            ;;

        -h|--help)
            print_usage
            exit 0
            ;;

        -*)
            echo " ====> ERROR: unknown option '$1'"
            print_usage
            exit 1
            ;;

        *)
            if [ -n "${time_now}" ]; then
                echo " ====> ERROR: more than one reference time was provided"
                exit 1
            fi

            time_now="$1"
            shift
            ;;
    esac
done

# -----------------------------------------------------------------------------------------
# Select execution mode
# -----------------------------------------------------------------------------------------
if [ -n "${time_start}" ] || [ -n "${time_end}" ]; then

    execution_mode="history"

    if [ -z "${time_start}" ] || [ -z "${time_end}" ]; then
        echo " ====> ERROR: both --time-start and --time-end are required"
        exit 1
    fi

    if [ -n "${time_now}" ]; then
        echo " ====> ERROR: --time cannot be used with --time-start and --time-end"
        exit 1
    fi

else

    execution_mode="realtime"

    if [ -z "${time_now}" ]; then
        time_now=$(date "+%Y-%m-%d %H:%M")
    fi
fi

if [ "${execution_mode}" = "history" ] && ${previous_day}; then
    echo " ====> ERROR: --previous_day can only be used in realtime mode"
    exit 1
fi

# -----------------------------------------------------------------------------------------
# Validate times
# -----------------------------------------------------------------------------------------
if [ "${execution_mode}" = "realtime" ]; then

    if ! date -d "${time_now}" >/dev/null 2>&1; then
        echo " ====> ERROR: invalid reference time '${time_now}'"
        exit 1
    fi

    selected_date=$(date -d "${time_now}" "+%Y-%m-%d")

    if ${previous_day}; then
        selected_date=$(date -d "${selected_date} -1 day" "+%Y-%m-%d")
    fi

    if [ -n "${time_hour}" ]; then
        if [[ "${time_hour}" =~ ^([01][0-9]|2[0-3])$ ]]; then
            time_hour="${time_hour}:00"
        elif [[ ! "${time_hour}" =~ ^([01][0-9]|2[0-3]):[0-5][0-9]$ ]]; then
            echo " ====> ERROR: invalid hour '${time_hour}'"
            echo " ====> EXPECTED FORMAT: HH or HH:MM"
            exit 1
        fi
    else
        time_hour=$(date -d "${time_now}" "+%H:%M")
    fi

    time_now="${selected_date} ${time_hour}"

else

    if [ -z "${time_hour}" ]; then
        time_hour="${time_hour_default}"
    fi

    if ! date -d "${time_start}" >/dev/null 2>&1; then
        echo " ====> ERROR: invalid start time '${time_start}'"
        exit 1
    fi

    if ! date -d "${time_end}" >/dev/null 2>&1; then
        echo " ====> ERROR: invalid end time '${time_end}'"
        exit 1
    fi

    if [[ "${time_hour}" =~ ^([01][0-9]|2[0-3])$ ]]; then
        time_hour="${time_hour}:00"
    elif [[ ! "${time_hour}" =~ ^([01][0-9]|2[0-3]):[0-5][0-9]$ ]]; then
        echo " ====> ERROR: invalid hour '${time_hour}'"
        echo " ====> EXPECTED FORMAT: HH or HH:MM"
        exit 1
    fi

    date_start=$(date -d "${time_start}" "+%Y-%m-%d")
    date_end=$(date -d "${time_end}" "+%Y-%m-%d")

    epoch_start=$(date -d "${date_start} 00:00" "+%s")
    epoch_end=$(date -d "${date_end} 00:00" "+%s")

    if [ "${epoch_start}" -gt "${epoch_end}" ]; then
        echo " ====> ERROR: start date is after end date"
        echo " ====> START DATE: ${date_start}"
        echo " ====> END DATE  : ${date_end}"
        exit 1
    fi

    # Count dates inclusively.
    number_of_days=0
    count_date="${date_start}"

    while [ "$(date -d "${count_date}" "+%s")" -le "${epoch_end}" ]; do
        number_of_days=$((number_of_days + 1))
        count_date=$(date -d "${count_date} +1 day" "+%Y-%m-%d")
    done
fi

# -----------------------------------------------------------------------------------------
# Start information
# -----------------------------------------------------------------------------------------
echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo " ====> EXECUTION MODE: ${execution_mode}"
echo " ====> FORCE RUN     : ${force_run}"

if [ "${execution_mode}" = "realtime" ]; then

    echo " ====> PREVIOUS DAY  : ${previous_day}"
    echo " ====> SELECTED HOUR : ${time_hour}"
    echo " ====> REFERENCE TIME: ${time_now}"

else

    echo " ====> TIME START    : ${date_start} ${time_hour}"
    echo " ====> TIME END      : ${date_end} ${time_hour}"
    echo " ====> DAILY HOUR    : ${time_hour}"
    echo " ====> NUMBER OF DAYS: ${number_of_days}"
fi

# -----------------------------------------------------------------------------------------
# Create lock folder
# -----------------------------------------------------------------------------------------
echo " ====> CREATE LOCK FOLDER ..."

mkdir -p "${folder_lock}"

echo " ====> CREATE LOCK FOLDER ... DONE: ${folder_lock}"

# -----------------------------------------------------------------------------------------
# Load Python environment
# -----------------------------------------------------------------------------------------
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
# Force mode: remove stale lock and kill lock owner, if available
# -----------------------------------------------------------------------------------------
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
# Run with lock
# -----------------------------------------------------------------------------------------
echo " ====> RUN POINT EXTRACTION ..."
echo " ====> PYTHON SCRIPT : ${fp_script}"
echo " ====> SETTINGS FILE: ${fp_settings}"
echo " ====> LOCK FILE     : ${fp_lock}"

(
    flock -n 9 || {
        echo " ====> RUN POINT EXTRACTION ... SKIPPED"
        echo " ====> ANOTHER INSTANCE IS ALREADY RUNNING"
        echo " ====> USE -f OR --force TO OVERRIDE"
        exit 1
    }

    if [ "${execution_mode}" = "realtime" ]; then

        run_time_step "${time_now}"

    else

        current_date="${date_start}"
        step_id=0

        while [ "$(date -d "${current_date}" "+%s")" -le "${epoch_end}" ]; do

            step_id=$((step_id + 1))
            current_time="${current_date} ${time_hour}"

            echo
            echo " ==================================================================================="
            echo " ====> HISTORY STEP ${step_id}/${number_of_days}"
            echo " ====> CURRENT TIME: ${current_time}"
            echo " ==================================================================================="

            run_time_step "${current_time}"

            current_date=$(date -d "${current_date} +1 day" "+%Y-%m-%d")
        done
    fi

) 9>"${fp_lock}"

exit_code=$?

if [ ${exit_code} -ne 0 ]; then
    echo " ====> RUN POINT EXTRACTION ... FAILED"
    echo " ====> EXIT CODE: ${exit_code}"
    exit ${exit_code}
fi

echo " ====> RUN POINT EXTRACTION ... DONE"
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="

