#!/usr/bin/env bash

# ==============================================================================
# SM MODEL POINTS 2 GRID - REALTIME
#
# Single-time examples:
#
#   ./run_points2grid.sh
#
#   ./run_points2grid.sh "2026-07-13 10:00"
#
#   ./run_points2grid.sh --hour "08:00"
#
#   ./run_points2grid.sh --previous-day --hour "08:00"
#
# Period examples:
#
#   ./run_points2grid.sh \
#       --time-start "2026-07-01 00:00" \
#       --time-end   "2026-07-13 00:00"
#
# Run once per day at 23:00:
#
#   ./run_points2grid.sh \
#       --time-start "2026-07-01" \
#       --time-end   "2026-07-13" \
#       --hour "23:00" \
#       --step-hours 24
#
# The period is processed backward:
#
#   time_end -> time_start
#
# If one execution fails, the script prints a warning and continues with the
# next time step.
#
# ==============================================================================

script_name='SM MODEL POINTS 2 GRID - HISTORY'
script_version="1.3.0"
script_date='2026/07/21'

export TZ="Europe/Rome"

# ------------------------------------------------------------------------------
# Paths

fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

fp_script="/hydro/library/fp_package_connectors/sm_model_tools/organizer_grid/sm_model_points2grid.py"
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/organizer_grid/sm_model_points2grid_history.json"

folder_lock="/hydro/lock/sm_cnr"
#fp_lock="${folder_lock}/run_points2grid_history.lock"

# ------------------------------------------------------------------------------
# Default options

time_arg=""
time_start=""
time_end=""

force_run=false
hour_override=""
previous_day=false
step_hours=24

period_mode=false

successful_runs=0
failed_runs=0
number_steps=1

# ------------------------------------------------------------------------------
# Usage

usage() {

    cat << EOF

Usage:

  $0 [REFERENCE_TIME] [OPTIONS]

Single-time mode:

  $0
  $0 "2026-07-13 10:00"
  $0 --hour "08:00"
  $0 --previous-day --hour "08:00"

Period mode:

  $0 --time-start "YYYY-MM-DD HH:MM" \\
     --time-end   "YYYY-MM-DD HH:MM" \\
     [--step-hours HOURS] \\
     [--hour "HH:MM"]

Options:

  -f, --force
      Remove the existing lock and terminate a process currently using it.

  --hour, -hour "HH:MM"
      Override the execution time.

      Valid format:
          HH:MM

      Valid range:
          "00:00" to "23:59"

      Examples:
          --hour "00:00"
          --hour "06:00"
          --hour "23:00"

      In period mode, the selected time is applied to both period bounds.

  --previous-day, --previous_day
      Shift the selected reference time or the complete period backward
      by one day.

  --time-start DATETIME
      Start of the processing period.

  --time-end DATETIME
      End of the processing period.

  --step-hours HOURS
      Period iteration interval in hours.

      Default:
          24

  -h, --help
      Show this help message.

Period direction:

  The period is always processed backward:

      time_end -> time_start

Failure handling:

  If one time step fails, a warning is printed and execution continues with
  the next time step.

EOF
}

# ------------------------------------------------------------------------------
# Parse command-line arguments

while [[ $# -gt 0 ]]; do

    case "$1" in

        -f|--force)
            force_run=true
            shift
            ;;

        --hour|-hour)
            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --hour requires a value in HH:MM format"
                echo "Example: --hour \"23:00\""
                exit 1
            fi

            hour_override="$2"
            shift 2
            ;;

        --previous-day|--previous_day)
            previous_day=true
            shift
            ;;

        --time-start|--time_start)
            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --time-start requires a datetime value"
                exit 1
            fi

            time_start="$2"
            shift 2
            ;;

        --time-end|--time_end)
            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --time-end requires a datetime value"
                exit 1
            fi

            time_end="$2"
            shift 2
            ;;

        --step-hours|--step_hours)
            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --step-hours requires a positive integer"
                exit 1
            fi

            step_hours="$2"
            shift 2
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
            if [[ -n "${time_arg}" ]]; then
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

# ------------------------------------------------------------------------------
# Validate hour override
#
# Accepted examples:
#
#   00:00
#   06:00
#   12:30
#   23:59

if [[ -n "${hour_override}" ]]; then

    if ! [[ "${hour_override}" =~ ^([01][0-9]|2[0-3]):([0-5][0-9])$ ]]; then
        echo "ERROR: Invalid --hour value '${hour_override}'"
        echo "Expected format: HH:MM"
        echo "Valid range    : 00:00 to 23:59"
        echo "Example        : --hour \"23:00\""
        exit 1
    fi

fi

# ------------------------------------------------------------------------------
# Validate step interval

if ! [[ "${step_hours}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Invalid --step-hours value '${step_hours}'"
    echo "Expected a positive integer"
    exit 1
fi

if [[ "${step_hours}" -le 0 ]]; then
    echo "ERROR: Invalid --step-hours value '${step_hours}'"
    echo "Expected a value greater than zero"
    exit 1
fi

# ------------------------------------------------------------------------------
# Validate period bounds

if [[ -n "${time_start}" && -z "${time_end}" ]]; then
    echo "ERROR: --time-start requires --time-end"
    exit 1
fi

if [[ -z "${time_start}" && -n "${time_end}" ]]; then
    echo "ERROR: --time-end requires --time-start"
    exit 1
fi

if [[ -n "${time_start}" && -n "${time_end}" ]]; then
    period_mode=true
fi

if ${period_mode} && [[ -n "${time_arg}" ]]; then

    echo "ERROR: A positional reference time cannot be used with period mode"
    echo
    echo "Use either:"
    echo
    echo "  $0 \"YYYY-MM-DD HH:MM\""
    echo
    echo "or:"
    echo
    echo "  $0 --time-start \"...\" --time-end \"...\""

    exit 1

fi

# ------------------------------------------------------------------------------
# Normalize time configuration

if ${period_mode}; then

    if ! time_start=$(date -d "${time_start}" "+%Y-%m-%d %H:%M" 2>/dev/null); then
        echo "ERROR: Invalid --time-start value"
        exit 1
    fi

    if ! time_end=$(date -d "${time_end}" "+%Y-%m-%d %H:%M" 2>/dev/null); then
        echo "ERROR: Invalid --time-end value"
        exit 1
    fi

    # Apply the HH:MM override to both period bounds.
    #
    # Example:
    #
    #   --time-start "2026-07-01"
    #   --time-end   "2026-07-13"
    #   --hour       "23:00"
    #
    # becomes:
    #
    #   2026-07-01 23:00
    #   2026-07-13 23:00

    if [[ -n "${hour_override}" ]]; then

        start_date=$(date -d "${time_start}" "+%Y-%m-%d")
        end_date=$(date -d "${time_end}" "+%Y-%m-%d")

        time_start="${start_date} ${hour_override}"
        time_end="${end_date} ${hour_override}"

    fi

    # Shift the complete period backward by one day.

    if ${previous_day}; then

        time_start=$(date -d "${time_start} -1 day" "+%Y-%m-%d %H:%M")
        time_end=$(date -d "${time_end} -1 day" "+%Y-%m-%d %H:%M")

    fi

    epoch_start=$(date -d "${time_start}" "+%s")
    epoch_end=$(date -d "${time_end}" "+%s")

    if [[ "${epoch_start}" -gt "${epoch_end}" ]]; then

        echo "ERROR: time_start must be earlier than or equal to time_end"
        echo "TIME START: ${time_start}"
        echo "TIME END  : ${time_end}"

        exit 1

    fi

    step_seconds=$((step_hours * 3600))
    period_seconds=$((epoch_end - epoch_start))

    number_steps=$((period_seconds / step_seconds + 1))

else

    # --------------------------------------------------------------------------
    # Single-time mode

    if [[ -z "${time_arg}" ]]; then

        time_now=$(date "+%Y-%m-%d %H:%M")

    else

        if ! time_now=$(date -d "${time_arg}" "+%Y-%m-%d %H:%M" 2>/dev/null); then
            echo "ERROR: Invalid reference time '${time_arg}'"
            exit 1
        fi

    fi

    base_date=$(date -d "${time_now}" "+%Y-%m-%d")
    original_time=$(date -d "${time_now}" "+%H:%M")

    if ${previous_day}; then
        base_date=$(date -d "${base_date} -1 day" "+%Y-%m-%d")
    fi

    if [[ -n "${hour_override}" ]]; then

        time_now="${base_date} ${hour_override}"

    else

        time_now="${base_date} ${original_time}"

    fi

    number_steps=1

fi

# ------------------------------------------------------------------------------
# Header

echo " ==================================================================================="
echo " ==> ${script_name}"
echo " ==> Version     : ${script_version}"
echo " ==> Release Date: ${script_date}"
echo " ==> START ..."
echo " ====> SYSTEM TIME    : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo " ====> PERIOD MODE    : ${period_mode}"

if ${period_mode}; then

    echo " ====> TIME START     : ${time_start}"
    echo " ====> TIME END       : ${time_end}"
    echo " ====> DIRECTION      : TIME END -> TIME START"
    echo " ====> STEP HOURS     : ${step_hours}"
    echo " ====> NUMBER OF RUNS : ${number_steps}"

else

    echo " ====> REFERENCE TIME : ${time_now}"

fi

echo " ====> HOUR OVERRIDE  : ${hour_override:-none}"
echo " ====> PREVIOUS DAY   : ${previous_day}"
echo " ====> FORCE RUN      : ${force_run}"

# ------------------------------------------------------------------------------
# Prepare lock folder

mkdir -p "${folder_lock}"

# Time tag (history mode)
if [[ -n "$time_start" && -n "$time_end" ]]; then
    time_tag="_$(date -d "$time_start" +"%Y%m%d%H%M")_$(date -d "$time_end" +"%Y%m%d%H%M")"
else
    time_tag=""
fi

# Lock file
fp_lock="${folder_lock}/run_sm_model_grid_history${time_tag}.lock"

# ------------------------------------------------------------------------------
# Load Python environment

echo " ====> LOAD PYTHON ENVIRONMENT ..."

if [[ -f "${fp_env_file}" ]]; then

    # shellcheck source=/dev/null
    source "${fp_env_file}"

else

    echo " ====> ENVIRONMENT FILE NOT FOUND: ${fp_env_file}"
    exit 1

fi

echo " ====> LOAD PYTHON ENVIRONMENT ... DONE"

# ------------------------------------------------------------------------------
# Validate required files

if [[ ! -f "${fp_script}" ]]; then

    echo "ERROR: Python script not found"
    echo "SCRIPT: ${fp_script}"

    exit 1

fi

if [[ ! -f "${fp_settings}" ]]; then

    echo "ERROR: Settings file not found"
    echo "SETTINGS: ${fp_settings}"

    exit 1

fi

# ------------------------------------------------------------------------------
# Force lock removal

if ${force_run}; then

    echo " ====> FORCE MODE ENABLED"

    if command -v fuser >/dev/null 2>&1; then
        fuser -k "${fp_lock}" >/dev/null 2>&1 || true
    fi

    rm -f "${fp_lock}"

fi

# ------------------------------------------------------------------------------
# Function to execute one time step

run_time_step() {

    local reference_time="$1"
    local step_index="$2"
    local total_steps="$3"

    local run_exit_code=0

    echo
    echo " -----------------------------------------------------------------------------------"
    echo " ====> RUN ${step_index}/${total_steps}"
    echo " ====> REFERENCE TIME: ${reference_time}"
    echo " ====> PYTHON SCRIPT : ${fp_script}"
    echo " ====> SETTINGS FILE : ${fp_settings}"

    local cmd=(
        python "${fp_script}"
        -settings_file "${fp_settings}"
        -time "${reference_time}"
    )

    echo " ====> COMMAND: ${cmd[*]}"

    if "${cmd[@]}"; then

        successful_runs=$((successful_runs + 1))

        echo " ====> RUN ${step_index}/${total_steps} ... DONE"

    else

        run_exit_code=$?
        failed_runs=$((failed_runs + 1))

        echo " ====> WARNING: RUN ${step_index}/${total_steps} ... FAILED"
        echo " ====> WARNING: REFERENCE TIME: ${reference_time}"
        echo " ====> WARNING: EXIT CODE     : ${run_exit_code}"
        echo " ====> WARNING: CONTINUE WITH NEXT TIME STEP"

    fi

    # Always return success so that set -e does not stop the complete period
    # when one Python execution fails.

    return 0
}

# ------------------------------------------------------------------------------
# Acquire lock

echo " ====> LOCK FILE: ${fp_lock}"

exec 9>"${fp_lock}"

if ! flock -n 9; then

    echo " ====> RUN POINTS2GRID ... SKIPPED"
    echo " ====> ANOTHER INSTANCE IS ALREADY RUNNING"
    echo " ====> USE -f OR --force TO OVERRIDE"

    exit 1

fi

echo " ====> LOCK ACQUIRED"

# ------------------------------------------------------------------------------
# Execute

if ${period_mode}; then

    current_time="${time_end}"
    step_index=0

    while true; do

        current_epoch=$(date -d "${current_time}" +%s)

        [[ ${current_epoch} -lt ${epoch_start} ]] && break

        step_index=$((step_index + 1))

        run_time_step \
            "${current_time}" \
            "${step_index}" \
            "${number_steps}"

        if (( step_hours % 24 == 0 )); then

            # Calendar-day iteration.
            # Keeps the selected HH:MM fixed across DST changes.
            days=$((step_hours / 24))

            current_date=$(date -d "${current_time}" +%F)
            current_hour=$(date -d "${current_time}" +%H:%M)

            current_time=$(date \
                -d "${current_date} -${days} day ${current_hour}" \
                "+%Y-%m-%d %H:%M")

        else

            # Sub-daily iteration.
            current_epoch=$((current_epoch - step_hours * 3600))
            current_time=$(date -d "@${current_epoch}" "+%Y-%m-%d %H:%M")

        fi

    done

else

    run_time_step \
        "${time_now}" \
        1 \
        1

fi

# ------------------------------------------------------------------------------
# Final summary

echo
echo " ==================================================================================="
echo " ====> EXECUTION SUMMARY"
echo " ====> EXPECTED RUNS  : ${number_steps}"
echo " ====> SUCCESSFUL RUNS: ${successful_runs}"
echo " ====> FAILED RUNS    : ${failed_runs}"

if [[ "${failed_runs}" -gt 0 ]]; then

    echo " ====> FINAL STATUS   : COMPLETED WITH WARNINGS"

else

    echo " ====> FINAL STATUS   : COMPLETED SUCCESSFULLY"

fi

echo " ==> ${script_name}"
echo " ==> Version     : ${script_version}"
echo " ==> Release Date: ${script_date}"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="

# The processing may contain failed time steps, but the complete period has been
# processed. Therefore, return success.

exit 0


