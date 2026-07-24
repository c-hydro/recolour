#!/usr/bin/env bash

# ==============================================================================
# RECOLOUR APPS - GRID METRICS
#
# The Python application manages the complete time period internally.
# This Bash runner executes the Python application once.
#
# Required:
#
#   --time-start "YYYY-MM-DD HH:MM"
#
# Optional:
#
#   --time-end "YYYY-MM-DD HH:MM"
#   --hour "HH:MM"
#   --previous-day
#   --reset-checkpoint
#   --force
#
# If --time-end is omitted:
#
#   without --previous-day:
#       time_end = current day at --hour
#
#   with --previous-day:
#       time_end = previous day at --hour
#
# If --hour is omitted:
#
#       hour = "23:00"
#
# Examples:
#
#   ./app_grid_metrics_history.sh \
#       --time-start "2024-01-01 23:00"
#
#   ./app_grid_metrics_history.sh \
#       --time-start "2024-01-01 23:00" \
#       --hour "06:00"
#
#   ./app_grid_metrics_history.sh \
#       --time-start "2024-01-01 23:00" \
#       --previous-day
#
#   ./app_grid_metrics_history.sh \
#       --time-start "2024-01-01 23:00" \
#       --previous-day \
#       --hour "06:00"
#
#   ./app_grid_metrics_history.sh \
#       --time-start "2024-01-01 23:00" \
#       --time-end "2025-12-31 23:00"
#
# ==============================================================================

set -u
set -o pipefail

# ------------------------------------------------------------------------------
# Script information
script_name="RECOLOUR APPS - GRID METRICS - UPDATE"
script_version="1.1.0"
script_date="2026/07/23"

export TZ="Europe/Rome"

# ------------------------------------------------------------------------------
# Paths
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"
fp_script="/hydro/library/fp_package_connectors/sm_model_tools/organizer_metrics/app_grid_metrics.py"
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/organizer_metrics/app_grid_metrics_update.json"

folder_lock="/hydro/lock/sm_cnr"
lock_name="run_app_grid_metrics_update"

# ------------------------------------------------------------------------------
# Default arguments

time_start=""
time_end=""
run_hour="23:00"

previous_day=false
force_run=false
reset_checkpoint=false

# ------------------------------------------------------------------------------
# Usage

usage() {

    cat << EOF

Usage:

  $0 --time-start "YYYY-MM-DD HH:MM" [OPTIONS]

Required arguments:

  --time-start, --time_start DATETIME

      Define the starting time of the analysis period.

      Expected format:

          "YYYY-MM-DD HH:MM"

      Example:

          --time-start "2024-01-01 23:00"


Optional arguments:

  --time-end, --time_end DATETIME

      Define the ending time of the analysis period.

      Expected format:

          "YYYY-MM-DD HH:MM"

      If omitted, time_end is automatically calculated using the current
      date and the hour selected with --hour.

      Without --previous-day:

          time_end = current day at --hour

      With --previous-day:

          time_end = previous day at --hour


  --hour, --run-hour, --run_hour TIME

      Define the hour used to calculate time_end when --time-end is omitted.

      Expected format:

          "HH:MM"

      Default:

          "23:00"

      Examples:

          --hour "00:00"
          --hour "06:00"
          --hour "23:00"


  --previous-day, --previous_day

      When --time-end is omitted, calculate time_end using the previous day.

      Example:

          current system date: 2026-07-23
          selected hour:       23:00

          resulting time_end:  2026-07-22 23:00


  --reset-checkpoint, --reset_checkpoint

      Pass --reset-checkpoint to the Python application.

      An existing checkpoint is ignored and overwritten.


  -f, --force

      Remove an existing Bash lock and terminate a process currently
      associated with the lock file, when possible.


  -h, --help

      Show this help message.


Examples:

  Use today at 23:00 as time_end:

      $0 \\
          --time-start "2024-01-01 23:00"


  Use today at 06:00 as time_end:

      $0 \\
          --time-start "2024-01-01 23:00" \\
          --hour "06:00"


  Use the previous day at 23:00 as time_end:

      $0 \\
          --time-start "2024-01-01 23:00" \\
          --previous-day


  Use the previous day at 06:00 as time_end:

      $0 \\
          --time-start "2024-01-01 23:00" \\
          --previous-day \\
          --hour "06:00"


  Explicitly define time_end:

      $0 \\
          --time-start "2024-01-01 23:00" \\
          --time-end   "2025-12-31 23:00"


  Explicitly define time_end and reset the checkpoint:

      $0 \\
          --time-start "2024-01-01 23:00" \\
          --time-end   "2025-12-31 23:00" \\
          --reset-checkpoint

EOF
}

# ------------------------------------------------------------------------------
# Parse command-line arguments

while [[ $# -gt 0 ]]; do

    case "$1" in

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

        --hour|--run-hour|--run_hour)

            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --hour requires a time value"
                exit 1
            fi

            run_hour="$2"
            shift 2
            ;;

        --previous-day|--previous_day)

            previous_day=true
            shift
            ;;

        --reset-checkpoint|--reset_checkpoint)

            reset_checkpoint=true
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
            echo
            usage
            exit 1
            ;;

        *)

            echo "ERROR: Unexpected positional argument '$1'"
            echo
            usage
            exit 1
            ;;

    esac

done

# ------------------------------------------------------------------------------
# Validate mandatory time_start

if [[ -z "${time_start}" ]]; then

    echo "ERROR: --time-start must be specified"
    echo
    echo "Expected format:"
    echo
    echo "  --time-start \"YYYY-MM-DD HH:MM\""
    echo
    echo "Example:"
    echo
    echo "  $0 --time-start \"2024-01-01 23:00\""

    exit 1

fi

# ------------------------------------------------------------------------------
# Validate run hour
#
# Accepted format:
#
#   HH:MM
#
# Examples:
#
#   00:00
#   06:00
#   23:00

if ! [[ "${run_hour}" =~ ^([01][0-9]|2[0-3]):[0-5][0-9]$ ]]; then

    echo "ERROR: Invalid --hour value"
    echo "VALUE          : ${run_hour}"
    echo "EXPECTED FORMAT: HH:MM"
    echo "EXAMPLE        : 23:00"

    exit 1

fi

# ------------------------------------------------------------------------------
# Normalize and validate time_start

time_start_input="${time_start}"

if ! time_start=$(date -d "${time_start_input}" "+%Y-%m-%d %H:%M" 2>/dev/null); then

    echo "ERROR: Invalid --time-start value"
    echo "VALUE          : ${time_start_input}"
    echo "EXPECTED FORMAT: YYYY-MM-DD HH:MM"

    exit 1

fi

# ------------------------------------------------------------------------------
# Set default time_end
#
# If time_end is omitted:
#
#   previous_day=false:
#       current date at run_hour
#
#   previous_day=true:
#       previous date at run_hour

time_end_automatic=false

if [[ -z "${time_end}" ]]; then

    time_end_automatic=true

    if ${previous_day}; then
        reference_date=$(date -d "1 day ago" "+%Y-%m-%d")
    else
        reference_date=$(date "+%Y-%m-%d")
    fi

    time_end="${reference_date} ${run_hour}"

fi

# ------------------------------------------------------------------------------
# Normalize and validate time_end

time_end_input="${time_end}"

if ! time_end=$(date -d "${time_end_input}" "+%Y-%m-%d %H:%M" 2>/dev/null); then

    echo "ERROR: Invalid --time-end value"
    echo "VALUE          : ${time_end_input}"
    echo "EXPECTED FORMAT: YYYY-MM-DD HH:MM"

    exit 1

fi

# ------------------------------------------------------------------------------
# Validate period order

epoch_start=$(date -d "${time_start}" "+%s")
epoch_end=$(date -d "${time_end}" "+%s")

if [[ "${epoch_start}" -gt "${epoch_end}" ]]; then

    echo "ERROR: time_start must be earlier than or equal to time_end"
    echo "TIME START: ${time_start}"
    echo "TIME END  : ${time_end}"

    exit 1

fi

# ------------------------------------------------------------------------------
# Define time-tagged lock file

mkdir -p "${folder_lock}"

time_start_tag=$(date -d "${time_start}" "+%Y%m%d%H%M")
time_end_tag=$(date -d "${time_end}" "+%Y%m%d%H%M")

time_tag="_${time_start_tag}_${time_end_tag}"

fp_lock="${folder_lock}/${lock_name}${time_tag}.lock"

# ------------------------------------------------------------------------------
# Header

echo " ==================================================================================="
echo " ==> ${script_name}"
echo " ==> Version         : ${script_version}"
echo " ==> Release Date    : ${script_date}"
echo " ==> START ..."
echo " ====> SYSTEM TIME      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo " ====> SETTINGS FILE    : ${fp_settings}"
echo " ====> PYTHON SCRIPT    : ${fp_script}"
echo " ====> TIME START       : ${time_start}"
echo " ====> TIME END         : ${time_end}"
echo " ====> TIME END AUTO    : ${time_end_automatic}"
echo " ====> RUN HOUR         : ${run_hour}"
echo " ====> PREVIOUS DAY     : ${previous_day}"
echo " ====> RESET CHECKPOINT : ${reset_checkpoint}"
echo " ====> FORCE RUN        : ${force_run}"
echo " ====> LOCK FILE        : ${fp_lock}"

# ------------------------------------------------------------------------------
# Load Python environment

echo " ====> LOAD PYTHON ENVIRONMENT ..."

if [[ -f "${fp_env_file}" ]]; then

    # shellcheck source=/dev/null
    source "${fp_env_file}"

else

    echo "ERROR: Environment file not found"
    echo "FILE: ${fp_env_file}"

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

    if [[ -f "${fp_lock}" ]]; then

        if command -v fuser >/dev/null 2>&1; then

            echo " ====> TERMINATE PROCESS USING LOCK FILE, IF ANY"
            fuser -k "${fp_lock}" >/dev/null 2>&1 || true

        fi

        echo " ====> REMOVE LOCK FILE: ${fp_lock}"
        rm -f "${fp_lock}"

    else

        echo " ====> LOCK FILE DOES NOT EXIST"

    fi

fi

# ------------------------------------------------------------------------------
# Acquire lock

exec 9>"${fp_lock}"

if ! flock -n 9; then

    echo " ====> GRID METRICS ... SKIPPED"
    echo " ====> ANOTHER INSTANCE IS ALREADY RUNNING"
    echo " ====> LOCK FILE: ${fp_lock}"
    echo " ====> USE -f OR --force TO OVERRIDE"

    exit 1

fi

echo " ====> LOCK ACQUIRED"

# Write the current Bash process ID into the lock file.

printf '%s\n' "$$" >&9

# ------------------------------------------------------------------------------
# Build Python command

cmd=(
    python "${fp_script}"
    -settings_file "${fp_settings}"
    -time_start "${time_start}"
    -time_end "${time_end}"
)

if ${reset_checkpoint}; then
    cmd+=(--reset-checkpoint)
fi

# ------------------------------------------------------------------------------
# Execute application

echo
echo " -----------------------------------------------------------------------------------"
echo " ====> RUN GRID METRICS ..."
echo " ====> COMMAND:"

printf "      "
printf "%q " "${cmd[@]}"
printf "\n"

echo " -----------------------------------------------------------------------------------"
echo

start_epoch=$(date "+%s")

"${cmd[@]}"
run_exit_code=$?

end_epoch=$(date "+%s")
elapsed_seconds=$((end_epoch - start_epoch))

# ------------------------------------------------------------------------------
# Final status

echo
echo " ==================================================================================="

if [[ "${run_exit_code}" -eq 0 ]]; then

    echo " ====> GRID METRICS ... DONE"
    echo " ====> FINAL STATUS   : COMPLETED SUCCESSFULLY"

else

    echo " ====> GRID METRICS ... FAILED"
    echo " ====> FINAL STATUS   : FAILED"
    echo " ====> EXIT CODE      : ${run_exit_code}"

fi

echo " ====> TIME START     : ${time_start}"
echo " ====> TIME END       : ${time_end}"
echo " ====> ELAPSED TIME   : ${elapsed_seconds} seconds"
echo " ==> ${script_name}"
echo " ==> Version         : ${script_version}"
echo " ==> Release Date    : ${script_date}"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="

exit "${run_exit_code}"
