#!/bin/bash -e

# ==============================================================================
# SM MODEL POINTS MERGER
#
# Modes:
#
#   HISTORY
#   -------
#   Iterate from --time-start to --time-end, using one reference time per day.
#
#   The reference time for each iteration is created using:
#
#       YYYY-MM-DD + --hour-current-day
#
#   Bash receives hours in HH:00 format:
#
#       --hour-current-day "23:00"
#       --hour-previous-day "00:00"
#
#   Python receives integer hours:
#
#       --hour-current-day 23
#       --hour-previous-day 0
#
#   Example:
#
#       ./sm_model_points_merger.sh \
#           --mode history \
#           --time-start "2024-01-01" \
#           --time-end "2024-12-31" \
#           --hour-current-day "23:00" \
#           --hour-previous-day "00:00" \
#           --force
#
#
#   NRT
#   ---
#   Run Python once using --time-now as reference time.
#
#   Example:
#
#       ./sm_model_points_merger.sh \
#           --mode nrt \
#           --time-now "2026-07-20 14:58" \
#           --hour-current-day "11:00" \
#           --hour-previous-day "00:00"
#
# ==============================================================================

set -uo pipefail

# ------------------------------------------------------------------------------
# Script information
script_name="SM MODEL POINTS MERGER - AIRT - HISTORY (DB)"
script_version="1.1.1"
script_date="2026/07/23"

# ------------------------------------------------------------------------------
# Timezone
export TZ="Europe/Rome"

# ------------------------------------------------------------------------------
# Environment
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python application
fp_script="/hydro/library/fp_package_connectors/sm_model_tools/organizer_ts/sm_model_points_merger.py"

# JSON settings
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/organizer_ts/sm_model_points_merger_airt_history.json"

# ------------------------------------------------------------------------------
# Logging
folder_log="/hydro/log/sm_cnr"
mkdir -p "${folder_log}"

time_system_tag="$(date '+%Y%m%d_%H%M')"
fp_log="${folder_log}/sm_model_points_merger_${time_system_tag}_airt_history.log"

# ------------------------------------------------------------------------------
# Lock
folder_lock="/hydro/lock/sm_cnr"
mkdir -p "${folder_lock}"

fp_lock="${folder_lock}/run_sm_model_points_merger_airt_history.lock"

# ------------------------------------------------------------------------------
# Defaults
script_mode="nrt"

time_now=""
time_start=""
time_end=""

# Bash receives hours using HH:00.
hour_current_day="23:00"
hour_previous_day="23:00"

force_run=false

# ------------------------------------------------------------------------------
# Functions
usage() {
    cat <<EOF

Usage:

  $0 --mode history [OPTIONS]
  $0 --mode nrt [OPTIONS]

Modes:

  --mode history
      Iterate from --time-start to --time-end, inclusive.
      Python is executed once for every day.

  --mode nrt
      Execute Python once using --time-now as reference time.

Common options:

  --hour H
  --hour-current-day "HH:00"
      Reference hour used for the current-day file.
      Valid range: "00:00" to "23:00".
      Default: "${hour_current_day}"

  --hour-previous-day "HH:00"
      Hour used for the previous-day and historical files.
      Valid range: "00:00" to "23:00".
      Default: "${hour_previous_day}"

  -f, --force
      Remove an existing lock file before execution.
      This option is handled only by Bash.

  -h, --help
      Show this help.

History options:

  --time-start "YYYY-MM-DD"
      First date to process.

      A datetime is also accepted, but only its date component is used.
      The reference hour is defined by --hour-current-day.

  --time-end "YYYY-MM-DD"
      Last date to process, inclusive.

      A datetime is also accepted, but only its date component is used.
      The reference hour is defined by --hour-current-day.

NRT options:

  --time-now "YYYY-MM-DD HH:MM"
      Reference time passed to Python.
      Default: current system time in Europe/Rome.

Examples:

  $0 \
      --mode history \
      --time-start "2024-01-01" \
      --time-end "2024-01-10" \
      --hour-current-day "23:00" \
      --hour-previous-day "00:00" \
      --force

  $0 \
      --mode nrt \
      --time-now "2026-07-20 14:58" \
      --hour-current-day "11:00" \
      --hour-previous-day "00:00"

  $0 \
      --mode nrt \
      --hour-current-day "06:00" \
      --hour-previous-day "23:00" \
      --force

EOF
}

log_message() {
    local message="$1"

    printf '%s\n' "${message}"
    printf '%s\n' "${message}" >> "${fp_log}"
}

error_exit() {
    local message="$1"

    log_message "ERROR: ${message}"
    exit 1
}

warning_message() {
    local message="$1"

    log_message "WARNING: ${message}"
}

validate_hour() {
    local hour_value="$1"
    local option_name="$2"

    if ! [[ "${hour_value}" =~ ^([01][0-9]|2[0-3]):00$ ]]; then
        error_exit \
            "${option_name} must use HH:00 format between 00:00 and 23:00."
    fi
}

hour_to_integer() {
    local hour_value="$1"
    local hour_component

    hour_component="${hour_value%%:*}"

    # Base 10 is explicitly used because 08 and 09 could otherwise be
    # interpreted as invalid octal numbers.
    printf '%d' "$((10#${hour_component}))"
}

validate_datetime() {
    local time_value="$1"
    local option_name="$2"

    if ! date -d "${time_value}" '+%Y-%m-%d %H:%M' > /dev/null 2>&1; then
        error_exit \
            "Invalid value '${time_value}' for ${option_name}. Expected a valid date or datetime."
    fi
}

cleanup_lock() {
    if [[ -f "${fp_lock}" ]]; then
        rm -f "${fp_lock}"
    fi
}

print_command() {
    local command_print

    printf -v command_print '%q ' "$@"
    log_message "==> COMMAND: ${command_print}"
}

run_python_command() {
    local reference_time="$1"
    local hour_current_int="$2"
    local hour_previous_int="$3"

    local command_args=(
        python
        "${fp_script}"
        -settings_file
        "${fp_settings}"
        --time
        "${reference_time}"
        --hour-current-day
        "${hour_current_int}"
        --hour-previous-day
        "${hour_previous_int}"
    )

    print_command "${command_args[@]}"

    "${command_args[@]}" 2>&1 | tee -a "${fp_log}"

    # Return the exit code of Python, not the exit code of tee.
    return "${PIPESTATUS[0]}"
}

# ------------------------------------------------------------------------------
# Parse arguments
while [[ $# -gt 0 ]]; do

    case "$1" in

        --mode)
            [[ $# -ge 2 ]] || error_exit "Missing value after --mode."
            script_mode="$2"
            shift 2
            ;;

        --time-now)
            [[ $# -ge 2 ]] || error_exit "Missing value after --time-now."
            time_now="$2"
            shift 2
            ;;

        --time-start)
            [[ $# -ge 2 ]] || error_exit "Missing value after --time-start."
            time_start="$2"
            shift 2
            ;;

        --time-end)
            [[ $# -ge 2 ]] || error_exit "Missing value after --time-end."
            time_end="$2"
            shift 2
            ;;

        --hour|--hour-current-day)
            [[ $# -ge 2 ]] || error_exit "Missing value after $1."
            hour_current_day="$2"
            shift 2
            ;;

        --hour-previous-day)
            [[ $# -ge 2 ]] || \
                error_exit "Missing value after --hour-previous-day."

            hour_previous_day="$2"
            shift 2
            ;;

        -f|--force)
            force_run=true
            shift
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            error_exit "Unknown argument: $1"
            ;;

    esac

done

# ------------------------------------------------------------------------------
# Normalize and validate mode
script_mode="$(
    printf '%s' "${script_mode}" |
        tr '[:upper:]' '[:lower:]'
)"

if [[ "${script_mode}" != "history" && "${script_mode}" != "nrt" ]]; then
    error_exit "--mode must be either 'history' or 'nrt'."
fi

# ------------------------------------------------------------------------------
# Validate common hour arguments
validate_hour "${hour_current_day}" "--hour-current-day"
validate_hour "${hour_previous_day}" "--hour-previous-day"

# Convert HH:00 values to integers for Python.
hour_current_day_int="$(
    hour_to_integer "${hour_current_day}"
)"

hour_previous_day_int="$(
    hour_to_integer "${hour_previous_day}"
)"

# ------------------------------------------------------------------------------
# Validate mode-specific arguments
if [[ "${script_mode}" == "history" ]]; then

    if [[ -z "${time_start}" ]]; then
        error_exit "--time-start is required in history mode."
    fi

    if [[ -z "${time_end}" ]]; then
        error_exit "--time-end is required in history mode."
    fi

    validate_datetime "${time_start}" "--time-start"
    validate_datetime "${time_end}" "--time-end"

    # History mode uses only the date component from CLI values.
    history_date_start="$(date -d "${time_start}" '+%Y-%m-%d')"
    history_date_end="$(date -d "${time_end}" '+%Y-%m-%d')"

    # Use UTC for date-only calculations so daylight-saving transitions
    # cannot produce 23-hour or 25-hour differences.
    history_epoch_start="$(
        TZ=UTC date -d "${history_date_start} 00:00:00" '+%s'
    )"

    history_epoch_end="$(
        TZ=UTC date -d "${history_date_end} 00:00:00" '+%s'
    )"

    if (( history_epoch_start > history_epoch_end )); then
        error_exit \
            "--time-start (${history_date_start}) must not be later than --time-end (${history_date_end})."
    fi

else

    if [[ -z "${time_now}" ]]; then
        time_now="$(date '+%Y-%m-%d %H:%M')"
    fi

    validate_datetime "${time_now}" "--time-now"

    # Normalize NRT reference time.
    time_now="$(date -d "${time_now}" '+%Y-%m-%d %H:%M')"

fi

# ------------------------------------------------------------------------------
# Check required files
[[ -f "${fp_env_file}" ]] || \
    error_exit "Environment file not found: ${fp_env_file}"

[[ -f "${fp_script}" ]] || \
    error_exit "Python script not found: ${fp_script}"

[[ -f "${fp_settings}" ]] || \
    error_exit "Settings file not found: ${fp_settings}"

# ------------------------------------------------------------------------------
# Manage lock
#
# The force option is used only by this Bash script and is not passed to Python.
if [[ -f "${fp_lock}" ]]; then

    if [[ "${force_run}" == true ]]; then
        warning_message "Removing existing lock file: ${fp_lock}"
        rm -f "${fp_lock}"
    else
        error_exit "Lock file already exists: ${fp_lock}"
    fi

fi

touch "${fp_lock}"
trap cleanup_lock EXIT INT TERM

# ------------------------------------------------------------------------------
# Load Python environment
# shellcheck disable=SC1090
source "${fp_env_file}"

# ------------------------------------------------------------------------------
# Header
log_message "================================================================================"
log_message "==> ${script_name}"
log_message "==> Version: ${script_version}"
log_message "==> Release date: ${script_date}"
log_message "==> Mode: ${script_mode}"
log_message "==> System time: $(date '+%Y-%m-%d %H:%M:%S %Z')"
log_message "==> Current-day hour: ${hour_current_day} (${hour_current_day_int})"
log_message "==> Previous-day hour: ${hour_previous_day} (${hour_previous_day_int})"
log_message "==> Log file: ${fp_log}"
log_message "================================================================================"

# ------------------------------------------------------------------------------
# Execute history mode
if [[ "${script_mode}" == "history" ]]; then

    # Calculate the inclusive number of calendar dates.
    history_days=$(( 
        (history_epoch_end - history_epoch_start) / 86400 + 1
    ))

    log_message "==> HISTORY DATE START:       ${history_date_start}"
    log_message "==> HISTORY DATE END:         ${history_date_end}"
    log_message "==> HISTORY NUMBER OF DAYS:   ${history_days}"
    log_message "==> HISTORY REFERENCE HOUR:   ${hour_current_day}"
    log_message "==> HISTORY PREVIOUS HOUR:    ${hour_previous_day}"
    log_message "================================================================================"

    current_date="${history_date_start}"

    history_index=0
    history_success=0
    history_failed=0

    while true; do

        history_index=$((history_index + 1))

        # Build the reference time using the current date and configured
        # current-day hour.
        reference_time="${current_date} ${hour_current_day}"

        log_message "--------------------------------------------------------------------------------"
        log_message "==> HISTORY STEP: ${history_index}/${history_days}"
        log_message "==> HISTORY DATE: ${current_date}"
        log_message "==> HISTORY REFERENCE TIME: ${reference_time}"
        log_message "--------------------------------------------------------------------------------"

        set +e

        run_python_command \
            "${reference_time}" \
            "${hour_current_day_int}" \
            "${hour_previous_day_int}"

        exit_code=$?

        set -e

        if (( exit_code != 0 )); then
            history_failed=$((history_failed + 1))

            warning_message \
                "Python merger failed for ${reference_time} with exit code ${exit_code}. Continuing with the next date."
        else
            history_success=$((history_success + 1))

            log_message \
                "==> HISTORY STEP COMPLETED SUCCESSFULLY: ${reference_time}"
        fi

        # The end date is included.
        if [[ "${current_date}" == "${history_date_end}" ]]; then
            break
        fi

        # Increment by one calendar date. Using noon avoids edge cases around
        # daylight-saving transitions when Europe/Rome changes UTC offset.
        current_date="$(
            date -d "${current_date} 12:00:00 + 1 day" '+%Y-%m-%d'
        )"

    done

    log_message "================================================================================"
    log_message "==> HISTORY EXECUTION SUMMARY"
    log_message "==> Requested dates: ${history_days}"
    log_message "==> Successful:      ${history_success}"
    log_message "==> Failed:          ${history_failed}"
    log_message "================================================================================"

    if (( history_failed > 0 )); then
        warning_message \
            "History completed with ${history_failed} failed date(s). Check the log file for details."
    fi

# ------------------------------------------------------------------------------
# Execute NRT mode
else

    current_date="$(date -d "${time_now}" '+%Y-%m-%d')"

    time_current_day="${current_date} ${hour_current_day}"
    time_previous_day="${current_date} ${hour_previous_day}"

    log_message "==> NRT REFERENCE TIME:       ${time_now}"
    log_message "==> NRT CURRENT DATE:         ${current_date}"
    log_message "==> NRT CURRENT-DAY HOUR:     ${hour_current_day}"
    log_message "==> NRT PREVIOUS-DAY HOUR:    ${hour_previous_day}"
    log_message "==> NRT CURRENT FILE TIME:    ${time_current_day}"
    log_message "==> NRT HISTORY FILE TIME:    ${time_previous_day}"
    log_message "================================================================================"

    set +e

    run_python_command \
        "${time_now}" \
        "${hour_current_day_int}" \
        "${hour_previous_day_int}"

    exit_code=$?

    set -e

    if (( exit_code != 0 )); then
        error_exit "Python merger failed with exit code ${exit_code}."
    fi

fi

# ------------------------------------------------------------------------------
# Completed
log_message "================================================================================"
log_message "==> ${script_name} completed"
log_message "==> End time: $(date '+%Y-%m-%d %H:%M:%S %Z')"
log_message "================================================================================"

exit 0

