#!/usr/bin/env bash
set -Eeuo pipefail

# ==============================================================================
# SM MODEL POINTS MERGER
#
# Modes:
#
#   HISTORY
#   -------
#   Create a period ending at --time-now and extending --n-days backward.
#
#   Example:
#       ./sm_model_points_merger.sh \
#           --mode history \
#           --time-now "2026-07-12 00:00" \
#           --n-days 10
#
#   Period passed to Python:
#       time_start = 2026-07-02 00:00
#       time_end   = 2026-07-12 00:00
#
#
#   NRT - CURRENT DAY
#   -----------------
#   Use --hour to select the reference hour in the current day.
#
#   Example:
#       ./sm_model_points_merger.sh \
#           --mode nrt \
#           --time-now "2026-07-12 14:35" \
#           --hour 12
#
#   Time passed to Python:
#       2026-07-12 12:00
#
#
#   NRT - PREVIOUS DAY
#   ------------------
#   Use --previous-day and --hour-previous-day.
#
#   Example:
#       ./sm_model_points_merger.sh \
#           --mode nrt \
#           --time-now "2026-07-12 14:35" \
#           --previous-day \
#           --hour-previous-day 23
#
#   Time passed to Python:
#       2026-07-11 23:00
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Script information
script_name="SM MODEL POINTS MERGER"
script_version="1.0.0"
script_date="2026/07/20"

# ------------------------------------------------------------------------------
# Timezone
export TZ="Europe/Rome"

# ------------------------------------------------------------------------------
# Environment
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python application
fp_script="/hydro/library/fp_package_connectors/sm_model_tools/organizer_ts/sm_model_points_merger.py"

# JSON settings
fp_settings="/hydro/fp_tools_postprocessing/analyzer_sm_cnr/model_grid/organizer_ts/sm_model_points_merger.json"

# ------------------------------------------------------------------------------
# Logging
folder_log="/hydro/log/sm_cnr"
mkdir -p "${folder_log}"

time_system_tag="$(date '+%Y%m%d_%H%M')"
fp_log="${folder_log}/sm_model_points_merger_${time_system_tag}.log"

# ------------------------------------------------------------------------------
# Lock
folder_lock="/hydro/lock/sm_cnr"
mkdir -p "${folder_lock}"

fp_lock="${folder_lock}/run_sm_model_points_merger.lock"

# ------------------------------------------------------------------------------
# Defaults
script_mode="nrt"

time_now=""
n_days=1

hour_current_day=23
hour_previous_day=23

previous_day=false
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
      Create a time period from TIME_NOW - N_DAYS to TIME_NOW.

  --mode nrt
      Select one reference hour in the current or previous day.

Common options:

  --time-now "YYYY-MM-DD HH:MM"
      Reference time.
      Default: current system time in Europe/Rome.

  -f, --force
      Remove an existing lock file and force execution.

  -h, --help
      Show this help.

History options:

  --n-days N
      Number of days before --time-now.
      Default: ${n_days}

NRT current-day options:

  --hour H
  --hour-current-day H
      Hour selected in the current day.
      Valid range: 0-23.
      Default: ${hour_current_day}

NRT previous-day options:

  --previous-day
      Select the previous calendar day.

  --hour-previous-day H
      Hour selected when --previous-day is active.
      Valid range: 0-23.
      Default: ${hour_previous_day}

Examples:

  $0 \\
      --mode history \\
      --time-now "2026-07-12 00:00" \\
      --n-days 10

  $0 \\
      --mode nrt \\
      --time-now "2026-07-12 14:35" \\
      --hour 12

  $0 \\
      --mode nrt \\
      --time-now "2026-07-12 14:35" \\
      --previous-day \\
      --hour-previous-day 23

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

validate_hour() {
    local hour_value="$1"
    local option_name="$2"

    if ! [[ "${hour_value}" =~ ^[0-9]+$ ]]; then
        error_exit "${option_name} must be an integer between 0 and 23."
    fi

    if (( hour_value < 0 || hour_value > 23 )); then
        error_exit "${option_name} must be between 0 and 23."
    fi
}

validate_positive_integer() {
    local value="$1"
    local option_name="$2"

    if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
        error_exit "${option_name} must be a non-negative integer."
    fi
}

validate_datetime() {
    local time_value="$1"

    if ! date -d "${time_value}" '+%Y-%m-%d %H:%M' > /dev/null 2>&1; then
        error_exit "Invalid datetime '${time_value}'. Expected YYYY-MM-DD HH:MM."
    fi
}

cleanup_lock() {
    if [[ -f "${fp_lock}" ]]; then
        rm -f "${fp_lock}"
    fi
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

        --n-days)
            [[ $# -ge 2 ]] || error_exit "Missing value after --n-days."
            n_days="$2"
            shift 2
            ;;

        --hour|--hour-current-day)
            [[ $# -ge 2 ]] || error_exit "Missing value after $1."
            hour_current_day="$2"
            shift 2
            ;;

        --previous-day)
            previous_day=true
            shift
            ;;

        --hour-previous-day)
            [[ $# -ge 2 ]] || error_exit "Missing value after --hour-previous-day."
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
# Normalize mode
script_mode="$(printf '%s' "${script_mode}" | tr '[:upper:]' '[:lower:]')"

if [[ "${script_mode}" != "history" && "${script_mode}" != "nrt" ]]; then
    error_exit "--mode must be either 'history' or 'nrt'."
fi

# ------------------------------------------------------------------------------
# Default reference time
if [[ -z "${time_now}" ]]; then
    time_now="$(date '+%Y-%m-%d %H:%M')"
fi

validate_datetime "${time_now}"
validate_positive_integer "${n_days}" "--n-days"
validate_hour "${hour_current_day}" "--hour"
validate_hour "${hour_previous_day}" "--hour-previous-day"

# Normalize reference time
time_now="$(date -d "${time_now}" '+%Y-%m-%d %H:%M')"

# ------------------------------------------------------------------------------
# Check required files
[[ -f "${fp_env_file}" ]] || error_exit \
    "Environment file not found: ${fp_env_file}"

[[ -f "${fp_script}" ]] || error_exit \
    "Python script not found: ${fp_script}"

[[ -f "${fp_settings}" ]] || error_exit \
    "Settings file not found: ${fp_settings}"

# ------------------------------------------------------------------------------
# Manage lock
if [[ -f "${fp_lock}" ]]; then

    if [[ "${force_run}" == true ]]; then
        log_message "WARNING: removing existing lock file: ${fp_lock}"
        rm -f "${fp_lock}"
    else
        error_exit "Lock file already exists: ${fp_lock}"
    fi

fi

touch "${fp_lock}"
trap cleanup_lock EXIT INT TERM

# ------------------------------------------------------------------------------
# Load environment
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
log_message "==> Input reference time: ${time_now}"
log_message "==> Log file: ${fp_log}"
log_message "================================================================================"

# ------------------------------------------------------------------------------
# Build and run command
command_args=(
    python
    "${fp_script}"
    -settings_file
    "${fp_settings}"
)

if [[ "${script_mode}" == "history" ]]; then

    # The period includes:
    #
    #   time_start = time_now - n_days
    #   time_end   = time_now
    #
    # Example:
    #   time_now = 2026-07-12 00:00
    #   n_days   = 10
    #
    # Result:
    #   time_start = 2026-07-02 00:00
    #   time_end   = 2026-07-12 00:00

    time_start="$(date -d "${time_now} - ${n_days} days" '+%Y-%m-%d %H:%M')"
    time_end="${time_now}"

    log_message "==> HISTORY TIME START: ${time_start}"
    log_message "==> HISTORY TIME END:   ${time_end}"
    log_message "==> HISTORY N DAYS:     ${n_days}"

    command_args+=(
        --time-start
        "${time_start}"
        --time-end
        "${time_end}"
    )

else

    current_date="$(date -d "${time_now}" '+%Y-%m-%d')"

    if [[ "${previous_day}" == true ]]; then

        selected_date="$(date -d "${current_date} - 1 day" '+%Y-%m-%d')"

        printf -v selected_hour '%02d' "${hour_previous_day}"

        time_reference="${selected_date} ${selected_hour}:00"

        log_message "==> NRT DAY MODE:       previous day"
        log_message "==> NRT SELECTED DATE:  ${selected_date}"
        log_message "==> NRT SELECTED HOUR:  ${selected_hour}:00"
        log_message "==> NRT REFERENCE TIME: ${time_reference}"

    else

        selected_date="${current_date}"

        printf -v selected_hour '%02d' "${hour_current_day}"

        time_reference="${selected_date} ${selected_hour}:00"

        log_message "==> NRT DAY MODE:       current day"
        log_message "==> NRT SELECTED DATE:  ${selected_date}"
        log_message "==> NRT SELECTED HOUR:  ${selected_hour}:00"
        log_message "==> NRT REFERENCE TIME: ${time_reference}"

    fi

    command_args+=(
        --time
        "${time_reference}"
    )

fi

# Add force flag to Python too, if supported
if [[ "${force_run}" == true ]]; then
    command_args+=(
        --force
    )
fi

# ------------------------------------------------------------------------------
# Print command safely
printf -v command_print '%q ' "${command_args[@]}"
log_message "==> COMMAND: ${command_print}"

# ------------------------------------------------------------------------------
# Run
set +e

"${command_args[@]}" 2>&1 | tee -a "${fp_log}"
exit_code=${PIPESTATUS[0]}

set -e

# ------------------------------------------------------------------------------
# Check result
if (( exit_code != 0 )); then
    error_exit "Python merger failed with exit code ${exit_code}."
fi

log_message "================================================================================"
log_message "==> ${script_name} completed successfully"
log_message "==> End time: $(date '+%Y-%m-%d %H:%M:%S %Z')"
log_message "================================================================================"

exit 0