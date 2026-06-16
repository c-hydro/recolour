#!/bin/bash -e

# -----------------------------------------------------------------------------------------
# Script information
script_name='SM MODEL POINTS EXTRACTOR - RAIN - REALTIME'
script_version="1.0.4"
script_date='2026/06/12'

# Environment file
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python script and settings
fp_script="/hydro/library/fp_package_connectors/sm_model_tools/sm_model_points_extractor.py"
fp_settings="/hydro/library/fp_package_connectors/sm_model_tools/sm_model_points_extractor_rain.json"

# Lock folder/file
folder_lock="/hydro/lock/sm_cnr"
fp_lock="${folder_lock}/run_extract_points_rain.lock"

# -----------------------------------------------------------------------------------------
# Parse arguments
# Accepted:
#   ./sm_model_points_extractor.sh
#   ./sm_model_points_extractor.sh "2026-06-10 12:00"
#   ./sm_model_points_extractor.sh -f
#   ./sm_model_points_extractor.sh --force
#   ./sm_model_points_extractor.sh "2026-06-10 12:00" -f
#   ./sm_model_points_extractor.sh -f "2026-06-10 12:00"
# -----------------------------------------------------------------------------------------
time_now=""
force_run=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)
            force_run=true
            shift
            ;;
        *)
            time_now="$1"
            shift
            ;;
    esac
done

if [ -z "${time_now}" ]; then
    time_now=$(date "+%Y-%m-%d %H:%M")
fi

# -----------------------------------------------------------------------------------------

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo " ====> REFERENCE TIME: ${time_now}"
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
echo " ====> LOCK FILE    : ${fp_lock}"

(
    flock -n 9 || {
        echo " ====> RUN POINT EXTRACTION ... SKIPPED"
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
    echo " ====> RUN POINT EXTRACTION ... FAILED"
    exit ${exit_code}
fi

echo " ====> RUN POINT EXTRACTION ... DONE"
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
