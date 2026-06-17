#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - COPERNICUS SWI DOWNLOADER - RUNNER'
script_version="1.0.0"
script_date='2026/06/11'

# Environment file
fp_env_file="/hydro/library/fp_libs_python_sm_cnr/sm_cnr_settings"

# Python script and settings
fp_script="/hydro/library/fp_package_recolour/tools/algorithm_downloader/copernicus/copernicus_downloader_swi.py"
fp_settings="/hydro/library/fp_package_recolour/tools/algorithm_downloader/copernicus/copernicus_downloader_swi.json"

# Machine/reference time
time_now=${1:-$(date "+%Y-%m-%d %H:%M")}

# get user and password from environment
export CDSE_USERNAME="${CDSE_USERNAME:-fabio.delogu@cimafoundation.org}"

if [ -z "${CDSE_PASSWORD}" ]; then
    echo " ====> ERROR: CDSE_PASSWORD is not set"
    echo " ====> Set it with: export CDSE_PASSWORD='your_password'"
    exit 1
fi
#-----------------------------------------------------------------------------------------

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."

echo " ====> REFERENCE TIME: ${time_now}"

echo " ====> LOAD PYTHON ENVIRONMENT ..."

if [ -f "${fp_env_file}" ]; then
    source "${fp_env_file}"
else
    echo " ====> LOAD PYTHON ENVIRONMENT ... FAILED"
    echo " ====> ENVIRONMENT FILE NOT FOUND: ${fp_env_file}"
    exit 1
fi

echo " ====> LOAD PYTHON ENVIRONMENT ... DONE"

echo " ====> RUN COPERNICUS SWI DOWNLOADER ..."
echo " ====> PYTHON SCRIPT: ${fp_script}"
echo " ====> SETTINGS FILE: ${fp_settings}"

python "${fp_script}" \
    -settings_file "${fp_settings}" \
    -time "${time_now}" \
    "${@:2}"

echo " ====> RUN COPERNICUS SWI DOWNLOADER ... DONE"

echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
