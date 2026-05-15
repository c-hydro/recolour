#!/bin/bash

################################################################################
# Script Name: run_validation_hsaf.sh
# Description: Run HSAF soil moisture validation workflow
################################################################################

set -euo pipefail

echo "=============================================================================="
echo ">>> START HSAF VALIDATION WORKFLOW"
echo "=============================================================================="

# --------------------------------------------------------------------------------
# USER VARIABLES
# --------------------------------------------------------------------------------

# Environment
ENV_FOLDER="/root/envs/env_validation_hsaf"
# Python executable
PYTHON_BIN="/root/envs/env_validation_hsaf/bin/python"

# Main application
APP_FILE="/root/projects/validation_sm_hsaf/algorithm_runner/app_validation_main.py"
# Settings file
SETTINGS_FILE="/root/projects/validation_sm_hsaf/algorithm_runner/2026_validation_global_ascat_gldas_cci_h16_v1.json"

# Log folder
LOG_FOLDER="/root/projects/validation_sm_hsaf/log/"
# Log filename
LOG_FILE="${LOG_FOLDER}/validation_$(date +"%Y%m%d_%H%M%S").log"

# --------------------------------------------------------------------------------
# CHECK ENVIRONMENT
# --------------------------------------------------------------------------------

echo ""
echo ">>> Checking environment..."

if [ ! -d "${ENV_FOLDER}" ]; then
    echo "ERROR: Environment folder not found!"
    echo "${ENV_FOLDER}"
    exit 1
fi

echo ">>> Activating environment..."
source "${ENV_FOLDER}/bin/activate"

echo ">>> Environment loaded"
echo ">>> Python : $(which python)"
echo ">>> Pip    : $(which pip)"

# --------------------------------------------------------------------------------
# CHECK FILES
# --------------------------------------------------------------------------------

echo ""
echo ">>> Checking application files..."

if [ ! -f "${APP_FILE}" ]; then
    echo "ERROR: Application file not found!"
    echo "${APP_FILE}"
    exit 1
fi

if [ ! -f "${SETTINGS_FILE}" ]; then
    echo "ERROR: Settings file not found!"
    echo "${SETTINGS_FILE}"
    exit 1
fi

echo ">>> Application file OK"
echo ">>> Settings file OK"

# --------------------------------------------------------------------------------
# CREATE LOG FOLDER
# --------------------------------------------------------------------------------

mkdir -p "${LOG_FOLDER}"

echo ""
echo ">>> Log folder:"
echo "${LOG_FOLDER}"

echo ">>> Log file:"
echo "${LOG_FILE}"

# --------------------------------------------------------------------------------
# PRINT CONFIGURATION
# --------------------------------------------------------------------------------

echo ""
echo "=============================================================================="
echo ">>> CONFIGURATION"
echo "=============================================================================="

echo "ENV_FOLDER   : ${ENV_FOLDER}"
echo "PYTHON_BIN   : ${PYTHON_BIN}"
echo "APP_FILE     : ${APP_FILE}"
echo "SETTINGS_FILE: ${SETTINGS_FILE}"
echo "LOG_FILE     : ${LOG_FILE}"

echo "=============================================================================="

# --------------------------------------------------------------------------------
# RUN APPLICATION
# --------------------------------------------------------------------------------

echo ""
echo ">>> Running validation workflow..."
echo ""

START_TIME=$(date)

echo ">>> Start time: ${START_TIME}"

COMMAND="${PYTHON_BIN} ${APP_FILE} -settings_file ${SETTINGS_FILE}"

echo ""
echo ">>> Command:"
echo "${COMMAND}"
echo ""

${PYTHON_BIN} "${APP_FILE}" \
    -settings_file "${SETTINGS_FILE}" \
    2>&1 | tee "${LOG_FILE}"

END_TIME=$(date)

# --------------------------------------------------------------------------------
# END
# --------------------------------------------------------------------------------

echo ""
echo "=============================================================================="
echo ">>> VALIDATION COMPLETED"
echo "=============================================================================="

echo ">>> Start time : ${START_TIME}"
echo ">>> End time   : ${END_TIME}"
echo ">>> Log file   : ${LOG_FILE}"

echo "=============================================================================="
