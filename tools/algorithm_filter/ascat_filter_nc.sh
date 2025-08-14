#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - ASCAT H122 FILTER [NETCDF 2 DOMAIN]'
script_version="1.0.0"
script_date='2025/08/14'

# === CONFIG ===
CONDA_ENV_NAME="/root/envs/conda_default/envs/recolour_default"
CONDA_BIN="/root/envs/conda_default/miniconda/bin/conda"
SCRIPT_PATH="/root/package/package_recolour/tools/algorithm_filters/ascat_filter_nc.py"
CONFIG_PATH="/root/projects/hsaf/algorithm_filters/ascat_filter_nc_realtime.json"

# === DATE RANGE (CLI args override env vars, leave empty for real-time mode) ===
START_DATE="${1:-${START_DATE:-}}"  # 1st CLI arg or env var or empty
END_DATE="${2:-${END_DATE:-}}"      # 2nd CLI arg or env var or empty
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Cycle(s) over date(s)
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."

# --- Run in history mode if both START_DATE and END_DATE are set ---
if [[ -n "$START_DATE" && -n "$END_DATE" ]]; then

	echo " ===> HISTORY MODE: $START_DATE → $END_DATE (exclusive) ... "

	current_date="$START_DATE"
	while [[ "$current_date" < "$END_DATE" ]]; do
		next_date=$(date -I -d "$current_date + 1 day")

		START_TIME="${current_date} 00:00"
		END_TIME="${next_date} 00:00"

		echo " ====> PROCESSING FROM $START_TIME TO $END_TIME ... "
		"$CONDA_BIN" run -p "$CONDA_ENV_NAME" python "$SCRIPT_PATH" \
		  -settings_file "$CONFIG_PATH"

		current_date="$next_date"

		echo " ====> PROCESSING FROM $START_TIME TO $END_TIME ... DONE"
	done

	echo " ===> HISTORY MODE: $START_DATE → $END_DATE (exclusive) ... DONE"

# --- Otherwise run in real-time mode ---
else

	echo " ===> REAL-TIME MODE ... "
	TIME_NOW=$(date -u +"%Y-%m-%d %H:%M")
	echo " ====> CURRENT UTC TIME: $TIME_NOW"

	# Optional: set process time to one cycle before now
	PROCESS_TIME=$(date -u -d "$TIME_NOW - 1 hour" +"%Y-%m-%d %H:%M")
	echo " ====> PROCESSING UP TO $PROCESS_TIME ... "

	"$CONDA_BIN" run -p "$CONDA_ENV_NAME" python "$SCRIPT_PATH" \
		-settings_file "$CONFIG_PATH" \
		-time "$PROCESS_TIME"

	echo " ====> PROCESSING UP TO $PROCESS_TIME ... DONE"

	echo " ===> REAL-TIME MODE ... DONE "
fi

# Final message
echo -e "\n ==> All dates are completed."
echo " " 

echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

