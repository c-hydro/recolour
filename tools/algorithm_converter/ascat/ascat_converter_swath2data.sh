#!/bin/bash

# ----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - ASCAT CONVERTER [SWATH 2 TIFF/CSV]'
script_version="1.0.0"
script_date='2025/08/09'

# configure algorithm
CONDA_ENV_NAME="/root/envs/conda_default/envs/recolour_default"
CONDA_ENV_PATH='/root/envs/conda_default/miniconda/bin/'
SCRIPT_PATH="/root/package/package_recolour/tools/algorithm_converter/ascat/ascat_converter_bufr2data.py"
CONFIG_PATH="/root/projects/hsaf/algorithm_converter/ascat_converter_bufr2data_h16_history.json"

# configure period
START_DATE="2020-07-01"
END_DATE="2020-08-01"  # exclusive upper bound (won’t run this date)
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# activate conda env
source $CONDA_ENV_PATH/activate $CONDA_ENV_NAME

# initialize date
current_date="$START_DATE"

# count total days
total_days=$(( ( $(date -d "$END_DATE" +%s) - $(date -d "$START_DATE" +%s) ) / 86400 ))
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Cycle(s) over date(s)
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."

# === LOOP OVER EACH DAY ===
day_counter=0
while [[ "$current_date" < "$END_DATE" ]]; do

	# info time start
	yyyymmdd=$(date -d "$current_date" +"%Y%m%d")
	echo " ===> PROCESS DATE: $yyyymmdd ... "

	# run the Python script for the current date
	python "$SCRIPT_PATH" "$CONFIG_PATH" "$yyyymmdd"

	# Increment the day counter
	day_counter=$((day_counter + 1))

	# progress bas
	progress=$(( (day_counter * 100) / total_days ))
	bar=$(printf "%${progress}s" | tr ' ' '#')
	printf " :::: PROCESS PROGRESS: [%-100s] %d%% Complete\n" "$bar" "$progress"

	# increment current_date by one day
	current_date=$(date -I -d "$current_date + 1 day")
	
	# info time end
	echo " ===> PROCESS DATE: $yyyymmdd ... DONE"
  
done

# Final message
echo -e "\n ==> ALL DATES ARE COMPLETED."
echo " " 

echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

