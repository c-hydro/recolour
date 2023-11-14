#!/bin/bash -e

# Activate virtual environment
PATH=/home/hsaf/recolour/conda//bin:$PATH
export PATH
source activate recolour_libraries

script_folder='/home/hsaf/recolour/tools/algorithm_grid2ts/ecmwf/'
script_file='app_grid2ts_ecmwf.py'
export PYTHONPATH="${PYTHONPATH}:$script_folder"

setting_file='app_grid2ts_ecmwf_ts_dr_h142.json'

time_now=$(date -u +"%Y%m%d%T00")

python $script_file -settings_file $setting_file -time $time_now
