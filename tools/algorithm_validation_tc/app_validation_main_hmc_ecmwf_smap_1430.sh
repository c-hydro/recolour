#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR RUNNER - APP VALIDATION - HMC/ECMWF/SMAP DR - HISTORY'
script_version="1.6.0"
script_date='2023/07/26'

# Script settings
script_folder_app='/home/hsaf/recolour/tools/algorithm_validation/'
script_file_app='app_validation_main.py'
script_folder_settings='/home/hsaf/recolour/tools/algorithm_validation/'

script_file_settings='app_validation_italy_hmc_ecmwf_smap_1430.json'

# Venv settings
virtual_env_folder='/home/hsaf/recolour/conda/bin/'
virtual_env_name='recolour_libraries'
# Time settings
time_now=$(date -u +"%Y%m%d%H00")
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Add script folder to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder_app"
# Add Venv path
export PATH=$virtual_env_folder:$PATH
source activate $virtual_env_name
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Script path(s)
script_path_app=${script_folder_app}${script_file_app}
script_path_settings=${script_folder_settings}${script_file_settings}
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."
echo " ==> Remember to have a look at the flags for saving/cleaning ancillary files ..."
echo " ==> COMMAND LINE: " python $script_path_app -script_path_settings $setting_file -time $time_now

python $script_path_app -settings_file $script_path_settings # -time $time_now

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------


