#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - TIME SERIES ECMWF LAYER 0-28 CM - DATA RECORD - REALTIME'
script_version="1.2.0"
script_date='2024/02/29'

# Script settings
script_folder_app='/share/home/idrologia/package/package_recolour/app/app_ts/join_ts_periods/'
script_file_app='app_sm_ts_join_periods.py'
script_folder_settings='/share/home/idrologia/project/sm_ts/algorithm_runner_dr/'
script_file_settings='app_sm_ts_join_periods_ecmwf_layer_0_28_cm_liguria_realtime.json'

# Venv settings
virtual_env_folder='/home/idrologia/library/conda_recolour/bin/'
virtual_env_name='recolour_runner_libraries'

# Time settings (-u to get gmt time)
time_now=$(date +"%Y-%m-%d %H:00")
#time_now="2024-02-29 00:00" # history case
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

# Run time
time_run=$(date -d "$time_now" +'%Y-%m-%d %H:00')

# Run python script (using setting and time)
echo -n " ===> COMMAND LINE: python $script_path_app -settings_file $script_path_settings -time "$time_run" "
python $script_path_app -settings_file $script_path_settings -time "$time_run"

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

