#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR RUNNER - HSAF ASCAT H16 [METOP-SM] - DATA RECORD'
script_version="1.1.0"
script_date='2024/02/18'

# Script settings
script_folder_app='/home/idrologia/package/package_recolour/app/app_cell/convert_swath2cell/'
script_file_app='app_swath2cell_ascat.py'
script_folder_settings='/share/home/idrologia/project/validation_hsaf/algorithm_runner/convert_swath2cell/'
script_file_settings='app_swath2cell_ascat_h16_dr.json'

# Venv settings
virtual_env_folder='/home/idrologia/library/conda_recolour/bin/'
virtual_env_name='recolour_adapter_4_hyde_libraries'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Time settings (-u to get gmt time)
time_now=$(date +"%Y-%m-%d %H:00")
time_now="2024-01-01 00:00" # complete period

# Time period execution
# "dr :: 2016-01-01 -- 2023-09-20 (2820 days)"
time_period_days=1 # fix on 0 to provide actual time only, if providing a time period always add 1
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

# Iterate over days
time_run=$(date -d "$time_now" +'%Y-%m-%d %H:00')
for time_period_step in $(seq 0 $time_period_days); do

    # Parse time information
    time_step=$(date -d "$time_run ${time_period_step} days ago" +'%Y-%m-%d %H:00')

	# Run python script (using setting and time)
	echo -n " ===> COMMAND LINE: python $script_path_app -settings_file $script_path_settings -time "$time_step" "

    python $script_path_app -settings_file $script_path_settings -time "$time_step"

    echo " ... DONE!"

done

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

