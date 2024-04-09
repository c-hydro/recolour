#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - CREATE MAP TC PRODUCT - NRT - HISTORY'
script_version="1.1.0"
script_date='2023/12/19'

# Script settings
script_folder_app='/home/idrologia/package/package_recolour/app/app_map/create_map_tc/'
script_file_app='app_map_grid_tc.py'
script_folder_settings='/home/idrologia/project/sm_tc/app_map/create_map_tc/'
script_file_settings='app_map_grid_tc_hmc_nrt_history.json'

# Venv settings
virtual_env_folder='/home/idrologia/library/conda_recolour/bin/'
virtual_env_name='recolour_runner_libraries'

# Time settings (-u to get gmt time)
time_now=$(date +"%Y-%m-%d %H:00")
time_now="2023-03-2 12:00" # history hmc ref
#time_now="2023-12-21 23:00" 


# Time period execution
# "dr :: 2016-01-01 -- 2023-09-20 (2820 days)"
# "dr :: 2008-01-01 -- 2015-12-31 (2922 days)"
time_period_days=750 # fix on 0 to provide actual time only, if providing a time period always add 1
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


