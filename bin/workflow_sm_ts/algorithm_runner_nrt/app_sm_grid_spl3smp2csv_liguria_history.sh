#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - CREATE SPL3SMP GRID2CSV TS - NRT - HISTORY'
script_version="1.2.0"
script_date='2024/02/29'

# Script settings
script_folder_app='/share/home/idrologia/package/package_recolour/app/app_ts/convert_grid_smap2csv/'
script_file_app='app_sm_grid_smap2csv.py'
script_folder_settings='/share/home/idrologia/project/sm_ts/algorithm_runner/'
script_file_settings='app_sm_grid_spl3smp2csv_liguria_history.json'

# Venv settings
virtual_env_folder='/home/idrologia/library/conda_recolour/bin/'
virtual_env_name='recolour_runner_libraries'

# Time settings (-u to get gmt time)
time_now=$(date +"%Y-%m-%d %H:00")
time_now="2024-01-01 00:00" # history case

# Time period execution
time_period_months=48 
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
for time_period_step in $(seq 0 $time_period_months); do

    # Parse time information
    time_step=$(date -d "$time_run ${time_period_step} month ago" +'%Y-%m-%d %H:00')

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

