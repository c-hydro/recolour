#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - DOWNLOADER - SOIL MOISTURE SMAP L2 - NRT - REALTIME'
script_version="1.0.0"
script_date='2023/09/20'

# Script settings
script_folder_app='/home/idrologia/package/package_recolour/tools/algorithm_downloader/smap/l2/'
script_file_app='smap_downloader_spl2smp_e.py'
script_folder_settings='/share/home/idrologia/project/sm_tc/tools/algorithm_downloader/smap/l2/'
script_file_settings='smap_downloader_spl2smp_e_nrt_realtime.json'

# Venv settings
virtual_env_folder='/home/idrologia/library/conda_recolour/bin/'
virtual_env_name='recolour_runner_libraries'

# Time settings (-u to get gmt time)
time=$(date +"%Y-%m-%d %H:00")
#time="2023-10-11 23:00"

# # Time period execution (-u to get gmt time)
time_period_days=2
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

# Iterate over hours
time_run=$(date -d "$time" +'%Y-%m-%d %H:00')
for time_period_step in $(seq 0 $time_period_days); do
    
    # Parse time information
    time_step=$(date -d "$time_run ${time_period_step} days ago" +'%Y-%m-%d %H:00')

	# Run python script (using setting and time)
	echo -n " ===> COMMAND LINE: " python $script_path_app -settings_file $script_path_settings -time "$time_step" 

	python $script_path_app -settings_file $script_path_settings -time "$time_step"
	
	echo " ... DONE!"

done

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------


