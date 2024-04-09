#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - DOWNLOADER - SOIL MOISTURE SMAP L3 - NRT - HISTORY (V005)'
script_version="3.1.1"
script_date='2024/01/22'

# Script settings
script_folder_app='/home/hsaf/recolour/tools/algorithm_downloader/smap/l3/'
script_file_app='smap_downloader_spl3smp_e.py'
script_folder_settings='/home/hsaf/recolour/bin/tools/algorithm_downloader/smap/l3/'
script_file_settings='smap_downloader_spl3smp_e_nrt_history_v005.json'

# Venv settings
virtual_env_folder='/home/hsaf/recolour/conda/bin/'
virtual_env_name='recolour_libraries'

# Time period execution
time_period_hour=0 # fixed on zero to provide actual time only

# Time settings (-u to get gmt time)
time_now=$(date +"%Y-%m-%d %H:00")
time_now="2023-11-05 00:00" # from 2015-03-31 to 2023-12-02
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
time_run=$(date -d "$time_now" +'%Y-%m-%d %H:00')
for time_period_step in $(seq 0 $time_period_hour); do
    
    # Parse time information
    time_step=$(date -d "$time_run ${time_period_step} hour ago" +'%Y-%m-%d %H:00')

	# Run python script (using setting and time)
	echo -n " ===> COMMAND LINE: "python $script_path_app -settings_file $script_path_settings -time "$time_step" 

	python $script_path_app -settings_file $script_path_settings -time "$time_step"
	
	echo " ... DONE!"

done

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------


