#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - COMPUTE METRICS - ITALY - ASCAT:GLDAS:CCI - PARALLEL NRT'
script_version="1.1.0"
script_date='2024/03/19'

# Get file information
script_folder_app='/home/idrologia/package/package_recolour/tools/algorithm_validation_hsaf/'
script_file_app='app_validation_main.py'
script_folder_settings='/home/idrologia/project/validation_hsaf/algorithm_runner/compute_metrics/'
script_file_settings='app_validation_ascat_gldas_cci_italy_nrt_parallel.json'

# Venv settings
virtual_env_folder='/home/idrologia/library/conda_recolour/bin/'
virtual_env_name='recolour_adapter_4_hyde_libraries'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Time settings (-u to get gmt time)
time_now=$(date +"%Y-%m-%d %H:00")
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

# compute time run
time_run=$(date -d "$time_now" +'%Y-%m-%d %H:00')

# run python script (using settings)
echo " ===> COMMAND LINE: python $script_path_app -settings_file $script_path_settings "

python $script_path_app -settings_file $script_path_settings

echo " ===> ... DONE!"

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------



