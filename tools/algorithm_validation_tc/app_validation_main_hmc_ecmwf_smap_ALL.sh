#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR RUNNER - APP VALIDATION - ECMWF/HMC/SMAP DR - HISTORY'
script_version="1.6.0"
script_date='2023/07/26'

# Script settings
script_folder_app='/home/hsaf/recolour/tools/algorithm_validation/'
script_file_app='app_validation_main.py'
script_folder_settings='/home/hsaf/recolour/tools/algorithm_validation/'

script_file_settings_list=(
#'app_validation_italy_hmc_ecmwf_smap_1357.json'
#'app_validation_italy_hmc_ecmwf_smap_1358.json'
#'app_validation_italy_hmc_ecmwf_smap_1359.json'
#'app_validation_italy_hmc_ecmwf_smap_1393.json'
#'app_validation_italy_hmc_ecmwf_smap_1394.json'
#'app_validation_italy_hmc_ecmwf_smap_1395.json'
'app_validation_italy_hmc_ecmwf_smap_1429.json'
'app_validation_italy_hmc_ecmwf_smap_1430.json'
'app_validation_italy_hmc_ecmwf_smap_1431.json'
)

# Time settings
time_now=$(date -u +"%Y%m%d%H00")

# Venv settings
virtual_env_folder='/home/hsaf/recolour/conda/bin/'
virtual_env_name='recolour_libraries'
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
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."
echo " ==> Remember to have a look at the flags for saving/cleaning ancillary files ..."


for script_file_settings in "${script_file_settings_list[@]}"; do

    # Set settings file
    script_path_settings=${script_folder_settings}${script_file_settings}

    echo " ==> COMMAND LINE: python "$script_path_app" -settings_file "$script_path_settings""

    python $script_path_app -settings_file $script_path_settings

    echo " ==> ... DONE! "

done

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------


