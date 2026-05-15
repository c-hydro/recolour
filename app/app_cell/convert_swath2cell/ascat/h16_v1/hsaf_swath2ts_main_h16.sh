#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='HSAF - Swath2TS Runner - H16'
script_version="1.1.0"
script_date='2026/05/10'

# Python virtual environment information
env_libs='run_recolour_system_env_validation_tools.sh'
source ${env_libs}
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get file information
script_file_1='/root/projects/validation_sm_hsaf/algorithm_converter/ascat_2018/h16_swath2idx_ts_main.py'
script_file_2='/root/projects/validation_sm_hsaf/algorithm_converter/ascat_2018/h16_idx2cont_ts_main.py'
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."

# Run python script (using setting and time)
python $script_file_1
python $script_file_2

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------
