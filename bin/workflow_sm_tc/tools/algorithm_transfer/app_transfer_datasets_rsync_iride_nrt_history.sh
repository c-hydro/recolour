#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR - TRANSFER DATASETS - IRIDE - NRT HISTORY'
script_version="1.0.0"
script_date='2021/11/18'

# Virtualenv default definition(s)
virtualenv_folder='/home/idrologia/library/conda_recolour/bin/'
virtualenv_name='recolour_runner_libraries'

# Default script folder(s)
script_folder='/home/idrologia/package/package_recolour/tools/algorithm_transfer/'
configuration_folder='/home/idrologia/project/sm_tc/tools/algorithm_transfer/'
package_folder='/home/idrologia/package/package_recolour/'

# Execution example:
# python3 hyde_tools_transfer_datasets.py -settings_algorithm configuration.json -time "2020-11-02 12:00"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get file information
script_file=${script_folder}'app_transfer_datasets.py'
settings_file=${configuration_folder}'app_transfer_datasets_rsync_iride_nrt_history.json'

# Get information (-u to get gmt time)
time_now=$(date -u +"%Y-%m-%d %H:00")
time_now="2024-03-10 12:00"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$virtualenv_folder:$PATH
source activate $virtualenv_name

# Add path to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder"
export PYTHONPATH="${PYTHONPATH}:$package_folder"
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."
echo " ==> COMMAND LINE: " python $script_file -settings_file $settings_file -time $time_now

# Run python script (using setting and time)
python $script_file -settings_file $settings_file -time "$time_now"

# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

