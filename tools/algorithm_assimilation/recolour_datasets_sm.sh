#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='DTE - PREPROCESSING SM'
script_version="1.0.0"
script_date='2021/05/14'

virtualenv_folder='/home/fabio/fp_virtualenv_python3/'
virtualenv_name='fp_virtualenv_python3_hyde_libraries'
script_folder='/home/fabio/Desktop/PyCharm_Workspace/project/dte/app_sm/'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get script information
script_file='dte_preprocessing_datasets_sm.py'
settings_file='dte_preprocessing_datasets_sm.json'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$virtualenv_folder/bin:$PATH
source activate $virtualenv_name

# Add path to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder"
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> COMMAND LINE: " python3 $script_file -settings_file $settings_file
echo " ===> EXECUTION ..."

# Execution pid
execution_pid=$$

# Run python script (using setting and time)
python3 $script_file -settings_file $settings_file


echo " ===> EXECUTION ... DONE"
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------	

