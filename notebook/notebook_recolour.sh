#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR NOTEBOOK - EXECUTION'
script_version="1.6.4"
script_date='2023/10/20'

# Argument(s) script definition
script_folder_root='.'
script_folder_library_base='.'
script_folder_library_extended='.'

# Argument(s) env definition(s)
fp_env_tag='recolour_notebook'

fp_env_folder_root='../conda/'
fp_env_folder_libraries='%ENV_TAG_libraries'

# Command-lines:
cmd_runner_jupyter="jupyter-lab"

# Jupyter lab interactive graph:
# https://towardsdatascience.com/how-to-produce-interactive-matplotlib-plots-in-jupyter-environment-1e4329d71651
# NodeJS=12
# https://github.com/jupyterlab/jupyterlab/issues/7526
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Set the arguments
fp_env_folder_root=${fp_env_folder_root/'%ENV_TAG'/$fp_env_tag}
fp_env_folder_libraries=${fp_env_folder_libraries/'%ENV_TAG'/$fp_env_tag}

echo ""
echo " ==> VIRTUAL ENVIRONMENT SELECTED:"
echo ""
echo " ==> Tag of virtual environmenf [string: path]-> ${fp_env_tag}"
echo " ==> Directory of libraries [string: path]-> ${fp_env_folder_root}"
echo " ==> Name of virtual environment [string: name] -> ${fp_env_folder_libraries}"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$fp_env_folder_root/bin:$PATH
source activate $fp_env_folder_libraries

# Add path to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder_root"
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."
echo " ===> EXECUTION ..."

time_now=$(date -d "$time_now" +'%Y-%m-%d %H:00')
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Check jupiter-nb env installation
echo " ====> CHECK JUPITER-NB ENV INSTALLATION ..."
if ! type "$cmd_runner_jupyter" > /dev/null; then
    echo " ====> CHECK JUPITER-NB ENV INSTALLATION ... FAILED. INSTALL THE JUPITER-NB IN THE SELECTED VIRTUAL ENVIRONMENT. EXIT WITH ERROR(S)"
    exit 1
else
    echo " ====> CHECK JUPITER-NB ENV INSTALLATION ... DONE"
fi
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Check jupiter-nb base library installation
echo " ====> CHECK JUPITER-NB PACKAGE BASE LIBRARY FOLDER ..."

if [ -d "$script_folder_library_base" ]; then
  echo " ====> CHECK JUPITER-NB BASE LIBRARY FOLDER ... DONE"
else
  echo " ====> CHECK JUPITER-NB BASE LIBRARY FOLDER ... FAILED IN $script_folder_library_base. EXIT WITH ERROR(S)"
  exit 2
fi
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Check jupyter-nb package library installation
echo " ====> CHECK JUPITER-NB PACKAGE EXTENDED LIBRARY FOLDER ..."

if [ -d "$script_folder_library_extended" ]; then
  echo " ====> CHECK JUPITER-NB PACKAGE EXTENDED LIBRARY FOLDER ... DONE"
else
  echo " ====> CHECK JUPYTER-NB PACKAGE EXTENDED LIBRARY FOLDER ... NOT AVAILABLE"
fi
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Run jupiter-nb
echo " ====> RUN JUPYTER-NB ... "
${cmd_runner_jupyter}
echo " ====> RUN JUPYTER-NB ... DONE"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script end
echo " ===> EXECUTION ... DONE"
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------



