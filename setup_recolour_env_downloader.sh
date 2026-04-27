#!/bin/bash -e

# -----------------------------------------------------------------------------------------
# Script information
script_name='RECOLOUR ENVIRONMENT - DOWNLOADER - CONDA'
script_version="2.1.0"
script_date='2026/04/23'

# Modern Miniconda installer
fp_env_file_miniconda='https://repo.anaconda.com/miniconda/Miniconda3-py311_26.1.1-1-Linux-x86_64.sh'

# Argument(s) default definition(s)
fp_env_tag_default='recolour_downloader'
fp_env_folder_root_default='./conda_p311/'
fp_env_file_reference_default='%ENV_TAG_settings'
fp_env_folder_libraries_default='%ENV_TAG_libraries'
fp_env_file_requirements_default='requirements_%ENV_TAG.yaml'

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo ""

script_args_n=$#
script_args_values=$@

echo " ==> Script arguments number: ${script_args_n}"
echo " ==> Script arguments values: ${script_args_values}"
echo ""

# -----------------------------------------------------------------------------------------
# Get arguments
if [ $# -eq 0 ]; then
    fp_env_tag=$fp_env_tag_default
    fp_env_folder_root=$fp_env_folder_root_default
    fp_env_file_reference=$fp_env_file_reference_default
    fp_env_folder_libraries=$fp_env_folder_libraries_default
    fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 1 ]; then
    fp_env_tag=$1
    fp_env_folder_root=$fp_env_folder_root_default
    fp_env_file_reference=$fp_env_file_reference_default
    fp_env_folder_libraries=$fp_env_folder_libraries_default
    fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 2 ]; then
    fp_env_tag=$1
    fp_env_folder_root=$2
    fp_env_file_reference=$fp_env_file_reference_default
    fp_env_folder_libraries=$fp_env_folder_libraries_default
    fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 3 ]; then
    fp_env_tag=$1
    fp_env_folder_root=$2
    fp_env_file_reference=$3
    fp_env_folder_libraries=$fp_env_folder_libraries_default
    fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 4 ]; then
    fp_env_tag=$1
    fp_env_folder_root=$2
    fp_env_file_reference=$3
    fp_env_folder_libraries=$4
    fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 5 ]; then
    fp_env_tag=$1
    fp_env_folder_root=$2
    fp_env_file_reference=$3
    fp_env_folder_libraries=$4
    fp_env_file_requirements=$5
else
    echo " ==> ERROR: too many arguments"
    exit 1
fi

# Replace tag placeholders
fp_env_folder_root=${fp_env_folder_root/'%ENV_TAG'/$fp_env_tag}
fp_env_file_reference=${fp_env_file_reference/'%ENV_TAG'/$fp_env_tag}
fp_env_folder_libraries=${fp_env_folder_libraries/'%ENV_TAG'/$fp_env_tag}
fp_env_file_requirements=${fp_env_file_requirements/'%ENV_TAG'/$fp_env_tag}

echo " ==> ENV TAG: ${fp_env_tag}"
echo " ==> ROOT: ${fp_env_folder_root}"
echo " ==> SETTINGS FILE: ${fp_env_file_reference}"
echo " ==> ENV NAME: ${fp_env_folder_libraries}"
echo " ==> YAML FILE: ${fp_env_file_requirements}"
echo ""

# -----------------------------------------------------------------------------------------
# Create root folder
if [ ! -d "$fp_env_folder_root" ]; then
    mkdir -p "$fp_env_folder_root"
fi

# -----------------------------------------------------------------------------------------
# Check local conda installation
echo " ====> CHECK PYTHON ENVIRONMENT ... "
if [ -d "${fp_env_folder_root}/bin" ]; then
    echo " ====> CHECK PYTHON ENVIRONMENT ... FOUND."
    fp_env_install=false
else
    echo " ====> CHECK PYTHON ENVIRONMENT ... NOT FOUND."
    fp_env_install=true
fi

# -----------------------------------------------------------------------------------------
# Install local Miniconda
echo " ====> INSTALL PYTHON ENVIRONMENT ... "
if $fp_env_install; then
    wget "$fp_env_file_miniconda" -O miniconda.sh

    if [ -d "$fp_env_folder_root" ]; then
        rm -rf "$fp_env_folder_root"
    fi

    bash miniconda.sh -b -p "$fp_env_folder_root"
    rm -f miniconda.sh
    echo " ====> INSTALL PYTHON ENVIRONMENT ... DONE!"
else
    echo " ====> INSTALL PYTHON ENVIRONMENT ... DONE. PREVIOUSLY INSTALLED"
fi

# -----------------------------------------------------------------------------------------
# Activate local conda base
export PATH="${fp_env_folder_root}/bin:${PATH}"
source "${fp_env_folder_root}/bin/activate"

# -----------------------------------------------------------------------------------------
# Configure conda
echo " ====> CONFIGURE CONDA ... "
conda config --set channel_priority strict
conda config --remove-key channels > /dev/null 2>&1 || true
conda config --add channels conda-forge
echo " ====> CONFIGURE CONDA ... DONE!"

# -----------------------------------------------------------------------------------------
# Create environment from YAML
echo " ====> INSTALL PYTHON LIBRARIES ... "

if [ ! -f "$fp_env_file_requirements" ] ; then
    echo " =====> ERROR: YAML requirements file not found: ${fp_env_file_requirements}"
    exit 1
fi

if conda env list | awk '{print $1}' | grep -Fxq "${fp_env_folder_libraries}"; then
    echo " =====> REMOVE PREVIOUS ENVIRONMENT: ${fp_env_folder_libraries}"
    conda env remove --yes --name "${fp_env_folder_libraries}"
fi

echo " =====> CREATE ENVIRONMENT USING CONDA ... "
conda env create --name "${fp_env_folder_libraries}" --file "${fp_env_file_requirements}"
echo " =====> CREATE ENVIRONMENT USING CONDA ... DONE"

echo " =====> ACTIVATE ENVIRONMENT: ${fp_env_folder_libraries}"
source activate "${fp_env_folder_libraries}"

echo " =====> INSTALL EXTRA PACKAGES ..."
python -m pip install -U "pip>=25.3" wheel
printf "setuptools<81\n" > build-constraints.txt
python -m pip install --build-constraint build-constraints.txt gldas
python -m pip install -U earthaccess
echo " =====> INSTALL EXTRA PACKAGES ... DONE"

# -----------------------------------------------------------------------------------------
# Checks
echo " =====> PYTHON VERSION ..."
python --version

echo " =====> RUN ECCODES CHECK ... "
python -m eccodes selfcheck

echo " =====> RUN CFGRIB CHECK ... "
python -m cfgrib selfcheck

echo " =====> RUN PYGRIB CHECK ... "
python -c "import pygrib; print('pygrib import OK')"

echo " =====> RUN XARRAY CHECK ... "
python -c "import xarray as xr; print('xarray version:', xr.__version__)"

echo " =====> RUN RASTERIO CHECK ... "
python -c "import rasterio; print('rasterio version:', rasterio.__version__)"

echo " =====> RUN GLDAS CHECK ... "
python -c "import gldas; print('gldas import OK')"

echo " =====> RUN EARTHACCESS CHECK ... "
python -c "import earthaccess; print('earthaccess version:', earthaccess.__version__)"

echo " ====> INSTALL PYTHON LIBRARIES ... DONE!"

# -----------------------------------------------------------------------------------------
# Create environmental file
echo " ====> CREATE ENVIRONMENTAL FILE ... "

cd "$fp_env_folder_root"

if [ -f "$fp_env_file_reference" ] ; then
    rm "$fp_env_file_reference"
fi

echo "PATH=${fp_env_folder_root}/bin:"'$PATH' >> "$fp_env_file_reference"
echo "export PATH" >> "$fp_env_file_reference"
echo "source ${fp_env_folder_root}/bin/activate" >> "$fp_env_file_reference"
echo "source activate ${fp_env_folder_libraries}" >> "$fp_env_file_reference"

echo " ====> CREATE ENVIRONMENTAL FILE ... DONE!"

# -----------------------------------------------------------------------------------------
# Info script end
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
