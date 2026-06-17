#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='SM CNR ENVIRONMENT - PYTHON3 LIBRARIES FOR SM-CNR - CONDA'
script_version="1.0.0"
script_date='2026/06/11'

# Miniconda installer
fp_env_file_miniconda='https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh'

# Default settings
fp_env_tag_default='sm_cnr'

fp_env_folder_root_default='/hydro/library/fp_libs_python_sm_cnr/'
fp_env_file_reference_default='%ENV_TAG_settings'
fp_env_folder_libraries_default='%ENV_TAG_libraries'

fp_env_file_requirements_default='requirements_%ENV_TAG.yaml'
#-----------------------------------------------------------------------------------------

echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."

script_args_n=$#
script_args_values=$@

echo ""
echo " ==> Script arguments number: ${script_args_n}"
echo " ==> Script arguments values: ${script_args_values}"
echo ""

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

fp_env_folder_root=${fp_env_folder_root/'%ENV_TAG'/$fp_env_tag}
fp_env_file_reference=${fp_env_file_reference/'%ENV_TAG'/$fp_env_tag}
fp_env_folder_libraries=${fp_env_folder_libraries/'%ENV_TAG'/$fp_env_tag}
fp_env_file_requirements=${fp_env_file_requirements/'%ENV_TAG'/$fp_env_tag}

echo " ==> ARGS SELECTED:"
echo ""
echo " ==> Env tag: ${fp_env_tag}"
echo " ==> Conda root folder: ${fp_env_folder_root}"
echo " ==> Environment settings file: ${fp_env_file_reference}"
echo " ==> Conda environment name: ${fp_env_folder_libraries}"
echo " ==> Requirements YAML file: ${fp_env_file_requirements}"
echo ""

#-----------------------------------------------------------------------------------------
# Install Miniconda
echo " ====> CHECK MINICONDA ENVIRONMENT ..."

if [ -d "${fp_env_folder_root}/bin" ]; then
    echo " ====> CHECK MINICONDA ENVIRONMENT ... FOUND"
    fp_env_install=false
else
    echo " ====> CHECK MINICONDA ENVIRONMENT ... NOT FOUND"
    fp_env_install=true
fi

if $fp_env_install; then
    echo " ====> INSTALL MINICONDA ENVIRONMENT ..."

    wget "${fp_env_file_miniconda}" -O miniconda.sh

    if [ -d "${fp_env_folder_root}" ]; then
        rm -rf "${fp_env_folder_root}"
    fi

    bash miniconda.sh -b -p "${fp_env_folder_root}"
    rm -f miniconda.sh

    echo " ====> INSTALL MINICONDA ENVIRONMENT ... DONE"
else
    echo " ====> INSTALL MINICONDA ENVIRONMENT ... PREVIOUSLY INSTALLED"
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Init conda
export PATH="${fp_env_folder_root}/bin:$PATH"
eval "$("${fp_env_folder_root}/bin/conda" shell.bash hook)"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Install python libraries
echo " ====> INSTALL PYTHON LIBRARIES ..."

if conda env list | awk '{print $1}' | grep -qx "${fp_env_folder_libraries}"; then
    echo " =====> CONDA ENVIRONMENT '${fp_env_folder_libraries}' ALREADY EXISTS"
else
    if [ -f "${fp_env_file_requirements}" ]; then
        echo " =====> USE OF CONDA REQUIREMENTS FILE YAML: ${fp_env_file_requirements}"
        conda env create --file "${fp_env_file_requirements}"
    else
        echo " =====> REQUIREMENTS YAML NOT FOUND"
        echo " =====> USE OF CONDA GENERIC COMMAND-LINE"

        conda create --yes \
            --name "${fp_env_folder_libraries}" \
            -c conda-forge \
            python=3.11 \
            pandas \
            requests \
            rasterio \
            numpy \
            xarray \
            netcdf4 \
            dask \
            bottleneck \
            pip
    fi
fi

echo " ====> INSTALL PYTHON LIBRARIES ... DONE"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Create environment settings file
echo " ====> CREATE ENVIRONMENTAL FILE ..."

cd "${fp_env_folder_root}"

if [ -f "${fp_env_file_reference}" ]; then
    rm "${fp_env_file_reference}"
fi

{
    echo "export PATH=${fp_env_folder_root}/bin:\$PATH"
    echo "eval \"\$(${fp_env_folder_root}/bin/conda shell.bash hook)\""
    echo "conda activate ${fp_env_folder_libraries}"
} >> "${fp_env_file_reference}"

echo " ====> CREATE ENVIRONMENTAL FILE ... DONE"
echo " ====> ENVIRONMENTAL FILE: ${fp_env_folder_root}/${fp_env_file_reference}"
#-----------------------------------------------------------------------------------------

echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
