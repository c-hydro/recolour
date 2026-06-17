#!/bin/bash

CONDA_DIR="./conda_soil/"
ENV_NAME="analyzer_soil_tools"

source "${CONDA_DIR}/etc/profile.d/conda.sh"

conda activate "${ENV_NAME}"

echo "Analyzer Soil Tools environment loaded"
echo "Conda Env: ${ENV_NAME}"
echo "Python:    $(which python)"
echo "Pip:       $(which pip)"
python --version
