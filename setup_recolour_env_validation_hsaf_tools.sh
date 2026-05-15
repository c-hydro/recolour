#!/bin/bash

# =============================================================================
# Setup environment for validation HSAF
# =============================================================================

ENV_NAME="env_validation_tools"

echo "Creating virtual environment: ${ENV_NAME}"

# check venv support
python3 -m venv --help >/dev/null 2>&1

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: python3-venv is missing"
    echo ""
    echo "Install with:"
    echo "sudo apt install python3-venv"
    echo ""
    exit 1
fi

# create environment
python3 -m venv "${ENV_NAME}"

# stop if creation failed
if [ ! -f "${ENV_NAME}/bin/activate" ]; then
    echo ""
    echo "ERROR: virtual environment creation failed"
    exit 1
fi

# activate environment
source "${ENV_NAME}/bin/activate"

# upgrade pip stack
python -m pip install --upgrade pip setuptools wheel

# install packages
python -m pip install \
    ascat \
    netCDF4 \
    xarray \
    pandas \
    numpy \
    scipy \
    h5py \
    rasterio \
    memon \
    pygeobase \
    psutil \
    progressbar

# =============================================================================
# Create loader script
# =============================================================================

cat << EOF > run_recolour_system_env_validation_tools.sh
#!/bin/bash

ENV_FOLDER="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )/${ENV_NAME}"

source "\${ENV_FOLDER}/bin/activate"

echo "Validation HSAF environment loaded"
echo "Python: \$(which python)"
echo "Pip:    \$(which pip)"
EOF

chmod +x run_recolour_system_env_validation_tools.sh

echo ""
echo "Environment created successfully"
echo ""
echo "To load the environment use:"
echo "source run_recolour_system_env_validation_tools.sh"
