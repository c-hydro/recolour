#!/bin/bash

# =============================================================================
# Setup environment
# =============================================================================

ENV_NAME="${1:-env_analyzer_soil_tools}"
LOADER_NAME="${2:-env_analyzer_soil_tools.sh}"

echo "Creating virtual environment: ${ENV_NAME}"

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

python3 -m venv "${ENV_NAME}"

if [ ! -f "${ENV_NAME}/bin/activate" ]; then
    echo ""
    echo "ERROR: virtual environment creation failed"
    exit 1
fi

source "${ENV_NAME}/bin/activate"

python -m pip install --upgrade pip setuptools wheel

python -m pip install \
    matplotlib \
    netCDF4 \
    xarray \
    pandas \
    numpy \
    scipy \
    rasterio \
    geopandas \
    shapely \
    owslib

# =============================================================================
# Create loader script
# =============================================================================

cat << EOF > "${LOADER_NAME}"
#!/bin/bash

ENV_FOLDER="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )/${ENV_NAME}"

source "\${ENV_FOLDER}/bin/activate"

echo "Validation HSAF environment loaded"
echo "Python: \$(which python)"
echo "Pip:    \$(which pip)"
EOF

chmod +x "${LOADER_NAME}"

echo ""
echo "Environment created successfully"
echo ""
echo "To load the environment use:"
echo "source ${LOADER_NAME}"
