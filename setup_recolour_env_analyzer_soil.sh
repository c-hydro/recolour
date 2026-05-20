#!/bin/bash

set -e

# =============================================================================
# Setup conda environment
# =============================================================================

ENV_NAME="${1:-analyzer_soil_tools}"
LOADER_NAME="${2:-env_analyzer_soil_tools.sh}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
CONDA_DIR="${CONDA_DIR:-./envs/}"

echo "========================================="
echo " Installing Conda Environment"
echo "========================================="
echo ""
echo "Environment Name : ${ENV_NAME}"
echo "Python Version   : ${PYTHON_VERSION}"
echo "Conda Directory  : ${CONDA_DIR}"
echo ""

# =============================================================================
# Install Miniconda if missing
# =============================================================================

if [ ! -d "${CONDA_DIR}" ]; then

    echo "Miniconda not found"
    echo "Installing Miniconda..."

    TMP_INSTALLER="/tmp/miniconda_installer.sh"

    wget -O "${TMP_INSTALLER}" \
        https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    bash "${TMP_INSTALLER}" -b -p "${CONDA_DIR}"

    rm -f "${TMP_INSTALLER}"

    echo "Miniconda installed successfully"

else

    echo "Miniconda already installed"

fi

# =============================================================================
# Load conda
# =============================================================================

if [ ! -f "${CONDA_DIR}/etc/profile.d/conda.sh" ]; then
    echo ""
    echo "ERROR: conda.sh not found"
    exit 1
fi

source "${CONDA_DIR}/etc/profile.d/conda.sh"

# =============================================================================
# Create environment
# =============================================================================

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then

    echo ""
    echo "Conda environment already exists: ${ENV_NAME}"

else

    echo ""
    echo "Creating conda environment: ${ENV_NAME}"

    conda create -y \
        -n "${ENV_NAME}" \
        python="${PYTHON_VERSION}"

fi

# =============================================================================
# Activate environment
# =============================================================================

conda activate "${ENV_NAME}"

echo ""
echo "Activated environment: ${ENV_NAME}"

PY_VER="$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")')"

echo "Python version: ${PY_VER}"

# =============================================================================
# Install packages
# =============================================================================

echo ""
echo "Installing packages from conda-forge..."

conda install -y -c conda-forge \
    matplotlib \
    netcdf4 \
    xarray \
    pandas \
    numpy \
    scipy \
    pyproj \
    gdal \
    rasterio \
    geopandas \
    shapely \
    owslib

# =============================================================================
# Create loader script
# =============================================================================

cat << EOF > "${LOADER_NAME}"
#!/bin/bash

CONDA_DIR="${CONDA_DIR}"
ENV_NAME="${ENV_NAME}"

source "\${CONDA_DIR}/etc/profile.d/conda.sh"

conda activate "\${ENV_NAME}"

echo "Analyzer Soil Tools environment loaded"
echo "Conda Env: \${ENV_NAME}"
echo "Python:    \$(which python)"
echo "Pip:       \$(which pip)"
python --version
EOF

chmod +x "${LOADER_NAME}"

# =============================================================================
# Final message
# =============================================================================

echo ""
echo "========================================="
echo " Environment created successfully"
echo "========================================="
echo ""
echo "To load the environment use:"
echo ""
echo "source ${LOADER_NAME}"
echo ""
