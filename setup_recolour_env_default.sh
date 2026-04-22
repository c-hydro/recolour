#!/usr/bin/env bash

# ============================================
# Recolour: Miniconda + Environment Installer
# ============================================

set -euo pipefail

# === USER CONFIG ===
ENV_NAME="recolour_default"
INSTALL_BASE="${1:-}"   # e.g. /opt/envs or /home/fabio/Desktop/Recolour_Workspace/conda_default
ENV_YML_PATH="${2:-}"   # e.g. /home/fabio/Desktop/Recolour_Workspace/recolour_default_environment.yml
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

# === USAGE / ARG CHECKS ===
if [[ -z "$INSTALL_BASE" || -z "$ENV_YML_PATH" ]]; then
  echo "Usage: $0 <install_base_path> <environment.yml>"
  echo "Example:"
  echo "  $0 /home/user/mytools ./environment.yml"
  exit 1
fi

if [[ ! -f "$ENV_YML_PATH" ]]; then
  echo "❌ Environment file not found: $ENV_YML_PATH"
  exit 1
fi

# === PATHS ===
MINICONDA_DIR="${INSTALL_BASE}/miniconda"
ENV_DIR="${INSTALL_BASE}/envs/${ENV_NAME}"
INFO_FILE="${ENV_DIR}/.env_info"

# Make sure base directories exist and are writable
mkdir -p "${INSTALL_BASE}" "${INSTALL_BASE}/envs" || {
  echo "❌ Cannot create or write to ${INSTALL_BASE}. Check permissions."
  exit 1
}

if ! touch "${INSTALL_BASE}/.write_test" 2>/dev/null; then
  echo "❌ Directory not writable: ${INSTALL_BASE}"
  echo "   Try using a path in your home directory or run with sudo (if appropriate)."
  exit 1
else
  rm -f "${INSTALL_BASE}/.write_test"
fi

# === STEP 1: Install Miniconda if needed ===
if [[ ! -d "$MINICONDA_DIR" ]]; then
  echo "🔽 Miniconda not found at: $MINICONDA_DIR"
  echo "   Downloading from: $MINICONDA_URL"
  INSTALLER_PATH="${INSTALL_BASE}/miniconda.sh"

  # Try wget first
  if ! wget -q "$MINICONDA_URL" -O "$INSTALLER_PATH"; then
    echo "⚠️  wget failed, trying curl..."
    if ! curl -L "$MINICONDA_URL" -o "$INSTALLER_PATH"; then
      echo "❌ Failed to download Miniconda using wget AND curl."
      echo "   URL: $MINICONDA_URL"
      echo "   Target: $INSTALLER_PATH"
      exit 2
    fi
  fi

  echo "📦 Installing Miniconda to $MINICONDA_DIR..."
  bash "$INSTALLER_PATH" -b -p "$MINICONDA_DIR"

  echo "🧹 Cleaning up installer..."
  rm -f "$INSTALLER_PATH"
else
  echo "✅ Miniconda already installed at $MINICONDA_DIR"
fi

# === STEP 2: Initialize Conda (temporary shell hook) ===
if [[ ! -x "${MINICONDA_DIR}/bin/conda" ]]; then
  echo "❌ conda executable not found at ${MINICONDA_DIR}/bin/conda"
  exit 1
fi

echo "🔧 Initializing conda for this shell..."
eval "$("${MINICONDA_DIR}/bin/conda" shell.bash hook)"
conda config --set auto_activate_base false

# === STEP 3: Create environment if not exists ===
if [[ -d "$ENV_DIR" ]]; then
  echo "⚠️  Environment already exists at:"
  echo "    $ENV_DIR"
  echo "    Skipping creation (no update performed)."
else
  echo "📁 Creating conda environment at $ENV_DIR..."
  conda env create -f "$ENV_YML_PATH" -p "$ENV_DIR"
  if [[ $? -ne 0 ]]; then
    echo "❌ Environment creation failed."
    exit 3
  fi
fi

# === STEP 4: Write activation info ===
echo "🔐 Writing activation info to $INFO_FILE"
mkdir -p "$ENV_DIR"
{
  echo "# To activate the environment, run:"
  echo "source \"${MINICONDA_DIR}/bin/activate\" \"$ENV_DIR\""
} > "$INFO_FILE"

# === DONE ===
echo
echo "✅ Environment setup complete!"
echo "📌 To activate your environment later, run:"
echo "    source \"${MINICONDA_DIR}/bin/activate\" \"$ENV_DIR\""
echo

