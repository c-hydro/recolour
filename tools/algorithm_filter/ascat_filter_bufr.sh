#!/bin/bash

# === CONFIG ===
CONDA_ENV_NAME="/root/envs/conda_default/envs/recolour_default"
CONDA_ENV_PATH='/root/envs/conda_default/miniconda/bin/'
SCRIPT_PATH="/root/projects/hsaf/algorithm_filters/ascat_filter_bufr.py"
CONFIG_PATH="/root/projects/hsaf/algorithm_filters/ascat_filter_bufr.json"


# === DATE RANGE ===
START_DATE="2025-06-01"
END_DATE="2025-06-15"  # exclusive upper bound (won’t run this date)

# === ACTIVATE CONDA ENV ===
source $CONDA_ENV_PATH/activate $CONDA_ENV_NAME

# === LOOP OVER EACH DAY ===
current_date="$START_DATE"

while [[ "$current_date" < "$END_DATE" ]]; do
  next_date=$(date -I -d "$current_date + 1 day")

  START_TIME="${current_date} 00:00"
  END_TIME="${next_date} 00:00"

  echo "🔁 Processing $START_TIME to $END_TIME"
  python "$SCRIPT_PATH" "$START_TIME" "$END_TIME" "$CONFIG_PATH"

  current_date="$next_date"
done

echo "✅ All done!"

