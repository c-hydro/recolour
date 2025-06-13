#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='ASCAT POINTS - RUN APP_POINTS_ASCAT'
script_version="1.2.0"
script_date='2025/06/13'

# Script settings
script_folder_app='./'
script_file_app='app_points_ascat.py'
script_folder_settings='./'
script_file_settings_default='app_points_ascat.json'

# Time period execution (days)
time_period_default=10

#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Parse inputs: [RUN_DATE] [SETTINGS_JSON] [TIME_PERIOD_DAYS] [PYTHON_SCRIPT]
#   RUN_DATE: override end date in YYYY-MM-DD format (defaults to today)
if [ $# -gt 4 ]; then
  echo "Usage: $0 [RUN_DATE] [SETTINGS_JSON] [TIME_PERIOD_DAYS] [PYTHON_SCRIPT]"
  exit 1
fi

run_date_arg="${1:-}"
settings_file="${2:-$script_folder_settings$script_file_settings_default}"
time_period_days="${3:-$time_period_default}"
python_script="${4:-$script_folder_app$script_file_app}"

# Validate RUN_DATE if provided
if [ -n "$run_date_arg" ]; then
  if ! date -d "$run_date_arg" >/dev/null 2>&1; then
    echo "Error: RUN_DATE '$run_date_arg' is not a valid date (YYYY-MM-DD)."
    exit 1
  fi
  END_DATE="$run_date_arg"
else
  END_DATE=$(date -I)
fi

# Validate settings file
if [ ! -r "$settings_file" ]; then
  echo "Error: Cannot read settings file '$settings_file'."
  exit 1
fi

# Validate time_period_days
if ! [[ "$time_period_days" =~ ^[0-9]+$ ]]; then
  echo "Error: TIME_PERIOD_DAYS must be a non-negative integer."
  exit 1
fi

# Validate python script
if [ ! -f "$python_script" ]; then
  echo "Error: Python script '$python_script' not found."
  exit 1
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Environment setup (optional)
# export PYTHONPATH="${PYTHONPATH}:$script_folder_app"
# virtual_env_folder='/path/to/venv/bin/'
# virtual_env_name='your_env_name'
# export PATH="$virtual_env_folder:$PATH"
# source activate $virtual_env_name
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Compute start date based on time_period_days
START_DATE=$(date -I -d "$time_period_days days ago")

echo " ==================================================================================="
echo " ==> $script_name (Version: $script_version Release_Date: $script_date)"
echo " ==> SETTINGS   : $settings_file"
echo " ==> RUN_DATE   : $END_DATE"
echo " ==> PERIOD_DAYS: Last $time_period_days days (from $START_DATE to $END_DATE)"
echo " ==> START ..."
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Loop over each day from START_DATE → END_DATE, at 00:00 and 12:00
current="$START_DATE"
while [[ "$current" < "$END_DATE" || "$current" == "$END_DATE" ]]; do
  for hour in "00:00" "12:00"; do
    time_step="${current} ${hour}"
    echo -n " ===> Running: python $python_script -settings_file $settings_file -time \"$time_step\" "
    python "$python_script" \
      -settings_file "$settings_file" \
      -time "$time_step"
    echo " ... DONE!"
  done
  current=$(date -I -d "$current + 1 day")
done

#-----------------------------------------------------------------------------------------
echo " ==> $script_name (Version: $script_version Release_Date: $script_date)"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
#-----------------------------------------------------------------------------------------

