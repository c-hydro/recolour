#!/bin/bash

# ----------------------------------------------------------------------------------------
# Script information
script_name='HMC FILE ORGANIZER'
script_version="2.2.0"
script_date='2025/07/21'
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- CONFIGURATION ---
DOMAIN_STRING=(
  'AvetoTrebbiaDomain' 'BormidaMDomain' 'BormidaSDomain' 
  'CentaDomain' 'CentroPonenteDomain' 'EntellaDomain' 'ErroDomain' 
  'FinaleseDomain' 'ImperieseDomain' 'LevanteGenoveseDomain' 
  'MagraDomain' 'OrbaSturaDomain' 'PonenteDomain' 'PonenteGenoveseDomain'
  'RoiaDomain' 'SavoneseDomain' 'ScriviaDomain' 'TanaroDomain')

START_DATE="2025-04-16"
END_DATE="2025-04-19"

SRC_FOLDER="/home/cfmi.arpal.org/satsuolo/MODELLO_CONTINUUM/archive"
DEST_FOLDER="/home/cfmi.arpal.org/satsuolo/hmc_grid_selected"

TIME_FORMAT_FOLDER='%Y/%m/%d'
TIME_FORMAT_FILE='%Y%m%d%H00'

LOG_FILE="/tmp/hmc_grid_copy_$(date +%Y%m%d%H%M%S).log"

# --- LOGGING FUNCTION ---
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# --- PATH RESOLVER FUNCTION ---
resolve_paths() {
  local domain=$1
  local time=$2

  folder=$(date -d "$time" +"$TIME_FORMAT_FOLDER")
  file=$(date -d "$time" +"$TIME_FORMAT_FILE")

  SRC_FOLDER_DEF="$SRC_FOLDER/$domain/weather_stations_realtime/$folder"
  DEST_FOLDER_DEF="$DEST_FOLDER/$domain/weather_stations_realtime/$folder"
  SRC_FILE_DEF="hmc.output-grid.${file}.nc.gz"
  SRC_FILE_PATH="$SRC_FOLDER_DEF/$SRC_FILE_DEF"
}
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- INITIAL LOG ---
echo " ===================================================================================" | tee -a "$LOG_FILE"
log "==> $script_name (Version: $script_version Release_Date: $script_date)"
log "==> START ..."
log "===> HMC Grid Copy Log - Started."

# --- CHECK DOMAIN ARRAY ---
if [ ${#DOMAIN_STRING[@]} -eq 0 ]; then
  log "ERROR: DOMAIN_STRING is empty. Aborting."
  exit 1
fi

# --- GENERATE TIME LIST: ALL HOURS FOR EACH DAY ---
TIME_LIST=()
current_day=$(date -d "$START_DATE" +%Y-%m-%d)
end_day=$(date -d "$END_DATE" +%Y-%m-%d)

while [ "$(date -d "$current_day" +%s)" -le "$(date -d "$end_day" +%s)" ]; do
  for hour in {0..23}; do
    TIME_LIST+=("$current_day $(printf "%02d:00" $hour)")
  done
  current_day=$(date -I -d "$current_day + 1 day")
done

# --- MAIN LOOP ---
for domain in "${DOMAIN_STRING[@]}"; do
  for time in "${TIME_LIST[@]}"; do

    resolve_paths "$domain" "$time"

    log "====> Processing: $SRC_FILE_PATH"
    mkdir -p "$DEST_FOLDER_DEF"

    if [ -f "$SRC_FILE_PATH" ]; then
      cp "$SRC_FILE_PATH" "$DEST_FOLDER_DEF/"
      log "====> SUCCESS: Copied to $DEST_FOLDER_DEF/"
    else
      log "====> WARNING: File not found: $SRC_FILE_PATH"
    fi

  done
done
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- FINAL LOG ---
log "===> Finished."
log "==> $script_name (Version: $script_version Release_Date: $script_date)"
log "==> ... END"
log "==> Bye, Bye"
echo " ===================================================================================" | tee -a "$LOG_FILE"

exit 0
# ----------------------------------------------------------------------------------------


