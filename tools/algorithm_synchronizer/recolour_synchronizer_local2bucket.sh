#!/usr/bin/env bash
#
# Sync /YYYY/MM trees between a local source and S3 using s3cmd.
# Modes:
#   - realtime: last N months including current month
#   - history : specific inclusive START..END month range
#
# Changelog:
#   1.0.0 (2025-09-24)
#     * FEAT: --n-months, --start-ym, --end-ym (retain backward-compat aliases)
#     * FIX: avoid resyncing unintended deeper levels; allows optional subpath filter
#     * IMP: safer span caps; clearer logs; preflight checks
#
# Example usage:
# recolour_synchronizer_local2bucket.sh --mode realtime --n-days 2
# recolour_synchronizer_local2bucket.sh --mode history --start-ym 2024-10 --end-ym 2025-02

set -euo pipefail

# ----------------------------------------------------------------------------------------
# script information
script_name='RECOLOUR - SYNCHRONISER - LOCAL TO S3 BUCKET'
script_version="1.0.0"
script_date='2025/09/24'

# timezone for date math (CET/CEST for CIMA)
export TZ="Europe/Rome"

# s3cmd config file (override with $S3CMD_CONFIG if set)
CONFIG_FILE="${S3CMD_CONFIG:-/root/cima-iride-s601}"

# base folders (no trailing slash) — can be overridden via env
SRC_BASE="${SOURCE_BASE:-/share/SM_TC/map/sm_v2}"
BUCKET_NAME="${S3_BUCKET:-cima-iride-s601}"
DEST_PREFIX_ENV="${DEST_PREFIX:-map_tc}"
DEST_PREFIX="${DEST_PREFIX_ENV%%/}"   # strip trailing slash if any

# Optional: limit sync to a subpath under each month (e.g., "tiles/")
SUBPATH_FILTER="${SUBPATH_FILTER:-}"   # no leading slash; empty = whole month

# Logging / locking (override with env if needed)
LOG_DIR="${LOG_DIR:-/var/log/s3-sync}"
LOCK_FILE="${LOCK_FILE:-/tmp/s3-sync-map_tc.lock}"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Defaults for runtime flags (month-based)
MODE=""                 # realtime | history
N_MONTHS=""            # integer >=1 (realtime)
START_YM=""            # YYYY-MM (history)
END_YM=""              # YYYY-MM (history)
DRY_RUN=false
DELETE_REMOVED=true
VERBOSE=false
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# method to print usage
print_usage() {
  cat <<'EOF'
Usage (month-based):
  sync_tree.sh --mode realtime --n-months N [--dry-run] [--keep] [--verbose]
  sync_tree.sh --mode history  --start-ym YYYY-MM --end-ym YYYY-MM [--dry-run] [--keep] [--verbose]

Options:
  --mode MODE         realtime | history
  --n-months N        number of months in realtime mode (>=1, includes current month)
  --start-ym YM       start month (YYYY-MM) for history mode (inclusive)
  --end-ym YM         end month (YYYY-MM) for history mode (inclusive)
  --dry-run           show actions only, make no changes
  --keep              do not --delete-removed on S3
  --verbose           verbose s3cmd output
  -h, --help          show this help

Backward-compat aliases (will be interpreted as months):
  --n-days => --n-months,  --start => --start-ym,  --end => --end-ym

Environment overrides:
  S3CMD_CONFIG, SOURCE_BASE, S3_BUCKET, DEST_PREFIX, LOG_DIR, LOCK_FILE, SUBPATH_FILTER
EOF
}

# parse args (order-agnostic)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="${2:-}"; shift 2 ;;
    --n-months) N_MONTHS="${2:-}"; shift 2 ;;
    --start-ym) START_YM="${2:-}"; shift 2 ;;
    --end-ym) END_YM="${2:-}"; shift 2 ;;
    # Backward-compat aliases
    --n-days) N_MONTHS="${2:-}"; shift 2 ;;
    --start) START_YM="${2:-}"; shift 2 ;;
    --end) END_YM="${2:-}"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --keep) DELETE_REMOVED=false; shift ;;
    --verbose) VERBOSE=true; shift ;;
    -h|--help) print_usage; exit 0 ;;
    *) echo "Unknown option: $1"; print_usage; exit 2 ;;
  esac
done
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Startup / Logging
mkdir -p "$LOG_DIR"
TIMESTAMP="$(date '+%Y-%m-%d_%H-%M-%S')"
LOG_FILE="${LOG_DIR}/map_tc_${TIMESTAMP}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# script start
echo "==================================================================================="
echo " ==> $script_name (Version: $script_version  Release_Date: $script_date)"
echo " ::: [$(date '+%F %T')] Starting sync"
echo " ::: Mode: $MODE"
echo " ::: Args: N_MONTHS=$N_MONTHS  START_YM=$START_YM  END_YM=$END_YM"
echo " ::: Config: CONFIG_FILE=$CONFIG_FILE  SRC_BASE=$SRC_BASE  BUCKET=$BUCKET_NAME  DEST_PREFIX=$DEST_PREFIX  SUBPATH_FILTER='${SUBPATH_FILTER}'"
echo " ::: Flags: dry_run=$DRY_RUN  delete_removed=$DELETE_REMOVED  verbose=$VERBOSE"
echo " ::: Log: $LOG_FILE"
echo " "

# check deps
require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo " ===> ERROR: '$1' not found in PATH"; exit 1; }; }
require_cmd s3cmd
require_cmd date

[[ -f "$CONFIG_FILE" ]] || { echo " ===> ERROR: s3cmd config '$CONFIG_FILE' not found"; exit 1; }
[[ -d "$SRC_BASE" ]] || { echo " ===> ERROR: source base '$SRC_BASE' not found"; exit 1; }

# DEST_PREFIX must not start with /
if [[ "$DEST_PREFIX" == /* ]]; then
  echo " ===> ERROR: DEST_PREFIX must not start with '/': '$DEST_PREFIX'"; exit 1
fi

# lightweight bucket reachability probe (non-fatal if 403 but exists)
if ! s3cmd -c "$CONFIG_FILE" ls "s3://$BUCKET_NAME/" >/dev/null 2>&1; then
  echo " ===> WARNING: unable to list s3://$BUCKET_NAME/ (permissions or bucket). Continuing…"
fi

# lock to avoid concurrent runs
exec 9> "$LOCK_FILE"
if ! flock -n 9; then
  echo " ===> WARNING: another sync is running (lock: $LOCK_FILE). Exiting."
  exit 0
fi
cleanup() { echo "[$(date '+%F %T')] Cleanup: releasing lock"; flock -u 9 || true; rm -f "$LOCK_FILE" || true; }
trap cleanup EXIT INT TERM
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# check date
is_valid_ym() { date -d "$1-01" '+%Y-%m' >/dev/null 2>&1; }
inc_month() { date -d "$1-01 +1 month" '+%Y-%m'; }
dec_month() { date -d "$1-01 -1 month" '+%Y-%m'; }
now_ym() { date '+%Y-%m'; }

# max month period
MAX_MONTHS_DEFAULT=240   # 20 years
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# s3cmd options
S3CMD_BASE_OPTS=(
  -c "$CONFIG_FILE"
  sync
  --recursive
  --progress
  --stats
  --acl-private
  --guess-mime-type
  --exclude '*~'
  --exclude '*.tmp'
  --exclude '.*.swp'
  --exclude '*.part'
)
$DRY_RUN && S3CMD_BASE_OPTS+=( --dry-run )
$DELETE_REMOVED && S3CMD_BASE_OPTS+=( --delete-removed )
$VERBOSE && S3CMD_BASE_OPTS+=( -d )
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Core: sync exactly one /YYYY/MM directory (optionally with SUBPATH_FILTER)
sync_one_month() {
  local y="$1" m="$2"
  local month_path="${y}/${m}"
  local src_dir="${SRC_BASE}/${month_path}/"
  local dest_uri_base="s3://${BUCKET_NAME}/${DEST_PREFIX}/${month_path}/"
  local src_path="$src_dir${SUBPATH_FILTER}"
  local dest_uri="$dest_uri_base${SUBPATH_FILTER}"

  echo " ====> [$(date '+%F %T')] Month: $month_path"
  echo " ====> Source       : $src_path"
  echo " ====> Destination  : $dest_uri"

  # If SUBPATH_FILTER is set, allow the month folder to exist while subpath may not
  if [[ ! -d "$src_dir" ]]; then
    echo " ===> WARNING: source month not found, skipping: $src_dir"
    return 0
  fi
  if [[ -n "$SUBPATH_FILTER" && ! -e "$src_path" ]]; then
    echo " ===> INFO: subpath not present under month, skipping: $src_path"
    return 0
  fi

  # Attempt to list dest prefix (not fatal if missing)
  if ! s3cmd -c "$CONFIG_FILE" ls "$dest_uri_base" >/dev/null 2>&1; then
    echo " ===> INFO: destination prefix not present yet (will be created if needed): $dest_uri_base"
  fi

  set -x
  s3cmd "${S3CMD_BASE_OPTS[@]}" "$src_path" "$dest_uri"
  local rc=$?
  set +x
  return $rc
}

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# sync data
errors=0

# info sync start
echo " ===> SYNC DATA ... "

# check mode
case "$MODE" in
  realtime)
    if [[ -z "$N_MONTHS" || ! "$N_MONTHS" =~ ^[0-9]+$ || "$N_MONTHS" -lt 1 ]]; then
      echo " ===> ERROR: --n-months must be an integer >= 1 for realtime mode."; exit 2
    fi
    if (( N_MONTHS > MAX_MONTHS_DEFAULT )); then
      echo " ===> ERROR: --n-months too large (>$MAX_MONTHS_DEFAULT). Refusing to run."; exit 2
    fi
    current="$(now_ym)"  # YYYY-MM
    for ((i=1; i<=N_MONTHS; i++)); do
      Y=$(date -d "$current-01" '+%Y'); M=$(date -d "$current-01" '+%m')
      if ! sync_one_month "$Y" "$M"; then errors=$((errors+1)); fi
      current="$(dec_month "$current")"
    done
    ;;

  history)
    if [[ -z "$START_YM" || -z "$END_YM" ]]; then
      echo " ===> ERROR: --start-ym and --end-ym are required for history mode."; exit 2
    fi
    if ! is_valid_ym "$START_YM" || ! is_valid_ym "$END_YM"; then
      echo " ===> ERROR: invalid month format; use YYYY-MM."; exit 2
    fi
    # Ensure START <= END
    if [[ "$(date -d "$START_YM-01" +%s)" -gt "$(date -d "$END_YM-01" +%s)" ]]; then
      echo " ===> ERROR: start month is after end month."; exit 2
    fi

    # Cap span to protect from runaway ranges
    span_count=0; tmp="$START_YM"
    while true; do
      span_count=$((span_count+1))
      [[ "$tmp" == "$END_YM" ]] && break
      tmp="$(inc_month "$tmp")"
      if (( span_count > MAX_MONTHS_DEFAULT )); then
        echo " ===> ERROR: month span exceeds $MAX_MONTHS_DEFAULT; refusing to run."; exit 2
      fi
    done

    current="$START_YM"
    while true; do
      Y=$(date -d "$current-01" '+%Y'); M=$(date -d "$current-01" '+%m')
      if ! sync_one_month "$Y" "$M"; then errors=$((errors+1)); fi
      [[ "$current" == "$END_YM" ]] && break
      current="$(inc_month "$current")"
    done
    ;;

  *)
    echo " ===> ERROR: --mode must be 'realtime' or 'history'."; print_usage; exit 2 ;;

esac

# info sync end
echo " ===> SYNC DATA … DONE"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# check errors
if [[ "$errors" -gt 0 ]]; then
  echo " ===> Completed with $errors error(s). Check the log: $LOG_FILE"
  exit 1
else
  echo " ===> Completed successfully. Log: $LOG_FILE"
fi
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# script end
echo " ==> $script_name (Version: $script_version  Release_Date: $script_date)"
echo " ==> … END"
echo " ==> Bye, Bye"
echo "==================================================================================="
# ----------------------------------------------------------------------------------------

