#!/bin/bash -e

# ----------------------------------------------------------------------------------------
# script generic information
script_name='RECOLOUR - DOWNLOADER - SOIL MOISTURE HMC - NRT'
script_version="3.2.3"
script_date='2025/10/23'

# SSH / machine args (key-based auth; inline in rsync -e)
machine_host=""
machine_port=""
machine_usr=""
machine_key=""   # private key path (paired with /root/.ssh/id_ecdsa.pub)

# optional: fail if key missing
if [ ! -f "${machine_key}" ]; then
  echo "ERROR: SSH private key not found at ${machine_key}"
  exit 1
fi

# Common SSH options for rsync (inline via -e)
SSH_E_OPTS="ssh -i ${machine_key} -p ${machine_port} -o BatchMode=yes -o StrictHostKeyChecking=no -o UserKnownHostsFile=/root/.ssh/known_hosts -o ConnectTimeout=20"

# rsync opts (avoid chown/chgrp/chmod on destination)
RSYNC_OPTS="-rlt --info=progress2 --no-owner --no-group --no-perms"

# script data condition(s)
data_reset=true # if true, reset destination file

# script data arg(s)
data_days_list=( 5 )

data_description_list=(
  "HMC  - Continuum Hydrological Model Soil Moisture 500m resolution - TIF"
)

data_period_list=( "2008-01-01 :: NOW" )

data_name_list=( "hmc" )

data_expected_list=( 1 )

data_active_list=( true )

data_file_src_list=( "SoilMoistureItaly_%YYYY%MM%DD230000.tif" )
data_file_dst_list=( "SoilMoistureItaly_%YYYY%MM%DD230000.tif" )

# case realtime
data_folder_src_list=(
  "/share/fp-ita/Publishfolder/grid/SoilMoisture/%YYYY/%MM/%DD/"
)

# case history (disabled)
# data_folder_src_list=(
#   "/home/silvestro/Flood_Proofs_Italia/Results_AllDomains/SoilMoisture_Historical/%YYYY/%MM/%DD/"
# )

data_folder_dst_list=(
  "/share/HMC/nrt/%YYYY/%MM/%DD/"
)

# Script time arg(s) (actual time)
time_now=$(date '+%Y-%m-%d %H:00')
# time_now='2023-04-23'
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# info script start
echo " ==================================================================================="
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> START ..."
echo " ===> INFO MACHINE -- HOST: ${machine_host} -- PORT: ${machine_port} -- USER: ${machine_usr}"
echo " ===> INFO SSH KEY -- ${machine_key}"

# parse and check time information
time_data_now=$(date -d "${time_now}" +'%Y%m%d%H%M')
echo " ===> INFO TIME -- TIME: ${time_data_now}"

# iterate over data name(s)
for ((i=0;i<${#data_description_list[@]};++i)); do

  # ----------------------------------------------------------------------------------------
  # get data information
  data_folder_src_tmp="${data_folder_src_list[i]}"
  data_file_src_tmp="${data_file_src_list[i]}"
  data_folder_dst_tmp="${data_folder_dst_list[i]}"
  data_file_dst_tmp="${data_file_dst_list[i]}"
  data_name_step="${data_name_list[i]}"
  data_description_step="${data_description_list[i]}"
  data_period_step="${data_period_list[i]}"
  data_days_step="${data_days_list[i]}"
  data_active_step="${data_active_list[i]}"

  # info description and name start
  echo " ====> PRODUCT NAME '${data_name_step}' ... "
  echo " ::: PRODUCT DESCRIPTION: ${data_description_step}"
  # ----------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------
  # parse data period to extract data start and data end
  data_period_nospace=$(echo "$data_period_step" | tr -d '[:space:]')

  # extract data start and data end
  IFS_DEFAULT="$IFS"
  IFS="::"; read data_period_start data_period_end <<< "$data_period_nospace"
  IFS="$IFS_DEFAULT"
  unset IFS_DEFAULT
  # adjust data format
  data_period_start=${data_period_start/':'/''}
  data_period_end=${data_period_end/':'/''}

  # time period now
  time_period_now=$(date -d "${time_now}" +'%Y%m%d')
  # time period start
  time_period_start=$(date -d "${data_period_start}" +'%Y%m%d')
  # time period end
  if [ "${data_period_end}" == "NOW" ] ; then
    time_period_end=$(date -d "${time_now}" +'%Y%m%d')
  else
    time_period_end=$(date -d "${data_period_end}" +'%Y%m%d')
  fi

  # info time(s)
  echo " ::: PRODUCT PERIOD: ${data_period_step}"
  echo " ::: PRODUCT NOW: ${time_period_now} -- PRODUCT START: ${time_period_start} PRODUCT END: ${time_period_end}"
  # ----------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------
  # check the time with the reference product period
  if [[ $time_period_now -ge $time_period_start ]] && [[ $time_period_now -le $time_period_end ]] ; then

    # ----------------------------------------------------------------------------------------
    # flag to activate datasets
    if [ "${data_active_step}" = true ] ; then

      # ----------------------------------------------------------------------------------------
      # iterate over days
      for day in $(seq ${data_days_step} -1 0); do

        # ----------------------------------------------------------------------------------------
        # get time step
        time_data_step=$(date -d "$time_now ${day} days ago" +'%Y%m%d%H%M')
        year_data_step=${time_data_step:0:4}
        month_data_step=${time_data_step:4:2}
        day_data_step=${time_data_step:6:2}
        hour_data_step='00'

        # info time download start
        echo " =====> TIME DOWNLOAD: ${time_data_step:0:8} ... "
        # ----------------------------------------------------------------------------------------

        # ----------------------------------------------------------------------------------------
        # Define dynamic folder(s)
        data_folder_src_step=${data_folder_src_tmp/'%YYYY'/$year_data_step}
        data_folder_src_step=${data_folder_src_step/'%MM'/$month_data_step}
        data_folder_src_step=${data_folder_src_step/'%DD'/$day_data_step}
        data_folder_src_step=${data_folder_src_step/'%HH'/$hour_data_step}
        data_folder_src_step="${data_folder_src_step%/}"  # strip trailing slash

        data_file_src_step=${data_file_src_tmp/'%YYYY'/$year_data_step}
        data_file_src_step=${data_file_src_step/'%MM'/$month_data_step}
        data_file_src_step=${data_file_src_step/'%DD'/$day_data_step}
        data_file_src_step=${data_file_src_step/'%HH'/$hour_data_step}

        data_folder_dst_step=${data_folder_dst_tmp/'%YYYY'/$year_data_step}
        data_folder_dst_step=${data_folder_dst_step/'%MM'/$month_data_step}
        data_folder_dst_step=${data_folder_dst_step/'%DD'/$day_data_step}
        data_folder_dst_step=${data_folder_dst_step/'%HH'/$hour_data_step}
        data_folder_dst_step="${data_folder_dst_step%/}"  # strip trailing slash

        data_file_dst_step=${data_file_dst_tmp/'%YYYY'/$year_data_step}
        data_file_dst_step=${data_file_dst_step/'%MM'/$month_data_step}
        data_file_dst_step=${data_file_dst_step/'%DD'/$day_data_step}
        data_file_dst_step=${data_file_dst_step/'%HH'/$hour_data_step}
        # ----------------------------------------------------------------------------------------

        # ----------------------------------------------------------------------------------------
        # Create folder(s)
        if [ ! -d "$data_folder_dst_step" ]; then
          mkdir -p "$data_folder_dst_step"
        fi
        # ----------------------------------------------------------------------------------------

        # ----------------------------------------------------------------------------------------
        # remove file (if flag_reset = true)
        if [ "${data_reset}" = true ] ; then
          if [ -e "${data_folder_dst_step}/${data_file_dst_step}" ]; then
            rm -f "${data_folder_dst_step}/${data_file_dst_step}"
          fi
        fi

        # check file exist or not in the destination folder
        if [ -e "${data_folder_dst_step}/${data_file_dst_step}" ]; then
          data_flag_download=false
        else
          data_flag_download=true
        fi
        # ----------------------------------------------------------------------------------------

        # ----------------------------------------------------------------------------------------
        # flag to activate download
        if [ "${data_flag_download}" = true ] ; then

          # Pre-flight: check remote file readability; if missing, skip quietly
          if ! ssh -e none -i "${machine_key}" -p "${machine_port}" \
            -o BatchMode=yes -o StrictHostKeyChecking=no -o UserKnownHostsFile=/root/.ssh/known_hosts -o ConnectTimeout=20 \
            "${machine_usr}@${machine_host}" \
            "test -r '${data_folder_src_step}/${data_file_src_step}'"; then
            echo " ..... SKIP: Remote file not found or not readable:"
            echo "             ${data_folder_src_step}/${data_file_src_step}"
            continue
          fi

          echo " ======> DOWNLOAD (rsync over ssh):"
          echo "         ${machine_usr}@${machine_host}:${data_folder_src_step}/${data_file_src_step}"
          echo "         -> ${data_folder_dst_step}/${data_file_dst_step}"

          # rsync command (inline -e), avoid owner/group/perms preservation (prevents chown errors)
          if rsync -e "${SSH_E_OPTS}" ${RSYNC_OPTS} \
            "${machine_usr}@${machine_host}:${data_folder_src_step}/${data_file_src_step}" \
            "${data_folder_dst_step}/${data_file_dst_step}"; then
            echo " ..... DONE"
          else
            rc=$?
            echo " ..... FAILED (exit code: ${rc})"
            case "${rc}" in
              23) echo "       Hint: Check remote path/read perms or destination disk space.";;
              24) echo "       Hint: Partial transfer; file may have changed during copy.";;
              255) echo "       Hint: SSH failed (key/known_hosts/port). Test: ssh -i ${machine_key} -p ${machine_port} ${machine_usr}@${machine_host} 'echo OK'";;
            esac
          fi

          echo " =====> TIME DOWNLOAD: ${time_data_step:0:8} ... DONE "
        else
          # info time download end
          echo " =====> TIME DOWNLOAD: ${time_data_step:0:8} ... SKIPPED. FILE ALREADY PRESENT."
        fi
        # ----------------------------------------------------------------------------------------

      done

      # info name end
      echo " ====> PRODUCT NAME '${data_name_step}' ... DONE"
      # ----------------------------------------------------------------------------------------

    else
      # info name end
      echo " ====> PRODUCT NAME '${data_name_step}' ... SKIPPED. DOWNLOAD IS NOT ACTIVATED."
    fi
    # ----------------------------------------------------------------------------------------

  else
    # info name end
    echo " ====> PRODUCT NAME '${data_name_step}' ... SKIPPED. TIME NOW NOT IN THE TIME PERIOD"
  fi
  # ----------------------------------------------------------------------------------------

done

# info script end
echo " ==> ${script_name} (Version: ${script_version} Release_Date: ${script_date})"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

