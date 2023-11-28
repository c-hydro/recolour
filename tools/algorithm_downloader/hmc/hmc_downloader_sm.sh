#!/bin/bash -e

# ----------------------------------------------------------------------------------------
# script generic information
script_name='RECOLOUR - DOWNLOADER - SOIL MOISTURE HMC'
script_version="3.1.0"
script_date='2023/11/29'

# script machine arg(s)
machine_reference="server_cima_67"
machine_url="130.251.104.67"
machine_usr=""
machine_pwd=""

# netrc scheme
# machine machine_string
# login login_string password pass_string

# scritp data condition(s)
data_reset=false # if true, reset destination file

# script data arg(s)
data_days_list=(
	5
)

data_description_list=(
	"HMC  - Continuum Hydrological Model Soil Moisture 500m resolution - TIF" 
)

data_period_list=(
	"2008-01-01 :: NOW"
)

data_name_list=(
	"hmc"
)

data_expected_list=(
	1
)

data_active_list=(
	true
)

data_file_src_list=(
	"SoilMoistureItaly_%YYYY%MM%DD230000.tif"
)

data_file_dst_list=(
	"SoilMoistureItaly_%YYYY%MM%DD230000.tif"
)

# case realtime
data_folder_src_list=(
	"/home/silvestro/Flood_Proofs_Italia/Results_AllDomains/SoilMoisture/%YYYY/%MM/%DD/"
	)
	
# case history
# data_folder_src_list=(
# 	"/home/silvestro/Flood_Proofs_Italia/Results_AllDomains/SoilMoisture_Historical/%YYYY/%MM/%DD/"
# 	)
	
data_folder_dst_list=(
	"/home/hsaf/share/recolour/hmc_grid/%YYYY/%MM/%DD/"
	)

# Script time arg(s) (actual time)
time_now=$(date '+%Y-%m-%d %H:00')
# time_now='2023-04-23'
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# info script start
echo " ==================================================================================="
echo " ==> "${script_name}" (Version: "${script_version}" Release_Date: "${script_date}")"
echo " ==> START ..."

# get credentials from .netrc (if not defined in the bash script)
if [[ -z ${machine_usr} || -z ${machine_pwd} ]]; then

	# check .netrc file availability
	netrc_file=~/.netrc
	if [ ! -f "${netrc_file}" ]; then
	  echo "${netrc_file} does not exist. Please create it to store login and password on your machine"
	  exit 0
	fi
	
	# get information from .netrc file
	machine_usr=$(awk '/'${machine_reference}'/{getline; print $2}' ~/.netrc)
	machine_pwd=$(awk '/'${machine_reference}'/{getline; print $4}' ~/.netrc)

fi
echo " ===> INFO MACHINE -- URL: ${machine_url} -- USER: ${machine_usr}"

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
	echo " ====> PRODUCT NAME '"${data_name_step}"' ... "
	echo " ::: PRODUCT DESCRIPTION: ${data_description_step}"
	# ----------------------------------------------------------------------------------------

	# ----------------------------------------------------------------------------------------
	# parse data period to extract data start and data end
	data_period_nospace=$(echo $data_period_step | tr -d '[:space:]')

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
				year_data_step=${time_data_step:0:4};
				month_data_step=${time_data_step:4:2}; day_data_step=${time_data_step:6:2}
				hour_data_step='00'

				# info time download start
				echo " =====> TIME DOWNLOAD: "${time_data_step:0:8}" ... "
				# ----------------------------------------------------------------------------------------

				# ----------------------------------------------------------------------------------------
				# Define dynamic folder(s)
				data_folder_src_step=${data_folder_src_tmp/'%YYYY'/$year_data_step}
				data_folder_src_step=${data_folder_src_step/'%MM'/$month_data_step}
				data_folder_src_step=${data_folder_src_step/'%DD'/$day_data_step}
				data_folder_src_step=${data_folder_src_step/'%HH'/$hour_data_step}

				data_file_src_step=${data_file_src_tmp/'%YYYY'/$year_data_step}
				data_file_src_step=${data_file_src_step/'%MM'/$month_data_step}
				data_file_src_step=${data_file_src_step/'%DD'/$day_data_step}
				data_file_src_step=${data_file_src_step/'%HH'/$hour_data_step}

				data_folder_dst_step=${data_folder_dst_tmp/'%YYYY'/$year_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%MM'/$month_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%DD'/$day_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%HH'/$hour_data_step}

				data_file_dst_step=${data_file_dst_tmp/'%YYYY'/$year_data_step}
				data_file_dst_step=${data_file_dst_step/'%MM'/$month_data_step}
				data_file_dst_step=${data_file_dst_step/'%DD'/$day_data_step}
				data_file_dst_step=${data_file_dst_step/'%HH'/$hour_data_step}
				# ----------------------------------------------------------------------------------------

				# ----------------------------------------------------------------------------------------
				# Create folder(s)
				if [ ! -d "$data_folder_dst_step" ]; then
					mkdir -p $data_folder_dst_step
				fi
				# ----------------------------------------------------------------------------------------

				# ----------------------------------------------------------------------------------------
				# remove file (if flag_reset = true)
				if [ "${data_reset}" = true ] ; then
					if [ -e ${data_folder_dst_step}/${data_file_dst_step} ]; then
						rm -rf ${data_folder_dst_step}/${data_file_dst_step}
					fi
				fi

				# check file exist or not in the destination folder
				if [ -e ${data_folder_dst_step}/${data_file_dst_step} ]; then
					data_flag_download=false
				else
					data_flag_download=true
				fi
				# ----------------------------------------------------------------------------------------

				# ----------------------------------------------------------------------------------------
				# flag to activate download
				if [ "${data_flag_download}" = true ] ; then

					# info download file name start
					echo -n " ======> DOWNLOAD FILE: ${data_folder_src_step}${data_file_src_step} IN ${data_folder_dst_step}"

					# get download file name
					rsync -ar ${machine_usr}@${machine_url}:${data_folder_src_step}/${data_file_src_step} ${data_folder_dst_step}/${data_file_dst_step} --info=progress2
					
					# info download file name end
					if [ $? -eq 0 ] > /dev/null 2>&1; then
				 		echo " ... DONE!"
					else
						echo " ... FAILED! DOWNLOAD ERROR!"
					fi

					# info time download end
					echo " =====> TIME DOWNLOAD: "${time_data_step:0:8}" ... DONE "

				else

					# info time download end
					echo " =====> TIME DOWNLOAD: "${time_data_step:0:8}" ... SKIPPED. FILE DOWNLOAD PREVIOUSLY."
				fi
				# ----------------------------------------------------------------------------------------

			done

			# info name end
			echo " ====> PRODUCT NAME '"${data_name_step}"' ... DONE"
			# ----------------------------------------------------------------------------------------

		else

			# ----------------------------------------------------------------------------------------
			# info name end
			echo " ====> PRODUCT NAME '"${data_name_step}"' ... SKIPPED. DOWNLOAD IS NOT ACTIVATED."
			# ----------------------------------------------------------------------------------------

		fi
		# ----------------------------------------------------------------------------------------

	else

		# ----------------------------------------------------------------------------------------
		# info name end
		echo " ====> PRODUCT NAME '"${data_name_step}"' ... SKIPPED. TIME NOW NOT IN THE TIME PERIOD"
		# ----------------------------------------------------------------------------------------

	fi
	# ----------------------------------------------------------------------------------------

done

# info script end
echo " ==> "${script_name}" (Version: "${script_version}" Release_Date: "${script_date}")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------
