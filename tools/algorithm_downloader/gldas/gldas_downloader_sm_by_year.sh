#!/bin/bash -e

# ----------------------------------------------------------------------------------------
# script generic information
script_name='RECOLOUR - DOWNLOADER - SOIL MOISTURE GLDAS'
script_version="2.0.0"
script_date='2024/03/07'

# script machine arg(s)
# netrc scheme -- machine machine_string login login_string password pass_string
machine_reference="urs.earthdata.nasa.gov"
machine_url="urs.earthdata.nasa.gov"
machine_usr="fabiodelogu"
machine_pwd="Cuba1978"

# script data condition(s)
data_reset=false # if true, reset destination file
# script data command
data_command=

# script data arg(s)
data_year_list=(2021 2022 2023)

data_description_list=(
	"GLDAS Noah Land Surface Model L4 3 hourly 0.25 x 0.25 degree V2.1 (GLDAS_NOAH025_3H) - NETCDF" 
)

data_period_list=(
	"2008-01-01 :: 2023-12-01"
)

data_name_list=(
	"gldas"
)

data_command_list=(
	"wget --load-cookies ~/.urs_cookies  --save-cookies ~/.urs_cookies --keep-session-cookies -r -c -nH -nd -np -A nc4 -P %FOLDER_DST --recursive --content-disposition"
)

data_active_list=(
	true
)

data_reference_src_list=(
	"https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH025_3H.2.1/%YYYY/%JDAY/"
)

# case realtime
data_reference_dst_list=(
	"/share/GLDAS/gldas_noah025_3h_v2.1_grid/%YYYY/%JDAY/"
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

# iterate over year
for ((i=0;i<${#data_year_list[@]};++i)); do 
	
	# ----------------------------------------------------------------------------------------
	# get time information
	data_year_step="${data_year_list[i]}"
	# start time informantion 
	echo " ====> PRODUCT YEAR: ${data_year_step} ... "
	
	# leap year or not
	if [ `expr $data_year_step % 400` -eq 0 ]
	then
		echo " ====> LEAP YEAR"
		data_jday_n=366
	elif [ `expr $data_year_step % 100` -eq 0 ]
	then
		echo " ====> NOT LEAP YEAR"
		echo data_jday_n=365
	elif [ `expr $data_year_step % 4` -eq 0 ]
	then
		echo " ====> LEAP YEAR"
		data_jday_n=366
	else
		echo " ====> NOT LEAP YEAR"
		data_jday_n=365
	fi
	
	# iterate over days
	for data_jday_step in $(seq -w 1 $data_jday_n); do
		
		# ----------------------------------------------------------------------------------------
		# start jday information
		echo " =====> PRODUCT JDAY: ${data_jday_step} ... "

		# iterate over data name(s)
		for ((i=0;i<${#data_description_list[@]};++i)); do

			# ----------------------------------------------------------------------------------------
			# get data information
			data_command_raw="${data_command_list[i]}"
			
			data_reference_src_raw="${data_reference_src_list[i]}"
			data_reference_dst_raw="${data_reference_dst_list[i]}"
			data_description_step="${data_description_list[i]}"
			data_name_step="${data_name_list[i]}"
			data_active_step="${data_active_list[i]}"
			data_period_step="${data_period_list[i]}"

			# info description and name start
			echo " ======> PRODUCT NAME '"${data_name_step}"' ... "
			echo " ::: PRODUCT DESCRIPTION: ${data_description_step}"
			
			# fill time value(s)
			data_reference_src_step=${data_reference_src_raw/'%YYYY'/$data_year_step}
			data_reference_src_step=${data_reference_src_step/'%JDAY'/$data_jday_step}
			
			data_reference_dst_step=${data_reference_dst_raw/'%YYYY'/$data_year_step}
			data_reference_dst_step=${data_reference_dst_step/'%JDAY'/$data_jday_step}
			
			data_command_step=${data_command_raw/'%FOLDER_DST'/$data_reference_dst_step}

			# flag to activate datasets
			if [ "${data_active_step}" = true ] ; then
				
				# create command-line
				exec_command_step=${data_command_step}' '${data_reference_src_step}

				# download data
				echo " =======> COMMAND LINE: $exec_command_step"
				echo -n " =======> DOWNLOAD ..."
			   
			    if ${exec_command_step} > /dev/null 2>&1; then
		            echo " DONE!"
		        else
		            echo " FAILED! Error in command execution!"
		            exit
		        fi

				# info name end
				echo " ======> PRODUCT NAME '"${data_name_step}"' ... DONE"

			
			else
				# info name end
				echo " ======> PRODUCT NAME '"${data_name_step}"' ... SKIPPED. DOWNLOAD IS NOT ACTIVATED."
			fi
			# ----------------------------------------------------------------------------------------
			
			# ----------------------------------------------------------------------------------------
			# remove tmp file(s)
			file_count=`ls -1 ${data_reference_dst_step}/*.tmp 2>/dev/null | wc -l`
			if [ $file_count != 0 ]
			then
				rm ${data_reference_dst_step}/*.tmp 
			fi
			# ----------------------------------------------------------------------------------------
			
		done
		
		# end jday information
		echo " =====> PRODUCT JDAY: ${data_jday_step} ... DONE"
		# ----------------------------------------------------------------------------------------
	
	done
	
	# end time informantion 
	echo " ====> DOWNLOAD DATA :: YEAR: ${data_year_step} ... DONE"
	# ----------------------------------------------------------------------------------------
	
done

# info script end
echo " ==> "${script_name}" (Version: "${script_version}" Release_Date: "${script_date}")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

