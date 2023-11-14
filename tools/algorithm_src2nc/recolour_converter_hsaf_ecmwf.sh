#!/bin/bash -e

# ----------------------------------------------------------------------------------------
# script generic information
script_name='RECOLOUR - ADAPTER - PRODUCT SOIL MOISTURE ECMWF-RZSM'
script_version="2.0.0"
script_date='2023/08/03'

# executable(s)
# exec_cdo="/home/fabio/Desktop/Apps/cdo-1.9.6_NC-4.6.0_HDF-1.8.17_ECCODES-2.12.5/bin/cdo" # local
exec_cdo="/home/fabio/fp_system_apps/cdo/bin/cdo" # local
# exec_cdo="/home/hsaf/library/cdo-1.7.2rc3_NC-4.1.2_HDF-1.8.17/bin/cdo" # server
exec_ncap2="/home/fabio/fp_system_apps/nco/bin/ncap2"
exec_ncks="/home/fabio/fp_system_apps/nco/bin/ncks"

# domain info
domain_name='europe'
domain_bb='-15,30,30,60'

# scritp data condition(s)
data_reset=false

# script data arg(s)
data_days_list=(
	10 
	10
	10
	10
	10
	10
	10
	10
	10
	10
)

data_description_list=(
	"H14  - Soil Wetness Profile Index in the roots region retrieved by Metop ASCAT surface wetness scatterometer assimilation method - 25 km - data"
	"H14  - Soil Wetness Profile Index in the roots region retrieved by Metop ASCAT surface wetness scatterometer assimilation method - 25 km - auxiliary"
	"H27  - Scatterometer Root Zone Soil Moisture (RZSM) Data Record 16km resolution"
	"H140 - Scatterometer Root Zone Soil Moisture (RZSM) Data Record 16km resolution - Extension"
	"H141 - Scatterometer Root Zone Soil Moisture Climate Data Record 10km resolutiong - GRIB"
	"H141 - Scatterometer Root Zone Soil Moisture Climate Data Record 10km resolutiong - NETCDF"
	"H142 - Scatterometer Root Zone Soil Moisture Climate Data Record 10km resolution - GRIB"
	"H142 - Scatterometer Root Zone Soil Moisture Climate Data Record 10km resolution - NETCDF"
	"H26  - Metop ASCAT NRT Root Zone Soil Moisture Profile Index 10km resolution - GRIB"
	"H26  - Metop ASCAT NRT Root Zone Soil Moisture Profile Index 10km resolution - NETCDF"
)

data_period_list=(
	"2012-07-01 :: NOW"
	"2018-10-04 :: NOW"
	"1992-01-01 :: 2014-12-31"
	"2015-01-01 :: 2016-12-31"
	"1992-01-01 :: 2018-12-31"
	"1992-01-01 :: 2018-12-31"
	"2019-01-01 :: 2021-12-31"
	"2019-01-01 :: 2021-12-31"
	"2021-11-04 :: NOW"
	"2021-11-04 :: NOW"
)

data_name_list=( 
	"h14" 
	"h14"
	"h27" 
	"h140" 
	"h141" 
	"h141"
	"h142" 
	"h142" 
	"h26"
	"h26"
)

data_active_list=(
	true 
	false
	false
	false
	false
	false
	false
	false
	false
	false
)

data_type_compressed_list=(
	"bz2" 
	"bz2"
	""
	""
	"tar.gz"
	""
	"tar.gz"
	""
	"bz2"
	""
)

data_type_format_list=(
	"grib" 
	"grib"
	"grib"
	"grib"
	""
	"netcdf"
	""
	"netcdf"
	"grib"
	"netcdf"
)

apply_fx_make_regular_grid_list=(
	true 
	false
	false
	false
	false
	false
	false
	false
	false
	false
)

apply_fx_make_conversion_list=(
	true 
	false
	false
	false
	false
	false
	false
	false
	false
	false
)

apply_fx_make_var_list=(
	true 
	false
	false
	false
	false
	false
	false
	false
	false
	false
)

apply_fx_make_domain_list=(
	true 
	false
	false
	false
	false
	false
	false
	false
	false
	false
)

apply_fx_make_file_grid_list=(
	true 
	false
	false
	false
	false
	false
	false
	false
	false
	false
)

data_file_src_list=(
	"h14_%YYYY%MM%DD_0000.grib.bz2"
	"t14_%YYYY%MM%DD_0000.grib.bz2"
	"h27_%YYYY%MM%DD00_T1279.grib"
	"h140_%YYYY%MM%DD00.grib"
	"%YYYY%MM%.tar.gz"
	"h141_%YYYY%MM%DD00_R01.nc"
	"%YYYY%MM%.tar.gz"
	"h142_%YYYY%MM%DD00_R01.nc"
	"h26_%YYYY%MM%DD00_TCO1279.grib.bz2"
	"h26_%YYYY%MM%DD00_R01.nc"
)

data_file_unzip_list=(
	"h14_%YYYY%MM%DD_0000_src.grib"
	"t14_%YYYY%MM%DD_0000_src.grib"
	""
	""
	""
	""
	""
	""
	"h26_%YYYY%MM%DD00_TCO1279_src.grib"
	""
)

data_file_grid_list=(
	"h14_%YYYY%MM%DD_0000_fdisk.grib"
	"t14_%YYYY%MM%DD_0000_fdisk.grib"
	""
	""
	""
	""
	""
	""
	"h26_%YYYY%MM%DD00_TCO1279_fdisk.grib"
	""
)
data_file_conversion_list=(
	"h14_%YYYY%MM%DD_0000_fdisk.nc"
	"t14_%YYYY%MM%DD_0000_fdisk.nc"
	""
	""
	""
	""
	""
	""
	"h26_%YYYY%MM%DD00_TCO1279_fdisk.nc"
	""
)

data_file_var_list=(
	"rzsm_h14_datasets_%YYYY%MM%DD_0000_fdisk.nc"
	"rzsm_h14_auxiliary_%YYYY%MM%DD_0000_fdisk.nc"
	""
	""
	""
	""
	""
	""
	"h26_%YYYY%MM%DD00_TCO1279_fdisk.nc"
	""
)

data_file_dst_list=(
	"rzsm_h14_datasets_%YYYY%MM%DD_0000_%DOMAIN.nc"
	"rzsm_h14_auxiliary_%YYYY%MM%DD_0000_%DOMAIN.nc"
	""
	""
	""
	""
	""
	""
	"rzsm_h26_datasets_%YYYY%MM%DD_0000_%DOMAIN.nc"
	""
)

grid_file_list=(
	"rzsm_h14_grid_%DOMAIN.nc"
	""
	""
	""
	""
	""
	""
	""
	"rzsm_h26_grid_%DOMAIN.nc"
	""
)

data_folder_src_list=(
	"$HOME/datasets/source/h14/%YYYY/%MM/%DD/" 
	"$HOME/datasets/source/h14/%YYYY/%MM/%DD/"
	"$HOME/datasets/source/h27_h140/%YYYY/%MM/%DD/" 
	"$HOME/datasets/source/h27_h140/%YYYY/%MM/%DD/" 
	"$HOME/datasets/source/h141_h142/%YYYY/%MM/"
	"$HOME/datasets/source/h141_h142/%YYYY/%MM/%DD/"
	"$HOME/datasets/source/h141_h142/%YYYY/%MM/"
	"$HOME/datasets/source/h141_h142/%YYYY/%MM/%DD/"
	"$HOME/datasets/source/h26/%YYYY/%MM/%DD/"
	"$HOME/datasets/source/h26/%YYYY/%MM/%DD/"
)

data_folder_dst_list=(
	"$HOME/datasets/destination/h14/%YYYY/%MM/%DD/" 
	"$HOME/datasets/destination/h14/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h27_h140/%YYYY/%MM/%DD/" 
	"$HOME/datasets/destination/h27_h140/%YYYY/%MM/%DD/" 
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/"
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/"
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h26/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h26/%YYYY/%MM/%DD/"
)

grid_folder_list=(
	"$HOME/datasets/destination/h14/%YYYY/%MM/%DD/" 
	"$HOME/datasets/destination/h14/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h27_h140/%YYYY/%MM/%DD/" 
	"$HOME/datasets/destination/h27_h140/%YYYY/%MM/%DD/" 
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/"
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/"
	"$HOME/datasets/destination/h141_h142/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h26/%YYYY/%MM/%DD/"
	"$HOME/datasets/destination/h26/%YYYY/%MM/%DD/"
)

# Script time arg(s)
time_now=$(date '+%Y-%m-%d %H:00')
time_now='2023-08-03 00:00' # DEBUG
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# info script start
echo " ==================================================================================="
echo " ==> "${script_name}" (Version: "${script_version}" Release_Date: "${script_date}")"
echo " ==> START ..."


echo " ===> INFO TIME -- TIME: ${time_data_now}"

# iterate over data name(s)
for ((i=0;i<${#data_description_list[@]};++i)); do
	
	# ----------------------------------------------------------------------------------------
	# get data information
    data_folder_src_tmp="${data_folder_src_list[i]}" 
    data_file_src_tmp="${data_file_src_list[i]}" 
    
    data_file_unzip_tmp="${data_file_unzip_list[i]}" 
    data_file_grid_tmp="${data_file_grid_list[i]}" 
    data_file_conversion_tmp="${data_file_conversion_list[i]}"
    data_file_var_tmp="${data_file_var_list[i]}"  
    
	data_folder_dst_tmp="${data_folder_dst_list[i]}" 
	data_file_dst_tmp="${data_file_dst_list[i]}" 
	
	grid_folder_tmp="${grid_folder_list[i]}"  
	grid_file_tmp="${grid_file_list[i]}" 
	
	data_name_step="${data_name_list[i]}"
	data_description_step="${data_description_list[i]}"
	data_period_step="${data_period_list[i]}" 
	data_days_step="${data_days_list[i]}"   
	data_active_step="${data_active_list[i]}"
	
	data_type_compressed_step="${data_type_compressed_list[i]}"
	data_type_format_step="${data_type_format_list[i]}"
	
	apply_fx_make_regular_grid_step="${apply_fx_make_regular_grid_list[i]}"
	apply_fx_make_conversion_step="${apply_fx_make_conversion_list[i]}"
	apply_fx_make_var_step="${apply_fx_make_var_list[i]}"
	apply_fx_make_domain_step="${apply_fx_make_domain_list[i]}"
	apply_fx_make_file_grid_step="${apply_fx_make_file_grid_list[i]}"

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
			for day in $(seq 0 ${data_days_step}); do
				
				# ----------------------------------------------------------------------------------------
				# get time step
				time_data_step=$(date -d "$time_now ${day} days ago" +'%Y%m%d%H%M')
				year_data_step=${time_data_step:0:4}; 
				month_data_step=${time_data_step:4:2}; day_data_step=${time_data_step:6:2}
				hour_data_step='00'

				# info time analysis start
				echo " =====> TIME ANALYSIS: "${time_data_step:0:8}" ... "
				# ----------------------------------------------------------------------------------------
				
				# ----------------------------------------------------------------------------------------
				# define folder(s) and filename(s)
				data_folder_src_step=${data_folder_src_tmp/'%YYYY'/$year_data_step}
				data_folder_src_step=${data_folder_src_step/'%MM'/$month_data_step}
				data_folder_src_step=${data_folder_src_step/'%DD'/$day_data_step}
				data_folder_src_step=${data_folder_src_step/'%HH'/$hour_data_step}
				data_folder_src_step=${data_folder_src_step/'%DOMAIN'/$domain_name}
				
				data_file_src_step=${data_file_src_tmp/'%YYYY'/$year_data_step}
				data_file_src_step=${data_file_src_step/'%MM'/$month_data_step}
				data_file_src_step=${data_file_src_step/'%DD'/$day_data_step}
				data_file_src_step=${data_file_src_step/'%HH'/$hour_data_step}
				data_file_src_step=${data_file_src_step/'%DOMAIN'/$domain_name}
				
				data_file_unzip_step=${data_file_unzip_tmp/'%YYYY'/$year_data_step}
				data_file_unzip_step=${data_file_unzip_step/'%MM'/$month_data_step}
				data_file_unzip_step=${data_file_unzip_step/'%DD'/$day_data_step}
				data_file_unzip_step=${data_file_unzip_step/'%HH'/$hour_data_step}
				data_file_unzip_step=${data_file_unzip_step/'%DOMAIN'/$domain_name}
				
				data_file_grid_step=${data_file_grid_tmp/'%YYYY'/$year_data_step}
				data_file_grid_step=${data_file_grid_step/'%MM'/$month_data_step}
				data_file_grid_step=${data_file_grid_step/'%DD'/$day_data_step}
				data_file_grid_step=${data_file_grid_step/'%HH'/$hour_data_step}
				data_file_grid_step=${data_file_grid_step/'%DOMAIN'/$domain_name}
				
				data_file_conversion_step=${data_file_conversion_tmp/'%YYYY'/$year_data_step}
				data_file_conversion_step=${data_file_conversion_step/'%MM'/$month_data_step}
				data_file_conversion_step=${data_file_conversion_step/'%DD'/$day_data_step}
				data_file_conversion_step=${data_file_conversion_step/'%HH'/$hour_data_step}
				data_file_conversion_step=${data_file_conversion_step/'%DOMAIN'/$domain_name}
				
				data_file_var_step=${data_file_var_tmp/'%YYYY'/$year_data_step}
				data_file_var_step=${data_file_var_step/'%MM'/$month_data_step}
				data_file_var_step=${data_file_var_step/'%DD'/$day_data_step}
				data_file_var_step=${data_file_var_step/'%HH'/$hour_data_step}
				data_file_var_step=${data_file_var_step/'%DOMAIN'/$domain_name}
				
				data_folder_dst_step=${data_folder_dst_tmp/'%YYYY'/$year_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%MM'/$month_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%DD'/$day_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%HH'/$hour_data_step}
				data_folder_dst_step=${data_folder_dst_step/'%DOMAIN'/$domain_name}
				
				data_file_dst_step=${data_file_dst_tmp/'%YYYY'/$year_data_step}
				data_file_dst_step=${data_file_dst_step/'%MM'/$month_data_step}
				data_file_dst_step=${data_file_dst_step/'%DD'/$day_data_step}
				data_file_dst_step=${data_file_dst_step/'%HH'/$hour_data_step}
				data_file_dst_step=${data_file_dst_step/'%DOMAIN'/$domain_name}
				
				grid_folder_step=${grid_folder_tmp/'%YYYY'/$year_data_step}
				grid_folder_step=${grid_folder_step/'%MM'/$month_data_step}
				grid_folder_step=${grid_folder_step/'%DD'/$day_data_step}
				grid_folder_step=${grid_folder_step/'%HH'/$hour_data_step}
				grid_folder_step=${grid_folder_step/'%DOMAIN'/$domain_name}
				
				grid_file_step=${grid_file_tmp/'%YYYY'/$year_data_step}
				grid_file_step=${grid_file_step/'%MM'/$month_data_step}
				grid_file_step=${grid_file_step/'%DD'/$day_data_step}
				grid_file_step=${grid_file_step/'%HH'/$hour_data_step}
				grid_file_step=${grid_file_step/'%DOMAIN'/$domain_name}
				
				# define file path(s)
				data_path_src_step=${data_folder_src_step}${data_file_src_step}
				data_path_unzip_step=${data_folder_dst_step}${data_file_unzip_step}
				data_path_grid_step=${data_folder_dst_step}${data_file_grid_step}
				data_path_conversion_step=${data_folder_dst_step}${data_file_conversion_step}
				data_path_var_step=${data_folder_dst_step}${data_file_var_step}
				data_path_dst_step=${data_folder_dst_step}${data_file_dst_step}
				grid_path_step=${grid_folder_step}${grid_file_step}
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
					if [ -e ${data_path_dst_step} ]; then
						rm -rf ${data_path_dst_step}
					fi
					if [ -e ${grid_path_step} ]; then
						rm -rf ${grid_path_step}
					fi
				fi	
				
				# check file exist or not in the destination folder
				if [ -e ${data_path_dst_step} ]; then
					data_flag_cmp=false
				else
					data_flag_cmp=true
				fi
				# ----------------------------------------------------------------------------------------
				
				# ----------------------------------------------------------------------------------------
				# active uncompress file
				echo -n " ======> UNZIP FILE: ${data_path_src_step} TO ${data_path_unzip_step} ..."
				if [ "${data_type_compressed_step}" == "bz2" ] ; then
					
					if [ -e ${data_path_unzip_step} ]; then
						rm -rf ${data_path_unzip_step}
					fi
					
					if ! [ -e ${data_path_unzip_step} ]; then
						
						if [ -e ${data_path_unzip_step} ]; then
							rm -rf ${data_path_unzip_step}
						fi
					
						if bzip2 -dc ${data_path_src_step} > ${data_path_unzip_step}; then
						    echo " DONE!"
						else
						    echo " FAILED! ERROR IN COMMAND EXECUTION"
						fi
					else
						echo " SKIPPED! FILE PREVIOUSLY UNZIPPED"
					fi
					
					data_path_tmp_step=${data_path_unzip_step}
					
				elif [ "${data_type_compressed_step}" == "" ] ; then
					echo " SKIPPED. UNZIP FILE IS NOT ACTIVATED"
					data_path_tmp_step=${data_path_src_step}
				else
					echo " FAILED. COMPRESSED TYPE IS NOT SUPPORTED"
					exit 0
				fi
				# ----------------------------------------------------------------------------------------
				
				# ----------------------------------------------------------------------------------------
				# active make regulare grid
				echo -n " ======> MAKE REGULAR GRID TO FILE: ${data_path_tmp_step} to ${data_path_grid_step} ..."
				if [ "${apply_fx_make_regular_grid_step}" = true ] ; then 
					
					if [ -e ${data_path_grid_step} ]; then
						rm -rf ${data_path_grid_step}
					fi

				    if $exec_cdo -R copy $data_path_tmp_step $data_path_grid_step > /dev/null 2>&1; then
				        echo " DONE!"
				    else
				        echo " FAILED! ERROR IN COMMAND EXECUTION"
				    fi
					
				else
					echo " SKIPPED. MAKE REGULAR GRID IS NOT ACTIVATED"
					data_path_grid_step=${data_path_tmp_step}
					cp ${data_path_tmp_step} ${data_path_grid_step}
				fi
				# ----------------------------------------------------------------------------------------

				# ----------------------------------------------------------------------------------------
				# active make conversion file from grib to netcdf fulldisk
				echo -n " ======> MAKE CONVERSION FILE: ${data_path_grid_step} to ${data_path_conversion_step} ..."
				if [ "${apply_fx_make_conversion_step}" = true ] ; then 
					
					if [ -e ${data_path_conversion_step} ]; then
						rm -rf ${data_path_conversion_step}
					fi
					
				    if $exec_cdo -f nc copy $data_path_grid_step $data_path_conversion_step > /dev/null 2>&1; then
				        echo " DONE!"
				    else
				        echo " FAILED! ERROR IN COMMAND EXECUTION"
				    fi

				else
					echo " SKIPPED. FLAG IS NOT ACTIVATED"
					data_path_conversion_step=${data_path_grid_step}
					cp ${data_path_grid_step} ${data_path_conversion_step}
				fi
				# ----------------------------------------------------------------------------------------	
				
				# ----------------------------------------------------------------------------------------	
				# active make variable(s)
				echo " ======> COMPUTE VARS FROM ${data_path_conversion_step} TO ${data_path_var_step} ..."
				if [ "${apply_fx_make_var_step}" = true ] ; then 
					
					if [ -e ${data_path_var_step} ]; then
						rm -rf ${data_path_var_step}
					fi
				
					# Compute VAR_0_7 ==> Equation: float(7)/float(7)*H14_SM_V_L1
					echo -n " =======> VAR SOIL MOISTURE VAR 0_7 ..."
					if [ -e ${data_path_conversion_step} ]; then
						
						if ${exec_ncap2} -s "var_0_7=1*var40" -v -A $data_path_conversion_step $data_path_var_step > /dev/null 2>&1; then
							echo " DONE!"
						else
							echo " FAILED! ERROR IN COMMAND EXECUTION"
						fi

					else
						echo " FAILED! FILE $data_path_conversion_step DOES NOT EXIST! VARIABLE IS NOT COMPUTED!"
					fi
					
					# Compute VAR_0_28 ==> Equation: float(7)/float(28)*H14_SM_V_L1 + float(28-7)/float(28)*H14_SM_V_L2
					echo -n " =======> VAR SOIL MOISTURE VAR 0_28 ..."
					if [ -e ${data_path_conversion_step} ]; then
						
						if ${exec_ncap2} -s "var_0_28=0.25*var40 + 0.75*var41" -v -A $data_path_conversion_step $data_path_var_step > /dev/null 2>&1; then
							echo " DONE!"
						else
							echo " FAILED! ERROR IN COMMAND EXECUTION"
						fi

					else
						echo " FAILED! FILE $data_path_conversion_step DOES NOT EXIST! VARIABLE IS NOT COMPUTED!"
					fi
					
					# Compute VAR_0_100 ==> Equation: float(7)/float(100)*H14_SM_V_L1 + float(28-7)/float(100)*H14_SM_V_L2 + (float(100-7-(28-7))/float(100))*H14_SM_V_L3
					echo -n " =======> VAR SOIL MOISTURE VAR 0_100 ..."
					if [ -e ${data_path_conversion_step} ]; then
				
						if ${exec_ncap2} -s "var_0_100=0.07*var40 + 0.21*var41 + 0.72*var42" -v -A $data_path_conversion_step $data_path_var_step > /dev/null 2>&1; then
							echo " DONE!"
						else
							echo " FAILED! ERROR IN COMMAND EXECUTION"
						fi

					else
						echo " FAILED! FILE $data_path_conversion_step DOES NOT EXIST! VARIABLE IS NOT COMPUTED!"
					fi
					
					echo " ======> COMPUTE VARS FROM ${data_path_conversion_step} TO ${data_path_var_step} ... DONE"
									
				else
					echo " ======> COMPUTE VARS FROM ${data_path_conversion_step} TO ${data_path_var_step} ... SKIPPED. FLAG IS NOT ACTIVATED"
					data_path_var_step=${data_path_conversion_step}
					cp ${data_path_conversion_step} ${data_path_var_step}
				fi
				# ----------------------------------------------------------------------------------------	

				# ----------------------------------------------------------------------------------------	
				# active make variable(s) over domain
				echo -n " ======> COMPUTE VARS OVER DOMAIN: ${data_path_var_step} to ${data_path_dst_step} ..."
				if [ "${apply_fx_make_domain_step}" = true ] ; then 
					if ! [ -e ${data_path_dst_step} ]; then
				
						if ${exec_cdo} sellonlatbox,${domain_bb} ${data_path_var_step} ${data_path_dst_step} > /dev/null 2>&1; then
							echo " DONE!"
						else
							echo " FAILED! ERROR IN COMMAND EXECUTION"
						fi

					else
						echo " SKIPPED! FILE PREVIOUSLY CREATED!"
					fi
	
				else
					echo " SKIPPED. COMPUTE VARIABLE OVER DOMAIN IS NOT ACTIVATED"
					data_path_dst_step=${data_path_conversion_step}
					cp ${data_path_var_step} ${data_path_dst_step}
				fi
				# ----------------------------------------------------------------------------------------	
				
				# ----------------------------------------------------------------------------------------
				# active compute grid file
				echo -n " =====> COMPUTE GRID FILE: ${data_path_dst_step} to ${grid_path_step} ..."
				if [ "${apply_fx_make_file_grid_step}" = true ] ; then 
					if [ -e ${data_path_dst_step} ]; then
						if ! [ -e ${grid_path_step} ]; then
							
							if ${exec_ncks} -v lon,lat ${data_path_dst_step} ${grid_path_step} > /dev/null 2>&1; then
								echo " DONE!"
							else
								echo " FAILED! ERROR IN COMMAND EXECUTION"
							fi	
						else
							echo " SKIPPED! FILE PREVIOUSLY CREATED!"		
						fi
					else
						echo " FAILED! ${data_path_dst_step} DOES NOT EXIST!"
					fi
				else
					echo " SKIPPED. COMPUTE FILE GRID IS NOT ACTIVATED"
				fi
				# ----------------------------------------------------------------------------------------
				
				# ----------------------------------------------------------------------------------------
				# remove tmp file(s)
				if [ -e ${data_path_unzip_step} ]; then
					if [[ "$data_path_src_step" != "$data_path_unzip_step" ]]; then
						rm ${data_path_unzip_step} 
					fi
				fi
				
				if [ -e ${data_path_grid_step} ]; then
					rm ${data_path_grid_step} 
				fi
				
				if [ -e ${data_path_conversion_step} ]; then
					rm ${data_path_conversion_step} 
				fi
				
				if [ -e ${data_path_var_step} ]; then
					rm ${data_path_var_step} 
				fi
				# ----------------------------------------------------------------------------------------
				
				# ----------------------------------------------------------------------------------------
				# info time analysis end
				if [ -e ${data_path_dst_step} ]; then
					echo " =====> TIME ANALYSIS: "${time_data_step:0:8}" ... DONE"
				else
					echo " =====> TIME ANALYSIS: "${time_data_step:0:8}" ... FAILED"
				fi
				# ----------------------------------------------------------------------------------------
				
			done
			
			# info name end
			echo " ====> PRODUCT NAME '"${data_name_step}"' ... DONE"
			# ----------------------------------------------------------------------------------------

		else
			
			# ----------------------------------------------------------------------------------------
			# info name end
			echo " ====> PRODUCT NAME '"${data_name_step}"' ... SKIPPED. PRODUCT IS NOT ACTIVATED."
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

