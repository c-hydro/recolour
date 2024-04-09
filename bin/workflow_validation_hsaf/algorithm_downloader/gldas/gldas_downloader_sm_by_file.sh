#!/bin/bash

# ----------------------------------------------------------------------------------------
# Script information
script_name='GLDAS FILES DOWNLOADER'
script_version="1.0.0"
script_date='2018/01/23'

# Time Informaation
date_start="2021-01-01"
date_end="2021-05-31"

# FTP information
ftp_file="GLDAS_NOAH025_3H.A%YYYY%MM%DD.*.021.nc4"
ftp_folder="https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH025_3H.2.1/"
ftp_site="urs.earthdata.nasa.gov"
ftp_user="fabiodelogu"
ftp_password="Cuba1978"

# Path(s) information
path_local_dataset="/share/HSAF_Data_Validation/datasets/gldas_noah025_3h_v2.1_grid/"
path_remote_dataset="hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH025_3H.2.1/"
path_temp_dataset="/share/HSAF_Data_Validation/datasets/temp/"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Prepare raw information
file_raw=$ftp_file
ftp_raw=$ftp_folder"%YYYY/%DOY/"

folder_remote_raw=$path_remote_dataset"%YYYY/%DOY/"
folder_temp_raw=$path_temp_dataset
folder_local_raw=$path_local_dataset"%YYYY/%DOY/"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Attach remote site
touch .netrc
echo "machine "$ftp_site" login "$ftp_user" password "$ftp_password >> .netrc
chmod 0600 .netrc
touch .urs_cookies
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Cycle(s) over date(s)
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."

date_step=$date_start
date_end=$(date -I -d "$date_end + 1 day")
while [ "$date_step" != "$date_end" ]; do 
	
	# ----------------------------------------------------------------------------------------
	# Define get and save information
	date_get=$(date -u -d $date_step +%Y%m%d)
	doy_get=$(date -u -d $date_step +%j)

	year_get=$(date -u -d "$date_step" +"%Y")
	month_get=$(date -u -d "$date_step" +"%m")
	day_get=$(date -u -d "$date_step" +"%d")
	
	file_get=${file_raw/'%YYYY'/$year_get}
	file_get=${file_get/'%MM'/$month_get}
	file_get=${file_get/'%DD'/$day_get}
	file_get=${file_get/'%DOY'/$doy_get}

	ftp_get=${ftp_raw/'%YYYY'/$year_get}
	ftp_get=${ftp_get/'%MM'/$month_get}
	ftp_get=${ftp_get/'%DD'/$day_get}
	ftp_get=${ftp_get/'%DOY'/$doy_get}

	folder_remote_get=${folder_remote_raw/'%YYYY'/$year_get}
	folder_remote_get=${folder_remote_get/'%MM'/$month_get}
	folder_remote_get=${folder_remote_get/'%DD'/$day_get}
	folder_remote_get=${folder_remote_get/'%DOY'/$doy_get}

	folder_temp_get=${folder_temp_raw/'%YYYY'/$year_get}
	folder_temp_get=${folder_temp_get/'%MM'/$month_get}
	folder_temp_get=${folder_temp_get/'%DD'/$day_get}
	folder_temp_get=${folder_temp_get/'%DOY'/$doy_get}
	
	folder_local_get=$folder_temp_get$folder_remote_get

	folder_local_save=${folder_local_raw/'%YYYY'/$year_get}
	folder_local_save=${folder_local_save/'%MM'/$month_get}
	folder_local_save=${folder_local_save/'%DD'/$day_get}
	folder_local_save=${folder_local_save/'%DOY'/$doy_get}
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------	
	# Create folder(s)
	if [ ! -d "$folder_local_get" ]; then
		mkdir -p $folder_local_get
	fi

	if [ ! -d "$folder_local_save" ]; then
		mkdir -p $folder_local_save
	fi
	# ----------------------------------------------------------------------------------------

	# ----------------------------------------------------------------------------------------
	# Execute wget command
	echo ""
	echo " ===> STEP :: DATE: "$date_step " DOY: "$doy_get " ===> START "
	echo " ====> FILE: "$file_get 
	echo " ====> FTP: "$ftp_get
	echo " ====> DATASET REMOTE FOLDER: "$folder_remote_get
	echo " ====> DATASET TEMP FOLDER: "$folder_local_get
	echo " ====> DATASET LOCAL FOLDER: "$folder_local_save
	
	wget -P $folder_temp_get --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off --no-parent -A $file_get $ftp_get
	
	mv $folder_local_get$file_get $folder_local_save

	rm -r $path_temp_dataset

	echo " ===> STEP :: DATE: "$date_step " DOY: "$doy_get " ===> END "
	echo ""
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------
	# Update date step
	date_step=$(date -I -u -d "$date_step + 1 day")
	# ----------------------------------------------------------------------------------------

done

echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------


