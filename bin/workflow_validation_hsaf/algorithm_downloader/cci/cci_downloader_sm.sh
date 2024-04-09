#!/bin/bash -e

# ----------------------------------------------------------------------------------------
# Example CCI product(s) downloader
# sftp esacci_sm_v042@ftp.geo.tuwien.ac.at:/_down/alldata_7zip_compressed/ESACCI-SOILMOISTURE-L3S-SSMV-PASSIVE_1978-2016-v04.2.zip /home/gabellani/HSAF_Data_Validation/datasets/esacci/

# NOTA BENE:
# La prima volta che ci colleghiamo ad un ftp bisogna comunque validare la sorgente dei dati e di conseguenza il downloader fallisce.
# Fare un download manuale dei dati ci consente di autenticare la authenticate key e di poter poi usare la procedura automatica.
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Script information
script_name='ESACCI FILES DOWNLOADER'
script_version="1.0.0"
script_date='2018/01/23'

# Time Informaation
date_start="2007-01-01"
date_end="2022-12-31"

# FTP information
site_file_raw="ESACCI-SOILMOISTURE-L3S-SSMV-PASSIVE-%YYYY%MM%DD000000-fv08.1.nc"
site_folder_raw="https://data.ceda.ac.uk/neodc/esacci/soil_moisture/data/daily_files/PASSIVE/v08.1/%YYYY/"

# Path(s) information
data_folder_raw="/share/CCI/esacci_soilmoisture_l3s_ssmv_passive_v08.1_grid/%YYYY/"
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
	date_get=$(date -u -d $date_step +'%Y%m%d')

	year_get=$(date -u -d "$date_step" +"%Y")
	month_get=$(date -u -d "$date_step" +"%m")
	day_get=$(date -u -d "$date_step" +"%d")
	
	site_file_get=${site_file_raw/'%YYYY'/$year_get}
	site_file_get=${site_file_get/'%MM'/$month_get}
	site_file_get=${site_file_get/'%DD'/$day_get}

	site_folder_get=${site_folder_raw/'%YYYY'/$year_get}
	site_folder_get=${site_folder_get/'%MM'/$month_get}
	site_folder_get=${site_folder_get/'%DD'/$day_get}

	data_folder_save=${data_folder_raw/'%YYYY'/$year_get}
	data_folder_save=${data_folder_save/'%MM'/$month_get}
	data_folder_save=${data_folder_save/'%DD'/$day_get}
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------	
	# Create folder(s)
	if [ ! -d "$data_folder_save" ]; then
		mkdir -p $data_folder_save
	fi
	# ----------------------------------------------------------------------------------------

	# ----------------------------------------------------------------------------------------
	# Execute wget command
	echo ""
	echo " ===> STEP :: DATE: "$date_step " ===> START "
	echo " ====> FILE: "$site_file_get 
	echo " ====> SITE FOLDER: "$site_folder_get
	echo " ====> DATA FOLDER: "$data_folder_save
	
	file_check=$data_folder_save$site_file_get

	if [ -f "$file_check" ]
	then
		echo " =====> DOWNLOAD SKIPPED :: $file_check found in data folder"
	else
		echo " =====> DOWNLOAD ACTIVATED :: $file_check not found in data folder"

		# create command-line
		exec_command_step="wget --load-cookies ~/.urs_cookies  --save-cookies ~/.urs_cookies --keep-session-cookies -r -c -nH -nd -np -A nc -P ${data_folder_save} --recursive --content-disposition ${site_folder_get}${site_file_get} " 

		# download data
		echo " ======> COMMAND LINE: $exec_command_step"
		echo -n " ======> DOWNLOAD ..."
	   
	    if ${exec_command_step} > /dev/null 2>&1; then
            echo " DONE!"
        else
            echo " FAILED! Error in command execution!"
            exit
        fi
		
	fi

	echo " ===> STEP :: DATE: "$date_step " ===> END "
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




