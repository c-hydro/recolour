#!/bin/bash -e

# ----------------------------------------------------------------------------------------
# Script information
script_name='ECMWF H27 FILES CONVERTER [GRIB 2 NC]'
script_version="1.0.0"
script_date='2019/03/12'

# Time Informaation
date_start="2007-01-01"
date_end="2016-12-31"

# Filename(s) information
file_in_raw="h27_%YYYY%MM%DD00_T1279.grib"
file_tmp_raw="h27_%YYYY%MM%DD_00_r.grib"
file_out_raw="rzsm_%YYYY%MM%DD00_dr.nc"

# Folder(s) information
folder_in_raw="/share/HSAF_Data_Validation/datasets/hsaf_mod/h27_grib/%YYYY/%MM/%DD/"
folder_tmp_raw="/share/HSAF_Data_Validation/datasets/hsaf_mod/h27_tmp/%YYYY/%MM/%DD/"
folder_out_raw="/share/HSAF_Data_Validation/datasets/hsaf_mod/h27_nc/%YYYY/%MM/%DD/"

# Cdo executable
cdo_exec="/share/HSAF_Data_Validation/apps/cdo-1.7.2rc3_NC-4.1.2_HDF-1.8.17/bin/cdo"
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
	
	file_in=${file_in_raw/'%YYYY'/$year_get}
	file_in=${file_in/'%MM'/$month_get}
	file_in=${file_in/'%DD'/$day_get}
	
	file_tmp=${file_tmp_raw/'%YYYY'/$year_get}
	file_tmp=${file_tmp/'%MM'/$month_get}
	file_tmp=${file_tmp/'%DD'/$day_get}
	
	file_out=${file_out_raw/'%YYYY'/$year_get}
	file_out=${file_out/'%MM'/$month_get}
	file_out=${file_out/'%DD'/$day_get}

	folder_in=${folder_in_raw/'%YYYY'/$year_get}
	folder_in=${folder_in/'%MM'/$month_get}
	folder_in=${folder_in/'%DD'/$day_get}
	
	folder_tmp=${folder_tmp_raw/'%YYYY'/$year_get}
	folder_tmp=${folder_tmp/'%MM'/$month_get}
	folder_tmp=${folder_tmp/'%DD'/$day_get}
	
	folder_out=${folder_out_raw/'%YYYY'/$year_get}
	folder_out=${folder_out/'%MM'/$month_get}
	folder_out=${folder_out/'%DD'/$day_get}
	# ----------------------------------------------------------------------------------------
    
	# ----------------------------------------------------------------------------------------	
	# Create folder(s)
	if [ ! -d "$folder_tmp" ]; then
		mkdir -p $folder_tmp
	fi
	if [ ! -d "$folder_out" ]; then
		mkdir -p $folder_out
	fi
	# ----------------------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------------------
    # Define path(s)
    path_file_in=$folder_in$file_in
    path_file_tmp=$folder_tmp$file_tmp
    path_file_out=$folder_out$file_out
    # ----------------------------------------------------------------------------------------
    
	# ----------------------------------------------------------------------------------------
	# Execute cdo command(s)
	echo ""
	echo " ===> STEP :: DATE: "$date_step " ===> START "
	echo " ====> FILE IN: "$path_file_in 
	echo " ====> FILE TMP: "$path_file_tmp
	echo " ====> FILE OUT: "$path_file_out
    
    if [ -e $path_file_in ]
    then
        echo " =====> FILE $path_file_in AVAILABLE"
        echo " ======> FILE GRID CONVERSION (SET REGULAR GRID) ... "
        $cdo_exec -R copy $path_file_in $path_file_tmp
        echo " ======> FILE GRID CONVERSION (SET REGULAR GRID) ... DONE"
        
        echo " ======> FILE FORMAT CONVERSION (GRIB 2 NC)... "
        $cdo_exec -f nc copy $path_file_tmp $path_file_out
        echo " ======> FILE FORMAT CONVERSION (GRIB 2 NC)... DONE"
    
    else
        echo " =====> FILE $path_file_in NOT AVAILABLE"
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
    
    
    
    
    
