#!/bin/bash

#-----------------------------------------------------------------------------------------
# Script information
script_name="APP SYSTEM - CLEANER DATASETS DEPRECATED - REALTIME"
script_version="1.0.0"
script_date="2024/03/19"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get script information
script_file="app_system_cleaner_datasets_deprecated.sh"

# Get time information (-u to get gmt time)
time_script_now=$(date -u +"%Y-%m-%d 00:00")
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Folder of remote and local machine(s)
group_datasets_name=(
	"DATA - HSAF SNOW - H10"
	"DATA - HSAF SNOW - H12"
	"DATA - HSAF SNOW - H13"
	"DATA - SMAP - SLP2SMP_E"
	"DATA - SMAP - SLP3SMP_E"
	"PRODUCT - SM TC - POINTS"
	"PRODUCT - SM TC - MAPS"
	"PRODUCT - PLANETEK"
)

group_folder_datasets=(
	"/share/HSAF_SNOW/ancillary/h10/" 
	"/share/HSAF_SNOW/ancillary/h12/" 
	"/share/HSAF_SNOW/ancillary/h13/"
	"/share/SMAP/ancillary/spl2smp_e/"
	"/share/SMAP/ancillary/spl3smp_e/"
	"/share/SM_TC/ancillary/points/"
    "/share/SM_TC/ancillary/maps/" 
    "/share/PLANETEK/"
)

group_file_datasets_clean=(
	true
	true
	true
	true
	true
	true
	true
	true
)

group_file_datasets_elapsed_days=(
	5
	5
	5
	10
	10
	15
	15
	15
)
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."
echo " ===> EXECUTION ..."

time_script_now=$(date -d "$time_script_now" +'%Y-%m-%d 00:00')
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Iterate over tags
for datasets_id in "${!group_datasets_name[@]}"; do

	# ----------------------------------------------------------------------------------------
	# Get values of tag(s) and folder(s)        
	datasets_name=${group_datasets_name[datasets_id]}
	
	folder_datasets=${group_folder_datasets[datasets_id]}

	file_datasets_clean=${group_file_datasets_clean[datasets_id]} 
	file_datasets_elapsed_days=${group_file_datasets_elapsed_days[datasets_id]}
	
	# Info datasets type start
	echo " ====> DATASETS TYPE ${datasets_name} ... "
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------
	# Check sync activation
	if ${file_datasets_clean} ; then
		
		# ----------------------------------------------------------------------------------------
		# Iterate over filename
		for file_datasets_name in $(find ${folder_datasets} -type f -mtime +${file_datasets_elapsed_days}); do
			echo " ====> DELETE FILENAME ${file_datasets_name} ... "
			
			if [ -f "$file_datasets_name" ] ; then
    			rm "$file_datasets_name"
    			echo " ====> DELETE FILENAME ${file_datasets_name} ... DONE"
			else
				echo " ====> DELETE FILENAME ${file_datasets_name} ... FAILED. FILE NOT FOUND"
			fi
			
		done
		# ----------------------------------------------------------------------------------------
		
		# ----------------------------------------------------------------------------------------
		# Find empty folders
		for folder_empty_name in $(find ${folder_datasets} -type d -empty); do
			
			echo " ====> DELETE EMPTY FOLDER ${folder_empty_name} ... "
			if [ -d "$folder_empty_name" ] ; then
				rmdir ${folder_empty_name} -vp --ignore-fail-on-non-empty {} 
				echo " ====> DELETE EMPTY FOLDER ${file_datasets_name} ... DONE"
			else
				echo " ====> DELETE EMPTY FOLDER ${file_datasets_name} ... FAILED. FOLDER NOT FOUND"
			fi
			
		done
		# ----------------------------------------------------------------------------------------
		
		# ----------------------------------------------------------------------------------------
		# Info datasets type end
		echo " ====> DATASETS TYPE ${datasets_name} ... DONE"
		# ----------------------------------------------------------------------------------------
		
	else
	
		# ----------------------------------------------------------------------------------------
		# Info tag end (not activated)
		echo " ====> DATASETS TYPE ${datasets_name} ... SKIPPED. SYNC NOT ACTIVATED"
		# ----------------------------------------------------------------------------------------
		
	fi
	# ----------------------------------------------------------------------------------------
	
done

# Info script end
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

