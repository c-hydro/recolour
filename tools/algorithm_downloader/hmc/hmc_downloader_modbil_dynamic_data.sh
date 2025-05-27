#!/bin/bash

#-----------------------------------------------------------------------------------------
# Script information
script_name='ARPAL - MODBIL - DOWNLOADER DYNAMIC DATA - HISTORY'
script_version="1.1.0"
script_date='2025/05/20'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Script flag(s)
flag_download=true
flag_update=false
flag_list_basin=1   # Choose 1, 2, 3, or 4

# Date range (inclusive)
start_date="2021-07-01"
end_date="2022-11-01"

# Unified domain list and corresponding group flags
domain_list_all=(
    "Tanaro" "Savonese"
    "Ponente" "Imperiese" "Finalese" "OrbaStura"
    "CentroPonente" "Centa" "LevanteGenovese" "Erro"
    "BormidaS" "BormidaM" "AvetoTrebbia"
    "PonenteGenovese" "Magra" "Scrivia"
    "Entella"
    "Roia"
)

domain_flag_list=(
    1 1 						# Group 1
    1 1 1 1 					# Group 1
    1 1 1 1 					# Group 1
    1 1 1   					# Group 1
    2 2 2                       # Group 2
    3                           # Group 3
    4                           # Group 4
)

#-----------------------------------------------------------------------------------------
# Extract domain_list and subfolder_sm based on selected flag
domain_list=()
subfolder=""

for i in "${!domain_list_all[@]}"; do
    if [[ "${domain_flag_list[$i]}" -eq "$flag_list_basin" ]]; then
        domain_list+=("${domain_list_all[$i]}")
    fi
done

case "$flag_list_basin" in
    1) subfolder="Maps" ;;
    2) subfolder="Maps" ;;        # Could also be "Maps" (2020-2022) "MapsPlants" (2013–2020)
    3) subfolder="Maps" ;;        # Could also be "Maps" (2021-2022) "Maps052020" (2013-2021)
    4) subfolder="Maps" ;;
    *) echo " ==> DOMAIN BAD CHOICE (1,2,3,4)"; exit 1 ;;
esac

# Template path and filename strings
template_remote_base="/home/cfmi.arpal.org/satsuolo/MODELLO_BILANCIO/%DOMAIN_NAME%Domain/Results/${subfolder}/"
template_local_path="/home/cfmi.arpal.org/satsuolo/QNAPDEV_SaturazioneSuolo/datasets/hmc/domains/modbil/%DOMAIN_NAME%Domain/%Y%/%m%/"
template_filename="%DOMAIN_NAME%DomainV_%YMD%2300.gz"
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
# Main logic
if [ "$flag_download" = true ]; then

    current_date="$start_date"
    while [[ "$current_date" != $(date -I -d "$end_date + 1 day") ]]; do

        year=$(date -d "$current_date" +%Y)
        month=$(date -d "$current_date" +%m)
        ymd=$(date -d "$current_date" +%Y%m%d)

        echo "[$script_name] Processing date: $ymd"

        for domain in "${domain_list[@]}"; do

            local_path="${template_local_path//%DOMAIN_NAME%/$domain}"
            local_path="${local_path//%Y%/$year}"
            local_path="${local_path//%m%/$month}"

            filename="${template_filename//%DOMAIN_NAME%/$domain}"
            filename="${filename//%YMD%/$ymd}"

            remote_path="${template_remote_base//%DOMAIN_NAME%/$domain}"
            remote_file="${remote_path}${filename}"
            local_file="${local_path}${filename}"

            echo "  Domain: $domain"
            echo "  Remote file: $remote_file"
            echo "  Local file:  $local_file"

            if [ ! -f "$remote_file" ]; then
                echo "  WARNING: Remote file not found: $remote_file"
                continue
            fi

            if [ -f "$local_file" ] && [ "$flag_update" = false ]; then
                echo "  SKIP: Local file already exists and update is disabled: $local_file"
                continue
            fi

            if [ ! -d "$local_path" ]; then
                echo "  Creating local directory: $local_path"
                mkdir -p "$local_path"
            fi

            echo "  Copying file ... "
            cp -f "$remote_file" "$local_file"
            echo "  Copying file ... DONE"
        done

        current_date=$(date -I -d "$current_date + 1 day")

    done

else
    echo "[$script_name] Download skipped. Set 'flag_download=true' to enable downloading."
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Info script end
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
#-----------------------------------------------------------------------------------------

