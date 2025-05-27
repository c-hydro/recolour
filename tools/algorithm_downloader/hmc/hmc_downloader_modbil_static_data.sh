#!/bin/bash

#-----------------------------------------------------------------------------------------
# Script information
script_name='ARPAL - MODBIL - DOWNLOADER STATIC DATA - HISTORY'
script_version="1.1.0"
script_date='2025/05/20'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Script flag
flag_download=true

# domain list
domain_list=(
    "Tanaro" "Savonese"
    "Ponente" "Imperiese" "Finalese" "OrbaStura" 
    "CentroPonente" "Centa"
    "LevanteGenovese" "Erro" 
    "BormidaS" "BormidaM" "AvetoTrebbia"
    "PonenteGenovese" "Magra" "Scrivia"
    "Entella"
    "Roia"
)

# Static path templates with placeholder
remote_template="/home/cfmi.arpal.org/satsuolo/MODELLO_BILANCIO/%DOMAIN_NAME%Domain/LandData/"
local_template="/home/cfmi.arpal.org/satsuolo/QNAPDEV_SaturazioneSuolo/datasets/hmc/geo/domains_modbil/%DOMAIN_NAME%Domain/"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Begin processing if flag is enabled
if [ "$flag_download" = true ]; then
    for domain in "${domain_list[@]}"; do
        # Replace placeholder with actual domain name
        remote_path="${remote_template//%DOMAIN_NAME%/$domain}"
        local_path="${local_template//%DOMAIN_NAME%/$domain}"

        echo "[$script_name] Processing domain: ${domain}"
        echo "Remote path: ${remote_path}"
        echo "Local path:  ${local_path}"

        # Check if remote path exists
        if [ ! -d "$remote_path" ]; then
            echo "[$script_name] WARNING: Remote path does not exist: $remote_path"
            continue
        fi

        # Create local directory only if it does not exist
        if [ ! -d "$local_path" ]; then
            echo "Creating local directory: $local_path"
            mkdir -p "$local_path"
        fi

        # Perform the rsync
        rsync -avh --progress "$remote_path" "$local_path"
    done
else
    echo "[$script_name] Download skipped. Set 'flag_download=true' to enable syncing."
fi
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Info script end
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
#-----------------------------------------------------------------------------------------

