#!/bin/bash

DOMAIN_STRING=(
	'AvetoTrebbiaDomain' 'BormidaMDomain' 'BormidaSDomain' 
	'CentaDomain' 'CentroPonenteDomain' 'EntellaDomain' 'ErroDomain' 
	'FinaleseDomain' 'ImperieseDomain' 'LevanteGenoveseDomain' 
	'MagraDomain' 'OrbaSturaDomain' 'PonenteDomain' 'PonenteGenoveseDomain'
	'RoiaDomain' 'SavoneseDomain' 'ScriviaDomain' 'TanaroDomain')
TIME_STRING=('2025/04/16 12:00')

TIME_FORMAT_FOLDER='%Y/%m/%d'
TIME_FORMAT_FILE='%Y%m%d%H00'

SRC_FOLDER="/home/cfmi.arpal.org/satsuolo/MODELLO_CONTINUUM/archive"
DEST_FOLDER="/home/cfmi.arpal.org/satsuolo/hmc_grid_selected"

for domain in "${DOMAIN_STRING[@]}"; do
  for time in "${TIME_STRING[@]}"; do
    folder=$(date -d "$time" +"$TIME_FORMAT_FOLDER")
    file=$(date -d "$time" +"$TIME_FORMAT_FILE")
    
    echo $folder $file
    
    SRC_FOLDER_DEF="$SRC_FOLDER/$domain/weather_stations_realtime/$folder"
    DEST_FOLDER_DEF="$DEST_FOLDER/$domain/weather_stations_realtime/$folder"
    
    SRC_FILE_DEF="hmc.output-grid.${file}.nc.gz"
    
    echo $SRC_FOLDER_DEF
    echo $SRC_FILE_DEF
    echo $DEST_FOLDER_DEF
    
    mkdir -p ${DEST_FOLDER_DEF}
    
    find $SRC_FOLDER_DEF -type f -name $SRC_FILE_DEF -exec cp {} $DEST_FOLDER_DEF \;
  done
done
