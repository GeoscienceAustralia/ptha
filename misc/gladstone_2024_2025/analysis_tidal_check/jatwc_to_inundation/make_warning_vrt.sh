#!/bin/bash
# ls -d Inundation_zones/*/
dirs=$(ls -d Inundation_zones/*/)

for dir in $dirs
do
    cd $dir
    gdalbuildvrt -resolution highest land_warning_arrival_time.vrt land_warning_arrival_time*.tif
    gdalbuildvrt -resolution highest land_warning_max_stage.vrt land_warning_max_stage*.tif
    gdalbuildvrt -resolution highest major_land_warning_arrival_time.vrt major_land_warning_arrival_time*.tif
    gdalbuildvrt -resolution highest major_land_warning_max_stage.vrt major_land_warning_max_stage*.tif
    gdalbuildvrt -resolution highest marine_warning_arrival_time.vrt marine_warning_arrival_time*.tif
    gdalbuildvrt -resolution highest marine_warning_max_stage.vrt marine_warning_max_stage*.tif
    gdalbuildvrt -resolution highest minor_land_warning_arrival_time.vrt minor_land_warning_arrival_time*.tif
    gdalbuildvrt -resolution highest minor_land_warning_max_stage.vrt minor_land_warning_max_stage*.tif
    gdalbuildvrt -resolution highest no_threat_warning_arrival_time.vrt no_threat_warning_arrival_time*.tif
    gdalbuildvrt -resolution highest no_threat_warning_max_stage.vrt no_threat_warning_max_stage*.tif
    cd ..
done
cd ..
