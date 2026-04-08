#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -lmem=190GB
#PBS -lncpus=48
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source R_431_NCI_modules.sh
#all_warning_zones=('Geraldton Coast' 'Gascoyne Coast' 'Ningaloo Coast' 'Pilbara Coast West')

## Geraldton already run
#all_warning_zones=('Gascoyne Coast' 'Ningaloo Coast' 'Pilbara Coast West')

#Rscript copy_elevation_rasters_locally.R

all_warning_zones=('Ningaloo Coast')

for warning_zone in "${all_warning_zones[@]}"; do 
    # Get JATWC H parameters
    Rscript compute_scenario_statistics_in_zone.R "$warning_zone"
    # Make temporary rasters
    Rscript map_threat_levels_in_zone.R "$warning_zone"
    # Convert rasters to polygons
    for warning_type in 'no_threat' 'marine_warning' 'land_warning' 'minor_land_warning' 'major_land_warning'; do
        Rscript convert_raster_zones_to_polygons_with_PTHA_limits.R "$warning_zone" $warning_type;
        done;
    # Make elevation contours
    Rscript make_elevation_contours.R "$warning_zone"
    done
