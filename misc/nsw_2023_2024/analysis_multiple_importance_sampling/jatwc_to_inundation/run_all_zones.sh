#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=05:00:00
#PBS -lmem=500GB
#PBS -lncpus=104
#PBS -ljobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85+gdata/fj6


# Convenience script to compute JATWC-style zones for all ATWS coastal zones in NSW,
# here including Lord Howe Island and also Norfolk Island.

source R_431_NCI_modules.sh

ATWS_ZONES_ALL=(\
"Eden Coast" \
"Batemans Coast" \
"Illawarra Coast" \
"Sydney Coast" \
"Hunter Coast" \
"Macquarie Coast" \
"Coffs Coast" \
"Byron Coast" \
"Lord Howe Island" \
"Norfolk Island" \
)

# Do the calculations for each zone.
for ATWS_ZONE in "${ATWS_ZONES_ALL[@]}"; do
    Rscript compute_scenario_statistics_in_zone.R "$ATWS_ZONE" ;
    Rscript map_threat_levels_in_zone.R "$ATWS_ZONE" ;
    # Make outputs for each warning type
    for warning_type in 'no_threat' 'marine_warning' 'land_warning' 'minor_land_warning' 'major_land_warning'; do
        Rscript convert_raster_zones_to_polygons_with_PTHA_limits.R "$ATWS_ZONE" "$warning_type" ;
        done ;
    done ;

