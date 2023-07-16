# Automate running in an interactive environment, assuming we've already got the elevation rasters
# and don't want to make contours

source ~/R_421_NCI_modules.sh

Rscript compute_scenario_statistics_in_zone.R 'Bunbury Geographe Coast'

Rscript map_threat_levels_in_zone.R 'Bunbury Geographe Coast'

for warning_type in 'no_threat' 'marine_warning' 'land_warning' 'minor_land_warning' \
    'major_land_warning'; do
    Rscript convert_raster_zones_to_polygons_with_PTHA_limits.R 'Bunbury Geographe Coast' $warning_type;
    done;

Rscript compute_max_depths_for_marine_warning_scenarios.R "Bunbury Geographe Coast"
