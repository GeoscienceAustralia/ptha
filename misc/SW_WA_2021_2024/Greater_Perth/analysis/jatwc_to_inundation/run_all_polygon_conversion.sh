for warning_type in 'no_threat' 'marine_warning' 'land_warning' 'minor_land_warning' 'major_land_warning'; do
    Rscript convert_raster_zones_to_polygons_with_PTHA_limits.R 'Perth Coast' $warning_type;
done

