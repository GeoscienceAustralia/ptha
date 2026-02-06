#
# The JATWC-style inundation zone polygons include areas where the 1/2500 84% tsunami was very small.
# Here we make new polygons which are the intersection of the JATWC-style polygons, and regions that have
# 1/2500 84% max-stage greater than 1.11 (=initial-sea-level + 0.01)
#

library(terra)

THRESHOLD_MAX_STAGE = 1.11 # (initial sea level + 1cm)

max_stage_1in2500_84pc = rast('../NSW_tsunami_modelling_project_final_outputs/Probabilistic_inundation_hazard/flow_variables_1in2500_at_84th_percentile/max_stage_with_exceedance_rate_1in2500_at_84th_percentile/max_stage_with_exceedance_rate_1in2500_84pc_masked_in_dry_areas.vrt')

land_warning_zone_shp = Sys.glob('../NSW_tsunami_modelling_project_final_outputs/JATWC_inundation_zones/Inundation_zones/*/*_land_warning_with-PTHA-exrate-limit_84pc_4e-04/*_land_warning_with-PTHA-exrate-limit_84pc_4e-04.shp')

land_warning_zone_names = basename(dirname(dirname(land_warning_zone_shp)))

for(i in 1:length(land_warning_zone_shp)){

    zone_name = land_warning_zone_names[i]
    zone_vect = vect(land_warning_zone_shp[i])
    zone_land_warning_folder_name = basename(dirname(land_warning_zone_shp[i]))

    zone_rast = crop(max_stage_1in2500_84pc, ext(zone_vect))

    thresh_rast = 1.0*(zone_rast > THRESHOLD_MAX_STAGE)

    thresh_rast[thresh_rast == 0] = NA

    thresh_rast_masked = mask(thresh_rast, zone_vect)

    thresh_rast_masked_vect = as.polygons(thresh_rast_masked)

    output_dir = paste0(zone_name, '/', zone_land_warning_folder_name, '_where_max_stage_exceeds_', THRESHOLD_MAX_STAGE)
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    writeVector(thresh_rast_masked_vect, 
        filename=paste0(output_dir, '/', basename(output_dir), '.shp'), 
        filetype='ESRI Shapefile', layer=basename(output_dir), 
        overwrite=TRUE)
}
