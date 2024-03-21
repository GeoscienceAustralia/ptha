#
# Collate files and documentation that will be needed by end-users. 
#
# The script below collates the data of interest and also enforces better
# folder naming than used in the current computational codes.
#
# To change the datasets that are delivered, all the paths need to be updated and checked.
#
# If you need to modify and rerun, then I suggest first deleting output_dir (to
# start with a 'clean slate').
#
output_dir = 'NSW_tsunami_modelling_project_Draft_outputs_First_batch_scenarios_ID710.5'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
#
# * PDF slides
#   * Use the slides to provide high-level documentation.
slides_to_copy = '../../../../doc/meetings/2024_03_08/NSW_Notes_on_modelling_2024_03_08.pdf'
if(file.exists(slides_to_copy)){
    # They are on my local machine [and a different location on the GA drives,
    # and the NSWSES MS Teams folder, and were also circulated by email]
    file.copy(slides_to_copy, 
        paste0(output_dir, '/NSW_Notes_on_modelling_2024_03_08.pdf'), 
        overwrite=TRUE)
}
#
# * Head directory README.md and README.pdf
#   * Use to provide high-level documentation -- where things are
file.copy('template_data_package/README.md', 
    paste0(output_dir, '/README.md'), 
    overwrite=TRUE)
file.copy('template_data_package/README.pdf', 
    paste0(output_dir, '/README.pdf'), 
    overwrite=TRUE)

#
# * Domain shapefiles
file.copy('jatwc_to_inundation/domains_shapefile', 
    output_dir, 
    recursive=TRUE, overwrite=TRUE)
#
# * elevation_in_model/
file.copy('jatwc_to_inundation/elevation_in_model', 
    output_dir,
    recursive=TRUE, overwrite=TRUE)
#
# * JATWC_inundation_zones/
dir.create(paste0(output_dir, '/JATWC_inundation_zones'), showWarnings=FALSE, recursive=TRUE)
#   * README.md and README.pdf
file.copy('template_data_package/JATWC_inundation_zones/README.md', 
    paste0(output_dir, '/JATWC_inundation_zones/README.md'), 
    overwrite=TRUE)
file.copy('template_data_package/JATWC_inundation_zones/README.pdf', 
    paste0(output_dir, '/JATWC_inundation_zones/README.pdf'), 
    overwrite=TRUE)
#   * ATWS_COASTAL_ZONE_POLYGONS/
dir.create(paste0(output_dir, '/JATWC_inundation_zones/ATWS_COASTAL_ZONE_POLYGONS'), recursive=TRUE, showWarnings=FALSE)
file.copy('jatwc_to_inundation/ATWS_ZONES/ATWS_Zones_V2_2_4/',
    paste0(output_dir, '/JATWC_inundation_zones/ATWS_COASTAL_ZONE_POLYGONS'),
    recursive=TRUE, overwrite=TRUE)
#   * Inundaton_polygons_per_coastal_zone/
coastal_zone_names = Sys.glob('jatwc_to_inundation/Inundation_zones/*')
for(coastal_zone_name in coastal_zone_names){
#     * Eden Coast/
#       * no_threat_...polygon
#       * marine_warning_...polygon
#       * land_warning_...polygon
#       * no_threat_max_stage ...
#       * marine_warning_max_stage ...
#       * no_threat_arrival_time ...
#       * marine_warning_arrival_time ...
#       * land_warning_arrival_time ...
#     * [Other zones ....]
    coastal_zone_output_dir = paste0(output_dir, 
        '/JATWC_inundation_zones/Inundation_zones/', basename(coastal_zone_name), '/')
    dir.create(coastal_zone_output_dir, recursive=TRUE, showWarnings=FALSE)

    for(threat_level in c('no_threat', 'marine_warning', 'land_warning')){

        # Copy the inundation footprints (limited to 1/2500 @ 84th percentile)
        poly_result = paste0(coastal_zone_name, '/', basename(coastal_zone_name), '_', 
            threat_level, '_with-PTHA-exrate-limit_84pc_4e-04')
        file.copy(poly_result,
            paste0(coastal_zone_output_dir, '/'), 
            recursive=TRUE, overwrite=TRUE)

        # Copy the max-stage over any scenario (excluding land-warning, since
        # we haven't enforced any limitation of 1/2500 @ 84%)
        if(threat_level != 'land_warning'){
            max_stage_tifs = Sys.glob(paste0(coastal_zone_name, '/', threat_level, '_max_stage*.tif'))
            max_stage_tifs_output_dir = paste0(coastal_zone_output_dir, '/', 
                basename(coastal_zone_name), '_', threat_level, '_max_stage_tifs')
            dir.create(max_stage_tifs_output_dir, recursive=TRUE, showWarnings=FALSE)
            file.copy(max_stage_tifs, max_stage_tifs_output_dir, overwrite=TRUE)
        }

        # Copy the minimum arrival time over any scenario
        arrival_time_tifs = Sys.glob(paste0(coastal_zone_name, '/', threat_level, '_arrival_time*.tif'))
        arrival_time_tifs_output_dir = paste0(coastal_zone_output_dir, '/', 
            basename(coastal_zone_name), '_', threat_level, '_arrival_time_tifs')
        dir.create(arrival_time_tifs_output_dir, recursive=TRUE, showWarnings=FALSE)
        file.copy(arrival_time_tifs, arrival_time_tifs_output_dir, overwrite=TRUE)
    
    }
}

#
# * Probabilistic_inundation_hazard
pi_output_dir = paste0(output_dir, '/Probabilistic_inundation_hazard/')
dir.create(pi_output_dir, recursive=TRUE, showWarnings=FALSE)

#   * README.md and README.pdf
file.copy('template_data_package/Probabilistic_inundation_hazard/README.md',
    paste0(pi_output_dir, '/README.md'), overwrite=TRUE)
file.copy('template_data_package/Probabilistic_inundation_hazard/README.pdf',
    paste0(pi_output_dir, '/README.pdf'), overwrite=TRUE)

#   * exceedance_rate_1mm_depth
dir.create(paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth'),
    recursive=TRUE, showWarnings=FALSE)

#     * logic_tree_mean/
lgtm_depth_folder = paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth/logic_tree_mean/')
dir.create(lgtm_depth_folder, recursive=TRUE, showWarnings=FALSE)
lgtm_depth_files = Sys.glob('probabilistic_inundation/ptha18-NSW2023b-ID710.5-sealevel110cm/highres_depth_with_variance/ptha18-NSW2023b-ID710.5-sealevel110cm-depth-LogicTreeMean-sum_of_source_zones/summed_HS_domain_*_depth_as_max_stage_minus_elevation0_domain__exceedance_rate_with_threshold_0.001.tif')
file.copy(lgtm_depth_files, lgtm_depth_folder, overwrite=TRUE)

#     * 84th_percentile
eu84_depth_folder = paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth/84th_percentile/')
dir.create(eu84_depth_folder, recursive=TRUE, showWarnings=FALSE)
eu84_depth_files = Sys.glob('probabilistic_inundation/ptha18-NSW2023b-ID710.5-sealevel110cm/highres_depth_epistemic_uncertainty/84pc/ptha18-NSW2023b-ID710.5-sealevel110cm-depth_exrate_0.001_0.84_sum_of_source_zones/random_sum_of_source_zones_depth_rast_threshold_0.001_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
file.copy(eu84_depth_files, eu84_depth_folder, overwrite=TRUE)

#   * depth_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL
depth_84pc_folder = paste0(pi_output_dir, '/depth_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL/')
dir.create(depth_84pc_folder, recursive=TRUE, showWarnings=FALSE)
depth_84pc_files = Sys.glob('probabilistic_inundation/highres_domains_depth_at_epistemic_uncertainty_84pc_masked_at_elevation_below_zero/sum_of_all_sources_depth_where_elevation_exceeds_zero_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
file.copy(depth_84pc_files, depth_84pc_folder, overwrite=TRUE)

#   * depth_above_initial_condition_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL
depth_ai_84pc_folder = paste0(pi_output_dir, '/depth_above_initial_condition_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL/')
dir.create(depth_ai_84pc_folder, recursive=TRUE, showWarnings=FALSE)
depth_ai_84pc_files = Sys.glob('probabilistic_inundation/highres_domains_depth_above_initial_condition_at_epistemic_uncertainty_84pc/sum_of_all_sources_depth_above_initial_condition_where_elevation_exceeds_0_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
file.copy(depth_ai_84pc_files, depth_ai_84pc_folder, overwrite=TRUE)

#   * maximum_stage_with_exceedance_rate_1in2500_at_84th_percentile
max_stage_84pc_drymasked_folder = paste0(pi_output_dir, '/max_stage_with_exceedance_rate_1in2500_at_84th_percentile_masked_in_dry_areas/')
dir.create(max_stage_84pc_drymasked_folder, recursive=TRUE, showWarnings=FALSE)
max_stage_84pc_drymasked_files = Sys.glob('probabilistic_inundation/highres_domains_max_stage_at_epistemic_uncertainty_84pc_masked_in_dry_areas/sum_of_all_sources_max_stage_masked_in_dry_regions_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
file.copy(max_stage_84pc_drymasked_files, max_stage_84pc_drymasked_folder)
#
