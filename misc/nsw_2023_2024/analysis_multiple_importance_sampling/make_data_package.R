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

# Function often needed for folders with tif files
copy_tifs_and_make_vrt<-function(tifs_to_copy, output_folder, output_vrt_basename){
    
    # Check we found some files
    failed = (sum(file.exists(tifs_to_copy)) == 0)
    if(failed){
        errmes = paste0('Could not find tifs ', tifs_to_copy)
        stop(errmes)
    }

    # Copy the tifs
    dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)
    copied = file.copy(tifs_to_copy, output_folder, overwrite=TRUE)
    if(any(!copied)) print(paste0('failed to copy ', tifs_to_copy[!copied]))

    # Make a vrt
    MYDIR = getwd()    
    on.exit(setwd(MYDIR)) # Ensure we end up in MYDIR
    setwd(output_folder)
    gdal_command = paste0('gdalbuildvrt -resolution highest ', output_vrt_basename, ' *.tif')
    errcode = system(gdal_command)
    return(errcode)
}

output_dir = 'NSW_tsunami_modelling_project_final_outputs_20241106'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
#
# * PDF slides
#   * Use the slides to provide high-level documentation.
slides_to_copy = '../../../../doc/meetings/2024_11_06/NSW_Notes_on_modelling_2024_11_06.pdf'
if(file.exists(slides_to_copy)){
    # They are on my local machine [and a different location on the GA drives,
    # and the NSWSES MS Teams folder, and were also circulated by email]
    file.copy(slides_to_copy, 
        paste0(output_dir, basename(slides_to_copy)), 
        overwrite=TRUE)
}
#
# * Head directory README.md and README.pdf and possibly other pdf files
#   * Use to provide high-level documentation -- where things are
file.copy('template_data_package/README.md',
    paste0(output_dir, '/README.md'),
    overwrite=TRUE)
pdf_files = Sys.glob('template_data_package/*.pdf')
file.copy(pdf_files,
    paste0(output_dir, '/', basename(pdf_files)),
    overwrite=TRUE)

#
# * Domain shapefiles -- same for every batch of scenarios, just use one
file.copy('../analysis_scenarios_ID710.5/jatwc_to_inundation/domains_shapefile', 
    output_dir, 
    recursive=TRUE, overwrite=TRUE)
#
# * elevation_in_model/ -- same for every batch of scenarios, just use one
file.copy('../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model', 
    output_dir,
    recursive=TRUE, overwrite=TRUE)

#
# * JATWC_inundation_zones/
#

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
            max_stage_tifs_vrt_name = paste0(basename(coastal_zone_name), '_all_', threat_level, '_max_stage.vrt')
            copy_tifs_and_make_vrt(max_stage_tifs, max_stage_tifs_output_dir, max_stage_tifs_vrt_name)
        }

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

#    * Image of flood hazard categories
file.copy('probabilistic_inundation/Flood_hazard_categories_FB03.png', 
    paste0(pi_output_dir, '/Flood_hazard_categories_FB03.png'), overwrite=TRUE)

#   * exceedance_rate_1mm_depth
dir.create(paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth'),
    recursive=TRUE, showWarnings=FALSE)

#     * logic_tree_mean/
lgtm_depth_folder = paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth/logic_tree_mean/')
lgtm_depth_files = Sys.glob('probabilistic_inundation/ptha18-NSW2023-MIS-sealevel110cm/highres_depth_with_variance/ptha18-NSW2023-MIS-sealevel110cm-depth-LogicTreeMean-sum_of_source_zones/summed_HS_domain_*_depth_as_max_stage_minus_elevation0_domain__exceedance_rate_with_threshold_0.001.tif')
copy_tifs_and_make_vrt(lgtm_depth_files, lgtm_depth_folder, "exceedance_rate_of_inundation_depth_1mm_logic_tree_mean.vrt")

#     * 84th_percentile
eu84_depth_folder = paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth/84th_percentile/')
eu84_depth_files = Sys.glob('probabilistic_inundation/ptha18-NSW2023-MIS-sealevel110cm/highres_depth_epistemic_uncertainty/84pc/random_sum_of_source_zones-depth_exrate_0.001_0.84_sum_of_source_zones/random_sum_of_source_zones_depth_rast_threshold_0.001_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(eu84_depth_files, eu84_depth_folder, "exceedance_rate_of_inundation_depth_1mm_84pc.vrt")

#     * 16th_percentile
eu16_depth_folder = paste0(pi_output_dir, '/exceedance_rate_of_inundation_1mm_depth/16th_percentile/')
eu16_depth_files = Sys.glob('probabilistic_inundation/ptha18-NSW2023-MIS-sealevel110cm/highres_depth_epistemic_uncertainty/16pc/random_sum_of_source_zones-depth_exrate_0.001_0.16_sum_of_source_zones/random_sum_of_source_zones_depth_rast_threshold_0.001_percentile_16_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(eu16_depth_files, eu16_depth_folder, "exceedance_rate_of_inundation_depth_1mm_16pc.vrt")

#
# FLOW VARIABLES, 1/2500 84TH PERCENTILE
#

#   * depth_above_initial_condition_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL
depth_84pc_folder = paste0(pi_output_dir, '/flow_variables_1in2500_at_84th_percentile/depth_above_initial_condition_with_exceedance_rate_1in2500_at_84th_percentile_masked_below_MSL/')
depth_84pc_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_depth_above_initial_condition_1in2500_84pc_at_sites_with_elevation_above_0/sum_of_all_sources_depth_above_initial_condition_where_elevation_exceeds_0_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(depth_84pc_files, depth_84pc_folder, "depth_above_initial_condition_with_exceedance_rate_1in2500_84pc_masked_below_MSL.vrt")

#   * maximum_stage_with_exceedance_rate_1in2500_at_84th_percentile
max_stage_84pc_drymasked_folder = paste0(pi_output_dir, '/flow_variables_1in2500_at_84th_percentile/max_stage_with_exceedance_rate_1in2500_at_84th_percentile/')
max_stage_84pc_drymasked_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_stage_1in2500_84pc/sum_of_all_sources_max_stage_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(max_stage_84pc_drymasked_files, max_stage_84pc_drymasked_folder, "max_stage_with_exceedance_rate_1in2500_84pc_masked_in_dry_areas.vrt")

#   * maximum_flux_with_exceedance_rate_1in2500_at_84th_percentile
max_flux_84pc_drymasked_folder = paste0(pi_output_dir, '/flow_variables_1in2500_at_84th_percentile/max_flux_with_exceedance_rate_1in2500_at_84th_percentile/')
max_flux_84pc_drymasked_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_flux_1in2500_84pc/sum_of_all_sources_max_flux_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(max_flux_84pc_drymasked_files, max_flux_84pc_drymasked_folder, "max_flux_with_exceedance_rate_1in2500_84pc_masked_in_dry_areas.vrt")

#   * maximum_speed_with_exceedance_rate_1in2500_at_84th_percentile
max_speed_84pc_drymasked_folder = paste0(pi_output_dir, '/flow_variables_1in2500_at_84th_percentile/max_speed_with_exceedance_rate_1in2500_at_84th_percentile/')
max_speed_84pc_drymasked_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_speed_1in2500_84pc/sum_of_all_sources_max_speed_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(max_speed_84pc_drymasked_files, max_speed_84pc_drymasked_folder, "max_speed_with_exceedance_rate_1in2500_84pc_masked_in_dry_areas.vrt")

# flood hazard categories 
flood_hazard_84pc_folder = paste0(pi_output_dir, '/flow_variables_1in2500_at_84th_percentile/flood_hazard_categories_1in2500_at_84th_percentile/')
flood_hazard_84pc_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_flood_hazard_categories_1in2500_84pc/sum_of_all_sources_hazard_categories_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(flood_hazard_84pc_files, flood_hazard_84pc_folder, "flood_hazard_categories_with_exceedance_rate_1in2500_84pc.vrt")

#
# FLOW VARIABLES, 1/250 50TH PERCENTILE
#


#   * depth_above_initial_condition_with_exceedance_rate_1in250_at_50th_percentile_masked_below_MSL
depth_50pc_folder = paste0(pi_output_dir, '/flow_variables_1in250_at_50th_percentile/depth_above_initial_condition_with_exceedance_rate_1in250_at_50th_percentile_masked_below_MSL/')
depth_50pc_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_depth_above_initial_condition_1in250_50pc_at_sites_with_elevation_above_0/sum_of_all_sources_depth_above_initial_condition_where_elevation_exceeds_0_rast_exrate_0.004_percentile_50_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(depth_50pc_files, depth_50pc_folder, "depth_above_initial_condition_with_exceedance_rate_1in250_50pc_masked_below_MSL.vrt")

#   * maximum_stage_with_exceedance_rate_1in250_at_50th_percentile
max_stage_50pc_drymasked_folder = paste0(pi_output_dir, '/flow_variables_1in250_at_50th_percentile/max_stage_with_exceedance_rate_1in250_at_50th_percentile/')
max_stage_50pc_drymasked_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_stage_1in250_50pc/sum_of_all_sources_max_stage_rast_exrate_0.004_percentile_50_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(max_stage_50pc_drymasked_files, max_stage_50pc_drymasked_folder, "max_stage_with_exceedance_rate_1in250_50pc_masked_in_dry_areas.vrt")

#   * maximum_flux_with_exceedance_rate_1in250_at_50th_percentile
max_flux_50pc_drymasked_folder = paste0(pi_output_dir, '/flow_variables_1in250_at_50th_percentile/max_flux_with_exceedance_rate_1in250_at_50th_percentile/')
max_flux_50pc_drymasked_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_flux_1in250_50pc/sum_of_all_sources_max_flux_rast_exrate_0.004_percentile_50_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(max_flux_50pc_drymasked_files, max_flux_50pc_drymasked_folder, "max_flux_with_exceedance_rate_1in250_50pc_masked_in_dry_areas.vrt")

#   * maximum_speed_with_exceedance_rate_1in250_at_50th_percentile
max_speed_50pc_drymasked_folder = paste0(pi_output_dir, '/flow_variables_1in250_at_50th_percentile/max_speed_with_exceedance_rate_1in250_at_50th_percentile/')
max_speed_50pc_drymasked_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_max_speed_1in250_50pc/sum_of_all_sources_max_speed_rast_exrate_0.004_percentile_50_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(max_speed_50pc_drymasked_files, max_speed_50pc_drymasked_folder, "max_speed_with_exceedance_rate_1in250_50pc_masked_in_dry_areas.vrt")

# flood hazard categories 
flood_hazard_50pc_folder = paste0(pi_output_dir, '/flow_variables_1in250_at_50th_percentile/flood_hazard_categories_1in250_at_50th_percentile/')
flood_hazard_50pc_files = Sys.glob('probabilistic_inundation/nsw_full_coast_MIS_flood_hazard_categories_1in250_50pc/sum_of_all_sources_hazard_categories_rast_exrate_0.004_percentile_50_subsam_1_Nrand_10000_seed_123_domain_index_*.tif')
copy_tifs_and_make_vrt(flood_hazard_50pc_files, flood_hazard_50pc_folder, "flood_hazard_categories_with_exceedance_rate_1in250_50pc.vrt")

# 
# ARRIVAL TIMES
#
at_output_dir = paste0(output_dir, '/Arrival_times/')
source_zone_ariv_dirs = Sys.glob('probabilistic_inundation/nsw_mis_arrival_time_min_and_scenario_average/*')
for(szad in source_zone_ariv_dirs){
    sz = basename(szad)

    # Minimum
    arrival_time_min_folder = paste0(at_output_dir, '/', sz, '/minimum_arrival_time/')
    arrival_time_min_tifs = Sys.glob(paste0(szad, '/minimum*.tif'))
    copy_tifs_and_make_vrt(arrival_time_min_tifs, arrival_time_min_folder, paste0('minimum_arrival_time_', sz, '.vrt'))

    # Mean 
    arrival_time_mean_folder = paste0(at_output_dir, '/', sz, '/mean_arrival_time/')
    arrival_time_mean_tifs = Sys.glob(paste0(szad, '/mean*.tif'))
    copy_tifs_and_make_vrt(arrival_time_mean_tifs, arrival_time_mean_folder, paste0('mean_arrival_time_', sz, '.vrt'))

}

