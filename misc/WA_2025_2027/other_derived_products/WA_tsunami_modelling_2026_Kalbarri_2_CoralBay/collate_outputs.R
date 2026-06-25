#
# Script to run on NCI that copies the definitive output products to folders here.
# Run on NCI from inside
#  /g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/other_derived_products/WA_tsunami_modelling_2026_Kalbarri_2_CoralBay
#

#
# INPUTS
#

# Only need one copy of the ATWS zones shapefile
ATWS_Zones = '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/ATWS_ZONES'

# All model outputs that should be provided -- for the first half there will be only one model.
model_output_basedir = 'model_outputs'
model_outputs = list(

    # Kalbarri 2 coral bay model
    kalbarri_2_coralbay = list(
        domains_shapefile = '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/domains_shapefile', 
        elevation_data_no_tidal_adjustment = '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/elevation_in_model_no_tidal_adjustment',
        elevation_data_with_tidal_adjustment = '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/elevation_in_model_with_tidal_adjustment',
        elevation_source_data = '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/elevation_source_file_index',
        tidal_adjustment = '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/tidal_adjustment',
        jatwc_inundation_zones = list(
            'Geraldton-Coast' = c(
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Geraldton-Coast/Geraldton-Coast_no_threat_with-PTHA-exrate-limit_84pc_4e-04',
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Geraldton-Coast/Geraldton-Coast_marine_warning_with-PTHA-exrate-limit_84pc_4e-04',
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Geraldton-Coast/Geraldton-Coast_land_warning_with-PTHA-exrate-limit_84pc_4e-04'),
            'Gascoyne-Coast' = c(
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Gascoyne-Coast/Gascoyne-Coast_no_threat_with-PTHA-exrate-limit_84pc_4e-04',
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Gascoyne-Coast/Gascoyne-Coast_marine_warning_with-PTHA-exrate-limit_84pc_4e-04',
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Gascoyne-Coast/Gascoyne-Coast_land_warning_with-PTHA-exrate-limit_84pc_4e-04'),
            'Ningaloo-Coast' = c(
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Ningaloo-Coast/Ningaloo-Coast_no_threat_with-PTHA-exrate-limit_84pc_4e-04',
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Ningaloo-Coast/Ningaloo-Coast_marine_warning_with-PTHA-exrate-limit_84pc_4e-04',
                '../../kalbarri_2_coralbay/analysis/jatwc_to_inundation/Inundation_zones/Ningaloo-Coast/Ningaloo-Coast_land_warning_with-PTHA-exrate-limit_84pc_4e-04')
            ),
        arrival_time = list(
            'outerrisesunda' = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/kalbarri2coralbay_arrival_time_min_and_scenario_average/outerrisesunda',
            'sunda2' = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/kalbarri2coralbay_arrival_time_min_and_scenario_average/sunda2'
            ),
        marine_warning_max_stage_AHD = list(
            'Geraldton-Coast' = '../../kalbarri_2_coralbay/analysis/convert_max_stage_outputs_to_AHD/marine_warning_max_stage_AHD/Geraldton-Coast/',
            'Gascoyne-Coast' = '../../kalbarri_2_coralbay/analysis/convert_max_stage_outputs_to_AHD/marine_warning_max_stage_AHD/Gascoyne-Coast/',
            'Ningaloo-Coast' = '../../kalbarri_2_coralbay/analysis/convert_max_stage_outputs_to_AHD/marine_warning_max_stage_AHD/Ningaloo-Coast/'
            ),
        max_stage_AHD_1in2500_84pc = '../../kalbarri_2_coralbay/analysis/convert_max_stage_outputs_to_AHD/kalbarri2coralbay_highres_domains_max_stage_AHD_percentile_0.84_exrate_0.0004_hazard/',
        max_depth_1in2500_84pc = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/kalbarri2coralbay_highres_domains_depth_percentile_0.84_exrate_0.0004_hazard',
        max_speed_1in2500_84pc = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/kalbarri2coralbay_highres_domains_max_speed_percentile_0.84_exrate_0.0004_hazard',
        max_flux_1in2500_84pc =  '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/kalbarri2coralbay_highres_domains_max_flux_percentile_0.84_exrate_0.0004_hazard',
        inundation_rate_logic_tree_mean = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/ptha18-kalbarri2coralbay-hazard/highres_depth_with_variance/ptha18-kalbarri2coralbay-hazard-depth-LogicTreeMean-sum_of_source_zones/',
        inundation_rate_84pc = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/ptha18-kalbarri2coralbay-hazard/highres_depth_epistemic_uncertainty/84pc/ptha18-kalbarri2coralbay-hazard-depth_exrate_0.001_0.84_sum_of_source_zones/',
        inundation_rate_16pc = '../../kalbarri_2_coralbay/analysis/probabilistic_inundation/ptha18-kalbarri2coralbay-hazard/highres_depth_epistemic_uncertainty/16pc/ptha18-kalbarri2coralbay-hazard-depth_exrate_0.001_0.16_sum_of_source_zones/'
        )
    )

#
# END INPUTS
#

# Check all the input files exist
files_exist = rapply(model_outputs, file.exists)
if(!all(files_exist)){
    print(files_exist)
    stop('Some files not found, deliberate halt')
}

stopifnot(file.exists(ATWS_Zones))

# Function often needed for folders with tif files
copy_tifs_and_make_vrt<-function(folder_with_tifs, output_folder, output_vrt_basename, tif_name_restriction_grep=""){
    
    # Find the tifs
    tifs_to_copy = Sys.glob(paste0(folder_with_tifs, '/*.tif'))
    # Optionally just get a subset of the tifs
    k = grep(tif_name_restriction_grep, basename(tifs_to_copy))
    tifs_to_copy = tifs_to_copy[k]
    # Check we found some files
    failed = (length(tifs_to_copy) == 0)
    if(failed){
        errmes = paste0('Could not find tifs in folder ', 
            folder_with_tifs, 
            ' that grep this (possibly empty) string: "', 
            tif_name_restriction_grep, '"')
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

# Each model's outputs are copied to a new location using this function
process_single_model<-function(model_files, model_output_dir, model_name){

    dir.create(model_output_dir, recursive=TRUE, showWarnings=FALSE)

    # Get the domain shapefile
    domains_shapefile_outputdir = paste0(model_output_dir, '/')
    dir.create(domains_shapefile_outputdir, showWarnings=FALSE)
    file.copy(model_files$domains_shapefile, domains_shapefile_outputdir, recursive=TRUE)

    # Get the elevation tifs -- no tidal adjustment
    copy_tifs_and_make_vrt(model_files$elevation_data_no_tidal_adjustment, 
        paste0(model_output_dir, '/elevation_in_model_no_tidal_adjustment'), 
        output_vrt_basename=paste0('all_elevation_no_tidal_adjustment_combined_', model_name, '.vrt'))

    # Get the elevation tifs -- with tidal adjustment
    copy_tifs_and_make_vrt(model_files$elevation_data_with_tidal_adjustment, 
        paste0(model_output_dir, '/elevation_in_model_with_tidal_adjustment'), 
        output_vrt_basename=paste0('all_elevation_with_tidal_adjustment_combined_', model_name, '.vrt'))

    # Get the raster with information on the elevation source
    elevation_source_data_output_dir = paste0(model_output_dir, '/')
    dir.create(elevation_source_data_output_dir, showWarnings=FALSE)
    file.copy(model_files$elevation_source_data, elevation_source_data_output_dir, recursive=TRUE)

    # Get the tidal adjustment raster
    tidal_adjustment_output_dir = paste0(model_output_dir, '/')
    dir.create(tidal_adjustment_output_dir, showWarnings=FALSE)
    file.copy(model_files$tidal_adjustment, tidal_adjustment_output_dir, recursive=TRUE)

    # Get the JATWC inundation zone info
    jatwc_output_dir = paste0(model_output_dir, '/jatwc_inundation_zones/')
    for(jatwc_zone in names(model_files$jatwc_inundation_zones)){
        # Polygon zones
        jatwc_zone_output_dir = paste0(jatwc_output_dir, '/', jatwc_zone)
        dir.create(jatwc_zone_output_dir, recursive=TRUE, showWarnings=FALSE)
        threat_folders = model_files$jatwc_inundation_zones[[jatwc_zone]]
        file.copy(threat_folders, jatwc_zone_output_dir, recursive=TRUE)

        # 'marine-warning-max-stage-AHD' tifs
        copy_tifs_and_make_vrt(model_files$marine_warning_max_stage_AHD[[jatwc_zone]],
            paste0(jatwc_zone_output_dir, '/', jatwc_zone, '_marine_warning_max_stage_AHD'),
            output_vrt_basename=paste0(jatwc_zone, '_marine_waring_max_stage_AHD.vrt'))
    }

    # Get the arrival times
    arrival_time_output_dir = paste0(model_output_dir, '/arrival_time')
    for(source_zone in names(model_files$arrival_time)){
        arrival_time_source_zone_dir = paste0(arrival_time_output_dir, '/', source_zone)
        # Minimum arrival time
        copy_tifs_and_make_vrt(model_files$arrival_time[[source_zone]], 
            output_folder = paste0(arrival_time_source_zone_dir, '/arrival_time_minimum'), 
            output_vrt_basename = paste0('all_arrival_time_minimum_', source_zone, '_', model_name, '.vrt'), 
            tif_name_restriction_grep='minimum_')
        # Scenario average arrival time
        copy_tifs_and_make_vrt(model_files$arrival_time[[source_zone]], 
            output_folder = paste0(arrival_time_source_zone_dir, '/arrival_time_scenario_average'), 
            output_vrt_basename = paste0('all_arrival_time_scenario_average_', source_zone, '_', model_name, '.vrt'), 
            tif_name_restriction_grep='mean_')
    }

    # max_depth
    copy_tifs_and_make_vrt(model_files$max_depth_1in2500_84pc,
        output_folder = paste0(model_output_dir, '/max_depth_1in2500_84pc'),
        output_vrt_basename = paste0('all_max_depth_1in2500_84pc_', model_name, '.vrt'))
    # max_stage (in AHD)
    copy_tifs_and_make_vrt(model_files$max_stage_AHD_1in2500_84pc,
        output_folder = paste0(model_output_dir, '/max_stage_AHD_1in2500_84pc'),
        output_vrt_basename = paste0('all_max_stage_AHD_1in2500_84pc_', model_name, '.vrt'))
    # max_speed
    copy_tifs_and_make_vrt(model_files$max_speed_1in2500_84pc,
        output_folder = paste0(model_output_dir, '/max_speed_1in2500_84pc'),
        output_vrt_basename = paste0('all_max_speed_1in2500_84pc_', model_name, '.vrt'))
    # max_flux
    copy_tifs_and_make_vrt(model_files$max_flux_1in2500_84pc,
        output_folder = paste0(model_output_dir, '/max_flux_1in2500_84pc'),
        output_vrt_basename = paste0('all_max_flux_1in2500_84pc_', model_name, '.vrt'))

    # rate of inundation (logic tree mean)
    copy_tifs_and_make_vrt(model_files$inundation_rate_logic_tree_mean,
        output_folder = paste0(model_output_dir, '/inundation_rate_logic_tree_mean'),
        output_vrt_basename = paste0('all_inundation_rate_logic_tree_mean_', model_name, '.vrt'),
        tif_name_restriction_grep = 'depth_as_max_stage_minus_elevation0_domain__exceedance_rate_with_threshold')

    # rate of inundation (84th percentile)
    copy_tifs_and_make_vrt(model_files$inundation_rate_84pc,
        output_folder = paste0(model_output_dir, '/inundation_rate_84pc'),
        output_vrt_basename = paste0('all_inundation_rate_84pc_', model_name, '.vrt'))

    # rate of inundation (16th percentile)
    copy_tifs_and_make_vrt(model_files$inundation_rate_16pc,
        output_folder = paste0(model_output_dir, '/inundation_rate_16pc'),
        output_vrt_basename = paste0('all_inundation_rate_16pc_', model_name, '.vrt'))

    return(0)
}


#
# Main function here
#
do_the_copy<-function(ATWS_Zones, model_outputs, model_output_basedir){
    file.copy(ATWS_Zones, '.', recursive=TRUE)
   
    for(model_name in names(model_outputs)){
        model_output_dir = paste0(model_output_basedir, '/', model_name)
        process_single_model(model_outputs[[model_name]], model_output_dir, model_name)
    }

}

do_the_copy(ATWS_Zones, model_outputs, model_output_basedir)
