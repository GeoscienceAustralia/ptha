#' Find the path to the "get_ptha_results.R" script
#' Use a function to make it work both on NCI and my home machine
find_ptha_results_file<-function(){
    get_ptha_results_file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access/get_PTHA_results.R'
    get_ptha_results_file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha_mm/ptha_access/get_PTHA_results.R'
    ifelse(file.exists(get_ptha_results_file_nci), get_ptha_results_file_nci, get_ptha_results_file_home)
}
# Path to the get_PTHA18_results.R script
get_ptha_results_script_file = find_ptha_results_file()

# Path to the get_detailed_PTHA18_source_zone_info.R script
get_detailed_ptha18_source_zone_info_script_file = paste0(
    dirname(get_ptha_results_script_file), 
    '/get_detailed_PTHA18_source_zone_info.R')

# Path to "plot.R" in the SWALS folder
swals_plots_script_file = paste0(dirname(dirname(get_ptha_results_script_file)), '/propagation/SWALS/plot.R')

# Default number of cores for parallel calculations
# Shared memory only. Too many cores may hit memory problems.
DEFAULT_MC_CORES = 48
DEFAULT_MC_CORES_SR = 104

# Path to an untarred multidomain directory with the same setup as the hazard runs.
# Used for occasional situations when we need to 
reference_multidomain_dir = normalizePath(
    '../../swals/OUTPUTS/extreme_source-full-ambient_sea_level_0.6/RUN_20231130_122017933/')

# Path to a single raster_output_file.tar (having the same model setup as all
# the random scenarios).
reference_raster_tar_file = Sys.glob(paste0(
    '../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random*/ptha18*/raster_output_files', c('.tar', '.tar.bz2')))[1]

# Named list with directories containing random scenario runs for each source-zone.
# Names must match ptha18 source-zone name
source_zone_modelled_tsunami_scenario_basedirs = list(
    'outerrisesunda' = '../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_outerrisesunda/',
    'sunda2'         = '../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/')


#' Get scenario metadata corresponding to a set of random model runs for a single source zone.
#'
#' Given a folder containing many model runs from a single source zone, this function should define
#' important variables.
#'
#' @param MD_BASE_DIR is a folder containing many random model runs for a single source zone, e.g.
#' '../../swals/OUTPUTS/ptha18-GreaterPerth-sealevel60cm/random_outerrisesunda/' 
#' @return The function environment, containing variables 
get_scenario_metadata_from_md_base_dir<-function(MD_BASE_DIR){

    source_info = basename(MD_BASE_DIR)
    stopifnot(startsWith(source_info, 'random_'))
    source_zone = gsub('random_', '', source_info) # This should name the unsegmented source in ptha18

    # Tarfiles with output rasters for all runs
    raster_tar_files = Sys.glob(paste0(MD_BASE_DIR, '/ptha18_*/raster_output_files', c('.tar', '.tar.bz2')))
    stopifnot(all(file.exists(raster_tar_files)))

    # Tarred multidomin directories
    tarred_multidomain_dirs = Sys.glob(paste0(MD_BASE_DIR, '/ptha18_*/RUN*', c('.tar', '.tar.bz2')))
    stopifnot(all(file.exists(raster_tar_files)) & 
        (length(raster_tar_files) == length(tarred_multidomain_dirs) ))

    if(source_zone == 'outerrisesunda'){
        # This source-zone doesn't have a segmented representation

        scenario_base = '../../sources/hazard/random_outerrisesunda/'
        all_source_names = 'outerrisesunda'
        all_source_samples = list(
            'outerrisesunda' = paste0(scenario_base, 'random_scenarios_outerrisesunda_unsegmented_HS.csv'))
        # Categorise the files above as unsegmented or segmented.
        # The source-representations are either unsegmented, or union(segments)
        UNSEGMENTED_INDEX = 1
        SEGMENTED_INDICES = c( )

    }else if(source_zone == 'sunda2'){
        # This source-zone has both an unsegmented and a "union-of-segments" representation

        scenario_base = '../../sources/hazard/random_sunda2/'

        # Files with random_scenarios for all logic-tree-branches for the sunda2 unsegmented, and the segments
        all_source_names = c('sunda2', 'sunda2_arakan', 'sunda2_andaman', 'sunda2_sumatra', 'sunda2_java')
        all_source_samples = list(
            'sunda2' = paste0(scenario_base, 'random_scenarios_sunda2_unsegmented_HS.csv'),
            'sunda2_arakan' = paste0(scenario_base, 'random_scenarios_sunda2_arakan_segment_HS.csv'),
            'sunda2_andaman' = paste0(scenario_base, 'random_scenarios_sunda2_andaman_segment_HS.csv'),
            'sunda2_sumatra' = paste0(scenario_base, 'random_scenarios_sunda2_sumatra_segment_HS.csv'),
            'sunda2_java' = paste0(scenario_base, 'random_scenarios_sunda2_java_segment_HS.csv')
        )
        # Categorise the files above as unsegmented or segmented.
        # The source-representations are either unsegmented, or union(segments)
        UNSEGMENTED_INDEX = 1
        SEGMENTED_INDICES = c(2,3,4,5)
    }else{
        stop(paste0('Unknown source zone: ', source_zone, ' from MD_BASE_DIR = ', MD_BASE_DIR))
    }

    unsegmented_wt = 1 - 0.5*(length(SEGMENTED_INDICES) > 0)
    segmented_wt   =     0.5*(length(SEGMENTED_INDICES) > 0)
    stopifnot(isTRUE(all.equal(unsegmented_wt + segmented_wt, 1.0)))

    stopifnot(UNSEGMENTED_INDEX == 1)
    if(length(SEGMENTED_INDICES) > 0){
        stopifnot(length(all_source_samples) == max(SEGMENTED_INDICES) & min(SEGMENTED_INDICES) == 2) # Beware missing a segment
    }

    stopifnot(all(all_source_names == names(all_source_samples)))
    stopifnot(all(unlist(lapply(all_source_samples, file.exists))))

    return(environment())
}

#' Given a vector of source-model row indices in the PTHA18 scenario database, find the
#' SWALS tarred multidomain_dir (or raster output files) that store the tsunami 
#' model results for each scenario.
#'
#' This depends on the naming convention of the SWALS model output files, so 
#' we make this function a user input. Likely one will just need to edit the
#' "matching_string" definition to conform to the model setup.
#'
#' @param row_indices PTHA18 row index (integer)
#' @param tarred_multidomain_data path of tarred multidomain directories OR tarred raster files. 
#' Typically files of the form
#'   ../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0024934_Mw_75_HS-full-ambient_sea_level_0.6/RUN_20230901_214821086.tar
#' OR
#'   ../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0024934_Mw_75_HS-full-ambient_sea_level_0.6/raster_output_files.tar
#' @param source_zone PTHA18 source zone, without any segment information (e.g. sunda2)
#' @param return_index If TRUE return the integer index of
#' tarred_multidomain_data corresponding to each value of row_indices. If FALSE
#' then return the corresponding entries of tarred_multidomain_data
#' @return The entry of tarred_multidomain_data with the scenario on the given source_zone and row_index.
find_matching_md_data<-function(row_indices, tarred_multidomain_data, source_zone, return_index=FALSE){

    # This test is true in my contexts (but not strictly needed)
    stopifnot(all(endsWith(tarred_multidomain_data, '.tar') | endsWith(tarred_multidomain_data, '.tar.bz2')))

    # Make a string with the start of the SWALS output folder name (beneath
    # ../../swals/OUTPUTS/random_sourcezone/...)
    matching_string = paste0('ptha18_random_scenarios_', source_zone, '_row_', 
        substring(as.character(1e+07 + row_indices), 2, 8), '_')

    # Match with the tarred_multidomain_data, with NA if we don't match or get multiple matches
    matching_ind = sapply(matching_string, function(x){
        p = grep(x, tarred_multidomain_data)
        if(length(p) != 1) p = NA 
        return(p)})
    if(any(is.na(matching_ind))){
        stop('Could not find simulation matching scenario')
    }

    if(return_index){
        return(matching_ind)
    }else{
        return(tarred_multidomain_data[matching_ind])
    }
}

#' Associated variable names with raster tif names via "raster_name_stub"
#'
#' @param VARIABLE_NAME
#' @return raster_name_stub. Inside the raster_output_files.tar archives, 
#' the rasters of interest should have names of the form
#' paste0(raster_name_stub, DOMAIN_INDEX, '.tif')
get_raster_name_stub_from_variable_name<-function(VARIABLE_NAME){
    if(VARIABLE_NAME == 'depth'){
        raster_name_stub =  'depth_as_max_stage_minus_elevation0_domain_'
    }else if(VARIABLE_NAME == 'max_stage'){
        raster_name_stub = 'max_stage_domain_'
    }else if(VARIABLE_NAME == 'max_flux'){
        raster_name_stub = 'max_flux_domain_'
    }else if(VARIABLE_NAME == 'max_speed'){
        raster_name_stub = 'max_speed_domain_'
    }else if(VARIABLE_NAME == 'arrival_time'){
        raster_name_stub = 'arrival_time_domain_'
    }else{
        stop(paste0('unsupported VARIABLE NAME', VARIABLE_NAME))
    }
    return(raster_name_stub)
}
