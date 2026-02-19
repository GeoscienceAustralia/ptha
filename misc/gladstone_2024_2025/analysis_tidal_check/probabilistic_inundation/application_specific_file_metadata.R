#' Find the path to the "get_ptha_results.R" script
#' Use a function to make it work both on NCI and my home machine
find_ptha_results_file<-function(){
    file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access/get_PTHA_results.R'
    file_home = ''
    ifelse(file.exists(file_nci), file_nci, file_home)
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

# Path to an untarred multidomain directory with the same setup as the hazard runs.
# Used for occasional situations when we need to 
# Paths are relative to the location of the script that sources this file.
reference_multidomain_dir = normalizePath(
    '../../../swals/OUTPUTS/v6/model_load_balanced-test_load_balance-ambient_sea_level_0/RUN_20240829_164918630/')

sea_level = "vary"

# Path to a single raster_output_file.tar (having the same model setup as all
# the random scenarios).
reference_raster_tar_file = Sys.glob(paste0('../../../swals/OUTPUTS/v6/ptha18_tidal_check/sea_level_', sea_level, '/*/*/raster_output_files.tar'))[1]

# Convenience for definition below
# .source_names = c('alaskaaleutians', 'newhebrides2', 'outerrisenewhebrides', 'puysegur2', 
#     'southamerica', 'kermadectonga2', 'outerrise_kermadectonga', 'outerrise_puysegur', 'solomon2')
.source_names = c('newhebrides2', 'southamerica', 'kermadectonga2', 'solomon2')

# Named list with directories containing random scenario runs for each source-zone.
# Names must match ptha18 source-zone name
source_zone_modelled_tsunami_scenario_basedirs = lapply(.source_names, function(x){
    paste0('../../../swals/OUTPUTS/v6/ptha18_tidal_check/sea_level_', sea_level, '/random_', x, '/')})
names(source_zone_modelled_tsunami_scenario_basedirs) = .source_names

#' Get scenario metadata corresponding to a set of random model runs for a single source zone.
#'
#' Given a folder containing many model runs from a single source zone, this function should define
#' important variables.
#'
#' @param MD_BASE_DIR is a folder containing many random model runs for a single source zone, e.g.
#' '../../../swals/OUTPUTS/ptha18-GreaterPerth-sealevel60cm/random_outerrisesunda/' 
#' @return The function environment, containing variables 
get_scenario_metadata_from_md_base_dir<-function(MD_BASE_DIR){

    source_info = basename(MD_BASE_DIR)
    stopifnot(startsWith(source_info, 'random_'))
    source_zone = gsub('random_', '', source_info) # This should name the unsegmented source in ptha18

    # Tarfiles with output rasters for all runs
    raster_tar_files = Sys.glob(paste0(MD_BASE_DIR, '/ptha18_*/raster_output_files.tar'))
    stopifnot(all(file.exists(raster_tar_files)))

    # Tarred multidomin directories
    tarred_multidomain_dirs = Sys.glob(paste0(MD_BASE_DIR, '/ptha18_*/RUN*', c('.tar', '.tar.bz2'))) # Match BOTH .tar and .tar.bz2
    stopifnot(all(file.exists(raster_tar_files)) & 
        (length(raster_tar_files) == length(tarred_multidomain_dirs) ))

    scenario_base = paste0('../../../sources/hazard/tide_check_50/random_', source_zone, '/')
    # FIXME: Consider renaming the following variables to avoid "all_" (since
    # in the current code, using importance sampling, they are of length=1)
    #   * Consider "all_source_names" --> vector with names to lookup ptha18_detailed$....source-details...
    #   * Consider "all_source_samples" --> source_samples_LTM
    all_source_names = gsub('random_', '', basename(scenario_base))
    stopifnot(length(all_source_names) == 1)
    # To reuse code from previous projects, it's convenient if
    # all_source_samples is a list, albeit of length 1 holding 1 filename
    all_source_samples = list() 
    all_source_samples$logic_tree_mean_curve_HS = paste0(scenario_base, 'random_scenarios_', all_source_names, 
        '_logic_tree_mean_HS.csv')
    stopifnot(length(all_source_samples) == 1)
    stopifnot(length(all_source_samples[[1]]) == 1)
    stopifnot(file.exists(all_source_samples[[1]]))

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
#'   ../../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0024934_Mw_75_HS-full-ambient_sea_level_0.6/RUN_20230901_214821086.tar
#' OR
#'   ../../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0024934_Mw_75_HS-full-ambient_sea_level_0.6/raster_output_files.tar
#' @param source_zone PTHA18 source zone, without any segment information (e.g. sunda2)
#' @param return_index If TRUE return the integer index of
#' tarred_multidomain_data corresponding to each value of row_indices. If FALSE
#' then return the corresponding entries of tarred_multidomain_data
#' @return The entry of tarred_multidomain_data with the scenario on the given source_zone and row_index.
find_matching_md_data<-function(row_indices, tarred_multidomain_data, source_zone, return_index=FALSE){

    # This test is true in my contexts (but not strictly needed)
    #stopifnot(all(endsWith(tarred_multidomain_data, '.tar'))) 
    stopifnot(all(endsWith(tarred_multidomain_data, '.tar') | endsWith(tarred_multidomain_data, '.tar.bz2'))) 

    # Make a string with the start of the SWALS output folder name (beneath
    # ../../../swals/OUTPUTS/random_sourcezone/...)
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
