#' Find the path to the "get_ptha_results.R" script
#' Use a function to make it work both on NCI and my home machine
find_ptha_results_file<-function(){
    get_ptha_results_file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
    get_ptha_results_file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
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
DEFAULT_MC_CORES = 104

# Path to an untarred multidomain directory with the same setup as the hazard runs.
# Used for occasional situations when we need to 
reference_multidomain_dir = normalizePath(
    '../../swals/OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/')

# Path to a single raster_output_file.tar (having the same model setup as all
# the random scenarios).
reference_raster_tar_file = Sys.glob('../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/*/*/raster_output_files.tar')[1]

# Convenience for definition below
.source_names = c('alaskaaleutians', 'newhebrides2', 'outerrisenewhebrides', 'puysegur2', 
    'southamerica', 'kermadectonga2', 'outerrise_kermadectonga', 'outerrise_puysegur', 'solomon2')

# Named list with directories containing random scenario runs for each source-zone.
# Names must match ptha18 source-zone name
source_zone_modelled_tsunami_scenario_basedirs = lapply(.source_names, function(x){
    c(
    paste0('../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_', x, '/'),
    paste0('../../swals/OUTPUTS/ptha18-NSW2023-ID1315.5-sealevel110cm/random_', x, '/'), 
    paste0('../../swals/OUTPUTS/ptha18-NSW2023-ID4186.3-sealevel110cm/random_', x, '/')) 
    })
names(source_zone_modelled_tsunami_scenario_basedirs) = .source_names


#' Get scenario metadata corresponding to a set of random model runs from multiple importance samples for a single source zone.
#'
#' Given folders containing many model runs from multiple importance samples and a single source zone, this function should define
#' important variables. It was edited from code designed for single importance samples (the whole code might warrent a refactor in future, but not possible now).
#'
#' @param MD_BASE_DIRS is a vector of folders for multiple importance sampling, each containing many random model runs for a single source zone, e.g.
#' '../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_alaskaaleutians', '../../swals/OUTPUTS/ptha18-NSW2023-ID1315.5-sealevel110cm/random_alaskaaleutians', '../../swals/OUTPUTS/ptha18-NSW2023-ID4186.3-sealevel110cm/random_alaskaaleutians'
#' @return The function environment, containing variables 
get_scenario_metadata_from_md_base_dir<-function(MD_BASE_DIRS){

    # Scenarios that match MD_BASE_DIRS in some order -- update as needed.
    scenario_sets = c('scenarios_ID710.5', 'scenarios_ID1315.5', 'scenarios_ID4186.3')

    # Find MD_BASE_DIR matching each set of scenarios
    md_base_dir_match = unlist(lapply(scenario_sets, function(scenarios_id1234){
            idstring = gsub('scenarios_', '', scenarios_id1234)
            matching_md_dir = grep(idstring, MD_BASE_DIRS)
            if(length(matching_md_dir) != 1){
                stop('Did not match sceenario_set to unique MD_BASE_DIR: ', scenarios_id1234, MD_BASE_DIRS)
            }
            return(MD_BASE_DIRS[matching_md_dir])
        }))

    # Get information on the source zone from the MD_BASE_DIRS
    source_info = basename(MD_BASE_DIRS)
    stopifnot(startsWith(source_info, 'random_'))
    stopifnot(all(source_info == source_info[1]))
    source_zone = gsub('random_', '', source_info[1]) # This should name the unsegmented source in ptha18
    source_info = source_info[1]

    # Get tarfiles with output rasters for all runs
    raster_tar_files = Sys.glob(paste0(MD_BASE_DIRS, '/ptha18_*/raster_output_files.tar'))
    stopifnot(all(file.exists(raster_tar_files)) & 
        (length(raster_tar_files) == length(unique(raster_tar_files))))

    # Tarred multidomin directories
    tarred_multidomain_dirs = unlist(lapply(MD_BASE_DIRS,
        function(x) Sys.glob(paste0(x, '/ptha18_*/RUN*', c('.tar', '.tar.bz2'))))) # Match BOTH .tar and .tar.bz2

    # Check that raster_tar_files and tarred_multidomain_dirs exist, and are aligned.
    stopifnot(
        all(file.exists(raster_tar_files)) & 
        (length(raster_tar_files) == length(tarred_multidomain_dirs) ) &
        all(dirname(raster_tar_files) == dirname(tarred_multidomain_dirs)))

    # Directory with importance samples for each set of scenarios
    scenario_bases = paste0('../../sources/hazard/', scenario_sets, '/random_', source_zone, '/')
    stopifnot(all(file.exists(scenario_bases)))

    # Source zones for each set of scenarios (should be the same for all cases)
    all_source_names = gsub('random_', '', basename(scenario_bases))
    stopifnot(length(unique(all_source_names)) == 1)
    all_source_names = all_source_names[1]

    # importance sample filenames for each set of scenarios
    all_sample_filenames = paste0(scenario_bases, 'random_scenarios_', all_source_names, 
        '_logic_tree_mean_HS.csv')
    stopifnot(all(file.exists(all_sample_filenames)))

    # Rasters showing the weight, arranged in the order of scenario_sets
    all_sample_weight_raster_files = paste0(
        '../../sources/hazard/multiple_importance_sampling_weights/importance_sample_weight_rasters/weight_raster_', 
        source_zone, '_r_1in10000_', scenario_sets, '.tif')
    stopifnot(all(file.exists(all_sample_weight_raster_files)))

    # Read all samples, appending useful data to facilitate later calculations.
    all_samples_tmp = lapply(all_sample_filenames, read.csv)
    all_samples_Nj = unlist(lapply(all_samples_tmp, nrow))
    for(i in 1:length(all_samples_tmp)){
        sample_name = scenario_sets[i]

        # For each scenario, find the matching raster_tar_file
        k = which(grepl(md_base_dir_match[i], raster_tar_files, fixed=TRUE))
        scenario_row_string = paste0('_', source_zone, '_row_', substring(format(1e+07 + all_samples_tmp[[i]]$ind), 2, 8), '_')
        matching_raster_tar_file = rep(NA, length(scenario_row_string))
        matching_raster_tar_file_index = rep(NA, length(scenario_row_string))
        for(j in 1:length(scenario_row_string)){
            file_match = grep(scenario_row_string[j], raster_tar_files[k], fixed=TRUE)
            if(length(file_match) != 1){
                print(raster_tar_files[k])
                print(file_match)
                stop(paste0('Did not find unique raster_tar file matching ', scenario_row_string[j], ' in raster tar files above')) 
            }
            matching_raster_tar_file[j] = raster_tar_files[k[file_match]]
            matching_raster_tar_file_index[j] = k[file_match]
        }

        # Append variables to the scenario table 
        all_samples_tmp[[i]] = cbind(all_samples_tmp[[i]], data.frame(
            # - sample name
            sample_name = rep(sample_name, all_samples_Nj[i]),
            # - index (arranged in the order of scenario_sets)
            sample_group_index = rep(i, all_samples_Nj[i]),
            # - fraction of scenarios in the parent sample, vs the pooled number of scenarios
            Nj_on_sumNj = rep(all_samples_Nj[i]/sum(all_samples_Nj), all_samples_Nj[i]),
            # - filename of the sample weight raster
            weight_raster = rep(all_sample_weight_raster_files[i], all_samples_Nj[i]),
            # - MD_BASE_DIR associated with the sample
            md_base_dir = rep(md_base_dir_match[i], all_samples_Nj[i]),
            # - raster_tar_file associated with the sample
            matching_raster_tar_file = matching_raster_tar_file,
            # - index into raster_tar_file associated with the sample
            matching_raster_tar_file_index = matching_raster_tar_file_index)
        )
    }

    # Store the combined samples in a data structure that facilitates code reuse
    all_samples = list()
    all_samples$logic_tree_mean_curve_HS = do.call(rbind, all_samples_tmp)
   
    # We don't expect double-ups in raster_tar_files (albeit some filenames will ultimately point to the same files) 
    # We do expect double-ups in the samples matching_raster_tar_file, since scenarios can be repeated, but they should
    # cover all raster_tar_files.
    stopifnot(length(raster_tar_files) == length(unique(all_samples[[1]]$matching_raster_tar_file)))

    return(environment())
}

##' Given a vector of source-model row indices in the PTHA18 scenario database, find the
##' SWALS tarred multidomain_dir (or raster output files) that store the tsunami 
##' model results for each scenario.
##'
##' This depends on the naming convention of the SWALS model output files, so 
##' we make this function a user input. Likely one will just need to edit the
##' "matching_string" definition to conform to the model setup.
##'
##' @param row_indices PTHA18 row index (integer)
##' @param tarred_multidomain_data path of tarred multidomain directories OR tarred raster files. 
##' Typically files of the form
##'   ../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0024934_Mw_75_HS-full-ambient_sea_level_0.6/RUN_20230901_214821086.tar
##' OR
##'   ../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0024934_Mw_75_HS-full-ambient_sea_level_0.6/raster_output_files.tar
##' @param source_zone PTHA18 source zone, without any segment information (e.g. sunda2)
##' @param return_index If TRUE return the integer index of
##' tarred_multidomain_data corresponding to each value of row_indices. If FALSE
##' then return the corresponding entries of tarred_multidomain_data
##' @return The entry of tarred_multidomain_data with the scenario on the given source_zone and row_index.
#find_matching_md_data<-function(row_indices, tarred_multidomain_data, source_zone, return_index=FALSE){
#
#    # FIXME: For multiple importance sampling this could map a row to the "wrong" sample (although the scenario will be the same in these cases).
#    # But we'll also need a way to store the weights for each sample.
#
#    # This test is true in my contexts (but not strictly needed)
#    #stopifnot(all(endsWith(tarred_multidomain_data, '.tar'))) 
#    stopifnot(all(endsWith(tarred_multidomain_data, '.tar') | endsWith(tarred_multidomain_data, '.tar.bz2'))) 
#
#    # Make a string with the start of the SWALS output folder name (beneath
#    # ../../swals/OUTPUTS/random_sourcezone/...)
#    matching_string = paste0('ptha18_random_scenarios_', source_zone, '_row_', 
#        substring(as.character(1e+07 + row_indices), 2, 8), '_')
#
#    # Match with the tarred_multidomain_data, with NA if we don't match or get multiple matches
#    matching_ind = sapply(matching_string, function(x){
#        p = grep(x, tarred_multidomain_data)
#        if(length(p) != 1) p = NA 
#        return(p)})
#    if(any(is.na(matching_ind))){
#        stop('Could not find simulation matching scenario')
#    }
#
#    if(return_index){
#        return(matching_ind)
#    }else{
#        return(tarred_multidomain_data[matching_ind])
#    }
#}

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
