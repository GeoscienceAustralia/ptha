#
# Compute rasters with rate (depth > threshold).
#
#library(terra)

library(utils)
suppressPackageStartupMessages(library(rptha))
suppressPackageStartupMessages(library(parallel))

#
# INPUTS
#

# To facilitate running multiple cases, it is easy to pass some input arguments
# via the commandline
input_args = commandArgs(trailingOnly=TRUE)
if(!(length(input_args) == 2 & is.finite(as.numeric(input_args[2])))){
    msg = paste0(
        'Usage is like this (e.g. to compute results for domain 4 using runs ', 
        'in ../../swals/OUTPUTS/ptha18_tonga_MSL0):\n',
        '    Rscript probabilistic_inundation.R ptha18_tonga_MSL0 4')
    stop(msg)
}

# The name of the folder [beneath swals/OUTPUTS/] where the runs for this
# series are located.
run_series_name = input_args[1] 
# The index of the domain to run
domain_index = as.numeric(input_args[2]) #4 

# Number of cores -- the exceedance-rate computation is distributed over these.
MC_CORES = 48 
if(MC_CORES > detectCores()){
    stop(' MC_CORES exceeds the number of cores on your machine. Bad idea.')
}

source_info = basename(run_series_name)
if(source_info == 'random_outerrisesunda'){
    # Key parameters for outerrisesunda source

    source_zone = 'outerrisesunda'
    scenario_base = '../../sources/hazard/random_outerrisesunda/'

    # Files with scenario row indices and other metadata
    scenario_data = paste0(scenario_base, 
        c('random_scenarios_outerrisesunda_unsegmented_HS.csv')
        )

    # Names for each source representation in scenario_data [unsegmented, and
    # various segments]
    names_scenario_data = gsub('.csv', '',
        gsub('random_scenarios_outerrisesunda_', '', basename(scenario_data),
             fixed=TRUE),
        fixed=TRUE)

    # The multidomain directories for all SWALS model runs, for every scenario.
    # These do NOT need to be ordered in the same was as the scenario_data.
    md_dirs = Sys.glob(paste0('../../swals/OUTPUTS/', run_series_name, 
        '/ptha18_random_scenarios_outerrisesunda_row_*/RUN*'))

}else if(source_info == 'random_sunda2'){
    # Key parameters for sunda2 source

    source_zone = 'sunda2'
    scenario_base = '../../sources/hazard/random_sunda2/'

    # Files with scenario row indices and other metadata
    scenario_data = paste0(scenario_base, 
        c('random_scenarios_sunda2_andaman_segment_HS.csv',
          'random_scenarios_sunda2_arakan_segment_HS.csv',
          'random_scenarios_sunda2_java_segment_HS.csv',
          'random_scenarios_sunda2_sumatra_segment_HS.csv',
          'random_scenarios_sunda2_unsegmented_HS.csv')
        )
    # Names for each source representation in scenario_data [unsegmented, and
    # various segments]
    names_scenario_data = gsub('.csv', '', 
        gsub('random_scenarios_sunda2_', '', basename(scenario_data), 
             fixed=TRUE), 
        fixed=TRUE)

    # The multidomain directories for all SWALS model runs, for every scenario.
    # These do NOT need to be ordered in the same was as the scenario_data.
    md_dirs = Sys.glob(paste0('../../swals/OUTPUTS/', run_series_name, 
        '/ptha18_random_scenarios_sunda2_row_*/RUN*'))

}else{
    stop(paste0('Unknown source_info: ', source_info))
}

# Output directory stub -- make sure it doesn't include "/"
run_series_name_nodir = gsub('/', '-', run_series_name, fixed=TRUE)

# Compute the exceedance-rates for each of these depths (output in raster format)
# SEPARATELY FOR THE UNSEGMENTED SOURCE AND EACH OF THE SEGMENTS. The combined 
# source would be 0.5*(sum of unsegmented + union-of-segments), since in PTHA18
# they are weighted 50-50. 
depth_thresholds_for_exceedance_rate_calculations = 0.001  # c(0.001, 0.1, 0.5, 1, 3, 5, 8, 10)

# Given a source-model row index in the PTHA18 scenario database, find the
# SWALS multidomain_dir that stores the tsunami model for that earthquake
# scenario. 
# This depends on the naming convention of the SWALS model output files, so 
# we make it a user input. Likely one will just need to edit the
# "matching_string" definition to conform to the model setup.
find_matching_md_dir<-function(row_indices, md_dirs, source_zone){
    # Make a string with the start of the SWALS output folder name (beneath
    # ../../swals/OUTPUTS/...)
    matching_string = paste0('ptha18_random_scenarios_', source_zone, '_row_', 
        substring(as.character(1e+07 + row_indices), 2, 8), '_')

    # Match with the md_dirs, with NA if we don't match or get multiple matches
    matching_ind = sapply(matching_string, f<-function(x){
        p = grep(x, md_dirs)
        if(length(p) != 1) p = NA 
        return(p)})
    if(any(is.na(matching_ind))){
        stop('Could not find simulation matching scenario')
    }

    return(md_dirs[matching_ind])
}

#
# END INPUTS
#


#
# READ THE SCENARIO DATABASES
#

# Get all the csv data in a list with good names
scenario_databases = lapply(scenario_data, read.csv)
names(scenario_databases) = names_scenario_data

# For each scenario in the scenario_database, append the associated md_dir that
# holds the SWALS model outputs.
for(i in 1:length(scenario_databases)){
    scenario_databases[[i]]$md_dir = find_matching_md_dir(
        scenario_databases[[i]]$inds, md_dirs, source_zone)
}



#' Given max-depth matrices and scenario rates, compute rate that a
#' depth_threshold is exceeded. 
#' 
#' NA values in the max-depth matrix will be treated as dry.
#'
#' @param included_indices a vector of non-repeated indices in
#' 1:length(scenario_rates) giving the rasters to include. This is used for
#' splitting the calculation in parallel -- a serial calculation would 
#' use 1:length(scenario_rates) -- whereas in parallel we split up the latter,
#' and subsequently sum them.
#' @param max_depth_files A list of rasters containing the max_depth (one for
#' each entry of md_dirs).
#' @param scenario_rates A vector with the individual scenario rates for each
#' entry of max_depth_files
#' @param depth_threshold The function will compute the exceedance rate of
#' (depth > depth_threshold).
#'
get_exceedance_rate_at_threshold_depth<-function(included_indices, 
    max_depth_files, scenario_rates, depth_threshold){

    stopifnot(length(scenario_rates) == length(max_depth_files))

    stopifnot(length(included_indices) == length(unique(included_indices)))

    stopifnot( (min(included_indices) >= 1) & 
               (max(included_indices) <= length(max_depth_files)) )

    local_max_depth_files = max_depth_files[included_indices]
    local_scenario_rates = scenario_rates[included_indices]

    # Sum [scenario_rate * above_depth_threshold] for all scenarios
    #r1 = as.matrix(terra::rast(local_max_depth_files[1]), wide=TRUE)
    r1 = as.matrix(raster(local_max_depth_files[1]))
    #target_dim = terra::dim(r1)
    target_dim = dim(r1)
    rm(r1)
    ex_rate = matrix(0, ncol=target_dim[2], nrow=target_dim[1])
    for(i in 1:length(local_scenario_rates)){
        # Read the raster
        #x_mat = as.matrix(terra::rast(local_max_depth_files[i]), wide=TRUE)
        x_mat = as.matrix(raster(local_max_depth_files[i]))
        # Make a raster that is 1 where we exceed the depth threshold, and 0 elsewhere
        # (including at NA sites)
        out = 1 - (x_mat <= depth_threshold)
        out[is.na(out)] = 0
        # Sum the exceedance-rates
        ex_rate = ex_rate + local_scenario_rates[i] * out
        rm(out, x_mat)
    }
    gc()
    return(ex_rate)
}

# Setup the parallel cluster. Note we will only do the exrate computation in
# parallel
local_cluster = makeCluster(MC_CORES)
#clusterCall(local_cluster, fun=function(){ library(terra); library(rptha) })
clusterCall(local_cluster, fun=function(){library(utils); suppressPackageStartupMessages(library(rptha)) })
clusterExport(local_cluster, varlist=ls(all=TRUE))


#
# Here we compute 'exceedance-rate' rasters at specified depth thresholds
#


# Make space for the outputs
dir.create(run_series_name_nodir, showWarnings=FALSE)

# For each scenario database, compute exceedance_rate rasters for a range of
# depth_thresholds
for(scenarios_name in names(scenario_databases)){

    # Useful to have one of the rasters ready as a template [to help with data
    # The raster depth file associated with each md_dir. There should only be one
    # per md_dir. All rasters must have the same extent and resolution.
    raster_files_one_domain = paste0('/vsitar/', dirname(md_dirs), 
        "/raster_output_files.tar/depth_as_max_stage_minus_elevation0_domain_", domain_index, ".tif")

    raster_template = raster(raster_files_one_domain[1])
    #raster_template = terra::rast(raster_files_one_domain[1])

    # This text will appear in the filename of all output rasters
    output_raster_name_tag = paste0('domain_', domain_index, '_', run_series_name_nodir)

    # Map the rows of the database to the rasters
    ind = match(dirname(scenario_databases[[scenarios_name]]$md_dir), 
                gsub('/vsitar/', '', dirname(dirname(raster_files_one_domain)), fixed=TRUE) )

    # Make rates for each raster.
    scenario_rates = rep(0, length(raster_files_one_domain))
    for(i in 1:length(ind)){
        # Here we loop over the scenario_database, and add the rate from the
        # table to scenario_rates. Notice this automatically treats double
        # counts, etc. Here we are only treating a single segment (or single
        # unsegmented) source at once
        scenario_rates[ind[i]] = scenario_rates[ind[i]] + 
            scenario_databases[[scenarios_name]]$importance_sampling_scenario_rates_basic[i]
    }

    # For each depth-threshold, make the exceedance-rate raster
    for(depth_threshold in depth_thresholds_for_exceedance_rate_calculations){

        # Compute the exceedance rates in parallel.  
        # Basic importance sampling weights
        exrates_parallel = parLapply(cl=local_cluster, 
            # Each process in the cluster operates on its own set of
            # the scenarios. We sum the results below 
            X = splitIndices(length(scenario_rates), MC_CORES),
            fun=get_exceedance_rate_at_threshold_depth,
            # The following arguments are not split -- the full vector
            # is passed to every process in the cluster
            max_depth_files=raster_files_one_domain, 
            scenario_rates=scenario_rates, 
            depth_threshold=depth_threshold)

        # Sum the exceedance rates from each cluster process
        combined_values = exrates_parallel[[1]]*0
        for(i in 1:length(exrates_parallel)){
            combined_values = combined_values + exrates_parallel[[i]]
        }

        # For the raster output, it is nice to set regions that are never
        # inundated to NA (genuinely NA regions that are not priority
        # domain will also be NA)
        combined_values[combined_values == 0] = NA

        # Convert to a raster and write to file
        #exrates_rast = terra::setValues(raster_template, combined_values)
        exrates_rast = setValues(raster_template, combined_values)

        raster_output_file = paste0(run_series_name_nodir, '/', 
            scenarios_name, '_', output_raster_name_tag, 
            '_exceedance_rate_with_threshold_', depth_threshold, 
            '.tif')

        #terra::writeRaster(exrates_rast, raster_output_file, 
        writeRaster(exrates_rast, raster_output_file, 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        rm(exrates_rast, combined_values, exrates_parallel)
        gc()
    }
}



