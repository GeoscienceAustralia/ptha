#
# Compute rasters with rate (depth > threshold).
#
library(rptha)
library(parallel)

#
# INPUTS
#

# To facilitate running multiple cases, it is easy to pass some input arguments
# via the commandline
input_args = commandArgs(trailingOnly=TRUE)
if(length(input_args) == 2 & is.finite(as.numeric(input_args[2]))){

    # The name of the folder [beneath swals/OUTPUTS/] where the runs for this
    # series are located.
    run_series_name = input_args[1] 
    # The index of the domain to run
    domain_index = as.numeric(input_args[2]) #4 

}else{
    msg = paste0(
        'Usage is like this (e.g. to compute results for domain 4 using runs ', 
        'in ../../swals/OUTPUTS/ptha18_tonga_MSL0):\n',
        '    Rscript probabilistic_inundation.R ptha18_tonga_MSL0 4')
    stop(msg)
}

# Number of cores -- the exceedance-rate computation is distributed over these.
MC_CORES = 48 
if(MC_CORES > detectCores()){
    stop(' MC_CORES exceeds the number of cores on your machine. Bad idea.')
}

# Get the files with scenario row indices and other metadata
scenario_data =  paste0('../../sources/random/',
    c('random_scenarios_kermadectonga2_hukurangi_segment_HS.csv',
      'random_scenarios_kermadectonga2_tonga_segment_HS.csv',
      'random_scenarios_kermadectonga2_kermadec_segment_HS.csv',
      'random_scenarios_kermadectonga2_unsegmented_HS.csv')
    )

# Names for each source representation in scenario_data [unsegmented, and
# various segments]
names_scenario_data = gsub('.csv', '', 
    gsub('random_scenarios_kermadectonga2_', '', basename(scenario_data), 
         fixed=TRUE), 
    fixed=TRUE)

# The multidomain directories for all SWALS model runs, for every scenario.
# These do NOT need to be ordered in the same was as the scenario_data.
md_dirs = Sys.glob(paste0('../../swals/OUTPUTS/', run_series_name, 
    '/ptha18_random_scenarios_kermadectonga2_row_*/RUN*'))
# The raster depth file associated with each md_dir. There should only be one
# per md_dir. All rasters must have the same extent and resolution.
raster_files_one_domain = paste0(md_dirs, 
    "/depth_as_max_stage_minus_elevation0_domain_", domain_index, ".tif")
# Useful to have one of the rasters ready as a template [to help with data
# export] 
raster_template = raster(raster_files_one_domain[1])

# This text will appear in the filename of all output rasters
output_raster_name_tag = paste0('domain_', domain_index, '_', run_series_name)

# Compute the exceedance-rates for each of these depths (output in raster format)
# SEPARATELY FOR THE UNSEGMENTED SOURCE AND EACH OF THE SEGMENTS. The combined 
# source would be 0.5*(sum of unsegmented + union-of-segments), since in PTHA18
# they are weighted 50-50. If you don't want to do these calculations, set it
# to an empty vector [ blahblah = c() ]
depth_thresholds_for_exceedance_rate_calculations = 
    c(0.001, 0.1, 0.5, 1, 3, 5, 8, 10)

# Compute depth-rasters corresponding to these exceedance-rates (output in raster format)
# for the COMBINED SEGMENTED + UNSEGEMENTED sources. If you don't want to do
# these calculations, set it to an empty vector [ blahblah = c() ]
exceedance_rate_thresholds_for_depth_calculations = c(
    # 2% in 50 years
    uniroot(f<-function(x){(1 - exp(-50*x)) - 0.02}, lower=0, upper=10, tol=1e-12)$root,
    # 10% in 50 years
    uniroot(f<-function(x){(1 - exp(-50*x)) - 0.1}, lower=0, upper=10, tol=1e-12)$root
    )


# Given a source-model row index in the PTHA18 scenario database, find the
# SWALS multidomain_dir that stores the tsunami model for that earthquake
# scenario. 
# This depends on the naming convention of the SWALS model output files, so 
# we make it a user input. Likely one will just need to edit the
# "matching_string" definition to conform to the model setup.
find_matching_md_dir<-function(row_indices, md_dirs){
    # Make a string with the start of the SWALS output folder name (beneath
    # Tonga_2020/swals/OUTPUTS/...)
    matching_string = paste0('ptha18_random_scenarios_kermadectonga2_row_', 
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

# We also do calculations with the 'not-self-normalised' importance sampling
# weights, and place them in their own directory with this name
alternate_run_series_name = paste0('alternate_', run_series_name)

# Get all the csv data in a list with good names
scenario_databases = lapply(scenario_data, read.csv)
names(scenario_databases) = names_scenario_data

# For each scenario in the scenario_database, append the associated md_dir that
# holds the SWALS model outputs.
for(i in 1:length(scenario_databases)){
    scenario_databases[[i]]$md_dir = find_matching_md_dir(
        scenario_databases[[i]]$scenario_row, md_dirs)
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

    # Read the rasters, and convert to a matrix with 0 or 1 [below or above
    # depth-threshold, with NA --> below] for each cell
    above_depth_threshold = lapply(local_max_depth_files, f<-function(x){
        x_mat = as.matrix(raster(x))
        out = 1 - (x_mat <= depth_threshold)
        out[is.na(out)] = 0
        rm(x_mat)
        return(out)
        })
    gc()

    # Sum [scenario_rate * above_depth_threshold] for all scenarios
    target_dim = dim(above_depth_threshold[[1]])
    ex_rate = matrix(0, ncol=target_dim[2], nrow=target_dim[1])
    for(i in 1:length(local_scenario_rates)){
        ex_rate = ex_rate + local_scenario_rates[i] * above_depth_threshold[[i]]
    }

    rm(above_depth_threshold)
    gc()
    return(ex_rate)
}

# Setup the parallel cluster. Note we will only do the exrate computation in
# parallel
local_cluster = makeCluster(MC_CORES)
clusterCall(local_cluster, fun=function(){ library(rptha) })
clusterExport(local_cluster, varlist=ls(all=TRUE))


#
# Here we compute 'exceedance-rate' rasters at specified depth thresholds
#


# Make space for the outputs
dir.create(run_series_name, showWarnings=FALSE)
dir.create(alternate_run_series_name, showWarnings=FALSE)

# For each scenario database, compute exceedance_rate rasters for a range of
# depth_thresholds
for(scenarios_name in names(scenario_databases)){

    # Map the rows of the database to the rasters
    ind = match(scenario_databases[[scenarios_name]]$md_dir, 
                dirname(raster_files_one_domain) )
    # Make rates for each raster.
    scenario_rates = rep(0, length(raster_files_one_domain))
    alternate_scenario_rates = rep(0, length(raster_files_one_domain))
    for(i in 1:length(ind)){
        # Here we loop over the scenario_database, and add the rate from the
        # table to scenario_rates. Notice this automatically treats double
        # counts, etc. Here we are only treating a single segment (or single
        # unsegmented) source at once
        scenario_rates[ind[i]] = scenario_rates[ind[i]] + 
            scenario_databases[[scenarios_name]]$scenario_rates[i]
        # Here we use the 'not-self-normalised' importance sampling based rates
        alternate_scenario_rates[ind[i]] = alternate_scenario_rates[ind[i]] + 
            scenario_databases[[scenarios_name]]$alternate_scenario_rates[i]
    }

    # For each depth-threshold, make the exceedance-rate raster
    for(depth_threshold in depth_thresholds_for_exceedance_rate_calculations){

        # Compute rasters with both 'self-normalised' importance sampling
        # rates, and also 'alternate' weights that are not self-normalised
        # (actually the 'alternate' weights are regular importance sampling
        # weights).
        for(exrates_type in c('selfNormalised', 'alternate')){
        
            if(exrates_type == 'selfNormalised'){ 

                # Compute the exceedance rates in parallel.  
                # Case with self-normalised importance sampling weights
                exrates_parallel = parLapply(cl=local_cluster, 
                    # Each process in the cluster operates on its own set of
                    # the scenarios. We sum the results below 
                    X = splitIndices(length(scenario_rates), MC_CORES),
                    fun=get_exceedance_rate_at_threshold_depth,
                    # The following arguments are not split -- the full vector
                    # is passed to every process in the cluster
                    max_depth_files=raster_files_one_domain, 
                    scenario_rates=scenario_rates, # Self-normalised
                    depth_threshold=depth_threshold)

            }else if(exrates_type == 'alternate'){

                # Compute the exceedance rates in parallel.  
                # Case with regular importance sampling weights (not self
                # normalised), which can have less bias, but has the property
                # that the weights of a group do not sum to 1 (so we don't use
                # it for Mw exceedance-rates). 
                exrates_parallel = parLapply(cl=local_cluster, 
                    # Each process in the cluster operates on its own set of
                    # the scenarios. We sum the results below 
                    X = splitIndices(length(scenario_rates), MC_CORES),
                    fun=get_exceedance_rate_at_threshold_depth,
                    # The following arguments are not split -- the full vector
                    # is passed to every process in the cluster
                    max_depth_files=raster_files_one_domain, 
                    scenario_rates=alternate_scenario_rates, # Regular IS weights
                    depth_threshold=depth_threshold)
            }

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
            exrates_rast = setValues(raster_template, combined_values)

            if(exrates_type == 'selfNormalised'){
                raster_output_file = paste0(run_series_name, '/', 
                    scenarios_name, '_', output_raster_name_tag, 
                    '_exceedance_rate_with_threshold_', depth_threshold, 
                    '.tif')
            }else if(exrates_type == 'alternate'){
                raster_output_file = paste0(alternate_run_series_name, '/', 
                    scenarios_name, '_', output_raster_name_tag, 
                    '_exceedance_rate_with_threshold_', depth_threshold, 
                    '.tif')
            }

            writeRaster(exrates_rast, raster_output_file, 
                options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
            rm(exrates_rast, combined_values, exrates_parallel)
            gc()
        }
    }
}


#
# Make a raster giving the depth at a specified target_exrate, for the COMBINED
# segmented + unsegmented sources
#

#' Store a subset of the rasters and associated scenario rates in the global
#' environment on the core
#'
#' @param core_inds indices of rasters in raster_files_one_domain that will be
#' stored on the core
#' @return nothing 
#'
make_persistent_data_on_core<-function(core_inds){

    if(length(core_inds) == 0){
        stop('The code is not setup to work with fewer indices than cores')
    }

    # The idea is to put relevant data from this function into the CORE_GLOBAL_ENV
    # Then it will be persistent on the core. Note that all lines that create
    # persistent data begin with CORE_GLOBAL_ENV... = 
    CORE_GLOBAL_ENV = globalenv()
    CORE_GLOBAL_ENV$core_inds = core_inds

    # Redundantly compute the combined 'segmented + unsegmented' rates
    # for all rasters, and then store ONLY those associated with
    # the raster_files_one_domain[core_inds]
    scenario_rates = rep(0, length(raster_files_one_domain))
    alternate_scenario_rates = rep(0, length(raster_files_one_domain))
    for(scenarios_name in names(scenario_databases)){
        # Map the rows of the database to the rasters
        ind = match(scenario_databases[[scenarios_name]]$md_dir, 
                    dirname(raster_files_one_domain) )
        for(i in 1:length(ind)){
            # 0.5* ( unsegmented + sum-of-segments)
            scenario_rates[ind[i]] = scenario_rates[ind[i]] + 
                0.5*scenario_databases[[scenarios_name]]$scenario_rates[i]
            alternate_scenario_rates[ind[i]] = alternate_scenario_rates[ind[i]] + 
                0.5*scenario_databases[[scenarios_name]]$alternate_scenario_rates[i]
        }
    }
    CORE_GLOBAL_ENV$core_segunseg_scenario_rates = scenario_rates[core_inds]
    CORE_GLOBAL_ENV$core_segunseg_alternate_scenario_rates = 
        alternate_scenario_rates[core_inds]

    # Read the rasters associated with raster_files_one_domain[core_inds]
    core_rasters = lapply(core_inds, 
        f<-function(x){
            tmp = as.matrix(raster(raster_files_one_domain[x]))
            tmp[is.na(tmp)] = 0
            return(tmp)
        })
    CORE_GLOBAL_ENV$core_rasters = core_rasters

    return(invisible())
}

#' Local-to-a-core exceedance-rate computation
#'
#' Given the threshold depth [a matrix with one value per raster cell], sum the
#' "local to a CORE" scenario rates where the depth exceeds the threshold. 
#'
#' @param dummy an unused argument that helps spread the calculation in
#' parallel -- each core works with data in its own CORE_GLOBAL_ENV. 
#' @param threshold a matrix (one value per raster cell) giving the depth for
#' which we will compute the exceedance rate.
#' @param exrate_type the type of importance sampling weight used to compute
#' the exceedance-rate.  Can be either 'selfNormalised' or 'alternate' (the
#' latter are regular importance-sampling weights).
#' @return a matrix giving the exceedance-rate for each raster cell
sum_rates_on_core_where_depth_exceeds_threshold<-function(dummy, threshold, 
    exrate_type){

    CORE_GLOBAL_ENV = globalenv()
    output = matrix(0, nrow=nrow(CORE_GLOBAL_ENV$core_rasters[[1]]), 
        ncol=ncol(CORE_GLOBAL_ENV$core_rasters[[1]]))

    for(i in 1:length(CORE_GLOBAL_ENV$core_rasters)){
        local_exceed = 1*(CORE_GLOBAL_ENV$core_rasters[[i]] > threshold)

        ## core_rasters do not have NA values [they were set to zero on read]
        #local_exceed[is.na(local_exceed)] = 0
        
        if(exrate_type == 'selfNormalised'){
            output = output + 
                local_exceed * CORE_GLOBAL_ENV$core_segunseg_scenario_rates[i]
        }else if(exrate_type == 'alternate'){
            output = output + 
                local_exceed * CORE_GLOBAL_ENV$core_segunseg_alternate_scenario_rates[i]
        }else{
            # No exrate value like this.
            # This should lead to errors in the overall calculation 
            output = output + NA
        }
    }

    return(output)
}

#' Parallel exceedance-rate computation
#'
#' This calls sum_rates_on_core_where_depth_exceeds_threshold on all cores in
#' the cluster, then does a parallel sum reduction.
#'
#' @param threshold matrix with the threshold (one value per raster cell)
#' @param exrate_type a character controlling the type of importance sampling
#' weights. Either 'selfNormalised' or 'alternate' (the latter use regular
#' importance sampling weights)
#'
get_threshold_exceedance_rates_with_reduction_over_cores<-function(threshold, 
    exrate_type){

    all_results = parLapply(local_cluster, 1:MC_CORES, 
        sum_rates_on_core_where_depth_exceeds_threshold, threshold=threshold, 
        exrate_type=exrate_type)

    # Sum the results from each core
    output = all_results[[1]] * 0
    for(i in 1:length(all_results)){
        output = output + all_results[[i]]
    }

    return(output)
}

#' Function to pass to the bisection search that finds the depth at a 
#' specified exceedance-rate
#'
#' @param threshold a matrix with one value per raster cell.
#' @param exrate_type either 'selfNormalised' or 'alternate' (the latter being
#' regular importance sampling)
#' @param target_exrate the exceedance rate of interest.
difference_between_target_exrate_and_exrate_at_threshold_depth<-function(
    threshold, exrate_type, target_exrate){

    # This has the property that if f(threshold) = 0, then threshold
    # corresponds to the depth with the target_exrate. 
    target_exrate - get_threshold_exceedance_rates_with_reduction_over_cores(
        threshold, exrate_type)
}


#' Local-to-a-core depth maxima and minima computation
#'
#' This is useful for setting the lower/upper bounds for the bisection search,
#' if used in conjunction with a reduction operation
depth_extrema_on_core<-function(dummy, bound_type){

    CORE_GLOBAL_ENV = globalenv()

    for(i in 1:length(CORE_GLOBAL_ENV$core_rasters)){

        local_depth = CORE_GLOBAL_ENV$core_rasters[[i]]
        ## Core rasters do not have NA values [set to zero on read]
        #local_depth[is.na(local_depth)] = 0
        if(i == 1){
            output = local_depth
        }

        if(bound_type == 'lower'){
            output = pmin(output, local_depth)
        }else if(bound_type == 'upper'){
            output = pmax(output, local_depth)
        }else{
            # Wrong input argument, this should lead to obvious errors.
            output = NA
        }
    }
    return(output)
}

#' Get the lower or upper bound of ALL the depth rasters over all cores.
#' 
#' @param bound_type, either 'lower' or 'upper'
#' 
depth_extrema_parallel<-function(bound_type='lower'){
    
    all_results = parLapply(local_cluster, 1:MC_CORES, depth_extrema_on_core, 
        bound_type=bound_type)

    for(i in 1:length(all_results)){
        if(i == 1){
            # Initialise the output
            output = all_results[[i]]
        }else{
            # Compute min or max
            if(bound_type == 'lower'){
                output = pmin(output, all_results[[i]])
            }else if(bound_type == 'upper'){
                output = pmax(output, all_results[[i]])
            }
        }
    }

    return(output)
}


# Read the rasters and make the scenario rates, distributed over cores, stored
# in their globalenv()
dummy = parLapply(local_cluster, 
    splitIndices(length(raster_files_one_domain), MC_CORES), 
    make_persistent_data_on_core)

# To compute depth-rasters corresponding to the exceedance rates, we use a
# bisection search. Bisection requires a 'lower-bound depth' and 'upper-bound
# depth' to guide the search, which MUST have opposite sign in the function
# being bisected to prevent spurious results
lower_bound_depth = depth_extrema_parallel(bound_type = 'lower') - 1e-05
upper_bound_depth = depth_extrema_parallel(bound_type = 'upper') + 1e-05

for(target_exrate in exceedance_rate_thresholds_for_depth_calculations){

    # Get the depth according to the importance sampling weights
    for(exrate_type in c('selfNormalised', 'alternate')){
        

        # First make sure the signs of the target function at the lower/upper
        # bounds are different. The bisection search does not check for the
        # case that sign sign( f(lower-bound) ) == sign( f(upper-bound) ), and
        # in that case it will converge toward the upper-bound. 
        lower_check = difference_between_target_exrate_and_exrate_at_threshold_depth(
            lower_bound_depth, exrate_type, target_exrate)
        upper_check = difference_between_target_exrate_and_exrate_at_threshold_depth(
            upper_bound_depth, exrate_type, target_exrate)
        if(any(sign(lower_check) == sign(upper_check))){
            msg = paste0('Error in bisection search: some function values at ',
                'lower_bound_depth and upper_bound_depth have the same sign, ',
                'so they do not bound the root')
            stop(msg)
        }

        # Bisection search, using the function that is internally run in
        # parallel
        TOL = 1e-06
        depth_at_target_exrate = bisection(
            difference_between_target_exrate_and_exrate_at_threshold_depth,
            lower=lower_bound_depth, upper=upper_bound_depth, 
            exrate_type=exrate_type, target_exrate=target_exrate,
            tolerance=TOL)
        like_zero = (depth_at_target_exrate$root < TOL)
        output_mat = depth_at_target_exrate$root
        output_mat[like_zero] = NA

        # Save it
        output_rast = setValues(raster_template, output_mat)
        if(exrate_type == 'alternate'){
            raster_output_file = paste0(alternate_run_series_name, 
                '/depth_at_exceedance_rate_1_in_', round(1/target_exrate), 
                 '_domain_',  domain_index, '.tif')
        }else{
            raster_output_file = paste0(run_series_name, 
                '/depth_at_exceedance_rate_1_in_', round(1/target_exrate), 
                 '_domain_',  domain_index, '.tif')
        }
        writeRaster(output_rast, raster_output_file, 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    }
}


# Clean up
stopCluster(local_cluster)
gc()
