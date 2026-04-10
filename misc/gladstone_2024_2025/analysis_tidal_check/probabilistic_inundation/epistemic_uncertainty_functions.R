asfm = new.env()
source('../application_specific_file_metadata.R', local=asfm)

isu = new.env()
# FIXME: Move the isu functionality into part of the PTHA repository.
# Need to edit documentation to make clear that these routines are for
# importance sampling without magnitude stratification (whereas the repo
# already includes many routines for stratified and stratified/importance
# sampling).
source('../../../sources/hazard/importance_sampling_utilities.R', local=isu)

#' Get key information from the source model (unsegmented, or one of the
#' segments) when scenarios have been sampled using importance_sampling
#' (without any magnitude stratification!)
#'
#' @param source_zone Name of source zone in PTHA18 (e.g. 'kermadectonga2')
#' @param source_segment_name Name of the source zone segment in 
#' ptha18_source_rate_env$crs_data$source_envs (e.g. 'kermadectonga2_tonga', or
#' 'kermadectonga' for the unsegmented source)
#' @param importance_sampling_scenarios_logic_tree_mean_on_source_zone
#' data.frame with random scenarios, sampled via importance sampling, with rates
#' corresponding to the "combined segmented/unsegmented logic-tree-mean". We won't
#' use the rates directly but there is sufficient information in the table to estimate
#' rates for any Mw-frequency curve. For the NSW study this was created somewhere 
#' in ../../sources/hazard/....
#' @param ptha18 an environment from sourcing the "get_PTHA_results.R" script
#' @param ptha18_source_rate_env an environment from sourcing the
#' "get_PTHA18_source_zone_info.R" script
#' @return list with useful scenario frequency information for the source and
#' segment (or unsegmented)
get_key_source_information_all_logic_tree_branches_IS<-function(
    source_zone, source_segment_name, 
    importance_sampling_scenarios_logic_tree_mean_on_source_zone, 
    ptha18, ptha18_source_rate_env){

    # Chance of sampling each sampled scenario. This was determined at the time of sampling
    sampling_prob_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean_on_source_zone$sampling_prob
    # Indices of sampled scenarios in the PTHA18 event tables - may have repeated values
    row_ind_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean_on_source_zone$inds 
    # Scenario magnitude
    mw_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean_on_source_zone$event_Mw

    # Get rate models with constant rigidity for all logic tree branches
    # (There are different posterior weights for models with variable rigidity,
    # here we will only use constant rigidity).
    rate_models_all_branches = ptha18_source_rate_env$crs_data$source_envs[[source_segment_name]]$mw_rate_function(
            NA, return_all_logic_tree_branches=TRUE)

    if(source_segment_name != source_zone){
        # Segmented source, e.g.
        #   source_segment_name = "kermadectonga2_tonga".
        # Then we want 
        #   source_zone_segment="tonga"
        source_zone_segment = substring(source_segment_name, 
            nchar(source_zone)+2, nchar(source_segment_name))
        stopifnot(paste0(source_zone, '_', source_zone_segment) == source_segment_name)
    }else{
        # Unsegmented source
        source_zone_segment = ""
    }

    # Get the PTHA18 scenario conditional probability (for ALL SCENARIOS)
    # according to the logic-tree-mean model over the segment. This also
    # includes information on PTHA18 LTM rates, but we don't need to use that.
    ptha18_conditional_probability_and_rates = 
        ptha18_source_rate_env$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
            source_zone=source_zone, segment=source_zone_segment)

    # The scenario conditional probability does not change between
    # logic-tree-branches on the source_zone and segment.
    # NOTE WE ASSUME CONSTANT RIGIDITY HS SCENARIOS
    ptha18_conditional_prob_given_mw_sampled_scenarios = 
        ptha18_conditional_probability_and_rates$HS_prob_given_Mw[row_ind_sampled_scenarios]

    # Range of Mw values in PTHA18 = 7.2, 7.3, 7.4, .... 9.6, 9.7, 9.8
    ptha18_unique_mw_values = sort(unique(ptha18_conditional_probability_and_rates$HS_mw))
    stopifnot(all(abs(diff(ptha18_unique_mw_values) - 0.1) < 1e-10))
    # Mw bin boundaries in PTHA18 = 7.15, 7.25, ... 9.65, 9.75, 9.85
    ptha18_unique_mw_bb = c(ptha18_unique_mw_values - 0.05, max(ptha18_unique_mw_values) + 0.05)

    # Mw bin corresponding to sampled scenarios. Here we use rounding to avoid the possibility of tiny
    # floating point differences affecting the result
    sampled_scenario_Mw_bin = match(round(mw_sampled_scenarios, 3), round(ptha18_unique_mw_values, 3))
    stopifnot(all(!is.na(sampled_scenario_Mw_bin)))

    # For all logic tree branches, get the PTHA18 rate for each PTHA18 Mw bin
    NLTB = nrow(rate_models_all_branches$all_par) # Number of logic tree branches
    mw_bin_rate_all_logic_trees = matrix(NA, nrow=length(ptha18_unique_mw_values), ncol=NLTB)
    for(i in 1:NLTB){
        # Get exceedance-rates at the boundaries of PTHA18 Mw bins
        exrates_tmp = approx(rate_models_all_branches$Mw_seq, 
            rate_models_all_branches$all_rate_matrix[i,],
            xout=ptha18_unique_mw_bb, rule=1)$y 
            # Do not use rule=2, since the upper limit might be non-zero which
            # will ruin the "diff" below.
        # The above interpolation will be NA if any Mw value exceeds
        # max(all_branches$Mw_seq). The latter has a varying range depending on
        # the source representation. But we know the exceedance-rate is zero in
        # this case
        k = is.na(exrates_tmp)
        if(any(k)) exrates_tmp[k] = 0

        # Store the rate of individual Mw bins
        mw_bin_rate_all_logic_trees[,i] = -diff(exrates_tmp)

    }

    return(list(
        source_zone=source_zone,
        source_segment_name=source_segment_name,
        source_zone_segment=source_zone_segment,
        rate_models_all_branches=rate_models_all_branches,
        sampling_prob_sampled_scenarios = sampling_prob_sampled_scenarios,
        row_ind_sampled_scenarios = row_ind_sampled_scenarios,
        mw_sampled_scenarios = mw_sampled_scenarios,
        ptha18_conditional_prob_given_mw_sampled_scenarios=ptha18_conditional_prob_given_mw_sampled_scenarios,
        ptha18_unique_mw_values=ptha18_unique_mw_values,
        ptha18_unique_mw_bb=ptha18_unique_mw_bb,
        sampled_scenario_Mw_bin=sampled_scenario_Mw_bin,
        mw_bin_rate_all_logic_trees = mw_bin_rate_all_logic_trees))
        
}

#' Convert the raster data into a format suited to pixel-by-pixel epistemic
#' uncertainty calculation
#'
#' @param scenario_row_index Vector with ptha18 scenario indices for the
#' sampled scenarios on the source zone. These integers will also appear 
#' within the path of the raster_tar_files
#' @param raster_tar_files tar files containing rasters from simulations of
#' each of the sampled scenarios. 
#' @param the source_zone (e.g. sunda2)
#' @param DOMAIN_INDEX index of the domain to make the site data for.
#' @param raster_name_stub Inside each raster_tar_file, the raster of interest
#' is named paste0(raster_name_stub, DOMAIN_INDEX, ".tif")
#' @param MC_CORES How many cores to use
#' @param MINIMUM_DEPTH value used to possibly constrain the threshold search later.
#' Could be a numeric value, or 'adaptive_minimum'
#' @param MAXIMUM_DEPTH value used to possibly constrain the threshold search later.
#' Could be a numeric value, or 'adaptive_maximum'
#' 
#' TEST MODE: If the environment variable TARGET_POINT_GLOBAL is defined, then
#' only the value at the target point will be extracted. This is useful for
#' testing.
make_all_pixel_data<-function(scenario_row_index, raster_tar_files, 
    source_zone, DOMAIN_INDEX, raster_name_stub, MC_CORES,
    MINIMUM_DEPTH='adaptive_minimum', MAXIMUM_DEPTH='adaptive_maximum'){

    character_minimum_depth = is.character(MINIMUM_DEPTH)
    if(character_minimum_depth){
        if(MINIMUM_DEPTH != 'adaptive_minimum'){
            stop(paste0('unknown character MINIMUM_DEPTH: ', MINIMUM_DEPTH))
        }
    }
    character_maximum_depth = is.character(MAXIMUM_DEPTH)
    if(character_maximum_depth){
        if(MAXIMUM_DEPTH != 'adaptive_maximum'){
            stop(paste0('unknown character MAXIMUM_DEPTH: ', MAXIMUM_DEPTH))
        }
    }

    # Make a mapping between the scenarios and the raster tar files
    raster_tar_file_row = asfm$find_matching_md_data(scenario_row_index, 
        raster_tar_files, source_zone, return_index=TRUE)
    scenarios_to_results_inds = raster_tar_file_row
    names(scenarios_to_results_inds) = NULL #""
    stopifnot(!any(is.na(scenarios_to_results_inds)))

    # Read all the rasters as matrices, in parallel for speed
    all_raster_files = paste0('/vsitar/', raster_tar_files, '/', 
        raster_name_stub, DOMAIN_INDEX, ".tif")
    if (nzchar(Sys.getenv('TARGET_POINT_GLOBAL'))){
        # extract only at the test point
        target_point = as.numeric(unlist(strsplit(Sys.getenv('TARGET_POINT_GLOBAL'), ',')))
        target_point = matrix(target_point, nrow=1)
        print(paste0(
            'TARGET_POINT_GLOBAL is defined -- only extracting at the target point',
            ' (', target_point[1], ',
            ', target_point[2], ')' ))
        all_depth_rasters_as_matrices = lapply(all_raster_files, function(x){
            r = raster(x)
            as.matrix(extract(r, target_point))
        })
    }
    else {
        library(parallel)
        local_cluster = makeCluster(MC_CORES)
        ignore = clusterCall(local_cluster, fun=function(){ 
            library(utils); suppressPackageStartupMessages(library(raster)) })
        all_depth_rasters_as_matrices = parLapplyLB(local_cluster,
            all_raster_files,
            function(x) as.matrix(raster(x)))
        stopCluster(local_cluster)
    }

    # The following are useful and will be returned to the main program
    output_matrix_dim = dim(all_depth_rasters_as_matrices[[1]])
    template_raster = raster(all_raster_files[1])

    # Pack into a big array, with fastest dimension containing the ith raster
    big_depth_array = array(NA,
        dim=c(length(all_depth_rasters_as_matrices),
              output_matrix_dim)
    )
    for(i in 1:length(all_depth_rasters_as_matrices)){
        big_depth_array[i,,] = all_depth_rasters_as_matrices[[i]]
    }
    rm(all_depth_rasters_as_matrices)
    gc(verbose=FALSE)
    # Now convert to a list (one entry per raster cell), with each entry a list
    # that contains the depths for all model runs, and a bit of bookkeeping
    # info.
    all_pixel_data = vector(mode='list', length=prod(dim(big_depth_array)[2:3]))
    for(j in 1:dim(big_depth_array)[3]){
        for(i in 1:dim(big_depth_array)[2]){
            counter = i + (j-1)*dim(big_depth_array)[2]
            all_pixel_data[[counter]] = list(
                model_runs_max_value=big_depth_array[,i,j], 
                i = i,
                j = j,
                counter=counter,
                min = ifelse(character_minimum_depth,
                            min(big_depth_array[,i,j], na.rm=TRUE),
                            MINIMUM_DEPTH),
                max = ifelse(character_maximum_depth,
                            max(big_depth_array[,i,j], na.rm=TRUE),
                            MAXIMUM_DEPTH)
            )
        }
    }
    rm(big_depth_array)
    gc(verbose=FALSE)

    outputs = list(all_pixel_data = all_pixel_data, 
                   scenarios_to_results_inds = scenarios_to_results_inds,
                   output_matrix_dim = output_matrix_dim,
                   template_raster = template_raster)

    # REMOVE ALL VARIABLES
    to_remove = setdiff(ls(all=TRUE), 'outputs')
    rm(list=to_remove)
    gc(verbose=FALSE)

    return(outputs)
}

#' Helper function for get_estimate below
#'
#' Get the depth exceedance rate at a single site at the given logic-tree 
#' uncertainty percentile for unsegmented and segmented models. Actually
#' this works for more than just "depth" -- other quantities can be passed in
#' its place, so long as the chosen_depth and input_scenario_depths are consistent.
#' Internally the code calls other functions that are written in terms of "stage"
#' rather than "depth", so some variable names use 'stage' instead of 'depth'.
#'
#' @param chosen_depth We compute the rate of (depth > chosen_depth)
#' @param input_scenario_depths vector with depths at one site for each scenario.
#' @param all_samples List of length 1 with a data.frame containing info on the
#' sampled scenarios. No rate information will be used. 
#' @param all_source_rate_info List with the result of
#' get_key_source_information_all_logic_tree_branches_IS for the unsegmented
#' and segmented source representations
#' @param unsegmented_list_index Index in all_source_rate_info[[ ]] holding the
#' unsegmented source scenarios
#' @param segment_list_indices vector of integers corresponding to indices in
#' all_source_rate_info[[ ]] holding the segmented source scenarios
#' @param ptha18 an environment from sourcing the "get_PTHA_results.R" script
#' @param percentile_probs The logic-tree uncertainty percentile at which we
#' evaluate the exceedance rate (in [0-1]).
#' @param Nrand Number of random samples used for percentile calculations 
#' @param use_numerical_probs See
#' ptha18$compute_exceedance_rate_percentiles_with_random_sampling
#' @param reproducible_random_seed a random seed to allow reproducible
#' calculations involving random sampling
#' @param unsegmented_wt the weight associated with the unsegmented source
#' representation
#' @param segmented_wt the weight associated with the union-of-segments source
#' representation
#' @param segments_copula_type See
#' ptha18$compute_exceedance_rate_percentiles_with_random_sampling
#' @return The rate at which depth > chosen_depth, at the given epistemic
#' uncertainty percentile 
#' 
get_exrate_at_depth_and_percentile_IS<-function(
    chosen_depth, 
    input_scenario_depths, 
    all_samples,
    all_source_rate_info,
    unsegmented_list_index,
    segment_list_indices,
    ptha18, 
    percentile_probs, 
    Nrand = 1e+06, 
    use_numerical_probs=TRUE,
    reproducible_random_seed = 1234,
    unsegmented_wt=0.5,
    union_of_segments_wt=0.5,
    segments_copula_type='comonotonic'
){

    sampled_scenarios_max_depth = input_scenario_depths

    # Store depth exrate curves for each logic tree and each source representation.
    all_depth_exrates_all_logic_trees = vector(mode='list', length=length(all_source_rate_info))
    names(all_depth_exrates_all_logic_trees) = names(all_source_rate_info)

    # We assume all source representations use the same scenarios,
    # Unlike in previous variants of this code, here the calculations
    # proceed with a single copy of the sampled scenarios. Details of
    # the rates in all_samples[[1]] will not be used -- all rate information
    # will be computed using logic-tree rates
    stopifnot(length(all_samples) == 1)

    # For each source representation (unsegmented, and each segment)
    for(nm_i in names(all_source_rate_info)){
        # Compute the depth-exceedance-rates for all logic-tree branches
       
        # The PTHA18 conditional probability (conditional on its magnitude bin) for each sampled scenario.
        # This changes with nm_i, but for a given nm_i it is identical for all logic-tree branches.
        ptha18_conditional_prob_given_mw_sampled_scenarios = 
            all_source_rate_info[[nm_i]]$ptha18_conditional_prob_given_mw_sampled_scenarios

        # Index of the PTHA18 Mw bin for each sampled scenario. This actually doesn't change with nm_i.
        sampled_scenario_Mw_bin = all_source_rate_info[[nm_i]]$sampled_scenario_Mw_bin

        # The probability of sampling each scenario used in the Monte Carlo method. This actually
        # doesn't change with nm_i.
        sampling_prob_sampled_scenarios = 
            all_source_rate_info[[nm_i]]$sampling_prob_sampled_scenarios

        stopifnot(all(sampled_scenario_Mw_bin == all_source_rate_info[[1]]$sampled_scenario_Mw_bin))
        stopifnot(all(sampling_prob_sampled_scenarios == all_source_rate_info[[1]]$sampling_prob_sampled_scenarios))

        # Double check there is exactly one input max depth for each sampled scenario
        N_scenarios = length(sampled_scenarios_max_depth)
        stopifnot(N_scenarios == length(sampled_scenario_Mw_bin))
        stopifnot(N_scenarios == length(sampling_prob_sampled_scenarios))

        NTS = length(chosen_depth) # Number of depth thresholds
        NLTB = nrow(all_source_rate_info[[nm_i]]$rate_models_all_branches$all_par) # Number of logic tree branches
        N_Mw_bins_PTHA18 = nrow(all_source_rate_info[[nm_i]]$mw_bin_rate_all_logic_trees) # Number of PTHA18 magnitude bins

        # According to PTHA18, the rate of event e in logic-tree-branch i is
        #   r_i(e) = mw_bin_rate_event_e_branch_i * conditional_prob_of_e_in_its_mw_bin 
        #
        # The exceedance-rate estimate we are interested in (using importance sampling) is
        #    R(i, thresh) = sum[  r_i(e) * 
        #        { I(e exceeds threshold) / ( chance_of_sampling_e_in_monte_carlo_method * Num_scenarios) } ]
        #
        # Combining the above two equations, we can write R(i, thresh) as
        #    R(i, thresh) = sum[  mw_bin_rate_event_e_branch_i * 
        #                 { (conditional_prob_of_e_in_its_mw_bin * I(e exceeds threshold)) / 
        #                   ( chance_of_sampling_e_in_monte_carlo_method * N_scenarios) } ]
        # Notice the term in curly braces above does not depend on the logic-tree branch.
        #
        # For computational efficiency we group elements of { } into magnitude
        # bins before the sum (as with a mw-bin, each scenario will have the
        # same value of mw_bin_rate_event_e_branch_i). This facilitates doing
        # the computation for all logic tree branches with a matrix
        # multiplication of relatively small dimension
        #
        # Compute the term in curly braces and then sum results within each magnitude bin
        terms_independent_of_logic_tree_branch = matrix(0, nrow=length(chosen_depth), ncol=N_Mw_bins_PTHA18)
        for(i in 1:length(chosen_depth)){
            # Term in curly braces above, arranged by sampled scenario
            terms_indep_ltb = ptha18_conditional_prob_given_mw_sampled_scenarios * 
                (sampled_scenarios_max_depth > chosen_depth[i])/
                (sampling_prob_sampled_scenarios*N_scenarios)

            # Aggregate the above by magnitude bin
            for(j in 1:N_Mw_bins_PTHA18){
                terms_independent_of_logic_tree_branch[i,j] = sum(terms_indep_ltb *
                    (sampled_scenario_Mw_bin == j))
            }

            # Sanity check to confirm I'm just rearranging the terms.
            stopifnot(isTRUE(all.equal(
                sum(terms_independent_of_logic_tree_branch[i,]),
                sum(terms_indep_ltb))))
        }

        # Calculation of depth exceedance-rate for all logic tree branches.
        exrates_by_depth_threshold = terms_independent_of_logic_tree_branch %*% 
            all_source_rate_info[[nm_i]]$mw_bin_rate_all_logic_trees
        # Notice
        # - terms_independent_of_logic_tree_branch has matrix dimensions [length(chosen_depth), N_Mw_bins_PTHA18] 
        # - mw_bin_rate_all_logic_trees has matrix dimensions [N_Mw_bins_PTHA18, NLTB]
        # - exrates_by_depth_threshold has matrix dimensions [length(chosen_depth), N_Mw_bins_PTHA18]
        
        all_depth_exrates_all_logic_trees[[nm_i]] = list(
            logic_tree_branch_exceedance_rates = exrates_by_depth_threshold, 
            threshold_stages = chosen_depth,
            logic_tree_branch_posterior_prob = all_source_rate_info[[nm_i]]$rate_models_all_branches$all_par_prob)

    }

    # Reproducible randomness (for sampling logic-trees) so we can do the optimization
    set.seed(reproducible_random_seed)

    # Mixed segmented, unsegmented, comonotonic segments
    mean_curve_env = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
        all_depth_exrates_all_logic_trees[[unsegmented_list_index]],
        all_depth_exrates_all_logic_trees[segment_list_indices], 
        N=Nrand,
        unsegmented_wt=unsegmented_wt,
        union_of_segments_wt=union_of_segments_wt,
        segments_copula_type=segments_copula_type,
        percentile_probs=percentile_probs,
        use_numerical_probs=use_numerical_probs)

    return(mean_curve_env$percentile_exrate)
}


#' Main function to compute the rate at which EXCEEDANCE_THRESHOLD is exceeded,
#' at a given epistemic uncertainty percentile. 
#' 
#' It will be executed in parallel (over pixels)
#'
#' @param pixel_data list containing the "depths" for all scenarios, at a single pixel with cell index i,j
#' @param all_samples list with ONE ENTRY containing a data.frame's with the sampled scenario info for the logic-tree-mean.
#' We won't use the logic-tree-mean rate info directly. Because length(all_samples)==1 it needn't be in a list, but
#' writing it this way facilitates reuse of older code.
#' @param all_source_rate_info list with scenario frequency information all logic tree branches (and other useful stuff). 
#' @param scenarios_to_results_inds vector of integers such that pixel_data$model_runs_max_value[scenarios_to_results_inds]
#' is ordered by scenario in the same way as all_samples[[1]].
#' @param EXCEEDANCE_THRESHOLD The threshold
#' @param PERCENTILE_TO_USE The percentile as a number in [0,1], e.g. 0.84 for the 84th percentile
#' @param NRAND Number of random samples used for percentile calculations 
#' @param SUB_SAMPLE Integer -- e.g. a value of 3 means that we only compute 1 pixel for every 3x3 block of cells.
#' @param NEEDS_INTERPOLATING If SUB_SAMPLE>1 then this value is assigned to 'skipped pixels' to denote that they should
#' be interpolated later.
#' @param ALWAYS_WET_EXRATE The solution for a cell that is always "wet" (i.e.
#' always above EXCEEDANCE_THRESHOLD) -- this is a common case so is provided for
#' optimization.
#' @param UNSEGMENTED_INDEX Index in all_samples[[ ]] holding the unsegmented source scenarios
#' @param SEGMENTED_INDICES vector of integers corresponding to indices in all_samples[[ ]] holdinig the segmented source scenarios
#' @param unsegmented_wt The weight assigned to the unsegmented model by the logic-tree
#' @param union_of_segments_wt The weight assigned to the union_of_segments model by the logic-tree
#' @param ptha18 
#' @param REPRODUCIBLE_SEED Seed used for reproducible randomness.
#' @return The rate at which EXCEEDANCE_THRESHOLD is exceeded at the PERCENTILE_TO_USE epistemic uncertainty.
#'
get_exrate_percentile_at_pixel<-function(
    pixel_data,
    all_samples,
    all_source_rate_info, 
    scenarios_to_results_inds,
    EXCEEDANCE_THRESHOLD,
    PERCENTILE_TO_USE,
    NRAND,
    SUB_SAMPLE,
    NEEDS_INTERPOLATING,
    ALWAYS_WET_EXRATE, 
    UNSEGMENTED_INDEX,
    SEGMENTED_INDICES, 
    unsegmented_wt,
    union_of_segments_wt, 
    ptha18,
    REPRODUCIBLE_SEED
    ){

    # Do not operate on every pixel -- aggregate
    if( ((pixel_data$i-1)%%SUB_SAMPLE != floor(SUB_SAMPLE/2)) |
        ((pixel_data$j-1)%%SUB_SAMPLE != floor(SUB_SAMPLE/2)) ){
        return(NEEDS_INTERPOLATING)
    }

    if(all(is.na(pixel_data$model_runs_max_value))) return(NA)

    stopifnot(length(PERCENTILE_TO_USE) == 1)

    # If not all pixel_data is NA, then remaining NA's represent dry values.
    # (This is how the rasters herein were created). In that case they won't count
    # and we need to convert NAs to values below the threshold
    if(any(is.na(pixel_data$model_runs_max_value))){
        k = which(is.na(pixel_data$model_runs_max_value))
        below_threshold = (EXCEEDANCE_THRESHOLD - 1)
        pixel_data$model_runs_max_value[k] = below_threshold
    }

    if(!is.null(ALWAYS_WET_EXRATE)){
        # If all depths are outside above the EXCEEDANCE_THRESHOLD, return an 'always-wet' value.
        # This is a common case so it's efficient to precompute it.
        if(all(pixel_data$model_runs_max_value > EXCEEDANCE_THRESHOLD)) return(ALWAYS_WET_EXRATE)


        # It is possible for ALWAYS_WET_EXRATE to be zero (e.g. low percentile curve with limited
        # samples) in which case it makes sense to return 0
        if(ALWAYS_WET_EXRATE == 0.0) return(ALWAYS_WET_EXRATE)   
    }

    # Never exceeding
    if(all(pixel_data$model_runs_max_value < EXCEEDANCE_THRESHOLD)) return(0.0)

    # Restructure the depth data so we have one entry per sample. This might
    # not be the same as pixel_data$model_runs_max_value, because some
    # scenarios might be sampled more than once (but the tsunami model is only
    # run once).
    random_scenario_depth = pixel_data$model_runs_max_value[scenarios_to_results_inds]

    exrate_at_percentile = get_exrate_at_depth_and_percentile_IS(
        chosen_depth=EXCEEDANCE_THRESHOLD,
        input_scenario_depths=random_scenario_depth,
        all_samples=all_samples,
        all_source_rate_info=all_source_rate_info,
        unsegmented_list_index = UNSEGMENTED_INDEX,
        segment_list_indices = SEGMENTED_INDICES,
        ptha18 = ptha18,
        percentile_probs = PERCENTILE_TO_USE,
        Nrand = NRAND,
        use_numerical_probs=TRUE,
        reproducible_random_seed= REPRODUCIBLE_SEED,
        unsegmented_wt=unsegmented_wt,
        union_of_segments_wt=union_of_segments_wt,
        segments_copula_type='comonotonic')

    return(exrate_at_percentile[1,1])
}
