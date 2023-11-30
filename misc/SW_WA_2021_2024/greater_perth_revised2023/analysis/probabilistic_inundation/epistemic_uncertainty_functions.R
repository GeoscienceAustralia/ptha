asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)

#' Get key information from the source model (unsegmented, or one of the
#' segments). 
#'
#' @param source_segment_name Name in 
#' ptha18_source_rate_env$crs_data$source_envs
#' @param random_scenarios data.frame with random scenarios for a given
#' source-zone and segment (or unsegmented)
#' @param ptha18 an environment from sourcing the "get_PTHA_results.R" script
#' @param ptha18_source_rate_env an environment from sourcing the
#' "get_PTHA18_source_zone_info.R" script
#' @return list with useful scenario frequency information for the source and
#' segment (or unsegmented)
get_logic_tree_branch_mw_bin_rates_and_posterior_probs<-function(
    source_segment_name, random_scenarios, 
    ptha18, ptha18_source_rate_env){

    # Unique magnitudes that we have random scenarios for, and the
    # magnitude-bin boundaries. Here we assume the bins are evenly spaced.
    unique_mw = ptha18$unique_sorted_with_check_for_even_spacing(random_scenarios$mw)
    dMw = (unique_mw[2] - unique_mw[1])
    unique_mw_bin_boundaries = c(unique_mw - (dMw/2), max(unique_mw) + dMw/2)

    # Get the logic-tree-mean rates within each magnitude bin (conveniently
    # stored in the random_scenarios)
    m_ind = match(unique_mw, random_scenarios$mw)
    stopifnot(!any(is.na(m_ind)))
    unique_mw_rates_source = random_scenarios$rate_with_this_mw[m_ind]

    inv_mw_bin_rates_logic_tree_mean = 1/unique_mw_rates_source
    # Avoid zero-division
    inv_mw_bin_rates_logic_tree_mean[unique_mw_rates_source == 0] = 0 

    # Get mw-exceedance-rate information for all logic tree branches in the
    # PTHA18
    all_branches = ptha18_source_rate_env$crs_data$source_envs[[source_segment_name]]$mw_rate_function(NA, 
        return_all_logic_tree_branches=TRUE)

    # For each logic-tree branch, make a matrix to store the rates in each mw
    # bin.
    # - columns correspond to different logic-tree branches,
    # - rows correspond to unique_mw
    num_logictree_branches = length(all_branches$all_par_prob)
    logic_tree_branch_mw_bin_rates = matrix(NA, 
        ncol=num_logictree_branches, nrow=length(unique_mw_rates_source))
    for(i in 1:num_logictree_branches){
        # Interpolate the exceedance-rate curve the magnitude-bin boundaries
        # (typically 7.15, 7.25, .... 9.55, 9.65)
        exrates_tmp = approx(all_branches$Mw_seq, 
            all_branches$all_rate_matrix[i,], 
            xout=unique_mw_bin_boundaries, rule=1)$y

        # The above interpolation will be NA if any Mw value exceeds
        # max(all_branches$Mw_seq). The latter has a varying range depending on
        # the source representation. But we know the exceedance-rate is zero in
        # this case
        k = is.na(exrates_tmp)
        if(any(k)) exrates_tmp[k] = 0
        # Store the individual bin rates
        logic_tree_branch_mw_bin_rates[,i] = -diff(exrates_tmp)
    }

    return(list(
        logic_tree_branch_mw_bin_rates = logic_tree_branch_mw_bin_rates,
        logic_tree_branch_posterior_prob = all_branches$all_par_prob,
        inv_mw_bin_rates_logic_tree_mean = inv_mw_bin_rates_logic_tree_mean,
        unique_mw = unique_mw))

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
make_all_pixel_data<-function(scenario_row_index, raster_tar_files, 
    source_zone, DOMAIN_INDEX, raster_name_stub, MC_CORES){

    # Make a mapping between the scenarios and the raster tar files
    raster_tar_file_row = asfm$find_matching_md_data(scenario_row_index, 
        raster_tar_files, source_zone, return_index=TRUE)

    scenarios_to_results_inds = raster_tar_file_row
    names(scenarios_to_results_inds) = ""
    stopifnot(!any(is.na(scenarios_to_results_inds)))

    # Read all the rasters as matrices, in parallel for speed
    library(parallel)
    local_cluster = makeCluster(MC_CORES)
    ignore = clusterCall(local_cluster, fun=function(){ 
        library(utils); suppressPackageStartupMessages(library(raster)) })
    all_raster_files = paste0('/vsitar/', raster_tar_files, '/', 
        raster_name_stub, DOMAIN_INDEX, ".tif")
    #all_depth_rasters_as_matrices = parLapply(local_cluster,
    all_depth_rasters_as_matrices = parLapplyLB(local_cluster,
        all_raster_files,
        function(x) as.matrix(raster(x)))
    stopCluster(local_cluster)

    # The following are useful and will be returned to the main program
    output_matrix_dim = dim(all_depth_rasters_as_matrices[[1]])
    template_raster = raster(all_raster_files[1])

    # Pack into a big array, with fastest dimension containing the ith raster
    big_depth_array = array(NA, 
        dim=c(length(all_depth_rasters_as_matrices), 
              output_matrix_dim))
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
            all_pixel_data[[counter]] = list(model_runs_max_value=big_depth_array[,i,j], 
                                            i = i, j = j, counter=counter)
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
#' BEWARE: This has length >= max(all_samples[[1]]$inds). It contains mostly NA
#' values except at indices in all_samples[[1]]$inds, where it contains the
#' corresponding depth. The format is a bit strange but corresponds to that in
#' ptha18$estimate_exrate_uncertainty
#' @param all_samples List with (for the unsegmented and segmented sources) the
#' sampled scenarios. We assume the scenarios are identical (but they have
#' different rates).
#' @param all_source_rate_info List with the result of
#' get_logic_tree_branch_mw_bin_rates_and_posterior_probs for each entry of
#' all_samples
#' @param unsegmented_list_index Index in all_samples[[ ]] holding the
#' unsegmented source scenarios
#' @param segment_list_indices vector of integers corresponding to indices in
#' all_samples[[ ]] holding the segmented source scenarios
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
get_exrate_at_depth_and_percentile<-function(
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

    all_depth_exrates_all_logic_trees = list()

    # Catch a common mistake regarding input_scenario_depths, which needs to 
    # have sufficient length to store non-sampled scenarios as well (even
    # though in practice they take on NA values).
    stopifnot(length(input_scenario_depths) >= max(all_samples[[1]]$inds))

    for(nm_i in names(all_source_rate_info)){
        # Compute the conditional probability of exceeding each
        # threshold_stage within each magnitude bin.
        #
        # IN PTHA18 THIS IS INDEPENDENT OF THE LOGIC-TREE-BRANCH.
        # The logic-tree branches give different rates to each magnitude bin,
        # but do not change the scenario conditional probabilities within a
        # magnitude bin.

        # We assume the same scenarios in all cases (but their rates differ)
        stopifnot(all(all_samples[[nm_i]]$inds == all_samples[[1]]$inds))

        conditional_prob_exceed_stage_mw = matrix(NA, nrow=length(chosen_depth),
            ncol=length(all_source_rate_info[[nm_i]]$unique_mw))

        for(i in 1:length(chosen_depth)){
            # Logic-tree mean exrates by magnitude bin
            exrate_by_mw_bin = ptha18$estimate_exrate_uncertainty(
                random_scenarios = all_samples[[nm_i]], 
                event_peak_stage = input_scenario_depths,
                threshold_stage = chosen_depth[i], 
                return_per_Mw_bin=TRUE)

            # This does not change with logic-tree branch in PTHA18. (The
            # conditional probabilities of scenarios within each mw bin are the
            # same for all logic-tree branches).
            conditional_prob_exceed_stage_mw[i,] = exrate_by_mw_bin$exrate*
                all_source_rate_info[[nm_i]]$inv_mw_bin_rates_logic_tree_mean

        }

        all_depth_exrates_all_logic_trees[[nm_i]] = list(
            source_segment_name = nm_i,
            unique_mw = all_source_rate_info[[nm_i]]$unique_mw,
            threshold_stages = chosen_depth,
            logic_tree_branch_exceedance_rates = (
                conditional_prob_exceed_stage_mw%*%
                all_source_rate_info[[nm_i]]$logic_tree_branch_mw_bin_rates),
            logic_tree_branch_mw_bin_rates = 
                all_source_rate_info[[nm_i]]$logic_tree_branch_mw_bin_rates,
            logic_tree_branch_posterior_prob = 
                all_source_rate_info[[nm_i]]$logic_tree_branch_posterior_prob,
            conditional_prob_exceed_stage_mw = conditional_prob_exceed_stage_mw
            )
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
#' @param all_samples list of data.frame's with the sampled scenario info for unsegmented/segments
#' @param all_source_rate_info list with mw bin rates and posterior probabilities for all logic tree branches. 
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
get_exrate_percentile_at_pixel<-function(pixel_data, all_samples, all_source_rate_info, 
    scenarios_to_results_inds,
    EXCEEDANCE_THRESHOLD, PERCENTILE_TO_USE,
    NRAND, SUB_SAMPLE, NEEDS_INTERPOLATING, ALWAYS_WET_EXRATE, 
    UNSEGMENTED_INDEX, SEGMENTED_INDICES, 
    unsegmented_wt, union_of_segments_wt, 
    ptha18,
    REPRODUCIBLE_SEED){

    if(all(is.na(pixel_data$model_runs_max_value))) return(NA)

    stopifnot(length(PERCENTILE_TO_USE) == 1)

    # Do not operate on every pixel -- aggregate
    if( ((pixel_data$i-1)%%SUB_SAMPLE != floor(SUB_SAMPLE/2)) |
        ((pixel_data$j-1)%%SUB_SAMPLE != floor(SUB_SAMPLE/2)) ){
        return(NEEDS_INTERPOLATING)
    }

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
    }

    # The functions below require the 'random_scenario_depth' data is packed
    # into a data structure with one entry per PTHA18 scenario (which aligns with the
    # PTHA18 scenario event data). Make this here -- it is only "not NA" for scenarios
    # that were actually simulated, and only those "not NA" values are accessed/used.
    random_scenario_depth = rep(NA, max(all_samples[[1]]$inds)) # Large enough to hold all not-NA values
    random_scenario_depth[all_samples[[1]]$inds] = 
        pixel_data$model_runs_max_value[scenarios_to_results_inds]

    exrate_at_percentile = get_exrate_at_depth_and_percentile(
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
