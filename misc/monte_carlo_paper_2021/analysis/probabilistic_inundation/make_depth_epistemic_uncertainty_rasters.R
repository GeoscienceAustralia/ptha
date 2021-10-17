#
# Compute the depth associated with a given exceedance-rate and
# epistemic-uncertainty percentile, as a raster.
#
library(rptha)

#
# Key input parameters
#

RUN_SERIES_NAME = commandArgs(trailingOnly=TRUE)[1] # 'ptha18_tonga_MSL0_meshrefine4'
DOMAIN_INDEX = as.numeric(commandArgs(trailingOnly=TRUE)[2]) # 5
PERCENTILE_TO_USE = as.numeric(commandArgs(trailingOnly=TRUE)[3]) # 0.84
TARGET_EXCEEDANCE_RATE = 1/2475
MC_CORES = 48

# Numerical parameters
#
# Do not operate on every pixel (speeds calculations)
SUB_SAMPLE = 3
# Results exceeding this value are set to MAX_DEPTH (speeds calculations)
MAX_DEPTH = 10 
# Results smaller than this value are set to MIN_DEPTH (speeds calculations)
MIN_DEPTH = 0.0001 
# How many random samples are used in the numerical percentile computation?
NRAND = 1e+04 
# Random seed used to make our 'random' function have deterministic results
# (required for root-finding). Try a different seed to check convergence
REPRODUCIBLE_SEED = 123  # 1234
# Put this value at sites that should later be interpolated (distinct from
# 'missing data' which is a deliberate gap). Don't make it too big, or
# something that will be coerced to integer (and thus not equal the numeric
# equivalent).
NEEDS_INTERPOLATING = -9999.1
# Convergence tol in "uniroot"
UNIROOT_TOL = 5.0e-03


# Get the ptha_access scripts
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
file_here = ifelse(file.exists(file_nci), file_nci, file_home)
source(file_here, local=ptha18, chdir=TRUE)

# Get the R session resulting from "compute_rates_all_sources.R", needed to 
# work with alternative logic-tree branches
ptha18_source_rate_env = new.env()
source(paste0(dirname(file_here), '/get_detailed_PTHA18_source_zone_info.R'),
       local=ptha18_source_rate_env, chdir=TRUE)

scenario_base = ifelse(grepl('_repeated', RUN_SERIES_NAME),
    '../../sources/random_repeated/',
    '../../sources/random/')

# Files with random_scenarios for all logic-tree-branches for the kermadectonga2 unsegmented, and the segments
all_source_names = c('kermadectonga2', 'kermadectonga2_tonga', 'kermadectonga2_kermadec', 'kermadectonga2_hikurangi')
all_source_samples = list(
    'kermadectonga2' = paste0(scenario_base, 'random_scenarios_kermadectonga2_unsegmented_HS.csv'),
    'kermadectonga2_tonga' = paste0(scenario_base, 'random_scenarios_kermadectonga2_tonga_segment_HS.csv'),
    'kermadectonga2_kermadec' = paste0(scenario_base, 'random_scenarios_kermadectonga2_kermadec_segment_HS.csv'),
    'kermadectonga2_hikurangi' = paste0(scenario_base, 'random_scenarios_kermadectonga2_hukurangi_segment_HS.csv'))
stopifnot(all(all_source_names == names(all_source_samples)))
stopifnot(all(unlist(lapply(all_source_samples, file.exists))))

# Categorise the files above as unsegmented or segmented.
# The source-representations are either unsegmented, or union(segments)
UNSEGMENTED_INDEX = 1
SEGMENTED_INDICES = c(2,3,4)


#
# End inputs
#


# Read the random scenarios (note: scenarios are the same in all cases, but the nominal rates differ)
all_samples = list()
for(nm_i in all_source_names){
    all_samples[[nm_i]] = read.csv(all_source_samples[[nm_i]])
}

#
# Get key information from the a source model (unsegmented, or one of the
# segments). This is the same for every site we might investigate.
#
get_logic_tree_branch_mw_bin_rates_and_posterior_probs<-function(source_segment_name, random_scenarios){

    # Get mw-exceedance-rate information for all logic tree branches in the
    # PTHA18
    all_branches = ptha18_source_rate_env$crs_data$source_envs[[source_segment_name]]$mw_rate_function(NA, 
        return_all_logic_tree_branches=TRUE)

    # Unique magnitudes that we have random scenarios for, and the
    # magnitude-bin boundaries
    unique_mw = ptha18$unique_sorted_with_check_for_even_spacing(random_scenarios$mw)
    dMw = (unique_mw[2] - unique_mw[1])
    unique_mw_bin_boundaries = c(unique_mw - (dMw/2), max(unique_mw) + dMw/2)

    # Get the logic-tree-mean rates within each magnitude bin (conveniently
    # stored in the random_scenarios)
    unique_mw_rates_source = 
        random_scenarios$rate_with_this_mw[match(unique_mw, random_scenarios$mw)]

    inv_unique_mw_rates_source = 1/unique_mw_rates_source
    inv_unique_mw_rates_source[unique_mw_rates_source == 0] = 0 # Avoid zero-division

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
        exrates_tmp = approx(all_branches$Mw_seq, all_branches$all_rate_matrix[i,], 
            xout=unique_mw_bin_boundaries, rule=1)$y

        # The above interpolation will be NA if any Mw value exceeds
        # all_branches$Mw_seq. The latter has a varying range depending on the
        # source representation. But we know the exceedance-rate is zero is
        # this case
        k = is.na(exrates_tmp)
        if(any(k)) exrates_tmp[k] = 0
        # Store the individual bin rates
        logic_tree_branch_mw_bin_rates[,i] = -diff(exrates_tmp)
    }

    return(list(
        logic_tree_branch_mw_bin_rates = logic_tree_branch_mw_bin_rates,
        logic_tree_branch_posterior_prob = all_branches$all_par_prob,
        inv_unique_mw_rates_source = inv_unique_mw_rates_source,
        unique_mw = unique_mw))

}

all_source_rate_info = list()
for(nm_i in all_source_names){
    all_source_rate_info[[nm_i]] = 
        get_logic_tree_branch_mw_bin_rates_and_posterior_probs(nm_i, all_samples[[nm_i]])
}

##
## Parallel broadcast of all_samples and all_source_rate_info
##


# Get the depth exceedance rate at the given percentile over all logic-tree
# branches, for unsegmented and segmented models. Later this can be combined
# with root-finding to get a combination of depth/percentile that meets a given
# exceedance-rate.
get_exrate_at_depth_and_percentile<-function(
    chosen_depth, 
    input_scenario_depths, 
    all_samples,
    all_source_rate_info,
    unsegmented_list_index,
    segment_list_indices,
    percentile_probs = PERCENTILE_TO_USE, 
    Nrand = 1e+06, 
    use_numerical_probs=TRUE,
    reproducible_random_seed = REPRODUCIBLE_SEED,
    unsegmented_wt=0.5,
    union_of_segments_wt=0.5,
    segments_copula_type='comonotonic'
    ){

    all_depth_exrates_all_logic_trees = list()

    for(nm_i in names(all_source_rate_info)){
        # Compute the conditional probability of exceeding each
        # threshold_stage, given the magnitude bin.
        # IN PTHA18 THIS IS INDEPENDENT OF THE LOGIC-TREE-BRANCH

        conditional_prob_exceed_stage_mw = matrix(NA, nrow=length(chosen_depth),
            ncol=length(all_source_rate_info[[nm_i]]$unique_mw))

        for(i in 1:length(chosen_depth)){

            exrate_by_mw_bin = ptha18$estimate_exrate_uncertainty(
                random_scenarios = all_samples[[nm_i]], 
                event_peak_stage = input_scenario_depths,
                threshold_stage = chosen_depth[i], 
                return_per_Mw_bin=TRUE)

            conditional_prob_exceed_stage_mw[i,] = 
                exrate_by_mw_bin$exrate*all_source_rate_info[[nm_i]]$inv_unique_mw_rates_source

        }

        all_depth_exrates_all_logic_trees[[nm_i]] = list(
            source_segment_name = nm_i,
            unique_mw = all_source_rate_info[[nm_i]]$unique_mw,
            threshold_stages = chosen_depth,
            logic_tree_branch_exceedance_rates = (
                conditional_prob_exceed_stage_mw%*%all_source_rate_info[[nm_i]]$logic_tree_branch_mw_bin_rates),
            logic_tree_branch_mw_bin_rates = all_source_rate_info[[nm_i]]$logic_tree_branch_mw_bin_rates,
            logic_tree_branch_posterior_prob = all_source_rate_info[[nm_i]]$logic_tree_branch_posterior_prob,
            conditional_prob_exceed_stage_mw = conditional_prob_exceed_stage_mw
            )
    }

    # Reproducible randomness so we can do the optimization
    set.seed(reproducible_random_seed)

    #my_numerical_probs = unique(sort(c(seq(0, 1, len=31), seq(0.78, 0.90, len=101))))
                        
    # Mixed segmented, unsegmented, comonotonic segments
    mean_curve_env = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
        all_depth_exrates_all_logic_trees[[unsegmented_list_index]],
        all_depth_exrates_all_logic_trees[segment_list_indices], 
        N=Nrand,
        unsegmented_wt=unsegmented_wt,
        union_of_segments_wt=union_of_segments_wt,
        segments_copula_type=segments_copula_type,
        percentile_probs=percentile_probs,
        #numerical_probs = my_numerical_probs,
        use_numerical_probs=use_numerical_probs)

    return(mean_curve_env$percentile_exrate)
}

# Function to use in combination with root-finding to compute the depth associated
# with a given exceedance-rate and percentile.
fun_to_zero<-function(
    chosen_depth, 
    scenario_depths,
    target_exceedance_rate = TARGET_EXCEEDANCE_RATE, 
    all_random_samples=all_samples,
    all_source_rate_representations=all_source_rate_info,
    use_numerical_probs=TRUE, 
    Nrand=1e+06, 
    percentile_probs=PERCENTILE_TO_USE, 
    reproducible_random_seed = REPRODUCIBLE_SEED,
    return_exrate_only=FALSE){

    if(length(chosen_depth) > 1) stop('chosen_depth must be of length=1')

    # Return log10( target_exrate ) - log10( exrate at chosen depth and percentile )
    exrate = get_exrate_at_depth_and_percentile(
        chosen_depth, 
        input_scenario_depths=scenario_depths, 
        all_samples=all_random_samples,
        all_source_rate_info=all_source_rate_representations,
        unsegmented_list_index = UNSEGMENTED_INDEX,
        segment_list_indices = SEGMENTED_INDICES,
        percentile_probs=percentile_probs,
        Nrand=Nrand,
        use_numerical_probs=use_numerical_probs,
        reproducible_random_seed = reproducible_random_seed)

    if(!return_exrate_only){
        # For root-finding, use difference of 'log10' exrates which is like optimising
        # a ratio -- also use protection from zero division
        zero_division_protection = (1e-20*target_exceedance_rate)
        output = log10 ( target_exceedance_rate ) - log10( zero_division_protection + exrate )
        return(c(output))
    }else{
        return(exrate)
    }
}


#
# Get the depths for every scenario, and pack into the correct data structure
#

make_all_site_data<-function(){
    # Get depth for every scenario, and pack into the right data structure
    rds_name = paste0(RUN_SERIES_NAME, '_depth_and_stage_exrate_curve_at_parliament.RDS')
    parliament_data = readRDS(rds_name)

    # Make arrays that nominally hold stage/elev for ALL SCENARIOS ON THIS
    # SOURCE-ZONE IN THE PTHA. Actually they will only be non-NA for randomly
    # sampled scenarios, however, the function we use them in
    # (estimate_exrate_uncertainty) requires this structure
    all_model_results = parliament_data$results_df
    random_scenarios = parliament_data$scenarios_databases$unsegmented_HS # Could alternatively use any segment, no change to results.
    scenarios_to_results_inds = match(random_scenarios$md_dir, all_model_results$md_dir)
    ## BROADCAST scenarios_to_results_inds

    # Read all the rasters as matrices, in parallel for speed
    library(parallel)
    local_cluster = makeCluster(MC_CORES)
    clusterCall(local_cluster, fun=function(){ library(raster) })
    all_raster_files = Sys.glob(paste0(all_model_results$md_dir, 
        "/depth_as_max_stage_minus_elevation0_domain_", DOMAIN_INDEX, ".tif"))
    all_depth_rasters_as_matrices = parLapply(local_cluster,
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
    gc()
    # Now convert to a list (one entry per raster cell), with each entry a list
    # that contains the depths for all model runs, and a bit of bookkeeping
    # info.
    all_site_data = vector(mode='list', length=prod(dim(big_depth_array)[2:3]))
    for(j in 1:dim(big_depth_array)[3]){
        for(i in 1:dim(big_depth_array)[2]){
            counter = i + (j-1)*dim(big_depth_array)[2]
            all_site_data[[counter]] = list(model_runs_max_depth=big_depth_array[,i,j], 
                                            i = i, j = j, counter=counter)
        }
    }
    rm(big_depth_array)
    gc()

    outputs = list(all_site_data = all_site_data, 
                   scenarios_to_results_inds = scenarios_to_results_inds,
                   output_matrix_dim = output_matrix_dim,
                   template_raster = template_raster)

    # REMOVE ALL VARIABLES
    to_remove = setdiff(ls(all=TRUE), 'outputs')
    rm(list=to_remove)
    gc()

    return(outputs)
}
outputs = make_all_site_data()
all_site_data = outputs$all_site_data
scenarios_to_results_inds = outputs$scenarios_to_results_inds
output_matrix_dim = outputs$output_matrix_dim
template_raster = outputs$template_raster
rm(outputs)
gc()

get_estimate<-function(site_data){
    # Main function that will be executed in parallel, with site_data holding
    # all_site_data[[j]] for the j'th pixel

    if(all(is.na(site_data$model_runs_max_depth))) return(NA)

    # Do not operate on every pixel -- aggregate
    if( ((site_data$i-1)%%SUB_SAMPLE != floor(SUB_SAMPLE/2)) |
        ((site_data$j-1)%%SUB_SAMPLE != floor(SUB_SAMPLE/2)) ){
        return(NEEDS_INTERPOLATING)
    }

    # If not all site_data is NA, then remaining NA's represent zero or near-zero values.
    # (This is how the rasters herein were created).
    # In that case we need to convert NAs to zero
    if(any(is.na(site_data$model_runs_max_depth))){
        k = which(is.na(site_data$model_runs_max_depth))
        site_data$model_runs_max_depth[k] = 0
    }

    # If all depths are outside this range, we can immediately return
    if(all(site_data$model_runs_max_depth > MAX_DEPTH)) return(MAX_DEPTH)
    if(all(site_data$model_runs_max_depth < MIN_DEPTH)) return(MIN_DEPTH)


    # The functions below require the 'random_scenario_depth' data is packed
    # into a data structure with one entry per PTHA18 scenario (which aligns with the
    # PTHA18 scenario event data). Make this here -- it is only "not NA" for scenarios
    # that were actually simulated
    random_scenario_depth = rep(NA, length(max(all_samples[[1]]$inds)))
    random_scenario_depth[all_samples[[1]]$inds] = 
        site_data$model_runs_max_depth[scenarios_to_results_inds]

    # Check lower depth limit -- if 'lower' is positive or zero, we can return MIN_DEPTH.
    lower = fun_to_zero(MIN_DEPTH, scenario_depths = random_scenario_depth, Nrand=NRAND)
    if(lower >= 0) return(MIN_DEPTH)

    # Check upper depth limit -- if upper is negative or zero, we can return MAX_DEPTH
    upper = fun_to_zero(MAX_DEPTH, scenario_depths = random_scenario_depth, Nrand=NRAND)
    if(upper <= 0) return(MAX_DEPTH)

    # Cheap estimate (fast)
    t0 = uniroot(
        fun_to_zero, 
        interval = c(MIN_DEPTH, MAX_DEPTH), 
        scenario_depths = random_scenario_depth, 
        Nrand=NRAND, 
        f.lower=lower, f.upper=upper, 
        tol=UNIROOT_TOL)

    t1 = t0

    # Refinement of estimate below -- was too slow herein 

    ## Find an root-finding interval around the cheap estimate, with new
    ## numerical parameters that will give a precise estimate
    #r1 = fun_to_zero(t0$root, 
    #    scenario_depths = random_scenario_depth,
    #    use_numerical_probs=FALSE, 
    #    Nrand=1e+06)
    #for(i in seq(-3,1)){
    #    val = min(max(t0$root - sign(r1)*10**i, MIN_DEPTH), MAX_DEPTH)
    #    #print(val)
    #    r2 = fun_to_zero(val, 
    #        scenario_depths=random_scenario_depth, 
    #        use_numerical_probs=FALSE, 
    #        Nrand=1e+06)
    #    if(r1*r2 < 0) break
    #}
    ## Store lower/upper intervals for root finding
    #if(val > t0$root){
    #    par_range = c(t0$root, val)
    #    f_range = c(r1, r2)
    #}else{
    #    par_range = c(val, t0$root)
    #    f_range = c(r2, r1)
    #}

    ## Proper root finding
    #t1 = uniroot(
    #     fun_to_zero, 
    #     interval = par_range,
    #     scenario_depths = random_scenario_depth, 
    #     Nrand=1e+06,
    #     use_numerical_probs=FALSE, 
    #     f.lower = f_range[1], 
    #     f.upper=f_range[2],
    #     tol=3.0e-04)
    
    return(t1$root)

}

try_get_estimate<-function(x) try(get_estimate(x))

#
# Clear unrequired memory before starting a parallel cluster
#
rm(ptha18_source_rate_env)
gc()

# Setup a cluster
library(parallel)
local_cluster = makeCluster(MC_CORES)
clusterCall(local_cluster, fun=function(){ library(rptha) })
# Copy over all variables except 'all_site_data', which will be too large
vars_to_export = setdiff(ls(all=TRUE), 'all_site_data')
clusterExport(local_cluster, varlist=vars_to_export)
# Check it worked
## clusterCall(local_cluster, fun=function() ls(all=TRUE, envir=.GlobalEnv))
all_pixel_results = parLapplyLB(cl = local_cluster, X = all_site_data, fun = try_get_estimate, chunk.size=10)

stopCluster(local_cluster)

# Pack to matrix
output_matrix = matrix(NA, ncol=output_matrix_dim[2], nrow=output_matrix_dim[1])
for(i in 1:length(all_pixel_results)){
    ind_i = all_site_data[[i]]$i
    ind_j = all_site_data[[i]]$j
    output_matrix[ind_i,ind_j] = all_pixel_results[[i]]
}

# Fill in the sub-sampled regions
for(j in 1:ncol(output_matrix)){
    for(i in 1:nrow(output_matrix)){
        if(is.na(output_matrix[i,j])) next
        if(output_matrix[i,j] == NEEDS_INTERPOLATING){
            # Find the 'nearest' cell that was given values
            i_ind = floor((i-1)/SUB_SAMPLE)*SUB_SAMPLE + floor(SUB_SAMPLE/2) + 1
            j_ind = floor((j-1)/SUB_SAMPLE)*SUB_SAMPLE + floor(SUB_SAMPLE/2) + 1
            output_matrix[i,j] = output_matrix[i_ind, j_ind]
        }
    }
}

# Save to file
output_raster = setValues(template_raster, c(t(output_matrix)))
dir.create(RUN_SERIES_NAME, showWarnings=FALSE)
output_raster_filename = paste0(RUN_SERIES_NAME, '/', 
    'depth_rast_invexrate_', round(1/TARGET_EXCEEDANCE_RATE), 
    '_percentile_', 100*PERCENTILE_TO_USE, 
    '_subsam_', SUB_SAMPLE,
    '_maxdpth_', MAX_DEPTH,
    '_mindpth_', MIN_DEPTH,
    '_Nrand_', NRAND,
    '_seed_', REPRODUCIBLE_SEED, 
    '_domain_index_', DOMAIN_INDEX, 
    '.tif') 
writeRaster(output_raster, output_raster_filename, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
