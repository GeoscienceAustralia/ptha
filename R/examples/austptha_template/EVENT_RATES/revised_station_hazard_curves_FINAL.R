#
# Workhorse function to do the revised stage-vs-exceedance-rate percentile
# calculations. Assumes we have already run the preprocessing script (i.e.
# revised_station_hazard_curves_PREPROCESSING.R)
#
#
compute_stage_exceedance_rate_curves<-function(source_name, point_inds_range,
    desired_max_stage_percentiles = c(0.025, 0.16, 0.5, 0.84, 0.975)){

    timer_start = Sys.time()

    mynode = paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')

    library(rptha)

    # Open files with tsunami max-stage information
    nc_tsunami = list()
    nc_tsunami$uniform = nc_open(paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_',
            source_name, '.nc'), readunlim=FALSE)
    nc_tsunami$stochastic = nc_open(paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_',
            source_name, '.nc'), readunlim=FALSE)
    nc_tsunami$variable_uniform = nc_open(paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_variable_uniform_slip_earthquake_events_tsunami_',
            source_name, '.nc'), readunlim=FALSE)


    npoints = diff(point_inds_range) + 1
    # Get max-stage information for the points of interest
    event_max_stage = list()
    for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){
        event_max_stage[[slip_type]] = ncvar_get(nc_tsunami[[slip_type]], 'max_stage', 
                                                 start=c( 1, point_inds_range[1]),
                                                 count=c(-1, npoints), collapse_degen=FALSE)
    }
    # Close netcdf files
    for(i in 1:length(nc_tsunami)) nc_close(nc_tsunami[[i]])

    # Build the output data structure
    output_data = list()
    output_data$point_inds_range = point_inds_range
    output_data$source_name = source_name
    for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){
        output_data[[slip_type]] = list()
        for(mu_type in c('fixed_mu', 'variable_mu')){
            output_data[[slip_type]][[mu_type]] = list()
            # Store 'logic-tree-mean' stage exceedance rate
            output_data[[slip_type]][[mu_type]]$stage_mean_exrate = matrix(NA, nrow=length(STAGE_SEQ), ncol=npoints)
            # Store 'logic-tree-percentile' stage exceedance rates
            output_data[[slip_type]][[mu_type]]$stage_percentile_exrates = 
                array(NA, dim=c(length(desired_max_stage_percentiles), length(STAGE_SEQ), npoints))
        }
    }


    for(point_ind in 1:npoints){
        for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){

            # Dry points will be NA. Causes the code trouble
            if(all(is.na(event_max_stage[[slip_type]][,point_ind]))){
                next
            }

            for(mu_type in c('fixed_mu', 'variable_mu')){

                ## Drop most cases
                #if(slip_type != 'stochastic' | mu_type != 'fixed_mu'){
                #    next
                #}

                # Get the entire family of stage-vs-exceedance-rate curves for each segment
                all_logic_tree_stage_exrates = vector(mode='list', length=length(SEG))
                for(i in 1:length(SEG)){
 
                    if(any(is.na(SEG[[i]]$event_conditional_probability[[slip_type]][[mu_type]]))){
                        stop("NA CONDITIONAL PROBABILITY")
                    }

                    if(any(is.na(SEG[[i]]$event_Mw[[slip_type]]))){
                        stop("NA Mw")
                    }

                    # Use sampled all_rate_matrix 
                    if(mu_type == 'fixed_mu'){
                        all_rate_matrix = SEG[[i]]$sampled_all_rate_matrix_fixedMu
                    }else{
                        all_rate_matrix = SEG[[i]]$sampled_all_rate_matrix_varyMu
                    }


                    all_logic_tree_stage_exrates[[i]] = convert_Mw_vs_exceedance_rates_2_stage_vs_exceedance_rates(
                        SEG[[i]]$Mw_seq,
                        all_rate_matrix,
                        SEG[[i]]$event_Mw[[slip_type]],
                        SEG[[i]]$event_conditional_probability[[slip_type]][[mu_type]],
                        event_max_stage[[slip_type]][,point_ind],
                        STAGE_SEQ, use_Fortran=TRUE, drop_small_events=TRUE)

                    # Quick sanity check
                    if(any( dim(all_logic_tree_stage_exrates[[i]]) != c(nrow(all_rate_matrix), length(STAGE_SEQ)))){
                        print(paste0('Problem with all_logic_tree_stage_exrates[[', i, ']]'))
                        print(paste0(point_ind, ' ', slip_type, ' ', mu_type))
                        print(paste0(point_inds_range, point_ind))
                        stop()
                    }

                }

                #
                # If there is more than 1 segment, we need to combine results carefully
                #
                if(length(SEG) > 1){
                    # SEGMENTED SOURCE 

                    #
                    # Combine segments
                    # Here we do the mean -- that's the simple case
                    #
                    stage_mean_exrate = STAGE_SEQ * 0
                    for(i in 1:length(SEG)){

                        if(mu_type == 'fixed_mu'){
                            seg_weights = SEG[[i]]$sampled_logic_tree_weights_fixedMu
                        }else{
                            seg_weights = SEG[[i]]$sampled_logic_tree_weights_varyMu
                        }

                        # stage exceedance-rates for this segment alone
                        local_mean_exrate = apply(all_logic_tree_stage_exrates[[i]], 2, 
                            f<-function(x) weighted.mean(x, w=seg_weights))

                        # Combine the possibilities that:
                        # A) all segments have their mean rate (given weight=0.5)
                        # B) Unsegmented version has the mean rate (given weight=0.5)
                        stage_mean_exrate = stage_mean_exrate + SEG[[i]]$row_weight * local_mean_exrate
                    }

                    #
                    # Now we do percentiles
                    #
                    # That's the complex case. 
                    #
                    # To do it, let's assume that there are 2 options
                    # A) The segmented models are all true simultaneously. In this
                    #    case we combine their rates assuming co-monotonicity. This has 50% weight
                    # B) The unsegmented model is true. This has 50% weight
                    # For a fixed stage, A) and B) both lead to a distribution of exceedance-rates.
                    #
                    # We then combine the distributions according to their weights to get the final result.
                    # The idea is that the exceedance-rate is like a random sample from a distribution
                    # that has a 50% chance of being A), and a 50% chance of being B)
                    #

                    # To do these calculations, we discretize the distributions at a bunch of percentiles
                    # Note these cover an even 
                    dp = 0.005
                    dense_percentiles = seq(dp/2, 1-dp/2, by=dp)

                    is_a_segment = unlist(lapply(SEG, f<-function(x) x$is_a_segment))
                    segments = which(is_a_segment)
                    unsegment = which(!is_a_segment)
                    if(length(unsegment) > 1) stop("More than one unsegmented branch! Impossible")
                    if(length(segments) < 2) stop("Less than 2 segments on a segmented source! Impossible")

                    # Compute case A. 
                    segmented_exrates_at_percentiles = matrix(0, nrow=length(dense_percentiles), ncol=length(STAGE_SEQ))
                    for(i in segments){

                        if(mu_type == 'fixed_mu'){
                            seg_weights = SEG[[i]]$sampled_logic_tree_weights_fixedMu
                        }else{
                            seg_weights = SEG[[i]]$sampled_logic_tree_weights_varyMu
                        }

                        local_exrates_at_percentiles = apply(all_logic_tree_stage_exrates[[i]], 2,
                            f<-function(x) weighted_percentile(x, weights=seg_weights, p=dense_percentiles))

                        # Co-monotonic -- the exceedance rates at each
                        # "dense_percentile" is true simultaneously on all
                        # logic-tree branches
                        segmented_exrates_at_percentiles = segmented_exrates_at_percentiles + local_exrates_at_percentiles
                    }

                    # Compute case B.
                    # For the unsegmented branch, also convert the logic-tree stage exceedance-rates into
                    # a matrix with the exceedance-rates at dense_percentiles
                    for(i in unsegment){
                        if(mu_type == 'fixed_mu'){
                            seg_weights = SEG[[i]]$sampled_logic_tree_weights_fixedMu
                        }else{
                            seg_weights = SEG[[i]]$sampled_logic_tree_weights_varyMu
                        }
                         
                        unsegmented_exrates_at_percentiles = apply(all_logic_tree_stage_exrates[[i]], 2,
                            f<-function(x) weighted_percentile(x, weights=seg_weights, p=dense_percentiles))
                    }

                    chance_of_unsegmented = SEG[[unsegment]]$row_weight
                    chance_of_segmented = 1 - chance_of_unsegmented

                    # Suppose "stage = STAGE_SEQ[j]"
                    # Then every entry of segmented_exrates_at_percentiles[,j] has the same chance of being correct
                    # (because dense_percentiles is even). This is also true of unsegmented_exrates_at_percentiles[,j]
                    # Thus it is straightforward to combine the distributions
                    combined_exrates_at_percentiles = rbind(segmented_exrates_at_percentiles, unsegmented_exrates_at_percentiles)
                    combined_weights = c(rep(chance_of_segmented, nrow(segmented_exrates_at_percentiles)),
                                       rep(chance_of_unsegmented, nrow(unsegmented_exrates_at_percentiles)))
                    combined_weights = combined_weights/sum(combined_weights)
                    # Finally we can do it
                    stage_percentile_exrates = apply(combined_exrates_at_percentiles, 2,
                        f<-function(x) weighted_percentile(x, weights=combined_weights, p=desired_max_stage_percentiles))

                }else{
                    # NO SEGMENTATION ON THIS SOURCE 
                    # EASY CASE! 

                    i = 1
                    if(mu_type == 'fixed_mu'){
                        seg_weights = SEG[[i]]$sampled_logic_tree_weights_fixedMu
                    }else{
                        seg_weights = SEG[[i]]$sampled_logic_tree_weights_varyMu
                    }

                    # Mean exceedance-rate for a given stage
                    stage_mean_exrate = apply(all_logic_tree_stage_exrates[[i]], 2, 
                        f<-function(x) weighted.mean(x, w=seg_weights))
                    # Percentiles for a given stage
                    stage_percentile_exrates = apply(all_logic_tree_stage_exrates[[i]], 2,
                        f<-function(x) weighted_percentile(x, weights=seg_weights, p=desired_max_stage_percentiles))
                }

                # Pack into our data structure
                output_data[[slip_type]][[mu_type]]$stage_mean_exrate[,point_ind] = stage_mean_exrate
                output_data[[slip_type]][[mu_type]]$stage_percentile_exrates[,,point_ind] = stage_percentile_exrates

            }
        }
    }
    timer_stop = Sys.time()
 
    output_data$time_taken = timer_stop - timer_start
    output_data$mynode = mynode
    return(output_data)
}

#
# Do the calculations for all sources
#

# Get the source names -- here we do it by reading key filesnames and stripping out the source-name
all_source_name_RDS_inputs = Sys.glob('./preprocessed_source_rate_revised_stage_percentiles/preprocessed_rate_info_*.RDS')
tmp = basename(all_source_name_RDS_inputs)
tmp = gsub("preprocessed_rate_info_", "", tmp)
all_source_names = gsub(".RDS", "", tmp, fixed=TRUE)

# Other key input variables
out_dir = 'preprocessed_source_rate_revised_stage_exrates_FULL/'
dir.create(out_dir, showWarnings=FALSE)
point_lower_index = as.numeric(commandArgs(trailingOnly=TRUE)[1]) # 1
point_upper_index = as.numeric(commandArgs(trailingOnly=TRUE)[2]) # 20185
MC_CORES = 16
SAMPLE_LOGIC_TREE = FALSE
NUM_LOGIC_TREE_SAMPLE = 3000 # Ignored if SAMPLE_LOGIC_TREE==FALSE
#

for(source_name in all_source_names){


    # Make a range of points for every core
    library(parallel)
    npts = (point_upper_index - point_lower_index+1)
    all_indices = splitIndices(npts, round(npts/5))
    for(i in 1:length(all_indices)){
        all_indices[[i]] = range(all_indices[[i]]) + point_lower_index - 1
    }

    #
    # Get the inputs
    #
    x = readRDS(paste0('preprocessed_source_rate_revised_stage_percentiles/preprocessed_rate_info_',
                       source_name, '.RDS'))
    # Get the mw-frequency info for all logic-tree curves, and event conditional probabilities, on all segments.
    SEG = x$seg 
    # Get the stages at which we evaluate the results
    STAGE_SEQ = x$stage_seq #2**seq(log(0.1, 2), log(20,2), len=40) #x$stage_seq
    rm(x)
    # SAMPLE FROM THE LOGIC-TREE TO REDUCE NUMBER OF BRANCHES
    for(i in 1:length(SEG)){
        if(SAMPLE_LOGIC_TREE){
        # Sample from the logic-tree
            logic_tree_sample_fixedMu = sample(1:nrow(SEG[[i]]$all_rate_matrix), size=NUM_LOGIC_TREE_SAMPLE,
                replace=TRUE, prob=SEG[[i]]$logic_tree_weights)
            logic_tree_sample_varyMu = sample(1:nrow(SEG[[i]]$all_rate_matrix), size=NUM_LOGIC_TREE_SAMPLE,
                replace=TRUE, prob=SEG[[i]]$logic_tree_weights_varyMu)

            SEG[[i]]$sampled_all_rate_matrix_fixedMu = SEG[[i]]$all_rate_matrix[logic_tree_sample_fixedMu,]
            SEG[[i]]$sampled_all_rate_matrix_varyMu = SEG[[i]]$all_rate_matrix[logic_tree_sample_varyMu,]
            SEG[[i]]$sampled_logic_tree_weights_fixedMu = rep(1/NUM_LOGIC_TREE_SAMPLE, NUM_LOGIC_TREE_SAMPLE)
            SEG[[i]]$sampled_logic_tree_weights_varyMu = rep(1/NUM_LOGIC_TREE_SAMPLE, NUM_LOGIC_TREE_SAMPLE)

            # Cleanup
            rm(logic_tree_sample_fixedMu, logic_tree_sample_varyMu)
        }else{
            SEG[[i]]$sampled_all_rate_matrix_fixedMu = SEG[[i]]$all_rate_matrix
            SEG[[i]]$sampled_all_rate_matrix_varyMu = SEG[[i]]$all_rate_matrix
            SEG[[i]]$sampled_logic_tree_weights_fixedMu = SEG[[i]]$logic_tree_weights
            SEG[[i]]$sampled_logic_tree_weights_varyMu = SEG[[i]]$logic_tree_weights_varyMu
        }
    }

    # Make a wrapper for the function we will run in parallel
    parallel_wrapper<-function(all_indices){
        compute_stage_exceedance_rate_curves(source_name=source_name, point_inds_range=all_indices)
    }

    gc()

    # Do the parallel calculation
    cl = makeForkCluster(nnodes=MC_CORES)
    clusterExport(cl, varlist=setdiff(ls(all=TRUE), 'cl'))
    all_outputs = parLapplyLB(cl, all_indices, parallel_wrapper, chunk.size=1)
    stopCluster(cl)
        
    # Store the main outputs
    output_file = paste0(out_dir, 'stage_exceedance_rate_percentiles_', source_name, '_p_', 
                         point_lower_index, '_', point_upper_index, '.RDS')
    saveRDS(all_outputs, output_file)

    # Also store the sub-sampled logic tree (for debugging)
    output_file2 = paste0(out_dir, 'stage_exceedance_rate_SEG_', source_name, '_p_', 
                         point_lower_index, '_', point_upper_index, '.RDS')
    saveRDS(SEG, output_file2)

    # Cleanup
    rm(all_outputs, SEG)
    gc()
}
