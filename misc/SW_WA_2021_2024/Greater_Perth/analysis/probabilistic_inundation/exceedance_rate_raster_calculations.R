library(raster)

# Get PTHA routines
ptha18 = new.env()
.file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
.file_home = '../../../../../../AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
source(ifelse(file.exists(.file_nci), .file_nci, .file_home), local=ptha18, chdir=TRUE)


#' Given max-depth matrices and scenario rates, compute rate that a
#' depth_threshold is exceeded. 
#' 
#' NA values in the max-depth matrix will be treated as dry.
#'
#' @param included_indices a vector of non-repeated indices in
#' 1:length(scenario_rates) giving the rasters to include. This CAN be used for
#' splitting the calculation in parallel (sum chunks in parallel, then sum the 
#' final result in serial). A serial calculation would use 1:length(scenario_rates).
#' It often may be more efficient to use coarser-grained parallelism (i.e. calling
#' this routine in serial, with different cores treating different domains).
#' @param max_depth_files A list of rasters containing the max_depth (one for
#' each entry of tarred_multidomain_dirs).
#' @param scenario_rates A vector with the individual scenario rates for each
#' entry of max_depth_files
#' @param depth_threshold The function will compute the exceedance rate of
#' (depth > depth_threshold), where 'depth' might be some other quantity (whatever is
#' stored in the max_depth_files).
#' @param print_progress Print index as each raster is processed
#'
get_exceedance_rate_at_threshold_depth<-function(included_indices, 
    max_depth_files, scenario_rates, depth_threshold, print_progress=FALSE){

    stopifnot(length(scenario_rates) == length(max_depth_files))

    stopifnot(length(included_indices) == length(unique(included_indices)))

    stopifnot( (min(included_indices) >= 1) & 
               (max(included_indices) <= length(max_depth_files)) )

    stopifnot(length(depth_threshold) == 1)

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
        if(print_progress) print(i)
        # Read the raster
        #x_mat = as.matrix(terra::rast(local_max_depth_files[i]), wide=TRUE)
        x_mat = as.matrix(raster(local_max_depth_files[i]))
        stopifnot(all(dim(x_mat) == dim(ex_rate)))

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

#' Exceedance-rates on one tile, for one source-zone representations.
#'
#' Compute 'exceedance-rate' raster for the given domain index and the
#' given 'scenarios_name' (which corresponds either to an unsegmented source-zone representation, 
#' or a segment) 
#'
#' @param input_domain_index_and_scenarios_name List of the form list(domain_index=3, scenarios_name='unsegmented_HS'),
#' where the domain_index is the integer in the domain tif file (inside a tar archive), and scenarios_name is one
#' entry of names(scenario_databases)
#' @param tarred_multidomain_dirs Vector of tarred multidomain directories
#' @param scenario_databases List with scenario databases
#' @param output_directory Directory in which tiffs will be saved, relative to working directory.
#' @param depth_thresholds_for_exceedance_rate_calculations Vector of depths for which we compute exceedance-rates.
#' @param raster_name_start leading characters of the rasters to be read in
#' each tarred raster dir, before the domain_index, e.g.'max_stage_domain_' 
#' @param print_progress_debug Print index as each raster is processed
#' @return Zero on successful completions -- but as a side-effect it creates the exceedance-rate raster.
compute_exceedance_rates_on_tile<-function(
    input_domain_index_and_scenarios_name,
    tarred_multidomain_dirs, 
    scenario_databases, 
    output_directory, 
    depth_thresholds_for_exceedance_rate_calculations, 
    raster_name_start,
    print_progress_debug=FALSE){

    # To ease parallel calculation we take the input parameters in a list, and
    # unpack here.
    domain_index = input_domain_index_and_scenarios_name$domain_index
    scenarios_name = input_domain_index_and_scenarios_name$scenarios_name

    # Useful to have one of the rasters ready as a template [to help with data
    # The raster depth file associated with each md_dir. There should only be one
    # per md_dir. All rasters must have the same extent and resolution.
    raster_files_one_domain = paste0('/vsitar/', dirname(tarred_multidomain_dirs), 
        "/raster_output_files.tar/", raster_name_start, 
        domain_index, ".tif")

    raster_template = raster(raster_files_one_domain[1])
    #raster_template = terra::rast(raster_files_one_domain[1])

    # This text will appear in the filename of all output rasters
    output_raster_name_tag = paste0('domain_', domain_index, '_', raster_name_start)

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

    t0 = sum(scenario_rates)
    t1 = sum(scenario_databases[[scenarios_name]]$importance_sampling_scenario_rates_basic)
    errmess = 'Scenario rates do not sum to the same value: Issue with bookkeeping or file naming'
    if(!isTRUE(all.equal(t0, t1))) stop(errmess)

    # For each depth-threshold, make the exceedance-rate raster. NOTE: Could
    # make calculations parallel over depth_threshold too.
    for(depth_threshold in depth_thresholds_for_exceedance_rate_calculations){

        tile_exceedance_rates = get_exceedance_rate_at_threshold_depth(
            included_indices = seq(1, length(raster_files_one_domain)),
            max_depth_files = raster_files_one_domain,
            scenario_rates = scenario_rates,
            depth_threshold=depth_threshold,
            print_progress=print_progress_debug)


        # For the raster output, it is nice to set regions that are never
        # inundated to NA (genuinely NA regions that are not priority
        # domain will also be NA)
        tile_exceedance_rates[tile_exceedance_rates == 0] = NA

        # Convert to a raster and write to file
        #exrates_rast = terra::setValues(raster_template, tile_exceedance_rates)
        exrates_rast = setValues(raster_template, tile_exceedance_rates)

        raster_output_file = paste0(output_directory, '/', 
            scenarios_name, '_', output_raster_name_tag, 
            '_exceedance_rate_with_threshold_', depth_threshold, 
            '.tif')

        #terra::writeRaster(exrates_rast, raster_output_file, 
        writeRaster(exrates_rast, raster_output_file, 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        rm(exrates_rast, tile_exceedance_rates)
        gc()
    }

    rm(raster_template); gc()
    return(0)
}


#' Given max-depth rasters and scenario rates, compute rate that a
#' depth_threshold is exceeded, and estimate the variance of the Monte-Carlo
#' error in that rate.
#' 
#' NA values in the max-depth matrix will be treated as dry. Note in this
#' routine the inputs MUST be arranged by SAMPLED SCENARIOS. So there may be
#' repeated scenario_rasters. This arrangement is suitable for computing the 
#' error variance.
#'
#' @param max_depth_files A list of rasters containing the max_depth (one for
#' each sampled scenario -- there can be repeats since we are sampling with replacement).
#' @param scenario_rates A vector with the individual scenario rates for each
#' entry of max_depth_files
#' @param scenario_magnitudes the magnitudes of each scenario.
#' @param depth_threshold The function will compute the exceedance rate of
#' (depth > depth_threshold), where 'depth' might be some other quantity (whatever is
#' stored in the max_depth_files).
#' @param print_progress Print index as each raster is processed
#'
get_exceedance_rate_and_error_variance_at_threshold_depth<-function(
    max_depth_files, scenario_rates, scenario_magnitudes, depth_threshold,
    print_progress=FALSE){

    # Logical checks on arguments
    stopifnot(length(scenario_rates) == length(max_depth_files))
    stopifnot(length(scenario_rates) == length(scenario_magnitudes))
    stopifnot(length(depth_threshold) == 1)

    # Get unique Mw values, confirming they are evenly spaced
    unique_Mw = ptha18$unique_sorted_with_check_for_even_spacing(scenario_magnitudes)

    local_max_depth_files = max_depth_files
    local_scenario_rates = scenario_rates

    # Make space to store exceedance-rate, and exceedance_rate variance
    r1 = as.matrix(raster(local_max_depth_files[1]))
    target_dim = dim(r1)
    rm(r1)
    ex_rate = matrix(0, ncol=target_dim[2], nrow=target_dim[1])
    ex_rate_variance = ex_rate

    for(mw in unique_Mw){
        k = which(scenario_magnitudes == mw)
        num_Mw = length(k) # Number of sceanrios in this magnitude bin

        # Within-bin values are needed for intermediate calculations
        ex_rate_in_bin = matrix(0, ncol=target_dim[2], nrow=target_dim[1])
        ex_rate_variance_in_bin = ex_rate_in_bin

        # Compute the within bin exceedance-rate, Eq 17 in Davies et al. 2022
        for(i in k){
            if(print_progress) print(i)
            # Read the raster
            x_mat = as.matrix(raster(local_max_depth_files[i]))
            stopifnot(all(dim(x_mat) == dim(ex_rate)))

            # Make a raster that is 1 where we exceed the depth threshold, and 0 elsewhere
            # (including at NA sites)
            out = 1 - (x_mat <= depth_threshold)
            out[is.na(out)] = 0
            # Sum the exceedance-rates
            ex_rate_in_bin = ex_rate_in_bin + local_scenario_rates[i] * out
        }
        # Add the rates in this magnitude bin to the total rates
        ex_rate = ex_rate + ex_rate_in_bin

        # Compute the within-bin variance of the error in the exceedance-rate
        # Eq 20 in Davies et al. 2022
        for(i in k){
            if(print_progress) print(i)
            # Read the raster
            x_mat = as.matrix(raster(local_max_depth_files[i]))
            stopifnot(all(dim(x_mat) == dim(ex_rate)))

            # Make a raster that is 1 where we exceed the depth threshold, and 0 elsewhere
            # (including at NA sites)
            out = 1 - (x_mat <= depth_threshold)
            out[is.na(out)] = 0

            # This is rearranged from Eq 20 in Davies et al., 2022
            ex_rate_variance_in_bin = ex_rate_variance_in_bin + 
                (local_scenario_rates[i] * out - ex_rate_in_bin/num_Mw)^2
        }
        # Add the variance in this bin to the total variance
        ex_rate_variance = ex_rate_variance + ex_rate_variance_in_bin
    
        gc()
    }

    return(list(exrate=ex_rate, exrate_var=ex_rate_variance))
}

#' Estimated Exceedance-rates and Monte-Carlo variance of the estimate, on one
#' tile, for one source-zone representations.
#'
#' Compute 'exceedance-rate' raster for the given domain index and the
#' given 'scenarios_name' (which corresponds either to an unsegmented source-zone representation, 
#' or a segment). This gives an estimate of the true exceedance rate. Also
#' compute an estimate of the variance of the exceedance-rate error. This corresponds
#' to Equations 17 and 20 from 
#'     Davies, G.; Weber, R.; Wilson, K. & Cummins, P. 
#'     From offshore to onshore probabilistic tsunami hazard assessment via efficient Monte-Carlo sampling, 2021
#'
#' @param input_domain_index_and_scenarios_name List of the form list(domain_index=3, scenarios_name='unsegmented_HS'),
#' where the domain_index is the integer in the domain tif file (inside a tar archive), and scenarios_name is one
#' entry of names(scenario_databases)
#' @param tarred_multidomain_dirs Vector of tarred multidomain directories
#' @param scenario_databases List with scenario databases
#' @param output_directory Directory in which tiffs will be saved, relative to working directory.
#' @param depth_thresholds_for_exceedance_rate_calculations Vector of depths for which we compute exceedance-rates.
#' @param raster_name_start leading characters of the rasters to be read in
#' each tarred raster dir, before the domain_index, e.g.'max_stage_domain_' 
#' @param print_progress_debug Print index as each raster is processed
#' @return Zero on successful completions -- but as a side-effect it creates the exceedance-rate raster.
compute_exceedance_rates_and_error_variance_on_tile<-function(
    input_domain_index_and_scenarios_name,
    tarred_multidomain_dirs, 
    scenario_databases, 
    output_directory, 
    depth_thresholds_for_exceedance_rate_calculations, 
    raster_name_start,
    print_progress_debug=FALSE){

    # To ease parallel calculation we take the input parameters in a list, and
    # unpack here.
    domain_index = input_domain_index_and_scenarios_name$domain_index
    scenarios_name = input_domain_index_and_scenarios_name$scenarios_name

    # Useful to have one of the rasters ready as a template [to help with data
    # The raster depth file associated with each md_dir. There should only be one
    # per md_dir. All rasters must have the same extent and resolution.
    raster_files_one_domain = paste0('/vsitar/', dirname(tarred_multidomain_dirs), 
        "/raster_output_files.tar/", raster_name_start, 
        domain_index, ".tif")

    raster_template = raster(raster_files_one_domain[1])
    #raster_template = terra::rast(raster_files_one_domain[1])

    # This text will appear in the filename of all output rasters
    output_raster_name_tag = paste0('domain_', domain_index, '_', raster_name_start)
    output_variance_raster_name_tag = paste0('domain_', domain_index, '_', raster_name_start, '_variance_of_')

    # Map the rows of the database to the rasters
    #ind = match(dirname(scenario_databases[[scenarios_name]]$md_dir), 
    #            gsub('/vsitar/', '', dirname(dirname(raster_files_one_domain)), fixed=TRUE) )

    # Map the rasters to rows of the database
    ind = match(dirname(scenario_databases[[scenarios_name]]$md_dir),
                gsub('/vsitar/', '', dirname(dirname(raster_files_one_domain)), fixed=TRUE))
    scenario_rasters = raster_files_one_domain[ind]
    scenario_rates = scenario_databases[[scenarios_name]]$importance_sampling_scenario_rates_basic
    scenario_mw    = scenario_databases[[scenarios_name]]$mw

    # For each depth-threshold, make the exceedance-rate raster. NOTE: Could
    # make calculations parallel over depth_threshold too.
    for(depth_threshold in depth_thresholds_for_exceedance_rate_calculations){

        tile_exceedance_rates_and_error_variance = get_exceedance_rate_and_error_variance_at_threshold_depth(
            max_depth_files = scenario_rasters,
            scenario_rates = scenario_rates,
            scenario_magnitudes = scenario_mw,
            depth_threshold=depth_threshold,
            print_progress=print_progress_debug)

        # For the raster output, it is nice to set regions that are never
        # inundated to NA (genuinely NA regions that are not priority
        # domain will also be NA)
        null_regions = (tile_exceedance_rates_and_error_variance$exrate == 0)
        tile_exceedance_rates_and_error_variance$exrate[null_regions] = NA
        tile_exceedance_rates_and_error_variance$exrate_var[null_regions] = NA
        rm(null_regions)

        # Convert to a raster and write to file

        #exrates_rast = terra::setValues(raster_template, tile_exceedance_rates_and_error_variance)
        exrates_rast = setValues(raster_template, tile_exceedance_rates_and_error_variance$exrate)
        raster_output_file = paste0(output_directory, '/', 
            scenarios_name, '_', output_raster_name_tag, 
            '_exceedance_rate_with_threshold_', depth_threshold, 
            '.tif')
        #terra::writeRaster(exrates_rast, raster_output_file, 
        writeRaster(exrates_rast, raster_output_file, 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

        exrates_var_rast = setValues(raster_template, tile_exceedance_rates_and_error_variance$exrate_var)
        raster_output_file = paste0(output_directory, '/', 
            scenarios_name, '_', output_variance_raster_name_tag, 
            '_exceedance_rate_with_threshold_', depth_threshold, 
            '.tif')
        #terra::writeRaster(exrates_rast, raster_output_file, 
        writeRaster(exrates_var_rast, raster_output_file,
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

        rm(exrates_rast, exrates_var_rast, tile_exceedance_rates_and_error_variance)
        gc()
    }

    rm(raster_template); gc()
    return(0)
}

# Quick tests of the exceedance-rate calculations
.test_exceedance_rate_raster_calculations<-function(){
    #
    # Compare the max-stage exceedance-rate against an independently
    # calculated value at a point, using code from '../max_stage_at_a_point'
    #
    # This test could be much improved -- it relies on a large set of existing simulated data.
    #
    input_domain_index_and_scenarios_name = list(
        domain_index = 39, # Compare against a site on domain 39
        scenarios_name = 'logic_tree_mean_curve_HS')

    scenario_databases = list()
    scenario_databases$logic_tree_mean_curve_HS = read.csv(
        '../../sources/hazard/random_sunda2/random_scenarios_sunda2_logic_tree_mean_curve_HS.csv')

    .local_find_matching_md_dir<-function(row_indices, tarred_multidomain_dirs){
        # Make a string with the start of the SWALS output folder name (beneath
        # ../../swals/OUTPUTS/...)
        matching_string = paste0('ptha18_random_scenarios_sunda2_row_', 
            substring(as.character(1e+07 + row_indices), 2, 8), '_')

        # Match with the tarred_multidomain_dirs, with NA if we don't match or get multiple matches
        matching_ind = sapply(matching_string, f<-function(x){
            p = grep(x, tarred_multidomain_dirs)
            if(length(p) != 1) p = NA 
            return(p)})
        if(any(is.na(matching_ind))) stop('Could not find simulation matching scenario')

        return(tarred_multidomain_dirs[matching_ind])
    }

    #
    tarred_multidomain_dirs = Sys.glob(
        '../../swals/OUTPUTS/ptha18-GreaterPerth-sealevel60cm/random_sunda2/ptha*/RUN*.tar.bz2')

    for(i in 1:length(scenario_databases)){
        scenario_databases[[i]]$md_dir = .local_find_matching_md_dir(
            scenario_databases[[i]]$inds, tarred_multidomain_dirs)
    }


    # Get "max-stage-exceedance' of "MSL + 1.0"
    MSL = 0.6
    depth_thresholds_for_exceedance_rate_calculations = MSL + 1.0

    output_directory = 'test_dir'
    raster_name_start = 'max_stage_domain_'
    dir.create(output_directory, showWarnings=FALSE)


    #
    # Test 1 -- basic exceedance-rate calculation
    #
    make_raster = compute_exceedance_rates_on_tile(
        input_domain_index_and_scenarios_name,
        tarred_multidomain_dirs,
        scenario_databases,
        output_directory,
        depth_thresholds_for_exceedance_rate_calculations,
        raster_name_start)

    x = raster('test_dir/logic_tree_mean_curve_HS_domain_39_max_stage_domain__exceedance_rate_with_threshold_1.6.tif')
    result = extract(x, matrix(c(113.0708085, -28.5679802), ncol=2))
    # Previously I computed the exceedance-rates at this site using very different code, see here:
    #    ../max_stage_at_point/
    # Visually the exceedance rate (for 1m above MSL) is very close to 1.0e-04.  
    # Actually it is 9.800535e-05
    if(abs(result - 9.800535e-05) < 1e-10){
        print('PASS')
    }else{
        print('FAIL')
    }

    #
    # Test 2 -- exceedance-rate with variance. Here the exceedance-rate should
    # be the same as the above up to floating-point differnces due to
    # reordering of a sum, and we also get the error variance
    output_directory = 'test_dir_2'
    dir.create(output_directory, showWarnings=FALSE)

    make_raster = compute_exceedance_rates_and_error_variance_on_tile(
        input_domain_index_and_scenarios_name,
        tarred_multidomain_dirs,
        scenario_databases,
        output_directory,
        depth_thresholds_for_exceedance_rate_calculations,
        raster_name_start)

    xnew = raster('test_dir_2/logic_tree_mean_curve_HS_domain_39_max_stage_domain__exceedance_rate_with_threshold_1.6.tif')
    xnew_var = raster('test_dir_2/logic_tree_mean_curve_HS_domain_39_max_stage_domain__variance_of__exceedance_rate_with_threshold_1.6.tif')
    result = extract(xnew, matrix(c(113.0708085, -28.5679802), ncol=2))
    # Previously I computed the exceedance-rates at this site using very different code, see here:
    #    ../max_stage_at_point/
    # Visually the exceedance rate (for 1m above MSL) is very close to 1.0e-04.  
    # Actually it is 9.800535e-05
    if(abs(result - 9.800535e-05) < 1e-10){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Check the exceedance-rates are just like before
    if(all(abs(as.matrix(xnew - x)) <= 1.0e-08*as.matrix(xnew), na.rm=TRUE)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # From a separate calculation, the variance should be
    # [1] 2.37206e-10
    result = extract(xnew_var, matrix(c(113.0708085, -28.5679802), ncol=2))
    if(abs(result - 2.37206e-10) < 1.0e-14){
        print('PASS')
    }else{
        print('FAIL')
    }

}

