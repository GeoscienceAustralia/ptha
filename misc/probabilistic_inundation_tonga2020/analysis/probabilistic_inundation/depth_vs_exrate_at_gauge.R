#
# Get the max-depth vs AEP at a given gauge, for each representation of the source-zone scenario rates.
# This is useful both to get results at known points of interest, and for testing. Actually there are
# two useful tests that can be done: 
# A) Cross-check with the raster exceedance-rate computations at known points. This is good because the
#    latter rely on the model max-depth rasters, whereas calculations here rely on the gauge points. Further
#    the probabilistic computations are quite separate [except they rely on the same scenarios], so it's a 
#    situation where we can 'do the same thing twice in two different ways', to sanity check the computations.
# B) Cross-check the max-stage-exceedance rate at a PTHA18 point computed with this model, vs the original PTHA18.
#    That's a good check because to get good agreement, a number of things need to be true:
#    - The 'random scenarios with importance sampling' for this study need to well represent the PTHA18 results
#    - The hydrodynamic modelling in this study needs to agree well with the PTHA18 results, notwithstanding they
#      are run on different grids, with different model resolutions, run-times
#      and solver options, and different input data.
#    While the above factors should lead to slight differences, we find good agreement, which
#    suggests everything is working as intended.
#

library(rptha)

file_nci = '../../../CODE/ptha/propagation/SWALS/plot.R'
file_home = '~/Code_Experiments/fortran/Structured_shallow_water/plot.R'
swals = new.env()
if(file.exists(file_nci)){
    source(file_nci, local=swals)
}else{
    source(file_home, local=swals)
}

#
# INPUT 
#

# Name of the folder inside ../../swals/OUTPUTS/ where the 'batch' of tsunami
# model runs is stored. 
run_series_name = commandArgs(trailingOnly=TRUE)[2] #'ptha18_tonga_MSL0'

# Vector with all model run directories [one per multidomain]
md_dirs = Sys.glob(paste0('../../swals/OUTPUTS/', run_series_name, '/ptha18_random_scenarios_kermadectonga2_row_*/RUN*'))

# Useful flag to control the cases to run
run_case = commandArgs(trailingOnly=TRUE)[1]
if(run_case == 'parliament'){
    # Case 1 -- Parliament house
    TARGET_GAUGE = matrix(c(184.8018, -21.1331), nrow=1)
    output_rds_name = paste0(run_series_name, '_depth_and_stage_exrate_curve_at_parliament.RDS')
}else if(run_case == 'ptha18_point_3458.3'){
    # Case 2 -- PTHA18 point
    TARGET_GAUGE = matrix(c(185.123916, -21.088881), nrow=1)
    output_rds_name = paste0(run_series_name, '_depth_and_stage_exrate_curve_at_ptha18_point_3458.3.RDS')
}else{
    stop(paste0('unknown value of run_case=', run_case))
}

# How many cores to use in parallel
MC_CORES = 48

# All files with randomly selected PTHA18 scenarios and their rates. These
# were created with the importance-sampling routine.
scenarios_data_files = c(
    '../../sources/random/random_scenarios_kermadectonga2_hukurangi_segment_HS.csv',
    '../../sources/random/random_scenarios_kermadectonga2_kermadec_segment_HS.csv',
    '../../sources/random/random_scenarios_kermadectonga2_tonga_segment_HS.csv',
    '../../sources/random/random_scenarios_kermadectonga2_unsegmented_HS.csv')
names_scenarios_databases = gsub('.csv', '', 
    gsub('random_scenarios_kermadectonga2_', '', basename(scenarios_data_files)))

# Give a vector of scenario_rows (i.e. corresponding to rows in the PTHA18
# events table for the source-zone), this function should return a
# corresponding vector of SWALS multidomain directories that contain the
# tsunami model run. 
# Because the matching might change from project to project, we make this an
# input argument
get_md_dirs_matching_scenario_rows<-function(scenario_rows, md_dirs){
    # For example, for scenario 12345, the multidomain_dir will contain the string '_0012345_',
    # and we match that. 
    row_string = paste0('_', substring(as.character(1e+07 + scenario_rows), 2, 8), '_Mw')
    match_md = sapply(row_string, f<-function(x){
        out = grep(x, md_dirs)
        if(length(out) != 1) out = NA 
        return(out)
        })
    if(any(is.na(match_md))) stop('problem matching md_dirs')
    return(md_dirs[match_md])
}


#
# END INPUTS
#

# Find the domain containing the gauge of interest.  Here we assume all md_dirs
# use the same model setup -- if not, the following lines need to go inside
# get_max_depth_at_gauge(md_dir)
domain_info = swals$find_domain_containing_point(TARGET_GAUGE, multidomain_dir=md_dirs[1]) 
DOMAIN_INDEX = domain_info$domain_index

# For a given model run [stored in md_dir], get the max-depth at the
# TARGET_GAUGE, using the Gauges_data_ID*.nc file.
get_max_depth_at_gauge<-function(md_dir){

    # Workaround for a case where, to make the model stable, I did fewer parallel partitions.
    troublesome_md_dir = 'ptha18_random_scenarios_kermadectonga2_row_0043831_Mw_95_HS-risetime_0-ambientsealevel_0.8-full-linear_with_manning-0.035-highres_tonga'
    if(basename(dirname(md_dir)) == troublesome_md_dir){
        # Need to compute new DOMAIN_INDEX for this case
        domain_info = swals$find_domain_containing_point(TARGET_GAUGE, multidomain_dir=md_dir) 
        local_domain_index = domain_info$domain_index
    }else{
        local_domain_index = DOMAIN_INDEX
    }
    all_domains = Sys.glob(paste0(md_dir, '/RUN*'))
    gauge_file = Sys.glob(paste0(all_domains[local_domain_index], '/Gauges_data_ID*.nc'))
    if(length(gauge_file) != 1) stop(paste0('did not find a unique gauge_file for ', md_dir))

    # The file will probably contain data for more than one gauge. Read it all here.
    fid = nc_open(gauge_file, readunlim=FALSE)
    time = ncvar_get(fid, 'time')
    lon = ncvar_get(fid, 'lon')
    lat = ncvar_get(fid, 'lat')
    elev = ncvar_get(fid, 'elevation0')
    stage = ncvar_get(fid, 'stage')
    nc_close(fid)

    # Figure out which gauge is nearest to the TARGET_GAUGE
    nearest_point = which.min(distHaversine(
        cbind(lon*0+TARGET_GAUGE[1], lat*0+TARGET_GAUGE[2]), 
        cbind(lon                  , lat                  )))

    # Store results at the nearest point
    max_stage_ind = which.max(stage[nearest_point,])
    output = data.frame(lon=lon[nearest_point], lat=lat[nearest_point], elev=elev[nearest_point], 
               max_stage=stage[nearest_point,max_stage_ind], time_of_max_stage=time[max_stage_ind])

    rm(all_domains, gauge_file, fid, time, lon, lat, elev, stage); gc()
    return(output)
}

# Start the cluster
library(parallel)
local_cluster = makeCluster(MC_CORES)
clusterCall(local_cluster, fun=function(){library(rptha)})
clusterExport(local_cluster, ls(all=TRUE))

# Read all the files and combine results into a data.frame
results_all_md_dirs = parLapply(local_cluster, md_dirs, get_max_depth_at_gauge)
results_df = do.call(rbind, results_all_md_dirs)

# Parallel part is finished.
stopCluster(local_cluster)

# Quick check that all the points are the same
stopifnot(all(results_df$lon == results_df$lon[1]))
stopifnot(all(results_df$lat == results_df$lat[1]))

# Read all the scenario metadata [unsegmented + segments]
scenarios_databases = lapply(scenarios_data_files, read.csv)
names(scenarios_databases) = names_scenarios_databases

# Find the multidomain_dir associated with each scenario
for(i in 1:length(scenarios_databases)){
    scenarios_databases[[i]]$md_dir = get_md_dirs_matching_scenario_rows(
        scenarios_databases[[i]]$scenario_row, md_dirs)
}

# Compute stage-vs-exceedanceRate and depth-vs-exceedanceRate curves
stage_exrate_curves = vector(mode='list', length=length(scenarios_databases))
names(stage_exrate_curves) = names_scenarios_databases
for(sd_names in names_scenarios_databases){
   
    # Populate scenario_rates based on the rates assigned to random scenarios
    # in the scenarios_databases 
    scenario_rates = rep(0, length=length(md_dirs))
    alternate_scenario_rates = rep(0, length=length(md_dirs))
    local_scenarios = scenarios_databases[[sd_names]]

    matching_inds = match(local_scenarios$md_dir, md_dirs)
    for(i in 1:length(matching_inds)){
        scenario_rates[matching_inds[i]] = scenario_rates[matching_inds[i]] + 
            local_scenarios$scenario_rates[i]
        alternate_scenario_rates[matching_inds[i]] = alternate_scenario_rates[matching_inds[i]] + 
            local_scenarios$alternate_scenario_rates[i]
    }

    # Depth vs exceedance-rate
    max_depths = results_df$max_stage - results_df$elev
    #chosen_depths = c(0.001, seq(0.05, 20, by=0.05))
    chosen_depths = seq(0.001, 20, by=0.001)
    exrate_of_chosen_depths = sapply(chosen_depths, f<-function(x) sum(scenario_rates * (max_depths > x)))
    alternate_exrate_of_chosen_depths = sapply(chosen_depths, f<-function(x) sum(alternate_scenario_rates * (max_depths > x)))

    # Stage vs exceedance-rate
    chosen_stages = chosen_depths
    exrate_of_chosen_stages = sapply(chosen_stages, f<-function(x) sum(scenario_rates * (results_df$max_stage > x)))
    alternate_exrate_of_chosen_stages = sapply(chosen_stages, f<-function(x) sum(alternate_scenario_rates * (results_df$max_stage > x)))

    stage_exrate_curves[[sd_names]] = list(
        depth=chosen_depths, exrate_depth=exrate_of_chosen_depths, alternate_exrate_depth=alternate_exrate_of_chosen_depths,
        stage=chosen_stages, exrate_stage=exrate_of_chosen_stages, alternate_exrate_stage=alternate_exrate_of_chosen_stages,
        source_name=sd_names, target_point=TARGET_GAUGE, 
        output_point = c(results_df$lon[1], results_df$lat[1]))
}

# Store these for later analysis
outputs = list(stage_exrate_curves=stage_exrate_curves, results_df = results_df)
saveRDS(outputs, output_rds_name)
