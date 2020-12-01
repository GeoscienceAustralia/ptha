#
# Get the PTHA18 access codes, needed for functions below
#
library(rptha)

file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
is_on_nci = file.exists(file_nci)
ptha18 = new.env()
source(ifelse(is_on_nci, file_nci, file_home), local=ptha18, chdir=TRUE)

# Get the full PTHA18 source-zone events table, unit-source data, etc
source_zone = 'kermadectonga2'
kt_HS_scenarios = ptha18$get_source_zone_events_data(source_zone, slip_type='stochastic')

# Get the random scenarios
sampled_scenarios = lapply(
    c("random_scenarios_kermadectonga2_hukurangi_segment_HS.csv", 
      "random_scenarios_kermadectonga2_kermadec_segment_HS.csv", 
      "random_scenarios_kermadectonga2_tonga_segment_HS.csv",
      "random_scenarios_kermadectonga2_unsegmented_HS.csv"),
    read.csv)

all_scenarios = do.call(rbind, sampled_scenarios)

if(is_on_nci){
    # The local caching can cause problems when running in parallel if we don't have all the rasters
    # created already [i.e. if multiple processes try to download the same file at once and corrupt each other].
    # Solve that by creating the rasters right here, with the same folder structure the caching would create.
    # That way no parallel downloads happen
    cache_dir = paste0('SOURCE_ZONES/', source_zone, '/EQ_SOURCE/Unit_source_data/', source_zone)
    dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
    copy_worked = file.copy(kt_HS_scenarios$unit_source_statistics$initial_condition_file, cache_dir)
    if(!all(copy_worked)) stop('Error when copying unit-source initial conditions')
}



#
# For each unique scenario, generate the initial condition with an appropriate
# filename
#
unique_scenario_ids = unique(all_scenarios$scenario_row)
unique_scenario_Mws = all_scenarios$mw[match(unique_scenario_ids, all_scenarios$scenario_row)]

# Store the initial conditions here
output_tif_dir = 'scenario_initial_conditions'
dir.create(output_tif_dir, showWarnings=FALSE)

# Create the raster for the ith unique_scenario_id
make_initial_condition_i<-function(i){
    target_row = unique_scenario_ids[i]
    expected_Mw = unique_scenario_Mws[i]

    local_Mw = kt_HS_scenarios$events$Mw[target_row]
    if(round(local_Mw, 1) != round(expected_Mw, 1)) stop('Logic error -- these magnitudes should be identical')

    # Make a row index character, padded with leading zeros as required, for consistency
    # of file names
    target_row_fixed_width = substring(as.character(1e+07 + target_row), 2, 8)

    output_filename = paste0(output_tif_dir, '/', source_zone, '_row_', 
                             target_row_fixed_width, 
                             '_Mw_', round(local_Mw*10), '_HS.tif')

    r1 = ptha18$get_initial_condition_for_event(kt_HS_scenarios, event_ID=target_row)

    writeRaster(r1, output_filename, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(r1); gc()
}

library(parallel)
MC_CORES=48
mclapply(1:length(unique_scenario_ids), make_initial_condition_i, mc.cores=MC_CORES)

# Eyeball the sources
plot_all_initial_conditions<-function(){

    all_ic = Sys.glob(paste0(output_tif_dir, '/*.tif'))
    pdf('all_initial_conditions_plot.pdf', width=6, height=9)
    for(i in 1:length(all_ic)){
        r1 = raster(all_ic[i])
        plot(r1)
        title(basename(all_ic[i]))
    }
    dev.off()

}
