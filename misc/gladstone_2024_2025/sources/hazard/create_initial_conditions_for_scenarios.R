sc = new.env()
source('sampling_config.R', local=sc)

ptha18 = new.env()
source(sc$get_ptha_results_script_file, ptha18, chdir=TRUE)

# Confirm we are on NCI -- necessary for when we copy unit-source tifs below.
is_on_NCI = startsWith(sc$get_ptha_results_script_file, '/g/data/w85/')
stopifnot(is_on_NCI)

# List of environments with scenarios and other useful things
source_zone_data = readRDS('source_zone_data_backup.RDS')

BASEDIR = normalizePath(getwd())

make_tifs_for_source<-function(szdi){
    on.exit(setwd(BASEDIR))

    # Move into the source-zone directory
    setwd(szdi$output_dir)

    source_zone = szdi$source_zone

    if(is_on_NCI){
        # The local caching can cause problems when running initial condition creation in parallel if we don't have all the rasters
        # created already [i.e. if multiple processes try to download the same file at once and corrupt each other].
        # Solve that by creating the rasters right here, with the same folder structure the caching would create.
        # That way no parallel downloads happen
        cache_dir = paste0('SOURCE_ZONES/', source_zone, '/EQ_SOURCE/Unit_source_data/', source_zone)
        dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
        source_file <- szdi$source_zone_scenarios$unit_source_statistics$initial_condition_file
        copy_worked = file.copy(source_file, cache_dir)
        if(!all(copy_worked)) stop('Error when copying unit-source initial conditions')
    }

    unique_scenario_ids = unique(szdi$output_sample$inds)
    unique_scenario_Mws = szdi$output_sample$event_Mw[match(unique_scenario_ids, szdi$output_sample$inds)]

    # Store the initial conditions here
    output_tif_dir = 'scenario_initial_conditions'
    dir.create(output_tif_dir, showWarnings=FALSE)

    # Create the raster for the ith unique_scenario_id
    make_initial_condition_i<-function(i){
        target_row = unique_scenario_ids[i]
        expected_Mw = unique_scenario_Mws[i]

        local_Mw = szdi$source_zone_scenarios$events$Mw[target_row]
        if(round(local_Mw, 1) != round(expected_Mw, 1)) stop('Logic error -- these magnitudes should be identical')

        # Make a row index character, padded with leading zeros as required, for consistency
        # of file names
        target_row_fixed_width = substring(as.character(1e+07 + target_row), 2, 8)

        output_filename = paste0(output_tif_dir, '/', source_zone, '_row_', 
                                 target_row_fixed_width, 
                                 '_Mw_', round(local_Mw*10), '_HS.tif')

        r1 = ptha18$get_initial_condition_for_event(szdi$source_zone_scenarios, event_ID=target_row)

        writeRaster(r1, output_filename, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

        rm(r1); gc()
        return(0)
    }

    library(parallel)
    tmp = mclapply(1:length(unique_scenario_ids), make_initial_condition_i, mc.cores=sc$MC_CORES)

}

making_all_tifs = lapply(source_zone_data, make_tifs_for_source)

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
