#
# Code to get 'good fitting' scenario for a bunch of historical events from the
# PTHA18 results {see GA Record for explanation}
#

# See ptha/ptha_access/README.md for a tutorial on using these functions

# Get the ptha_access functions (from the ptha repository)
source('../../get_PTHA_results.R', chdir=TRUE)

# These are hard-coded lists with our "best" scenarios, created with find_desired_event_rows.R
# There can be repetition [especially for VAUS, as it's easy to randomly repeat VAUS rupture]
source('best_fitting_HS.R')
source('best_fitting_VAUS.R')

for(slip_type in c('stochastic', 'variable_uniform')){

    if(slip_type == 'stochastic'){
        scenarios_like_events = scenarios_like_events_HS
    }else if(slip_type == 'variable_uniform'){
        scenarios_like_events = scenarios_like_events_VAUS
    }

    source_zone_data = vector(mode='list', length=length(scenarios_like_events))
    for(i in 1:length(scenarios_like_events)){
        print(paste0('Working on ', scenarios_like_events[[i]]$stats_file))

        # Tag for the filenames
        if(grepl('varyMu', scenarios_like_events[[i]]$stats_file)){
            name_info = paste0(slip_type, '_varyMu')
        }else{
            name_info = slip_type
        }

        # Get the scenario metadata for the desired events
        source_zone = scenarios_like_events[[i]]$source_zone
        desired_events = scenarios_like_events[[i]]$desired_event_rows

        # Keep trying the download if it timesout
        source_zone_data[[i]] = try(log('a'), silent=TRUE)
        counter = 0
        while(is(source_zone_data[[i]], 'try-error')){
            source_zone_data[[i]] = try(get_source_zone_events_data(source_zone, 
                slip_type=slip_type,
                desired_event_rows = desired_events))
            counter = counter + 1
            if(counter >= 10) stop('Too many failed downloads, abort')
        }

        deformation_rasters = vector(mode='list', length=length(desired_events))

        # Get the deformation rasters for the desired_events
        for(j in 1:length(desired_events)){
            deformation_rasters[[j]] = get_initial_condition_for_event(
                source_zone_data[[i]], j)
        }

        # For convenience, store the rasters inside the source_zone_data
        source_zone_data[[i]]$deformation_rasters = deformation_rasters
        source_zone_data[[i]]$name_info = name_info
    }

    # Save the image so we can easily investigate again later
    save.image(paste0('scenarios_like_historic_events_', slip_type, '.Rdata'))

    # Make plots
    pdf(paste0('scenario_initial_conditions_', slip_type, '.pdf'), width=10, height=10)
    for(i in 1:length(scenarios_like_events)){

        event_title = basename(scenarios_like_events[[i]]$stats_file)

        for(j in 1:length(source_zone_data[[i]]$deformation_rasters)){
            par(oma=c(0,0, 2, 0))
            plot(source_zone_data[[i]]$deformation_rasters[[j]])
            title(event_title, outer=TRUE, line=-3)
            title(paste0(source_zone_data[[i]]$desired_event_rows[j], ";  ", 
                         source_zone_data[[i]]$name_info), line=-4)
        }
    }
    dev.off()

    # Rasters to a file
    out_dir = paste0('ptha18_scenario_initial_conditions_', slip_type)
    dir.create(out_dir, showWarnings=FALSE)
    for(i in 1:length(scenarios_like_events)){
        for(j in 1:length(source_zone_data[[i]]$desired_event_rows)){
            out_name = paste0(out_dir, '/', scenarios_like_events[[i]]$source_zone, '_', 
                              scenarios_like_events[[i]]$desired_event_rows[j], '_', 
                              source_zone_data[[i]]$name_info, '_',
                              basename(scenarios_like_events[[i]]$stats_file), '.tif')
            writeRaster(source_zone_data[[i]]$deformation_rasters[[j]], out_name)
        }
    }

}
