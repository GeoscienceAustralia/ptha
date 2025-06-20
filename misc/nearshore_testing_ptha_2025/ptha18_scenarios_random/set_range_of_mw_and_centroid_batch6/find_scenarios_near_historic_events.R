#
# Record information on the historic events for which we seek "similar location
# and magnitude". 
# In this case we specify a magnitude range, as well as 2 points, which constrain 
# the scenario centroid location. We allow any down-dip location, but
# use the 2 points to constrain the along-strike location.
#
source('compare_with_data_environment.R')

target_events = list(

    'java2006-batch6' = list(
        Mw_lower = 7.7 - 0.15,
        Mw_upper = 7.7 + 0.15,
        event_point_1 = c(109.67, -10.),  # Inferred from Fujii and Satake 2006, cover region slightly wider than their inverted sources
        event_point_2 = c(106.9, -9.26),  # "                                                                                          "
        event_time = strptime('2006-07-17 08:19:27', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'),
        tsunami_event_dir = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/'
     ),

    'sumatra2007-batch6' = list(
        Mw_lower = 8.5 - 0.15,
        Mw_upper = 8.5 + 0.15,
        event_point_1 = c(101.25, -5.4),  # Inferred from Fujii and Satake 2008, cover region slightly wider than their inverted sources
        event_point_2 = c(99.12, -2.55),  # "                                                                                          "
        event_time = strptime('2007-09-12 11:10:27', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'),
        tsunami_event_dir = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/'
     )


    )

#
# Get the corresponding events for each source-zone
#
target_scenarios = vector(mode ='list', length=length(target_events))
names(target_scenarios) = names(target_events)
for(i in 1:length(target_events)){
    
    source_zone_env = make_scenario_selection_environment(target_events[[i]]$tsunami_event_dir)

    target_scenarios[[i]] = list('uniform_slip' = NA, 
                                 'variable_area_uniform_slip' = NA, 
                                 'heterogeneous_slip' = NA)
    target_scenarios[[i]]$uniform_slip = source_zone_env$find_events_with_magnitude_and_alongstrike_centroid_constraints(
        target_events[[i]]$Mw_lower,
        target_events[[i]]$Mw_upper,
        target_events[[i]]$event_point_1,
        target_events[[i]]$event_point_2)
    target_scenarios[[i]]$variable_area_uniform_slip = source_zone_env$find_events_with_magnitude_and_alongstrike_centroid_constraints(
        target_events[[i]]$Mw_lower,
        target_events[[i]]$Mw_upper,
        target_events[[i]]$event_point_1,
        target_events[[i]]$event_point_2,
        use_variable_uniform_slip_runs=TRUE)
    target_scenarios[[i]]$heterogeneous_slip = source_zone_env$find_events_with_magnitude_and_alongstrike_centroid_constraints(
        target_events[[i]]$Mw_lower,
        target_events[[i]]$Mw_upper,
        target_events[[i]]$event_point_1,
        target_events[[i]]$event_point_2,
        use_stochastic_slip_runs=TRUE)

    rm(source_zone_env)

}

save.image('find_scenarios_near_historical_events_R_image.RData')


#
# Sample 60 scenarios at random from each
#
NSCENARIOS = 60
sample_scenarios<-function(event_table){
    nr = nrow(event_table)
    if(nr == 1) stop('Only one scenario in event table -- beware the behaviour of "sample"!')
    k = sample(1:nr, size=NSCENARIOS, replace=TRUE, prob = event_table$rate_annual)
    return(k)
}

target_scenarios_random_indices = vector(mode='list', length=length(target_scenarios))
names(target_scenarios_random_indices) = names(target_scenarios)
for(i in 1:length(target_scenarios)){

    target_scenarios_random_indices[[i]] = list('uniform_slip' = list(), 
                                                'variable_area_uniform_slip' = list(), 
                                                'heterogeneous_slip' = list())
    for(st in names(target_scenarios_random_indices[[i]])){

        # Record the names of all scenario rows
        target_scenarios_random_indices[[i]][[st]]$scenario_rows = 
            as.numeric(rownames(target_scenarios[[i]][[st]]))

        # Make the random seed reproducible -- that way if I need to increase the number
        # of sceenarios, it is easy to do
        target_scenarios_random_indices[[i]][[st]]$random_seed = get_random_seed()

        # Get a random sample with 20 scenarios per event set
        target_scenarios_random_indices[[i]][[st]]$random_inds = sample_scenarios(
            target_scenarios[[i]][[st]])
    }
}

# Some random scenarios are be repeated 
save.image('find_scenarios_near_historical_events_R_image_B.RData')

## Check double-ups
lapply(target_scenarios_random_indices, 
      f<-function(x) lapply(x, f<-function(y) length(unique(y$random_inds))))

lapply(target_scenarios_random_indices, 
      f<-function(x) lapply(x, f<-function(y) (table(y$random_inds))))

# Get the ptha18-reading codes
ptha18 = new.env()
source('/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R',
       local=ptha18, chdir=TRUE)

#
# Get the event data from the PTHA18 database
#
random_scenarios = vector(mode='list', length=length(target_events))
names(random_scenarios) = names(target_events)
for(i in 1:length(target_scenarios_random_indices)){

    sourcename = basename(dirname(target_events[[i]]$tsunami_event_dir))
    print(sourcename)

    random_scenarios_data = list()

    # We have slightly different slip-type naming conventions in this script, vs the
    # ptha18 databases. 
    slip_types = c('uniform_slip', 'variable_area_uniform_slip', 'heterogeneous_slip')
    slip_types_ptha18 = c('uniform', 'variable_uniform', 'stochastic')

    # Loop over slip types and get the randomly chosen scenarios on this source
    for(j in 1:length(slip_types)){
        st = slip_types[j]
        st2 = slip_types_ptha18[j]

        scenario_rows = target_scenarios_random_indices[[i]][[st]]$scenario_rows
        random_inds = unique(target_scenarios_random_indices[[i]][[st]]$random_inds)
        target_rows = scenario_rows[random_inds]
        nonunique_target_rows = scenario_rows[target_scenarios_random_indices[[i]][[st]]$random_inds]

        # This can fail due to server issues -- keep trying
        random_scenarios_data[[st]] = try(log('a'), silent=TRUE); error_count = 0
        while(is(random_scenarios_data[[st]], 'try-error')){
            error_count = error_count + 1
            if(error_count == 11) stop('Too many download errors, aborting')
            random_scenarios_data[[st]] = try(ptha18$get_source_zone_events_data(sourcename,
                slip_type=st2, 
                desired_event_rows = target_rows))
        }
        # Some random scenarios will double up -- we store that info
        # Beware these indices refer to "indices in the subset of events with
        # specified magnitude + centroid range"
        random_scenarios_data[[st]]$original_random_scenario_counts = 
            table(target_scenarios_random_indices[[i]][[st]]$random_inds)
        random_scenarios_data[[st]]$original_random_scenario_inds = 
            (target_scenarios_random_indices[[i]][[st]]$random_inds)
        # 
        # For the desired event rows, count how many times each event is represented in the original
        # sample, so we can simply weight their contribution at the end. This can be derived from
        # the previous info, but easier to store it here.
        random_scenarios_data[[st]]$desired_event_rows_count = 
            sapply(random_scenarios_data[[st]]$desired_event_rows, 
                   f<-function(x) sum(x == nonunique_target_rows))

    }

    random_scenarios[[i]] = random_scenarios_data
}

# 
save.image('find_scenarios_near_historical_events_R_image_C.RData')

#
# Make the initial condition rasters for each
#
for(i in 1:length(random_scenarios)){

    random_scenarios_i = random_scenarios[[i]]
    output_dir = names(random_scenarios)[i]
    dir.create(output_dir, showWarnings=FALSE)

    # Loop over slip types and create the Okada-deformation raster
    slip_type = c('uniform_slip', 'variable_area_uniform_slip', 'heterogeneous_slip')
    for(st in slip_type){
        # For each event in the table of random events, get the raster
        nr = nrow(random_scenarios_i[[st]]$events)
        for(j in 1:nr){
            rast = ptha18$get_initial_condition_for_event(random_scenarios_i[[st]], event_ID=j)
            rast_file = paste0(output_dir, '/', output_dir, '_', st, '_', 
                               random_scenarios_i[[st]]$desired_event_rows[j], '_count_',
                               random_scenarios_i[[st]]$desired_event_rows_count[j],  '.tif')
            writeRaster(rast, file=rast_file, options=c('COMPRESS=DEFLATE'))
        }
    }
}

#
# 1. Eyeball -- maybe they are not sufficiently "in the right location" ?
# 2. Compute energy -- guess what kind of waves we will get.
# 3. 

global_elev = raster('../../elevation/derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif')
zc = readOGR('../../sources/figures/zero_contour/', layer='zero_contour')

raster_energies = vector(mode='list', length=length(random_scenarios))
names(raster_energies) = names(random_scenarios)
for(i in 1:length(random_scenarios)){
    print(i)
    # Find the rasters
    event_dir = names(random_scenarios)[i]
    event_rasts = Sys.glob(paste0(event_dir, '/*.tif'))


    # Make an elevation raster corresponding to these rasters
    r1 = raster(event_rasts[1])
    r1_p = rasterToPoints(r1)
    r1_elev = extract(global_elev, r1_p[,1:2])
    r1_p[,3] = r1_elev
    r1_elev = rasterFromXYZ(r1_p, crs=CRS(proj4string(r1)))

    r1_ext = extent(r1_elev)
    ref_lat = 0.5*(extent(r1_ext)@ymin + extent(r1_ext)@ymax)

    # Compute 'rho * g / 2 * dA' for each cell in m^2, and set this to zero on dry cells
    rho = 1024
    grav = 9.8
    rho_g_dA_on_2 = (r1_elev < 0) * rho * grav / 2 * area(r1_elev) * 1e+06

    # Make a plot of all the initial conditions
    pdf(paste0(event_dir, '/scenario_plots.pdf'), width=6, height = 6)

    pe = rep(NA, length(event_rasts))
    for(j in 1:length(event_rasts)){
        up = raster(event_rasts[j])
        # Potential energy calculation = integral { rho * g /2 * stage^2 }
        # Cross-checked against the SWALS potential energy calculation for some cases
        pe[j] = sum(as.matrix(rho_g_dA_on_2) * as.matrix(up)**2)

        plot(up, asp=1/cos(ref_lat/180*pi))
        title(paste0('Potential energy: ', signif(pe[j], 3), '\n ', basename(event_rasts[j])))
        points(target_events[[i]]$event_hypocentre[1], 
               target_events[[i]]$event_hypocentre[2], col='black')
        plot(zc, col='black', add=TRUE)

        points(rbind(target_events[[i]]$event_point_1, target_events[[i]]$event_point_2),
               col='purple', pch=1)
    }

    dev.off()

    raster_energies[[i]] = data.frame(rasts=event_rasts, potential_energy=pe)
}

save.image('find_scenarios_near_historical_events_R_image_D.RData')

#


