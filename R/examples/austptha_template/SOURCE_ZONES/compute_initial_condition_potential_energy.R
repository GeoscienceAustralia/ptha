#
# Script to compute the potential energy for each initial condition
#
# This calculation was implemented in late 2020, well after publication
# of the initial PTHA18.
#
# Here we zero out dry areas by interpolating from a "wet-or-dry" DEM, which is
# 1 where the PTHA18-input-DEM is below MSL, 0 above. The interpolation is
# bilinear. Finally we assume areas with a value less than 0.5 are dry. 
#
# Note this will not in general give an identical energy to SWALS, because the
# latter is using a different grid-size and a different interpolation to
# determine wet or dry areas. In general, re-gridding and choice of interpolation
# can have a noticable [yet small] effect on the energy.
#

library(rptha)
library(parallel)

MC_CORES = 48 # For Gadi


ptha18 = new.env()
ptha18_access_script = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
source(ptha18_access_script, local=ptha18, chdir=TRUE)

# Get the names of the source-zone directories by searching for
# unit-source-statistics netcdf files.
source_zone_unit_source_stats = Sys.glob('/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/*/TSUNAMI_EVENTS/unit_source*.nc')
source_zones = basename(dirname(dirname(source_zone_unit_source_stats)))

# For each source-zone
# - Read all the unit-sources, compute the area for every cell, set dry regions to NA, and convert to matrices
#     For each scenario_type (FAUS, VAUS, HS)
#     - Read the scenarios
#     - Compute the deformation as a matrix
#     - Call the potential energy calculation
#     - Save the energy in a vector with length = number of scenarios
#     - Write a netcdf file that contains the initial potential energy for each scenario, for this slip_type.
#     - We can optionally add some other things [e.g. scenario metadata ]

get_all_unit_sources_with_NA_dry_regions<-function(source_zone, unit_source_stats_file){

    # Get the tif names from the unit source statistics file.
    # This ensures the tifs are ordered in the same way as the statistics
    unit_source_statistics = read_table_from_netcdf(unit_source_stats_file)
    all_unit_sources = unit_source_statistics$initial_condition_file
   
    # Get a raster derived from the Australian PTHA DEM [which the flow solver
    # used to interpolate the bathymetry]. This has 1 in wet regions, and zero
    # in dry regions
    wet_or_dry_rast = ptha18$get_wet_or_dry_DEM()

    # Resample the wet_or_dry_raster to match up with the elevation.
    # Note on a single source-zone the unit-source rasters will always have the same
    # extent and resolution, so we only need to do this once
    wet_or_dry_on_uss = resample(wet_or_dry_rast, raster(all_unit_sources[1]))

    # Cleanup
    rm(wet_or_dry_rast)
    rm(unit_source_statistics)
    gc()

    all_unit_source_matrices = mclapply(all_unit_sources, f<-function(unit_source_file){
        r = raster(unit_source_file)
        # Here we determine whether cells are 'wet' or 'dry'. Note the treatment here
        # will affect the result. Initially I tried with "==0" to define dry cells. That
        # generally leads to greater potential energy than this approach. In general the
        # current approach seemed closer to the SWALS examples I looked at, but they were
        # different again [as expected, given calculations happen on a finer grid].
        r[wet_or_dry_on_uss < 0.5] = NA
        output = as.matrix(r)
        rm(r)
        gc()
        return(output)
    }, mc.cores=MC_CORES, mc.preschedule=TRUE)

    # Get the cell areas in m^2
    # Note these are computed based on the ellipsoid, and are similar but not
    # identical to the values computed assuming the earth is a sphere.
    cell_areas_on_uss = as.matrix(area(wet_or_dry_on_uss) * 1000 * 1000)
    ## Here we make a polygon and illustrate the differences in area using sphere/ellipsoid
    ## The idea is the polygon could represent one pixel in a raster. 
    # m1 = matrix(c(0      , -1.98, 
    #               0.03333, -1.98, 
    #               0.03333, -1.98-0.03333, 
    #               0      , -1.98 - 0.03333), 
    #             ncol=2, byrow=TRUE)
    ## Differences in the area arise due to ellipticity -- f=0 is a sphere, by
    ## default the wgs84 ellipsoid is used.
    #> areaPolygon(m1)
    #[1] 13665932
    #> sqrt(areaPolygon(m1))
    #[1] 3696.746
    #> areaPolygon(m1, f=0)
    #[1] 13757810
    #> sqrt(areaPolygon(m1, f=0))
    #[1] 3709.152
    gc()

    return(list(all_unit_sources = all_unit_sources, 
                all_unit_source_matrices = all_unit_source_matrices,
                cell_areas_m2 = cell_areas_on_uss))

}

# For a given source zone, slip type, and unit-source data, compute the initial potential energy
# for each scenario.
get_initial_potential_energy_for_scenarios<-function(source_zone, unit_sources_data, slip_type){

    allowed_slip_types = c('uniform', 'stochastic', 'variable_uniform')

    if(!(slip_type %in% allowed_slip_types)){
        print('The given slip_type:')
        print(slip_type)
        print('is not within the allowed values:')
        print(allowed_slip_types)
        stop('Unknown value of slip_type')
    }

    # Get the event data
    szed = ptha18$get_source_zone_events_data(source_zone, slip_type=slip_type)
    # Remove much of the data -- try to save memory in parallel
    szed2 = list()
    if('event_slip_string' %in% names(szed$events)){
        szed2$events = data.frame(szed$events[c('event_slip_string', 'event_index_string')])
    }else{
        szed2$events = data.frame(szed$events[c('slip', 'event_index_string')])
    }
    szed = szed2
    rm(szed2)
    gc()

    # This computes the potential energy based on the row index.
    compute_energy<-function(row_index, unit_sources_data, szed){

        event_inds = as.numeric(strsplit(szed$events$event_index_string[row_index], '-')[[1]])

        if('event_slip_string' %in% names(szed$events)){
            # Different slip for each unit-source
            event_slip = as.numeric(strsplit(szed$events$event_slip_string[row_index], '_')[[1]])
            stopifnot(length(event_inds) == length(event_slip))
        }else{
            # Same slip for each unit-source
            event_slip = szed$events$slip[row_index] + 0*event_inds
        }

        # Compute the deformation
        deformation = 0*unit_sources_data$all_unit_source_matrices[[1]]
        for(i in 1:length(event_inds)){
            deformation = deformation + event_slip[i] * 
                unit_sources_data$all_unit_source_matrices[[event_inds[i]]]
        }

        potential_energy = sea_surface_available_potential_energy(deformation, 
            unit_sources_data$cell_areas_m2,
            gravity = 9.8,
            seawater_density = 1024,
            MSL = 0)

        # Try to avoid accumulation of memory
        #to_remove = setdiff(ls(all=TRUE), c('potential_energy', 'row_index', 'unit_sources_data', 'szed'))
        to_remove = setdiff(ls(all=TRUE), c('potential_energy'))
        rm(to_remove)
        gc()
        return(potential_energy)
    }

    # Fail neatly
    try_compute_energy<-function(row_index, unit_sources_data=unit_sources_data, szed=szed){
        output = try(compute_energy(row_index, unit_sources_data, szed))
        return(output)
    }

    all_energies = mclapply(1:nrow(szed$events), try_compute_energy, 
                            unit_sources_data=unit_sources_data, szed=szed,
                            mc.cores=MC_CORES, mc.preschedule=TRUE)

    return(all_energies)

}


# Compute the potential energy for the initial condition of every scenario, on
# every source zone
for(si in 1:length(source_zones)){

    # Process the unit source tifs
    unit_source_data = get_all_unit_sources_with_NA_dry_regions(source_zones[si], source_zone_unit_source_stats[si])

    # Compute the potential energy for each slip type
    ipe_U = get_initial_potential_energy_for_scenarios(source_zones[si], unit_source_data, 'uniform')
    ipe_S = get_initial_potential_energy_for_scenarios(source_zones[si], unit_source_data, 'stochastic')
    ipe_V = get_initial_potential_energy_for_scenarios(source_zones[si], unit_source_data, 'variable_uniform')

    # Save the result to a temporary file -- can convert to netcdf for distribution later.
    output_file = paste0('./', source_zones[si], '/TSUNAMI_EVENTS/potential_energy_', source_zones[si], '.RDS')
    output = list(ipe_U = ipe_U, ipe_S = ipe_S, ipe_V=ipe_V, source_zone=source_zones[si])
    saveRDS(output, output_file)

    # Clear the memory
    rm(unit_source_data, ipe_S, ipe_U, ipe_V, output)
    gc()

}

#
# Write the potential energy to a netcdf file.
# We need a separate one for each source zone and slip type.
#
all_RDS_potential_energy_files = Sys.glob('./*/TSUNAMI_EVENTS/potential_energy*.RDS')
for(i_file in 1:length(all_RDS_potential_energy_files)){

    pe = readRDS(all_RDS_potential_energy_files[i_file])
   
    # For the 'defunct' source-zones, the RDS file is nearly empty and will contain try-errors
    # Skip those ones 
    has_try_errors = (
        any(unlist(lapply(pe$ipe_U, f<-function(x) is(x, 'try-error')))) |
        any(unlist(lapply(pe$ipe_S, f<-function(x) is(x, 'try-error')))) |
        any(unlist(lapply(pe$ipe_V, f<-function(x) is(x, 'try-error')))) )

    if(has_try_errors){
        print(paste0('skipping ', pe$source_zone))
        next
    }

    source_zone = pe$source_zone

    # Write one file for each slip type. 
    for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){

        # Put useful attributes in the netcdf file
        global_attributes = list()
        global_attributes$source_zone = source_zone
        global_attributes$slip_type = slip_type
        global_attributes$parent_script_name = parent_script_name()

        # Get the appropriate initial potential energy
        if(slip_type == 'uniform'){
            ipe = unlist(pe$ipe_U)
        }else if(slip_type == 'stochastic'){
            ipe = unlist(pe$ipe_S)
        }else if(slip_type == 'variable_uniform'){
            ipe = unlist(pe$ipe_V)
        }else{
            stop(paste0('unrecognized slip_type ', slip_type, 
                        ' on source_zone ', source_zone))
        }

        # Read the event metadata
        earthquake_events_file = paste0('./', source_zone, '/TSUNAMI_EVENTS/all_', 
            slip_type, '_slip_earthquake_events_', source_zone, '.nc')
        if(!file.exists(earthquake_events_file)){
            stop(paste0('could not find earthquake_events_file ', earthquake_events_file))
        }
        global_attributes$earthquake_events_file = normalizePath(earthquake_events_file)

        global_attributes$Note_on_potential_energy = 'The initial potential energy (units joules) has been computed by spatially integrating the potential energy of the tsunami initial condition. Land-areas were zeroed-out by interpolating from wet-dry mask, which was itself derived from the DEM used in PTHA18. Note differences with other potential energy calculations can occur -- typically the effect is small, but differences of a few percent are not unusual. Larger differences might be expected if most of the deformation occurs near wet-dry boundaries. The following factors can cause that: A) The cell-area computation will vary depending on whether a sphere or an ellipsoid is assumed, and; B) our hydrodynamic models often interpolate the initial condition and the elevation to a different resolution prior to computation, and; C) the designation of wet-or-dry areas used herein, derived by interpolating from a logical raster, is not identical to what would be obtained by interpolating DEM values, even if the resolution is otherwise the same.'
       global_attributes$Note_on_earthquake_scenarios = 'The potential energies herein correspond to model scenarios in the earthquake_events_file, so the number of rows in that file should equal the number of rows here. Herein we also provide the rate_annual, the weight_with_nonzero_rate, and Mw, mainly to cross-reference with the same variables in the earthquake_events_file (and thus check that we have correctly associated potential energies with their scenarios). Note that not all of the earthquake scenarios are considered possible according to PTHA18 (i.e. some have a rate of zero).'

        # Get the annual rate from the earthquake events file
        # This is mainly included to enable cross-checking with the earthquake
        # events file -- since the annual rate will vary often between scenarios.
        events = read_table_from_netcdf(earthquake_events_file)
        rate_annual = events$rate_annual
        weight_with_nonzero_rate = events$weight_with_nonzero_rate
        Mw = events$Mw
        rm(events)
        gc()

        # Make a new file for potential energy, without breaking the code
        potential_energy_output_nc_file = paste0(dirname(normalizePath(earthquake_events_file)),
             '/POTENTIAL_ENERGY_all_', slip_type, '_slip_earthquake_events_', source_zone, '_POTENTIAL_ENERGY.nc')

        output_data = data.frame(initial_potential_energy = ipe, 
                                 rate_annual=rate_annual, 
                                 weight_with_nonzero_rate=weight_with_nonzero_rate,
                                 Mw=Mw)

        write_table_to_netcdf(
            output_data, 
            filename=potential_energy_output_nc_file, 
            global_attributes_list=global_attributes,
            units=c('joules', 'events_per_year', '', 'moment_magnitude'),
            long_names=c('potential_energy_integrated_from_tsunami_initial_condition', 
                         'logic_tree_mean_annual_rate_of_events', 
                         'weighted_fraction_of_logic_tree_branches_having_nonzero_rate',
                         'moment_magnitude_with_constant_rigidity'),
            var_prec=c('double', 'double', 'double', 'double'),
            add_session_info_attribute=TRUE)

    }
    
}

