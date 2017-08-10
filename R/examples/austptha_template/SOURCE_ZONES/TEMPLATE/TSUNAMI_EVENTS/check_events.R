#
# Code which performs some basic logical checks on the main tsunami event
# output netcdf files.
#
# Checks include:
#     - mw is distributed as expected
#     - mw is consistent with integral of (slip x area)
#     - uniform slip area is equal to sum of area over unit sources
#     - stochastic slip 'number of unit sources' is equal to stochastic slip 'number of slip values'
#     - max_stage is non-negative at gauges with elevation < 0, except for a tolerably small fraction of 
#       'trapped' sites with initial deformation < 0 [e.g. lakes in the Okada subsidence zone]
#

library(rptha)
config_env = new.env()
source('config.R', local=config_env)

# Netcdf filenames with tsunami and all events
all_uniform_eq_events_file = Sys.glob(
    'all_uniform_slip_earthquake_events_tsunami_*.nc')
all_stochastic_eq_events_file = Sys.glob(
    'all_stochastic_slip_earthquake_events_tsunami_*.nc')
all_variable_uniform_eq_events_file = Sys.glob(
    'all_variable_uniform_slip_earthquake_events_tsunami_*.nc')
unit_source_statistics_file = Sys.glob('unit_source_statistics_*.nc')

# 
unit_source_statistics = read_table_from_netcdf(unit_source_statistics_file)

# Utility for testing
assert<-function(istrue, msg=''){

    if(istrue){
        print('PASS')
    }else{
        print('FAIL')
        stop(msg)
    }

}

#' Check that Mw-min, Mw-max / dMw range is represented
#'
#' @param fid_local netcdf file handle [i.e. result of nc_open('filename.nc')]
#' @param slip_type character. 'uniform' or 'stochastic' or 'variable_uniform'
#' @return the function environment. Should also print lots of 'PASS'
#'   statements, and will make a png plot of earthquake scaling for
#'   stochastic/variable_uniform slip. 
#'
run_checks<-function(fid_local, slip_type){

    # Check mw values are as desired
    mws = round(ncvar_get(fid_local, 'event_Mw'), 3)
    assert(min(mws) == config_env$Mw_min, 'mw_min problem')
    assert(max(mws) == config_env$Mw_max, 'mw_max problem')
    assert(all(abs(diff(sort(unique(mws))) - config_env$dMw) <= 1.0e-13), 
        'mw spacing problem')

    # Check mw is consistent with area and slip 
    moment = M0_2_Mw(mws, inverse=TRUE)
    if(slip_type == 'uniform'){

        event_area = ncvar_get(fid_local, 'event_area')    
        event_slip = ncvar_get(fid_local, 'event_slip')    
        moment_A = 3e+10 * event_area * 1e+06 * event_slip
        assert(all(abs(moment - moment_A) < 1.0e-06*moment), 
            'moment inconsistency')

    }else{
   
        # Indices in event 
        event_index_string = ncvar_get(fid_local, 'event_index_string')
        unit_sources_in = lapply(as.list(event_index_string), f<-function(x){
            get_unit_source_indices_in_event(data.frame(event_index_string = x))
            }
        )
    
        # Slip in event
        event_slip_string = ncvar_get(fid_local, 'event_slip_string')
        event_slips = lapply(as.list(event_slip_string), f<-function(x){
            as.numeric(strsplit(x, '_')[[1]])
            }
        )

        # Check number of unit sources and number of slips are equal, for each event
        l1 = unlist(lapply(unit_sources_in, length))
        l2 = unlist(lapply(event_slips, length))
        assert(all(l1 == l2), 
            'inconsistent number of slips and unit sources, stochastic')

        moment_A = unlist(lapply(1:length(event_index_string), f<-function(x){
            sum(event_slips[[x]] * 
                unit_source_statistics$width[unit_sources_in[[x]]] * 
                unit_source_statistics$length[unit_sources_in[[x]]] * 
                1e+06 * 3e+10)
        }))
        # Noting we store slip to a few significant figures, need some
        # tolerance here
        assert(all(abs(moment - moment_A) < 1.0e-03*moment), 
            'moment inconsistency stochastic')
    }

    # Check that sum of unit source areas is the same as reported area 
    if(slip_type == 'uniform'){

        event_index_string = ncvar_get(fid_local, 'event_index_string')
        unit_sources_in = lapply(as.list(event_index_string), f<-function(x){
            get_unit_source_indices_in_event(data.frame(event_index_string = x))
            }
        )
        event_areas_B = unlist(lapply(unit_sources_in, 
            f<-function(x){
                sum(unit_source_statistics$width[x]*
                    unit_source_statistics$length[x]) }
        ))
        assert(all(abs(event_area - event_areas_B) < 1.0e-06*event_area), 
            'area inconsistency uniform')

    }

    elev = ncvar_get(fid_local, 'elev')
    lon = ncvar_get(fid_local, 'lon')
    lat = ncvar_get(fid_local, 'lat')

    # Check max stage, in chunks for memory efficiency
    ngauges = length(lat)
    chunksize = ceiling(1000/(ceiling(length(moment)/10000)))
    gauge_chunks = parallel::splitIndices(ngauges, ceiling(ngauges/chunksize))

    # Should either be NA, or a positive number. Negative numbers were used for
    # 'incomplete' files.
    # Note that if you have a lake overlapping the 'depression' part of the
    # source-zone, then max_stage could be negative
    for(i in 1:length(gauge_chunks)){
        start = gauge_chunks[[i]][1]
        count = length(gauge_chunks[[i]])
        
        # Check that EITHER: max_stage>=0, OR is.na(max_stage), OR initial
        # stage is negative [the latter allows for trapped lakes in the Okada
        # subsidence zone to have negative max_stage].
        max_stage = ncvar_get(fid_local, 'max_stage', start=c(1, start), 
            count=c(-1, count))
        initial_stage = ncvar_get(fid_local, 'initial_stage', start=c(1, start), 
            count=c(-1, count))
        stage_na_or_positive = ( is.na(max_stage) | (max_stage >= 0))
        initial_stage_negative_but_not_missing = ((initial_stage < 0) & 
            (initial_stage > (config_env$null_double + 1)))
        assert(all(stage_na_or_positive | initial_stage_negative_but_not_missing), 
            paste0('negative peak stage error ', start, ' ', count))

        # It should be rare to have sites with max_stage negative.
        max_stage_negative_fraction = mean(max_stage < 0, na.rm=TRUE)
        assert(max_stage_negative_fraction < 1/100, 
            paste0('too many negative peak stages', start, ' ', count))
    
       
        # Find gauges that have NA -- these should be gauges with elevation >
        # 0, or gauges in the boundary condition exclusion zone
        na_gauges = gauge_chunks[[i]][which(is.na(max_stage[1,]))]
        if(length(na_gauges) > 0){
            assert(all(
                elev[na_gauges] >= 0 | 
                lat[na_gauges] >= config_env$lat_range[2] | 
                lat[na_gauges] <= config_env$lat_range[1]), 
                paste0('na gauges with elev < 0, not in boundary regions ', 
                    start, ' ', count)
                )
        }
    }

    if(slip_type != 'uniform'){
        #
        # Make a length/width scaling law plot
        #

        alongstrike_range = lapply(unit_sources_in, 
            f<-function(x) range(unit_source_statistics$alongstrike_number[x]))
        downdip_range = lapply(unit_sources_in, 
            f<-function(x) range(unit_source_statistics$downdip_number[x]))

        event_length = unlist(mapply(
            f<-function(strk_ind, dip_ind){
                # Compute the length, based on up-dip unit-sources
                r1 = strk_ind[1]:strk_ind[2]
                r2 = dip_ind[1]:dip_ind[2]
                tokeep = which(
                    (unit_source_statistics$downdip_number == r2[1]) & 
                    (unit_source_statistics$alongstrike_number %in% r1)
                )
                sum(unit_source_statistics$length[tokeep])
            }, 
            alongstrike_range,
            downdip_range 
            ))

        event_width = unlist(mapply(
            f<-function(strk_ind, dip_ind){
                # Compute the width, based on least-along-strike unit-sources
                r1 = strk_ind[1]:strk_ind[2]
                r2 = dip_ind[1]:dip_ind[2]
                tokeep = which(
                    (unit_source_statistics$downdip_number %in% r2) & 
                    (unit_source_statistics$alongstrike_number == r1[1])
                )
                sum(unit_source_statistics$width[tokeep])
            }, 
            alongstrike_range,
            downdip_range 
            ))


        # Area of 'non-zero' patches of the rectangular region of slip        
        event_area_nonzero = unlist(lapply(unit_sources_in, 
            f<-function(x) sum(unit_source_statistics$length[x] * 
                unit_source_statistics$width[x])))
        # Mean slip on non-zero patches of the rectangular region of slip
        event_mean_slip_nonzero = unlist(lapply(event_slips, 
            f<-function(x) mean(x) ))

        # Area of full rectangular region of slip, including effect of zero
        # patches
        event_area_all = unlist(mapply(
            f<-function(strk_ind, dip_ind){
                # Compute the width, based on least-along-strike unit-sources
                r1 = strk_ind[1]:strk_ind[2]
                r2 = dip_ind[1]:dip_ind[2]
                tokeep = which(
                    (unit_source_statistics$downdip_number %in% r2) & 
                    (unit_source_statistics$alongstrike_number %in% r1)
                )
                sum(unit_source_statistics$width[tokeep] *
                    unit_source_statistics$length[tokeep])
            }, 
            alongstrike_range,
            downdip_range 
            ))

        # Mean slip, including 'zero areas
        event_mean_slip_all = event_mean_slip_nonzero * event_area_nonzero / 
            event_area_all

        # Maximum slip
        event_max_slip = unlist(lapply(event_slips, 
            f<-function(x) max(x) ))

        # Compute theoretical results
        unique_mws = unique(mws)
        scaling_law_dim = lapply(as.list(unique_mws), 
            f<-function(x) Mw_2_rupture_size(x, detailed=TRUE, CI_sd=2))

        # Convenience function for the plot
        local_scaling_law_results<-function(var='area'){ 

            areas_scaling = unlist(lapply(scaling_law_dim, 
                f<-function(x) x$values[var]))
            areas_scaling_pci = unlist(lapply(scaling_law_dim, 
                f<-function(x) x$plus_CI[var]))
            areas_scaling_mci = unlist(lapply(scaling_law_dim, 
                f<-function(x) x$minus_CI[var]))

            points(unique_mws, areas_scaling, t='l', col='red')
            points(unique_mws, areas_scaling_pci, t='l', col='red', 
                lty='dashed')
            points(unique_mws, areas_scaling_mci, t='l', col='red', 
                lty='dashed')
        }

        # Convenience function to add the median value of a model variable to a plot
        add_quartiles<-function(mws, var){
            lowerq_var = aggregate(var, list(mws), f<-function(x) quantile(x, p=0.25, type=6))
            points(lowerq_var[,1], lowerq_var[,2], col='orange')
            upperq_var = aggregate(var, list(mws), f<-function(x) quantile(x, p=0.75, type=6))
            points(upperq_var[,1], upperq_var[,2], col='orange')
            median_var = aggregate(var, list(mws), f<-function(x) median(x))
            points(median_var[,1], median_var[,2], col='red')
        }

        # Plotting
        png_filename = paste0('event_size_scaling_', slip_type, '_', 
                basename(dirname(getwd())), '.png')
        png(png_filename, width=13, height=10, units='in', res=200)

        par(mfrow=c(3,3))

        # Area
        plot(mws, event_area_nonzero, log='y', 
            main='Event area with non-zero slip', 
            xlab='Mw', ylab='km^2')
        grid()
        local_scaling_law_results('area')
        add_quartiles(mws, event_area_nonzero)

        # Area
        plot(mws, event_area_all, log='y', 
            main='Event area (including zero slip cells)', 
            xlab='Mw', ylab='km^2')
        grid()
        local_scaling_law_results('area')
        add_quartiles(mws, event_area_all)

        # Length
        plot(mws, event_length, log='y', main='Event length', 
            xlab='Mw', ylab='km')
        grid()
        local_scaling_law_results('length')
        add_quartiles(mws, event_length)

        # Width 
        plot(mws, event_width, log='y', main='Event width', 
            xlab='Mw', ylab='km')
        grid()
        local_scaling_law_results('width')
        add_quartiles(mws, event_width)

        # Mean slip
        plot(mws, event_mean_slip_nonzero, log='y', 
            main='Event mean_slip on non-zero slip patches', 
            xlab='Mw', ylab='km')
        points(unique_mws, slip_from_Mw(unique_mws), t='l', col='red', 
            lty='dashed')
        points(unique_mws, slip_from_Mw(unique_mws)*3, t='l', col='red', 
            lty='dashed')
        add_quartiles(mws, event_mean_slip_nonzero)
        grid()

        # Mean slip
        plot(mws, event_mean_slip_all, log='y', 
            main='Event mean_slip including zero-slip patches', 
            xlab='Mw', ylab='km')
        points(unique_mws, slip_from_Mw(unique_mws), t='l', col='red', 
            lty='dashed')
        points(unique_mws, slip_from_Mw(unique_mws)*3, t='l', col='red', 
            lty='dashed')
        add_quartiles(mws, event_mean_slip_all)
        grid()


        # Peak slip
        plot(mws, event_max_slip, log='y', main='Event max_slip', 
            xlab='Mw', ylab='km')
        points(unique_mws, slip_from_Mw(unique_mws), t='l', col='red', 
            lty='dashed')
        points(unique_mws, slip_from_Mw(unique_mws)*3, t='l', col='red', 
            lty='dashed')
        add_quartiles(mws, event_max_slip)
        grid()

        # Peak slip / mean slip
        plot(mws, event_max_slip/event_mean_slip_nonzero, log='y', 
            main='Event max_slip/mean_non-zero_slip', 
            xlab='Mw', ylab='')
        grid()
        add_quartiles(mws, event_max_slip/event_mean_slip_nonzero)
        abline(h=3, col='red')

        plot(mws, event_max_slip/event_mean_slip_all, log='y', 
            main='Event max_slip/mean_slip_including_zeros', 
            xlab='Mw', ylab='')
        grid()
        add_quartiles(mws, event_max_slip/event_mean_slip_all)
        abline(h=3, col='red')

        dev.off()

    }

    return(environment())

}

# Check uniform slip
fid_uniform = nc_open(all_uniform_eq_events_file, readunlim=FALSE)
uniform_env = run_checks(fid_uniform, slip_type = 'uniform')
nc_close(fid_uniform)

# Check stochastic slip
fid_stochastic = nc_open(all_stochastic_eq_events_file, readunlim=FALSE)
stochastic_env = run_checks(fid_stochastic, slip_type = 'stochastic')
nc_close(fid_stochastic)

# Check variable_uniform slip
fid_variable_uniform = nc_open(all_variable_uniform_eq_events_file, readunlim=FALSE)
variable_uniform_env = run_checks(fid_variable_uniform, slip_type = 'variable_uniform')
nc_close(fid_variable_uniform)

