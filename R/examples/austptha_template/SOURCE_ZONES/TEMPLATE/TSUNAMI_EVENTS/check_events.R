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


all_uniform_eq_events_file = Sys.glob('all_uniform_slip_earthquake_events_tsunami_*.nc')
all_stochastic_eq_events_file = Sys.glob('all_stochastic_slip_earthquake_events_tsunami_*.nc')
unit_source_statistics_file = Sys.glob('unit_source_statistics_*.nc')

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

# Check that Mw-min, Mw-max / dMw range is represented
run_checks<-function(fid_local){

    is_uniform_slip = ('event_area' %in% names(fid_local$var))

    # Check mw values are as desired
    mws = round(ncvar_get(fid_local, 'event_Mw'), 3)
    assert(min(mws) == config_env$Mw_min, 'mw_min problem')
    assert(max(mws) == config_env$Mw_max, 'mw_max problem')
    assert(all(abs(diff(sort(unique(mws))) - config_env$dMw) <= 1.0e-13), 'mw spacing problem')

    # Check mw is consistent with area and slip 
    moment = M0_2_Mw(mws, inverse=TRUE)
    if(is_uniform_slip){

        event_area = ncvar_get(fid_local, 'event_area')    
        event_slip = ncvar_get(fid_local, 'event_slip')    
        moment_A = 3e+10 * event_area * 1e+06 * event_slip
        assert(all(abs(moment - moment_A) < 1.0e-06*moment), 'moment inconsistency')

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
        assert(all(l1 == l2), 'inconsistent number of slips and unit sources, stochastic')

        moment_A = unlist(lapply(1:length(event_index_string), f<-function(x){
            sum(event_slips[[x]] * 
                unit_source_statistics$width[unit_sources_in[[x]]] * 
                unit_source_statistics$length[unit_sources_in[[x]]] * 
                1e+06 * 3e+10)
        }))
        # Noting we store slip to a few significant figures, need some tolerance here
        assert(all(abs(moment - moment_A) < 1.0e-03*moment), 'moment inconsistency stochastic')
    }

    # Check that sum of unit source areas is the same as reported area 
    if(is_uniform_slip){

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
        assert(all(abs(event_area - event_areas_B) < 1.0e-06*event_area), 'area inconsistency uniform')

    }

    elev = ncvar_get(fid_local, 'elev')
    lon = ncvar_get(fid_local, 'lon')
    lat = ncvar_get(fid_local, 'lat')

    # Check max stage, in chunks for memory efficiency
    ngauges = length(lat)
    chunksize = 1000
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
        initial_stage_negative_but_not_missing = ((initial_stage < 0) & (initial_stage > (config_env$null_double + 1)))
        assert(all(stage_na_or_positive | initial_stage_negative_but_not_missing), 'negative peak stage error')

        # It should be rare to have sites with max_stage negative.
        max_stage_negative_fraction = mean(max_stage < 0, na.rm=TRUE)
        assert(max_stage_negative_fraction < 1/100, 'too many negative peak stages')
    
       
        # Find gauges that have NA -- these should be gauges with elevation >
        # 0, or gauges in the boundary condition exclusion zone
        na_gauges = gauge_chunks[[i]][which(is.na(max_stage[1,]))]
        if(length(na_gauges) > 0){
            assert(all(
                elev[na_gauges] >= 0 | 
                lat[na_gauges] >= config_env$lat_range[2] | 
                lat[na_gauges] <= config_env$lat_range[1]), 
                'na gauges with elev < 0, not in boundary regions'
                )
        }
    }

}

# Check uniform slip
fid_uniform = nc_open(all_uniform_eq_events_file, readunlim=FALSE)
run_checks(fid_uniform)
nc_close(fid_uniform)

# Check stochastic slip
fid_stochastic = nc_open(all_stochastic_eq_events_file, readunlim=FALSE)
run_checks(fid_stochastic)
nc_close(fid_stochastic)

