
library(rptha)
config_env = new.env()
source('config.R', local=config_env)


all_uniform_eq_events_file = Sys.glob('all_uniform_slip_earthquake_events_tsunami_*.nc')
all_stochastic_eq_events_file = Sys.glob('all_stochastic_slip_earthquake_events_tsunami_*.nc')
unit_source_statistics_file = Sys.glob('unit_source_statistics_*.nc')

unit_source_statistics = read_table_from_netcdf(unit_source_statistics_file)

assert<-function(istrue){
    if(istrue){
        print('PASS')
    }else{
        print('FAIL')
        stop()
    }
}

# Check that Mw-min, Mw-max / dMw range is represented
run_checks<-function(fid_local){

    is_uniform_slip = ('event_area' %in% names(fid_local$var))

    # Check mw values are as desired
    mws = round(ncvar_get(fid_local, 'event_Mw'), 3)
    assert(min(mws) == config_env$Mw_min)
    assert(max(mws) == config_env$Mw_max)
    assert(all(abs(diff(sort(unique(mws))) - config_env$dMw) <= 1.0e-13))

    # Check mw is consistent with area and slip 
    moment = M0_2_Mw(mws, inverse=TRUE)
    if(is_uniform_slip){

        event_area = ncvar_get(fid_local, 'event_area')    
        event_slip = ncvar_get(fid_local, 'event_slip')    
        moment_A = 3e+10 * event_area * 1e+06 * event_slip
        assert(all(abs(moment - moment_A) < 1.0e-06*moment))

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

        moment_A = unlist(lapply(1:length(event_index_string), f<-function(x){
            sum(event_slips[[x]] * 
                unit_source_statistics$width[unit_sources_in[[x]]] * 
                unit_source_statistics$length[unit_sources_in[[x]]] * 
                1e+06 * 3e+10)
        }))
        # Noting we store slip to a few significant figures, need some tolerance here
        assert(all(abs(moment - moment_A) < 1.0e-03*moment))
    }

    # Check that sum of unit source areas is the same as reported area 
    if(is_uniform_slip){
        event_index_string = ncvar_get(fid_local, 'event_index_string')
        unit_sources_in = lapply(as.list(event_index_string), f<-function(x){
            get_unit_source_indices_in_event(data.frame(event_index_string = x))
            }
        )
        event_areas_B = unlist(lapply(unit_sources_in, 
            f<-function(x) sum(unit_source_statistics$width[x]*unit_source_statistics$length[x])
        ))
        assert(all(abs(event_area - event_areas_B) < 1.0e-06*event_area))

    }
}

fid_uniform = nc_open(all_uniform_eq_events_file, readunlim=FALSE)
run_checks(fid_uniform)
nc_close(fid_uniform)

fid_stochastic = nc_open(all_stochastic_eq_events_file, readunlim=FALSE)
run_checks(fid_stochastic)
nc_close(fid_stochastic)

