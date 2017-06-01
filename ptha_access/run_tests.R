source('./R/sum_tsunami_unit_sources.R', local=TRUE)
source('./R/config.R', local=TRUE)

test_file = "unit_source_tsunami/RUN_20161121104520_puysegur_1_1/RUN_ID100001_20161123_082248.005/Gauges_data_ID100001.nc"
source_zone = 'puysegur'

gauge_netcdf_file = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
     source_zone, '/TSUNAMI_UNIT_SOURCE/', test_file) 

test_sum_tsunami_unit_sources(gauge_netcdf_file)


test_puysegur<-function(){
    source('./get_PTHA_results.R', local=TRUE)
   
    puysegur = get_source_zone_events_data('puysegur')
	model_240 = get_flow_time_series_at_hazard_point(puysegur, 240, c(55015.4, 55042.4))

    max_55015 = max(model_240$flow[['55015.4']][1,,1])
    max_55042 = max(model_240$flow[['55042.4']][1,,1])

    # Values known from previous checks
    er1 = abs(max_55015 - 0.077164)
    er2 = abs(max_55042 - 0.052797)

    m1 = which.max(model_240$flow[['55015.4']][1,,1])
    m2 = which.max(model_240$flow[['55042.4']][1,,1])

    if((er1 < 1.0e-05) & (er2 < 1.0e-05) & (m1 == 62) & (m2 == 76)){
        print('PASS')
    }else{
        print('FAIL: test_puysegur giving different results')

    }
}

test_puysegur()

