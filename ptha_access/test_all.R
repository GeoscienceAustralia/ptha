# Read key codes
source('./R/sum_tsunami_unit_sources.R', local=TRUE)
source('./R/config.R', local=TRUE)

# Make a file on NCI to use to test
test_file = "unit_source_tsunami/RUN_20161121104520_puysegur_1_1/RUN_ID100001_20161123_082248.005/Gauges_data_ID100001.nc"
source_zone = 'puysegur'
gauge_netcdf_file = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
     source_zone, '/TSUNAMI_UNIT_SOURCE/', test_file) 


# Run the unit tests in sum_tsunami_unit_sources.R
test_sum_tsunami_unit_sources(gauge_netcdf_file)


# Make a regression test to check that we get the same results
# for Puysegur. Obviously this will have to be updated each time
# the database changes.
test_puysegur<-function(){
    source('./get_PTHA_results.R', local=TRUE)

    puysegur = get_source_zone_events_data('puysegur')
    model_240 = get_flow_time_series_at_hazard_point(puysegur, 240, 
        hazard_point_ID = c(1.1, 10.1, 22.1, 55015.4, 55042.4))

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


    # Check location info is right -- this is based on 'dput' of an earlier result
    test_data = structure(list(lon = structure(c(129.089294433594, 129.18830871582, 
        127.597557067871, 160.256103515625, 161.841110229492), .Dim = 5L), 
        lat = structure(c(-8.40988159179688, -8.29211044311523, -8.10336685180664, 
        -46.8297233581543, -44.897777557373), .Dim = 5L), elev = structure(c(-775.012756347656, 
        -716.009460449219, -456.208404541016, -5052.998046875, -4819.00732421875
        ), .Dim = 5L), gaugeID = structure(c(1.1, 10.1, 22.1, 55015.4, 
        55042.4), .Dim = 5L)), .Names = c("lon", "lat", "elev", "gaugeID"
        ), row.names = c(NA, -5L), class = "data.frame")

    m1 = max(abs(model_240$location - test_data))
    print(m1)
    if(m1 > 1.0e-010){
        print('FAIL -- location information haz changed')
    }else{
        print('PASS')
    }
 
}
# Run the puysegur regression test
t1 = system.time(test_puysegur())

