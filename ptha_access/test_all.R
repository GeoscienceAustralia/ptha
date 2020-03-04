# Read key codes
source('R/sum_tsunami_unit_sources.R', local=TRUE)
source('R/config.R', local=TRUE)

# Make a file on NCI to use to test
#test_file = "unit_source_tsunami/RUN_20161121104520_puysegur_1_1/RUN_ID100001_20161123_082248.005/Gauges_data_ID100001.nc"
#source_zone = 'puysegur'
#gauge_netcdf_file = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
#     source_zone, '/TSUNAMI_UNIT_SOURCE/', test_file) 

# Find a file that contains hazard points. Easiest way is to read them from a tide gauge file
unit_source_stats_puysegur = paste0(.GDATA_OPENDAP_BASE_LOCATION, 
    'SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/unit_source_statistics_puysegur2.nc')
fid = nc_open(unit_source_stats_puysegur)
gauge_netcdf_file = ncvar_get(fid, 'tide_gauge_file', start=c(1, 1), count=c(fid$dim$max_nchar$len,1))[1]
nc_close(fid)
gauge_netcdf_file = adjust_path_to_gdata_base_location(gauge_netcdf_file)

# Run the unit tests in sum_tsunami_unit_sources.R
test_sum_tsunami_unit_sources(gauge_netcdf_file)


# Make a regression test to check that we get the same results
# for Puysegur. Obviously this will have to be updated each time
# the database changes.
test_puysegur2<-function(){
    source('./get_PTHA_results.R', local=TRUE)

    puysegur = get_source_zone_events_data('puysegur2')
    model_3051 = get_flow_time_series_at_hazard_point(puysegur, 3051, 
        hazard_point_ID = c(1.1, 10.1, 22.1, 55015.4, 55042.4))

    max_55015 = max(model_3051$flow[['55015.4']][1,,1])
    max_55042 = max(model_3051$flow[['55042.4']][1,,1])

    # Values known from previous checks
    er1 = abs(max_55015 - 0.1771731) #0.050279)
    er2 = abs(max_55042 - 0.3976239) #0.138759)

    m1 = which.max(model_3051$flow[['55015.4']][1,,1])
    m2 = which.max(model_3051$flow[['55042.4']][1,,1])

    if((er1 < 1.0e-05) & (er2 < 1.0e-05) & (m1 == 59) & (m2 == 61)){
        print('PASS')
    }else{
        print('FAIL: test_puysegur giving different results')

    }


    # Check location info is right -- this is based on 'dput' of an earlier result
    test_data = structure(list(lon = structure(c(129.089294433594, 129.18830871582, 
        127.597557067871, 160.256103515625, 161.841110229492), .Dim = 5L), 
        lat = structure(c(-8.40988159179688, -8.29211044311523, -8.10336685180664, 
        -46.8297233581543, -44.897777557373), .Dim = 5L), elev = structure(c(-775.006286621094, 
        -716.009460449219, -455.951629638672, -5052.998046875, -4819.00732421875
        ), .Dim = 5L), gaugeID = structure(c(1.1, 10.1, 22.1, 55015.4, 
        55042.4), .Dim = 5L)), class = "data.frame", row.names = c(NA, 
        -5L))

    m1 = max(abs(model_3051$location - test_data))
    if(m1 > 1.0e-010){
        print('FAIL -- location information haz changed')
    }else{
        print('PASS')
    }

    # Check if we unpack to gauge-based list, it still works
    model_3051b = get_flow_time_series_at_hazard_point(puysegur, 3051, 
        hazard_point_ID = c(1.1, 10.1, 22.1, 55015.4, 55042.4),
        store_by_gauge=FALSE)


    if(length(model_3051b$flow) == dim(model_3051$flow[[1]])[1]){
        print('PASS')
    }else{
        print('FAIL: store_by_gauge not working as desired')
    }

    FAILED=FALSE
    for(i in 1:length(model_3051b$flow)){
        for(j in 1:length(model_3051$flow)){
            if(! all(model_3051b$flow[[i]][j,,] == model_3051$flow[[j]][i,,])){
                FAILED=TRUE
            }
        }
    }
    if(FAILED){
        print('FAIL: store_by_gauge not ordered as desired')
    }else{
        print('PASS')
    }
 
}
# Run the puysegur regression test
t1 = system.time(test_puysegur2())


# Test that we can read a scenario with many unit sources -- i.e. one that
# would fail if the netcdf install is not sufficiently up to date.

test_large_event_read<-function(){

    source('./get_PTHA_results.R', local=TRUE)

    expected_event_index_string = "324-325-329-330-333-334-338-344-345-348-349-350-351-352-354-355-358-359-361-362-363-365-366-367-369-370-372-373-376-377-381-384-385-386-387-388-389-390-391-392-393-394-395-396-397-398-399-400-401-402-403-404-405-406-407-408-409-410-411-412-413-414-415-416-417-418-419-420-421-422-423-424-425-426-427-428-430-431-434-435-436-438-439-440-443-454-461-462-"

    expected_event_slip_string = "2.428_0.4631_6.739_1.026_6.742_8.005_7.856_2.007_2.8_7.307_1.599_19.78_15.27_4.94_21.14_29.31_25.13_18.47_10.37_19.73_20.71_26.1_4.098_3.605_35.83_6.34_6.517_34.34_3.984_30.11_41.04_0.7955_57.45_32.88_25.11_35.7_77.66_84.67_64.68_56.23_87.42_117_91.92_68.23_114.9_161.2_142.1_109.2_162.3_196.5_198.3_153.7_186_196.9_206.4_182.5_169.6_163_187.9_181.4_127.7_133.5_155.5_144_85.38_102.2_105.1_73.65_30.45_56.9_59.6_36.4_8.488_38.66_80.98_16.74_46.75_70.94_51.39_50.19_24.28_20.31_24.2_14.35_7.434_1.672_2.327_2.264_"

    x = get_source_zone_events_data('southamerica', slip_type='stochastic', desired_event_rows=150000)

    test_result = (x$events$event_slip_string == expected_event_slip_string &
        x$events$event_index_string == expected_event_index_string)

    if(test_result){
        print('PASS')
    }else{
        print('FAIL -- incorrect read of a long string of event metadata. This suggests your netcdf install is not sufficiently up to date (there was a bug in the remote reading of long character strings in netcdf-c, whcih was fixed in the bleeding edge source in late 2017)')
    }

}

test_large_event_read()
