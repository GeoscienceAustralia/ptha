source('./R/sum_tsunami_unit_sources.R', local=TRUE)
source('./R/config.R', local=TRUE)

test_file = "unit_source_tsunami/RUN_20161121104520_puysegur_1_1/RUN_ID100001_20161123_082248.005/Gauges_data_ID100001.nc"
source_zone = 'puysegur'

gauge_netcdf_file = paste0(.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/', 
     source_zone, '/TSUNAMI_UNIT_SOURCE/', test_file) 

test_sum_tsunami_unit_sources(gauge_netcdf_file)
