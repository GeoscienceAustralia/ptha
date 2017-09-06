#
# Some code for checking that models have finished [i.e. not been killed before finishing],
# and that the tide gauge outputs exist and are correctly ordered
#

# Read config info
config_env = new.env()
source('config.R', local=config_env)

unit_sources_output_basedir = config_env$all_runs_output_base_dir #'/g/data/.....

# Useful utility for 'fast' tail
tail_file<-function(filename, n=1){
    x = system(paste0('tail -n',n, ' ', filename), intern=TRUE)
    return(x)
}


#'
#' Check that model runs have logfiles, which appear to have finished
#'
check_models_have_been_run<-function(sourcename=config_env$site_name, verbose=FALSE){

    # Note: These file paths might need to be changed in other applications.
    all_logs = Sys.glob(paste0(unit_sources_output_basedir, '/unit_source_tsunami/RUN_*/*/log*'))
    all_unit_sources = Sys.glob(paste0(unit_sources_output_basedir, '/unit_source_tsunami/RUN_*'))

    if(verbose){
        print(paste0('Number of logfiles: ', length(all_logs)))
        print(paste0('Number of unit_sources: ', length(all_unit_sources)))
    }
    
    if((length(all_logs) == 0) | (length(all_unit_sources) == 0)){
        stop('Error: Could not find model files')
    }
    
    # Check that there are logs for every unit source
    all_unit_sources_with_logs = unique(dirname(dirname(all_logs)))
    
    us_match = match(all_unit_sources, all_unit_sources_with_logs)
    if(length(all_unit_sources) != length(all_unit_sources_with_logs)){
        stop('Error: Not all unit sources have logfiles (a)')
    }
    if(length(us_match) != length(unique(us_match))){
        stop('Error: Not all unit sources have logfiles (b)')
    }

    file_not_finished = FALSE
    for(i in 1:length(all_logs)){
        log_end = tail_file(all_logs[i])
        if(!grepl('#####', log_end)){
            file_not_finshed = TRUE
            print(paste0('Unfinished run: ', all_logs[i]))
        }
    }

    if(file_not_finished){ 
        stop('Error: Some runs incomplete')
    }else{
        #print('PASS')
    }
}


#
# Check tide gauge outputs exist, and that the coordinates are all ordered in the same way
#
check_model_gauge_integrity<-function(sourcename=config_env$site_name, verbose=FALSE){

    library(ncdf4)
    all_netcdf_gauges = Sys.glob(paste0(unit_sources_output_basedir,
        '/unit_source_tsunami/RUN_*/*/Gauge*.nc'))

    if(verbose){
        print(paste0('Number of netcdf files: ', length(all_netcdf_gauges)))
    }

    if(length(all_netcdf_gauges) == 0) stop('No netcdf gauges found')
    if(length(all_netcdf_gauges) == 1) stop('Only one netcdf gauge found')

    for(i in 1:length(all_netcdf_gauges)){
        f1 = all_netcdf_gauges[i]
        fid = nc_open(f1, readunlim=FALSE)
        lon = try(ncvar_get(fid, varid='lon', count=-1))
        if(class(lon) == 'try-error'){
            stop(paste0('Error in reading lon in file: ', f1))
        }
        lat = try(ncvar_get(fid, varid='lat', count=-1))
        if(class(lat) == 'try-error'){
            stop(paste0('Error in reading lat in file: ', f1))
        }
        time = try(ncvar_get(fid, varid='time', count=-1))
        if(class(time) == 'try-error'){
            stop(paste0('Error in reading time in file: ', f1))
        }

        nc_close(fid)

        if(i == 1){
            lon0 = lon
            lat0 = lat
            time0 = time
        }
        
        l1 = length(lon) 
        l2 = length(lat) 
        if(l1 != l2){
            stop(paste0('Error: Lon and Lat have different lengths, file:', f1))
        }

        if(i == 1) next

        if(l1 != length(lon0)){
            stop(paste0('Error: Lon and Lat have different lengths to first gauge file @file:', f1))
        }

        if(any((lon - lon0) != 0)){
            stop(paste0('Error: Lons do not agree with the first gauge file @file:', f1))
        }

        if(any((lat - lat0) != 0)){
            stop(paste0('Error: Lats do not agree with the first gauge file @file:', f1))
        }
    
        if(length(time) != length(time0)){
            stop(paste0('Time length does not agree with first gauge file @file: ', f1)) 
        }

        if(any(time - time0 != 0)){
            stop(paste0('Times do not agree with first gauge @file: ', f1))
        }

    }

    #print('PASS')
}
