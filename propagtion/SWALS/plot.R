#library(ncdf4)

#' Get name of folder most recently written to OUTPUTS (according to timestamp in name)
get_recent_output_folder<-function(){
    outputs = dir('OUTPUTS', pattern='RUN_')
    output_folder = sort(outputs, decreasing=TRUE)[1]
    output_folder = paste0('OUTPUTS/', output_folder)
    return(output_folder)
}

#' Get times at which output is written
get_time<-function(output_folder){
    model_time = try(scan(Sys.glob(paste0(output_folder, '/', 'Time_*'))[1]), silent=TRUE)
    return(model_time)
}

#' Get the model dimensions
get_model_dim<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model dimensions
    model_dim_line = grep('nx :', file_metadata)
    model_dim = as.numeric(scan(text=file_metadata[model_dim_line], what='character')[3:4])
    return(model_dim)
}

#' Get the model cell size
get_dx<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model cell size
    model_dim_line = grep('dx :', file_metadata)
    model_dx = as.numeric(scan(text=file_metadata[model_dim_line], what='character')[3:4])
    return(model_dx)
}

get_lower_left_corner<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model dimensions
    model_dim_line = grep('lower_left_corner:', file_metadata)
    model_dim = as.numeric(scan(text=file_metadata[model_dim_line], what='character')[2:3])
    return(model_dim)
}

get_model_output_precision<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])
    model_dim_line = grep('output_precision', file_metadata)
    dp = as.numeric(scan(text=file_metadata[model_dim_line], what='character')[2])
    return(dp)
}

get_gauges_mixed_binary_format<-function(output_folder){
    # Read the gauge data, which is in a mixture of ascii and binary

    real_precision = get_model_output_precision(output_folder)

    gauge_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Gauges_*nc_metadata')))
    times = get_time(output_folder)

    # If we didn't save times, guess a large number
    #if(length(times) == 0) times = rep(NA, 10000)

    # Get static variables (coordinates, elevation, ids, and maybe other things)
    static_var_line = grep('static_var:', gauge_metadata)
    nstatic_var = length(scan(text = gauge_metadata[static_var_line], what='character')[-1]) + 3
    static_var_ids = as.numeric(scan(text = gauge_metadata[static_var_line], what='character')[-1])
    static_var_possible_names = c('stage0', 'uh0', 'vh0', 'elevation0')

    static_var_file = Sys.glob(paste0(output_folder, '/', 'Gauges*nc_static'))[1]
    lt = 1e+05
    static_var = readBin(static_var_file, size=real_precision, n = nstatic_var * lt,
        what='numeric')
    static_var = matrix(static_var, ncol=nstatic_var, byrow=TRUE)

    # Get time-series
    time_series_var_line = grep('time_series_var:', gauge_metadata)
    ntime_series_var = length(scan(text = gauge_metadata[time_series_var_line], what='character')[-1])
    time_series_var_ids = as.numeric(scan(text = gauge_metadata[time_series_var_line], what='character')[-1])
    time_series_var_possible_names = c('stage', 'uh', 'vh', 'elevation')

    time_series_var_file = Sys.glob(paste0(output_folder, '/', 'Gauges*nc_timeseries'))[1]
    time_series_var = readBin(time_series_var_file, size = real_precision, 
        n = ntime_series_var * lt * length(static_var[,1]), what='numeric')
    lt = length(time_series_var)/(ntime_series_var * length(static_var[,1]))
    dim(time_series_var) = c(length(static_var[,1]), ntime_series_var, lt)

    # Read static variables. First 2 columns are lon/lat, so ignore those
    static_var_list=list()
    counter=2
    for(i in 1:length(static_var_possible_names)){
        nm = static_var_possible_names[i]
        if(i %in% static_var_ids){
            counter=counter+1
            static_var_list[[nm]] = static_var[,counter]
        }else{
            static_var_list[[nm]] = NA
        }
    }

    # Read time-series variables
    time_series_var_list=list()
    counter=0
    for(i in 1:length(time_series_var_possible_names)){
        nm = time_series_var_possible_names[i]
        if(i %in% time_series_var_ids){
            counter=counter+1
            time_series_var_list[[nm]] = time_series_var[,counter,]
        }else{
            time_series_var_list[[nm]] = NA
        }
    }

    # Make sure this has the same structure as for netcdf gauge output
    output = list(lon = static_var[,1], lat=static_var[,2], 
        gaugeID=static_var[,ncol(static_var)],
        static_var=static_var_list, 
        time_var = time_series_var_list)

    return(output)        
}

#' Read gauges from netcdf file.
get_gauges_netcdf_format<-function(output_folder){
    library(ncdf4)

    gauge_fid = nc_open(Sys.glob(paste0(output_folder, '/', 'Gauges_data_*.nc'))[1])
    lon = ncvar_get(gauge_fid, 'lon')
    lat = ncvar_get(gauge_fid, 'lat')
    time = ncvar_get(gauge_fid, 'time')
    gaugeID = ncvar_get(gauge_fid, 'gaugeID')

    gauge_output_names = names(gauge_fid)    

    time_series_names = c('stage', 'uh', 'vh', 'elevation')
    static_names = paste0(time_series_names, '0')

    static_var = list()
    for(i in 1:length(static_names)){
        tmp = try(ncvar_get(gauge_fid, static_names[i]), silent=TRUE)
        if(class(tmp) == 'try-error'){
            static_var[[static_names[i]]] = NA
        }else{
            static_var[[static_names[i]]] = tmp
        }
    }

    time_series_var = list()
    for(i in 1:length(time_series_names)){
        tmp = try(ncvar_get(gauge_fid, time_series_names[i]), silent=TRUE)
        if(class(tmp) == 'try-error'){
            time_series_var[[time_series_names[i]]] = NA
        }else{
            time_series_var[[time_series_names[i]]] = tmp
        }
    }

    outputs = list(lon=lon, lat=lat, time=time, gaugeID=gaugeID, 
        static_var=static_var, time_var=time_series_var)
    return(outputs)
}

# Read either the netcdf gauges OR the mixed binary format
get_gauges<-function(output_folder){

    outputs = try(get_gauges_mixed_binary_format(output_folder), silent=TRUE)
    if(class(outputs) == 'try-error'){
        outputs = try(get_gauges_netcdf_format(output_folder), silent=TRUE)
    }

    return(outputs)
}

#' Quickly read model outputs
#' @param var Start of name of file to read (containing a flow variable)
#' @param binary 
#' @param drop_walls If TRUE, then set values on the boundaries to values just inside the
#' boundaries. Useful when we use an extreme elevation for a reflective condition
#' @param output_folder folder containing output files. If NULL, read the most recent one
#' @return The variable as a 3d array 
get_gridded_variable<-function(var = 'Var_1', drop_walls=FALSE, output_folder = NULL){
   
    if(is.null(output_folder)){ 
        output_folder = get_recent_output_folder()
    }
    
    # FIXME: Update for parallel
    file_to_read = Sys.glob(paste0(output_folder, '/', var, '*'))[1]

    model_dim = get_model_dim(output_folder)

    # Read the binary file, using 'while' to get all the data
    file_connection = file(file_to_read, open='rb')

    real_precision = get_model_output_precision(output_folder)
    output_stage_scan = c()
    new_output_stage_scan = readBin(file_connection, what='numeric', n = prod(model_dim), size=real_precision)
    while(length(new_output_stage_scan)>0){
        output_stage_scan = c(output_stage_scan, new_output_stage_scan)
        new_output_stage_scan = readBin(file_connection, what='numeric', n = prod(model_dim), size=real_precision)
    }
    close(file_connection)

    num_times = length(output_stage_scan)/prod(model_dim)

    dim(output_stage_scan) = c(model_dim, num_times)

    if(drop_walls){
        output_stage_scan[1,,] = output_stage_scan[2,,]
        output_stage_scan[,1,] = output_stage_scan[,2,]
        output_stage_scan[model_dim[1],,] = output_stage_scan[model_dim[1]-1,,]
        output_stage_scan[,model_dim[2],] = output_stage_scan[, model_dim[2]-1,]
    }
   
    return(output_stage_scan) 
}


get_all_recent_results<-function(output_folder=NULL){
    if(is.null(output_folder)){
        output_folder = get_recent_output_folder()
    }
    stage = try(get_gridded_variable(var='Var_1', output_folder=output_folder), silent=TRUE)
    ud = try(get_gridded_variable(var='Var_2', output_folder=output_folder), silent=TRUE)
    vd = try(get_gridded_variable(var='Var_3', output_folder=output_folder), silent=TRUE)
    elev = try(get_gridded_variable(var='Var_4', output_folder=output_folder), silent=TRUE)
    # maxQ might not exist if the simulation has not finished
    maxQ = try(get_gridded_variable(var='Max_quantities', output_folder=output_folder), silent=TRUE)

    # time is often not written until the sim is finished
    time = try(get_time(output_folder), silent=TRUE) 
    nx = try(get_model_dim(output_folder), silent=TRUE)
    dx = try(get_dx(output_folder), silent=TRUE)
    lower_left_corner = try(get_lower_left_corner(output_folder), silent=TRUE)

    xs = lower_left_corner[1] + ((1:nx[1]) - 0.5) * dx[1]
    ys = lower_left_corner[2] + ((1:nx[2]) - 0.5) * dx[2]

    gauges = try(get_gauges(output_folder), silent=TRUE)

    return(list(stage=stage, ud=ud, vd=vd, elev=elev, maxQ=maxQ, time=time, 
        output_folder = output_folder, nx=nx, gauges = gauges, xs = xs, ys=ys, 
        lower_left_corner=lower_left_corner, dx=dx))
}

## Nice plot example
#persp(seq(0,200, len=700), seq(0,200, len=700), X[,,11], phi=40, border=NA,
#    col='green', shade=1, axes=FALSE, box=FALSE, asp=1, scale=FALSE)


# X = get_all_recent_results(c(1100, 800))
# library(raster)
# E = raster(t(X$elev[,,1]), xmn=495000, xmx=(495000+11000), ymn=1610000, ymx=1610000+8000, 
#     crs=CRS("+init=epsg:3123"))
# E = flip(E, direction='y')
# library(rasterVis)
# plot3D(E)

# Convert peak stage output to raster
#' @param swals_out result of get_all_recent_results
#' @param proj4string
#' @param na_above_zero logical. If TRUE, then when elevation is > 0, set all stages to NA
make_max_stage_raster<-function(swals_out, proj4string='+init=epsg:4326', na_above_zero=FALSE){
    library(raster)

    stg = swals_out$maxQ[,,1]
    if(na_above_zero){
        stg[swals_out$maxQ[,,2] > 0] = NA
    }

    dx = swals_out$xs[2] - swals_out$xs[1]
    dy = swals_out$ys[2] - swals_out$ys[1]
    xs = swals_out$xs
    ys = swals_out$ys
    max_stage = raster(t(stg), xmn=min(xs)-dx/2,
        xmx = max(xs)+dx/2, ymn=min(ys)-dy/2, ymx=max(ys)+dy/2)
    max_stage = flip(max_stage, direction='y')
    proj4string(max_stage) = proj4string
    return(max_stage)
}

