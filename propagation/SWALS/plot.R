#
# Scripts to work with SWALS output files. 
#
# It is not essential to use R to run SWALS. However, some kind of scripting
# language is useful to work with the outputs, especially if using a
# multidomain and distributed-memory parallel jobs. 
#
# In general any single domain will contain "halo regions" (i.e. where the flow
# state is received from other domains), and other regions where it is the
# "priority domain". In the halo regions the values that are stored in output
# files should not be used (because some other domain is the priority domain in
# those regions). In sum, care is required to combine results from multiple domains
# while only using "priority domain" cells. 
#


#' Get name of folder most recently written to OUTPUTS (according to timestamp in name)
get_recent_output_folder<-function(){
    outputs = dir('OUTPUTS', pattern='RUN_')
    output_folder = sort(outputs, decreasing=TRUE)[1]
    output_folder = paste0('OUTPUTS/', output_folder)
    return(output_folder)
}

#' Get times at which grid output is written
get_time<-function(output_folder){
    model_time = try(scan(Sys.glob(paste0(output_folder, '/', 'Time_*'))[1], quiet=TRUE), silent=TRUE)
    return(model_time)
}

#' Get the model dimensions
get_model_dim<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model dimensions
    model_dim_line = grep('nx :', file_metadata)
    model_dim = as.numeric(scan(text=file_metadata[model_dim_line], what='character', quiet=TRUE)[3:4])
    return(model_dim)
}

#' Get the model cell size
get_dx<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model cell size
    model_dim_line = grep('dx :', file_metadata)
    model_dx = as.numeric(scan(text=file_metadata[model_dim_line], what='character', quiet=TRUE)[3:4])
    return(model_dx)
}

#' Get the domain lower-left corner
get_lower_left_corner<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model dimensions
    model_dim_line = grep('lower_left_corner:', file_metadata)
    model_dim = as.numeric(scan(text=file_metadata[model_dim_line], what='character', quiet=TRUE)[2:3])
    return(model_dim)
}

#' Get the domain output precision
get_model_output_precision<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])
    model_dim_line = grep('output_precision', file_metadata)
    dp = as.numeric(scan(text=file_metadata[model_dim_line], what='character', quiet=TRUE)[2])
    return(dp)
}

#' Get the domain outputs in a DEFUNCT mixed-binary format
#'
#' This has been superceeded by netcdf -- but is occasionally useful if one has
#' netcdf compilation problems and wants to test something out. However not all
#' the functionality is supported with this file format, so it isn't a general
#' solution to netcdf problems.
#'
get_gauges_mixed_binary_format<-function(output_folder){
    # Read the gauge data, which is in a mixture of ascii and binary

    real_precision = get_model_output_precision(output_folder)

    gauge_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Gauges_*nc_metadata')))
    times = get_time(output_folder)

    # If we didn't save times, guess a large number
    #if(length(times) == 0) times = rep(NA, 10000)

    # Get static variables (coordinates, elevation, ids, and maybe other things)
    static_var_line = grep('static_var:', gauge_metadata)
    nstatic_var = length(scan(text = gauge_metadata[static_var_line], what='character', quiet=TRUE)[-1]) + 3
    static_var_ids = as.numeric(scan(text = gauge_metadata[static_var_line], what='character', quiet=TRUE)[-1])
    static_var_possible_names = c('stage0', 'uh0', 'vh0', 'elevation0')

    static_var_file = Sys.glob(paste0(output_folder, '/', 'Gauges*nc_static'))[1]
    lt = 1e+05
    static_var = readBin(static_var_file, size=real_precision, n = nstatic_var * lt,
        what='numeric')
    static_var = matrix(static_var, ncol=nstatic_var, byrow=TRUE)

    # Get time-series
    time_series_var_line = grep('time_series_var:', gauge_metadata)
    ntime_series_var = length(scan(text = gauge_metadata[time_series_var_line], what='character', quiet=TRUE)[-1])
    time_series_var_ids = as.numeric(scan(text = gauge_metadata[time_series_var_line], what='character', quiet=TRUE)[-1])
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
        time_var = time_series_var_list,
        priority_gauges_only=FALSE)

    return(output)        
}

#' Read domain gauges from netcdf file. 
#'
#' In general this is the output format we actually use -- the alternative 'home-brew binary format'
#' is not comprehensively supported.
#'
#' @param output_folder The domain's output folder
#' 
get_gauges_netcdf_format<-function(output_folder){
    library('ncdf4')

    gauge_fid = try(
        nc_open(Sys.glob(paste0(output_folder, '/', 'Gauges_data_*.nc'))[1], readunlim=FALSE),
        silent=TRUE)
    if(is(gauge_fid,'try-error')){
        return(gauge_fid)        
    }

    # Gauge coordinates
    lon = ncvar_get(gauge_fid, 'lon')
    lat = ncvar_get(gauge_fid, 'lat')
    time = ncvar_get(gauge_fid, 'time')
    gaugeID = ncvar_get(gauge_fid, 'gaugeID')

    gauge_output_names = names(gauge_fid)    

    time_series_names = c('stage', 'uh', 'vh', 'elevation')
    # The variables that are only written once [e.g. often chosen for
    # elevation] have a zero appended to the name
    static_names = paste0(time_series_names, '0')

    # Read static variables (i.e. things that are not written every time-step
    # -- sometime elevation0)
    static_var = list()
    for(i in 1:length(static_names)){
        if(static_names[i] %in% names(gauge_fid$var)){
            tmp = try(ncvar_get(gauge_fid, static_names[i]), silent=TRUE)
            if(is(tmp, 'try-error')) tmp = NA
        }else{
            tmp = NA
        }
        static_var[[static_names[i]]] = tmp
    }

    # Read teh time-series variables (i.e. typically stage, uh, vh)
    time_series_var = list()
    for(i in 1:length(time_series_names)){
        if(time_series_names[i] %in% names(gauge_fid$var)){
            tmp = try(ncvar_get(gauge_fid, time_series_names[i]), silent=TRUE)
            if(is(tmp, 'try-error')){
                tmp = NA
            }
            # Force it to be a matrix even if there is only 1 gauge
            if(length(dim(tmp)) < 2) dim(tmp) = c(1, length(tmp))
        }else{
            tmp = NA
        }
        time_series_var[[time_series_names[i]]] = tmp
    }

    # Check for the 'priority_gauges_only' variable -- this tells us if all the
    # stored gauges are inside the priority domain. If TRUE, that's convenient, because
    # we don't need to worry about throwing away non-priority-domain gauges.
    x = ncatt_get(gauge_fid, varid=0, 'priority_gauges_only')
    priority_gauges_only = FALSE
    if(x$hasatt){
        if(x$value == 'true') priority_gauges_only=TRUE
    }

    nc_close(gauge_fid)

    outputs = list(lon=lon, lat=lat, time=time, gaugeID=gaugeID, 
        static_var=static_var, time_var=time_series_var, 
        priority_gauges_only=priority_gauges_only)
    return(outputs)
}

#' Read EITHER the netcdf gauges OR the DEFUNCT mixed binary format
#'
get_gauges<-function(output_folder){

    outputs = try(get_gauges_netcdf_format(output_folder), silent=TRUE)
    if(is(outputs, 'try-error')){
        outputs = try(get_gauges_mixed_binary_format(output_folder), silent=TRUE)
    }

    return(outputs)
}


#' Get a global attribute from a netcdf file
#'
#' The netcdf file can be provided as a filename (in which case the function
#' opens,reads, and closes the file) or an existing file-handle (in which case
#' the function only reads the attribute)
#'
#' @param nc_file netcdf file name, or file object created by a call to 'ncdf4::nc_open'
#' @param attname name of global attribute in the file
#
get_nc_global_attribute<-function(nc_file, attname){
    if(is(nc_file, 'character')){
        fid = nc_open(nc_file, readunlim=FALSE)
        bbox = ncatt_get(fid, varid=0, attname=attname)
        nc_close(fid)
    }else{
        bbox = ncatt_get(fid, varid=0, attname=attname)
    }
    return(bbox)
}

#' Quickly read domain gridded outputs
#'
#' This function probably should not be called directly -- but it is called by
#' get_all_recent_results etc. FIXME: It should probably be re-written to use more sensible names. 
#'
#' @param var A variable name -- beware the names are archaic (derived from the
#'        old home-brew binary file format). The accepted variable names are denoted 
#'        below [left-hand-side] along with the mapping to more sensible names 
#'        in the netcdf file:
#'            Var_1 = 'stage'
#'            Var_2 = 'uh'
#'            Var_3 = 'vh'
#'            Var_4 = 'elev'
#'            x = 'x'
#'            y = 'y'
#'            Max_quantities = 'max_stage'
#'            elev0 = 'elevation0'
#'            is_priority_domain='is_priority_domain'
#' @param output_folder folder containing the domain output files. If NULL, read the most recent one
#' @return The variable as a 3d array 
#'
get_gridded_variable<-function(var = 'Var_1', output_folder = NULL){
   
    if(is.null(output_folder)){ 
        output_folder = get_recent_output_folder()
    }

    nc_file = Sys.glob(paste0(output_folder, '/Grid_output_*.nc'))

    # If netcdf output does not exist, then use home-brew binary (DEFUNCT)
    if(length(nc_file) != 1){
        # Home-brew binary -- not really well supported, avoid this.
    
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

    }else{

        # Netcdf output -- this is the usual approach. Here the naming convention is inherited from
        # the hacky old binary format.
        var_list = list(Var_1 = 'stage', Var_2 = 'uh', Var_3 = 'vh', Var_4 = 'elev', x = 'x', y = 'y',
            Max_quantities = 'max_stage', elev0 = 'elevation0', is_priority_domain='is_priority_domain')
        library('ncdf4')
        fid = nc_open(nc_file, readunlim=FALSE)
        var_name = var_list[[var]]
        if(var_name %in% c(names(fid$var), names(fid$dim))){
            output_stage_scan = try(ncvar_get(fid, var_name), silent=TRUE)
        }
        nc_close(fid)
    }
   
    return(output_stage_scan) 
}

#' Convenient function to read in pretty much all of the results in a SINGLE DOMAIN.
#'
#' For multidomains, see 'get_multidomain' (or just manually loop over all domains).
#' The default options in this function can sometimes be prohibitive due to memory issues, so
#' there are options to skip reading grids and gauges (i.e. to obtain minimal domain metadata,
#' which can later be used for more efficient reading of data subsets).
#' 
#' @param output_folder the domain output folder
#' @param read_grids Should we read the stage/uh/vh/elev grids? These can require large memory.
#' @param read_gauges Should we read the gauge outputs?
#' @param read_time_and_geometry Should we read the grid output time and geometry information?
#' @param always_read_priority_domain Only matters if read_grids=FALSE. In that case, setting this
#'        to TRUE ensures the priority domain information is read. By default, it is always read.
#' @param always_read_xy Only matters if read_grids=FALSE. In that case, setting this to TRUE ensures
#'        the x/y grid coordinates are read. Otherwise they are always read. 
#' @param quiet Try to minimise printed output. Note we cannot completely control this due to the ncdf4
#'        package (which reports errors in a manner that we can't seem to suppress).
#' @param always_read_max_grids Only matters if read_grids=FALSE. In that case, setting this to TRUE ensures
#'        the max_stage and elevation are read. Otherwise they are always read. 
#' @return A list with all the domain data.
get_all_recent_results<-function(output_folder=NULL, read_grids=TRUE, read_gauges=TRUE, 
                                 read_time_and_geometry = TRUE, 
                                 always_read_priority_domain=FALSE, 
                                 always_read_xy = read_time_and_geometry,
                                 quiet=FALSE, 
                                 always_read_max_grids=FALSE){

    if(quiet) sink(tempfile())

    if(is.null(output_folder)){
        output_folder = get_recent_output_folder()
    }

    if(read_grids){
        stage = try(get_gridded_variable(var='Var_1', output_folder=output_folder), silent=TRUE)
        ud = try(get_gridded_variable(var='Var_2', output_folder=output_folder), silent=TRUE)
        vd = try(get_gridded_variable(var='Var_3', output_folder=output_folder), silent=TRUE)
        elev = try(get_gridded_variable(var='Var_4', output_folder=output_folder), silent=TRUE)
        # maxQ might not exist if the simulation has not finished
        maxQ = try(get_gridded_variable(var='Max_quantities', output_folder=output_folder), silent=TRUE)
        elev0 = try(get_gridded_variable(var='elev0', output_folder=output_folder), silent=TRUE)
        is_priority_domain = try(get_gridded_variable(var='is_priority_domain', output_folder=output_folder), silent=TRUE)
    }else{
        stage = NULL
        ud = NULL
        vd = NULL
        elev = NULL
        maxQ = NULL
        elev0 = NULL
        is_priority_domain = NULL
        if(always_read_max_grids){
            maxQ = try(get_gridded_variable(var='Max_quantities', output_folder=output_folder), silent=TRUE)
            elev0 = try(get_gridded_variable(var='elev0', output_folder=output_folder), silent=TRUE)
        }
    }

    if(always_read_priority_domain & !read_grids){
        is_priority_domain = try(get_gridded_variable(var='is_priority_domain', output_folder=output_folder), silent=TRUE)
    }

    if(read_time_and_geometry){
        # time is often not written until the sim is finished
        time = try(get_time(output_folder), silent=TRUE) 
        nx = try(get_model_dim(output_folder), silent=TRUE)
        dx = try(get_dx(output_folder), silent=TRUE)
        lower_left_corner = try(get_lower_left_corner(output_folder), silent=TRUE)
    }else{
        time = NULL
        nx = NULL
        dx = NULL
        lower_left_corner = NULL
    }

    if(always_read_xy){
        # x and y --  get from netcdf if we can
        xs = try(get_gridded_variable(var='x', output_folder=output_folder), silent=TRUE)
        if(!is(xs, 'try_error')){
            ys = try(get_gridded_variable(var='y', output_folder=output_folder), silent=TRUE)
        }else{
            # Old-style, beware this does not work with decimated output
            xs = lower_left_corner[1] + ((1:nx[1]) - 0.5) * dx[1]
            ys = lower_left_corner[2] + ((1:nx[2]) - 0.5) * dx[2]
        }
    }else{
        xs = NULL
        ys = NULL
    }

    if(read_gauges){
        gauges = try(get_gauges(output_folder), silent=TRUE)
    }else{
        gauges = NULL
    }

    if(quiet) sink()

    return(list(stage=stage, ud=ud, vd=vd, elev=elev, maxQ=maxQ, elev0=elev0, 
        is_priority_domain = is_priority_domain, time=time, 
        output_folder = output_folder, nx=nx, gauges = gauges, xs = xs, ys=ys, 
        lower_left_corner=lower_left_corner, dx=dx, lw=(nx*dx)))
}

#' Read all domains in a multidomain directory into a list.
#'
#' This just applies get_all_recent_results to all domain output folders in a multidomain.
#'
#' @param multidomain_dir the directory with the results
#' @param ... Further arguments to get_all_recent_results
#' @return A list with one entry per domain, containing the result of get_all_recent_results
get_multidomain<-function(multidomain_dir, ...){

    all_domain_files = Sys.glob(paste0(multidomain_dir, '/RUN_*'))
    md = lapply(all_domain_files, f<-function(x) get_all_recent_results(x, ...))
    return(md)
}

#' Get the times at which grids or gauges are output in the multidomain
#'
#' Sometimes it is useful to know how many times there are
#'
#' @param multidomain_dir Directory containing the multidomain
#' @param output_type either 'grids' (default) or 'gauges'
#' @return a vector of times at which the 'grids' or 'gauges' are output
#'
get_multidomain_output_times<-function(multidomain_dir, output_type='grids'){
    library(ncdf4)

    if(!file.exists(multidomain_dir)){
        stop(paste0('Could not find multidomain_dir ', multidomain_dir))
    }

    # Find the netcdf grid or gauges files, depending 
    if(output_type =='grids'){
        target_files = Sys.glob(paste0(multidomain_dir, '/RUN_*/Grid_output_*.nc'))
    }else if(output_type == 'gauges'){
        target_files = Sys.glob(paste0(multidomain_dir, '/RUN_*/Gauges_data_*.nc'))
    }else{
        stop(paste0('output_type should be "grids" or "gauges", but I got: ', output_type))
    }


    if(length(target_files) < 1){
        stop(paste0('Could not find any files holding ', output_type, ' in multidomain ', multidomain_dir))
    }

    fid = nc_open(target_files[1], readunlim=FALSE)
    time = as.numeric(ncvar_get(fid, 'time'))
    nc_close(fid)

    return(time)
}

#' Convert peak stage output to raster using the domain object.
#'
#' Instead consider using 'merge_domains_nc_grids' which will combine partitioned domains, and
#' can output rasters. OR if you just want to make a plot, consider 'multidomain_image', which
#' works for the full multidomain.
#'
#' @param swals_out result of get_all_recent_results
#' @param proj4string
#' @param na_above_zero logical. If TRUE, then when elevation is > 0, set all stages to NA
#' @param return_elevation logical. If TRUE, return a list with max_stage AND elevation as rasters
#' @param na_outside_priority_domain If TRUE, set the raster to NA in non-priority-domain regions. 
#'        This is a good idea to simplify interpretation.
#' @return Either the max_stage_raster (if return_elevation=FALSE), or a list with both the max_stage
#'        and elevation rasters.
make_max_stage_raster<-function(swals_out, proj4string='+init=epsg:4326', na_above_zero=FALSE, 
                                return_elevation=FALSE, na_outside_priority_domain=FALSE){
    library('raster')

    if(length(dim(swals_out$maxQ)) == 3){
        # Using old binary format
        stg = swals_out$maxQ[,,1]
        elev = swals_out$maxQ[,,2]
    }else{
        # Using netcdf format 
        stg = swals_out$maxQ
        elev = swals_out$elev0
    }

    if(na_above_zero){
        stg[elev > 0] = NA
    }

    if(na_outside_priority_domain){
        priority_region = swals_out$is_priority_domain
        stg[!priority_region] = NA
        elev[!priority_region] = NA
    }

    # Make max-stage raster
    dx = swals_out$xs[2] - swals_out$xs[1]
    dy = swals_out$ys[2] - swals_out$ys[1]
    xs = swals_out$xs
    ys = swals_out$ys
    max_stage = raster(t(stg), xmn=min(xs)-dx/2,
        xmx = max(xs)+dx/2, ymn=min(ys)-dy/2, ymx=max(ys)+dy/2)
    max_stage = flip(max_stage, direction='y')
    proj4string(max_stage) = proj4string

    if(!return_elevation){
        return(max_stage)
    }else{
        # Make elevation raster
        elevation = raster(t(elev), xmn=min(xs)-dx/2,
            xmx = max(xs)+dx/2, ymn=min(ys)-dy/2, ymx=max(ys)+dy/2)
        elevation = flip(elevation, direction='y')
        proj4string(elevation) = proj4string

        output = list(max_stage = max_stage, elevation=elevation)
        return(output)
    }
    
        
}

#' Make an image plot for a multidomain
#'
#' This function will read the required output data from all domains in the multidomain,
#' and combine into a single plot. It generally doesn't need to use much memory, and is
#' particularly useful in cases where reading all the multidomain data would be prohibitive.
#'
#' @param multidomain_dir directory where all multidomain outputs live
#' @param variable 'max_stage' or 'stage' or any other variable in the netcdf file
#' @param time_index if the variable has a time dimension, the index which we plot
#' @param xlim xlimit for image
#' @param ylim ylimit for image
#' @param zlim z limits for variable in image plot
#' @param col colourscheme for variable in image plot
#' @param add if TRUE, add to an existing plot
#' @param var_transform_function if not NULL, a function that is used to transform the variable before plotting
#' @param NA_if_stage_not_above_elev logical. If TRUE, the set regions with stage <= (elev + 1e-03) to NA
#' @param use_fields logical. If TRUE, use image.plot from the fields package. Otherwise use graphics::image
#' @param clip_to_zlim logical. If TRUE, clip the variable limits to be within zlim before plotting.
#' @param buffer_is_priority_domain. If TRUE, replace "is_priority_domain" with a version that is TRUE
#' even if the neighbouring cell is a priority domain. This is an attempt to remove gaps in the image that
#' can occur if we use a spatial_stride != 1 
#' @param asp plot x-y aspect ratio
#' @param fields_axis_args If use_fields=TRUE, then this list is passed to image.plot(..., axis.args=fields_axis_args).
#' It can be used to control the labels on the colourbar.
#' @return Nothing, but make the plot.
multidomain_image<-function(multidomain_dir, variable, time_index, xlim, ylim, zlim, cols, add=FALSE,
    var_transform_function = NULL, NA_if_stage_not_above_elev = FALSE, use_fields=FALSE, clip_to_zlim=FALSE,
    buffer_is_priority_domain=FALSE, asp=1, fields_axis_args=list()){

    library('ncdf4')
    library(fields)
    # Find all netcdf
    all_nc = Sys.glob(paste0(multidomain_dir, '/*/Grid*.nc')) 

    # Start a new plot
    if(use_fields){
        if(!add) image.plot(matrix(0, ncol=2, nrow=2), asp=asp, xlim=xlim, ylim=ylim, zlim=zlim, 
                            col=cols, nlevel=length(cols)+1, axis.args=fields_axis_args)
    }else{
        if(!add) image(matrix(0, ncol=2, nrow=2), asp=asp, col='white', xlim=xlim, ylim=ylim, zlim=zlim)
    }

    # Loop over all domains, and add them to the image
    for(i in 1:length(all_nc)){

        # Open data
        fid = nc_open(all_nc[i], readunlim=FALSE)
      
        # Get x/y coords 
        xs = ncvar_get(fid, 'x')
        ys = ncvar_get(fid, 'y')
      
        # Get the variable of interest 
        if(fid$var[[variable]]$ndim > 2){
            var = ncvar_get(fid, variable, start=c(1,1, time_index), count=c(-1,-1,1))
        }else{
            var = ncvar_get(fid, variable)
        }

        # Find out where the priority domain is
        is_priority = ncvar_get(fid, 'is_priority_domain')

        if(buffer_is_priority_domain){
            # In some situations this can remove 'gaps' in the image that result
            # from using (stride > 1) in the output grids
            nx = length(xs)
            ny = length(ys)
            for(inner_i in (-1):1){
                for(inner_j in (-1):1){
                    is_priority[2:(nx-1), 2:(ny-1)] = pmax(
                        is_priority[2:(nx-1), 2:(ny-1)],
                        is_priority[(2:(nx-1)) + inner_i, (2:(ny-1)) + inner_j])
                }
            }
        }

        # Set areas outside of priority domain to NA
        var[is_priority != 1] = NA

        if(NA_if_stage_not_above_elev){
            # Read stage and elevation, and set var to NA there
            if(variable %in% c('max_stage', 'stage')){
                stage = var
            }else{
                stage = ncvar_get(fid, 'stage', start=c(1,1, time_index), count=c(-1,-1,1))
            }
            if('elev' %in% names(fid$var)){
                elevation = ncvar_get(fid, 'elev', start=c(1,1, time_index), count=c(-1,-1,1))
            }else{
                elevation = ncvar_get(fid, 'elevation0', start=c(1,1), count=c(-1,-1))
            }

            var[stage < elevation + 1.0e-03] = NA
        }
        nc_close(fid)

        if(!is.null(var_transform_function)) var = var_transform_function(var)
        
        #image(xs2, ys2, var, zlim=zlim, col=cols, add=TRUE)

        # Try giving the 'boundary' points to image
        dx = (xs[length(xs)] - xs[1])/(length(xs)-1)
        dy = (ys[length(ys)] - ys[1])/(length(ys)-1)
        #xs2 = c(xs - dx/2, xs[length(xs)] +dx/2)
        #ys2 = c(ys - dy/2, ys[length(ys)] +dy/2)
        xs2 = seq(xs[1] - dx/2, xs[length(xs)] + dx/2, len=length(xs)+1)
        ys2 = seq(ys[1] - dy/2, ys[length(ys)] + dy/2, len=length(ys)+1)

        if(clip_to_zlim){
            var = pmax(var, zlim[1])
            var = pmin(var, zlim[2])
        }

        #if(use_fields){
        #    image.plot(xs2, ys2, var, zlim=zlim, col=cols, add=TRUE, useRaster=TRUE)
        #}else{
            image(xs2, ys2, var, zlim=zlim, col=cols, add=TRUE, useRaster=TRUE)
        #}
    }

}

#
# Determine whether coordinates x,y are on the priority domain of "domain"
#
is_on_priority_domain<-function(x, y, domain){

    # Get indices of xi/yi on the domain (if x/y are not on the domain, this
    # may lead to indices outside the domain)
    x_int = ceiling( (x-domain$lower_left_corner[1])/(domain$lw[1] ) * length(domain$xs))
    y_int = ceiling( (y-domain$lower_left_corner[2])/(domain$lw[1] ) * length(domain$ys))

    is_on = rep(FALSE, length(x_int))
    nx = dim(domain$is_priority_domain) # NOTE: If the netcdf output is decimated, then possibly nx != domain$nx
    for(i in 1:length(is_on)){
        xi = x_int[i]
        yi = y_int[i]
        if( (xi > 0) & (xi <= nx[1]) & (yi > 0) & (yi <= nx[2])){
            # We are on the domain
            # Check if we are in the priority region
            if(domain$is_priority_domain[xi,yi] == 1) is_on[i] = TRUE
        }
    }

    return(is_on)

}


# Given a domain object (i.e. output from 'get_all_recent_results'), find the index and
# distance of the point (grid-cell) nearest to 'x', 'y'
nearest_point_in_domain<-function(x, y, domain){

    kx = which.min(abs(x-domain$xs))
    ky = which.min(abs(y-domain$ys))

    d_x = (domain$xs[kx]-x)
    d_y = (domain$ys[ky]-y)

    out = list(xind=NA, yind=NA, dist=NA)

    if( (kx <= length(domain$is_priority_domain[,1])) &
        (ky <= length(domain$is_priority_domain[1,]))){
        
        if(domain$is_priority_domain[kx,ky] == 1){
            out = list(xind=kx, yind=ky, dist=sqrt(d_x^2 + d_y^2))
        }
    }
    return(out)
}

# Find the x/y range of a multidomain. Here 'md' is a list of domain objects
# (e.g. output of "get_all_recent_results"), which collectively make the
# multidomain
multidomain_range<-function(md){

    all_ranges = lapply(md, f<-function(x) rbind(range(x$xs), range(x$ys)))

    min_x = min(unlist(lapply(all_ranges, f<-function(x) x[1,1])))
    min_y = min(unlist(lapply(all_ranges, f<-function(x) x[2,1])))
    max_x = min(unlist(lapply(all_ranges, f<-function(x) x[1,2])))
    max_y = min(unlist(lapply(all_ranges, f<-function(x) x[2,2])))

    out = rbind(c(min_x, max_x), c(min_y, max_y))
    return(out)
}

# Find the cell in a multidomain that is nearest to 'x,y'. It will return
# the xindex, yindex, distance, and index of the domain in md.
# Here 'md' is a list of domain objects (e.g. output of
# "get_all_recent_results"), which collectively make the multidomain
nearest_point_in_multidomain<-function(x, y, md){

    md_nearest = lapply(md, f<-function(z) nearest_point_in_domain(x, y, z))

    md_nearest_distances = unlist(lapply(md_nearest, f<-function(x) x$dist))

    if(all(is.na(md_nearest_distances))){
        # The point is not in any domain. This can happen e.g. when we ask for
        # points at model boundaries
        closest_domain = -Inf
        out = list(xind=NA, yind=NA, dist=NA)
    }else{
        closest_domain = which.min(md_nearest_distances)
        out = md_nearest[[closest_domain]]
    }
    out$closest_domain=closest_domain
    return(out)
}

# Given a collection of netcdf grid files which partition A SINGLE DOMAIN,
# merge them. This can provide a method for working with
# distributed-memory-parallel multidomains, by first combining each domain into
# a single output.
#
# There are 2 ways to select the netcdf files to merge: 
#    1) Provide 'nc_grid_files', OR 
#    2) Provide both 'multidomain_dir' and 'domain_index'
#
# @param nc_grid_files vector with all nc_grid files that together make up the domain
# @param multidomain_dir directory with multidomain
# @param domain_index index of the domain of interest
# @param desired_var name of variable in the netcdf file
# @param desired_time_index time slice of variable in the netcdf file (or NA for non-time variables)
# @param return_raster If FALSE, return a list with xs, ys, grid. Otherwise return a raster
# @param proj4string projection info for the raster if return_raster=TRUE
#
# ## Example 1 -- provide vector of file names
# all_nc = Sys.glob('RUN_ID000000000*00001_0*/Grid*.nc')
# p1 = merge_domains_nc_grids(nc_grid_files = all_nc)
# 
# ## Example 2 -- nicer -- provide multidomain directory, and domain index
# p1 = merge_domains_nc_grids(multidomain_dir='.', domain_index=3)
#
merge_domains_nc_grids<-function(nc_grid_files = NULL,  multidomain_dir=NA, domain_index = NA,
                                 desired_var = 'max_stage', desired_time_index = NA,
                                 return_raster=FALSE, proj4string="+init=epsg:4326"){
    library('ncdf4')
    library('raster')
    library('sp')
    library('rgdal')

    # Check input makes sense
    if(all(is.null(nc_grid_files)) & (is.na(multidomain_dir) | is.na(domain_index))){
        stop('Must provide EITHER nc_grid_files OR domain_index and multidomain_dir')
    }

    # Find the matching domain files
    if(all(is.null(nc_grid_files))){
        # NOTE: The number of '0' ahead of "domain_index" below protects us
        # against accidently matching other domains, but will eventually fail
        # if we have enough domains. For instance, if we have 10001 domains,
        # then domains '1' and '10001' would both match together, which
        # would be wrong. Some way off however!
        nc_grid_files = Sys.glob(paste0(multidomain_dir, '/RUN_ID0*000', domain_index, 
                                        '_*/Grid*000', domain_index, '.nc'))
    }
    
    # Open all the files 
    fids = sapply(nc_grid_files, f<-function(x){nc_open(x, readunlim=FALSE)}, simplify=FALSE)

    # Get the 'x' and 'y' and 'time' dimension variables
    xs = lapply(fids, f<-function(x) ncvar_get(x, 'x'))
    ys = lapply(fids, f<-function(x) ncvar_get(x, 'y'))
    ts = lapply(fids, f<-function(x) ncvar_get(x, 'time'))

    # Check that times are compatible in all files
    if(length(ts) > 1){
        # We will accept "numerically negligable" differences in times
        dts = c(0, diff(ts[[1]]))
        for(i in 2:length(ts)){
            if(!all(abs(ts[[i]] - ts[[1]]) <= dts/1000)){
                print(paste0('Times in file ', i, ': ', nc_grid_files[i],
                            ' are incompatible with times in file 1: ', nc_grid_files[1]))
                print(cbind(ts[[i]], ts[[1]]))
                stop('Times are incompatible')
            }
        }
    }
  
    # Check that the x/y spacings are compatible  

    dxs = unlist(lapply(xs, f<-function(x) mean(diff(x))))
    dys = unlist(lapply(ys, f<-function(x) mean(diff(x))))

    mean_dx = mean(dxs)
    mean_dy = mean(dys)

    dx_error_tol = 1.0e-03
    if(any(abs(dxs - mean_dx) > dx_error_tol*mean_dx)){
        stop('dx differences exceed tolerance')
    }
    if(any(abs(dys - mean_dy) > dx_error_tol*mean_dy)){
        stop('dy differences exceed tolerance')
    }

    # Make a 'full grid' xs and ys 

    x_range = range(unlist(lapply(xs, range)))
    y_range = range(unlist(lapply(ys, range)))

    nx = round(diff(x_range)/mean_dx + 1 )
    full_grid_xs = seq(x_range[1], x_range[2], length=nx)
    ny = round(diff(y_range)/mean_dy + 1 )
    full_grid_ys = seq(y_range[1], y_range[2], length=ny)

    # We also need 'full grid' times

    full_grid_times = ts[[1]]

    # The output variable
    output_full = matrix(NA, nrow=nx, ncol=ny)

    for(i in 1:length(fids)){
        ipd = ncvar_get(fids[[i]], 'is_priority_domain')
        if(is.na(desired_time_index)){
            local_var = ncvar_get(fids[[i]], desired_var)
        }else{
            local_var = ncvar_get(fids[[i]], desired_var, start=c(1, 1, desired_time_index), 
                                  count=c(-1, -1, 1))
        }
        # Set 'non-priority-domain' regions to NA
        if(any(ipd == 0)){
            local_var[ipd == 0] = NA
        }

        xi = which.min(abs(full_grid_xs - xs[[i]][1]))
        yi = which.min(abs(full_grid_ys - ys[[i]][1]))
        # Set only regions where ipd == 1
        output_full[xi:(xi+dim(local_var)[1]-1), yi:(yi+dim(local_var)[2]-1)][ipd == 1] = local_var[ipd == 1]
    }

    # Close the netcdf files
    lapply(fids, nc_close)

    if(!return_raster){
        # Return the grid, with x/y dimensions (e.g. for passing to 'image')
        output = list(xs = full_grid_xs, ys = full_grid_ys)
        output[[desired_var]] = output_full
        return(output)
    }else{
        # Return a raster
        nx = length(full_grid_xs)
        ny = length(full_grid_ys)

        dx = mean_dx
        dy = mean_dy

        xmn = full_grid_xs[1] - dx/2
        xmx = full_grid_xs[length(full_grid_xs)] + dx/2
        ymn = full_grid_ys[1] - dy/2
        ymx = full_grid_ys[length(full_grid_ys)] + dy/2

        r1 = raster(t(output_full), xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx,
                    crs=CRS(proj4string))
        r1 = flip(r1, direction='y')
        return(r1)
    }
}

#' Merge gauges from a multidomain, discarding those that are not in priority domains
#'
#' @param md list with the multidomain info.
#' @param multidomain_dir the directory containing all the multdomain outputs. 
#' If this is provided then md should be NA, and the code will read the necessary information
#' @param assume_all_gauges_are_priority If TRUE this assumes 
#'
merge_multidomain_gauges<-function(md = NA, multidomain_dir=NA, assume_all_gauges_are_priority=TRUE){
   
    # For all gauges, figure out if they are in the priority domain
    clean_gauges = list()
    counter = 0

    if(!is.na(multidomain_dir)){
        # In this case we passed "multidomain_dir" to the function
        #
        # Need to read the multidomain information

        if(!( (length(md) == 1))){
            stop('Cannot provide both multidomain_dir and md')
        }
        if(!is.na(md)){
            stop('Cannot provide both multidomain_dir and a non-NA value of md')
        }

        # Redefine md
        md_domains = Sys.glob(paste0(multidomain_dir, '/RUN_ID*'))
        if(length(md_domains) == 0) stop(paste0('Could not find any domains in multidomain_dir ', multidomain_dir))

        # Read the multidomain info
        # We can avoid most variables except the gauges
        md = lapply(md_domains, 
                    f<-function(x) get_all_recent_results(x, read_grids=FALSE, always_read_priority_domain=FALSE,
                                                          read_time_and_geometry=(!assume_all_gauges_are_priority)) )

        # Check if the gauges are all in their priority domain. If so, we can
        # avoid an expensive computation to derive this info
        all_domain_gauges_are_priority <-function(x){
            out = FALSE
            if(length(x) <= 1){
                # No md output -- something is wrong
                stop('Error in merge_multidomain_gauges: There is an empty NULL domain in the md list')
            }else{
                if( (length(x$gauges) == 0 ) | (is(x$gauges, 'try-error'))){
                    # No gauges -- no problem
                    out = TRUE
                }else{
                    # This will be FALSE for older files
                    if(x$gauges$priority_gauges_only) out = TRUE
                }
            }
            return(out)
        }
        all_gauges_are_priority = all(unlist(lapply(md, all_domain_gauges_are_priority)))

        if(!all_gauges_are_priority){
            # In this case, we will need to use the is_priority_domain data to
            # search the gauges.
            # This is necessary for older SWALS gauge netcdf files which did
            # not record whether all gauges were priority (and in general they
            # were not).
            for(j in 1:length(md)){
                md[[j]]$is_priority_domain = try(
                    get_gridded_variable(var='is_priority_domain', output_folder=md[[j]]$output_folder), 
                    silent=TRUE)
            }
        }

    }else{
        # In this case we passed an existing md list to the function
        #
        if(!is(md, 'list')){
            stop('Must provide named arguments with either "multidomain_dir" giving the multidomain directory, or a list "md" containing the output from get_all_recent_results for each domain in the multidomain')
        }

        # Useful to have these variables
        md_domains = lapply(md, f<-function(x) x$output_folder)
        multidomain_dir = dirname(md[[1]]$output_folder)
    }

    for(j in 1:length(md)){

        # Some domains might have no gauges -- skip them
        if( (length(md[[j]]$gauges) == 0) ) next
        if( (is(md[[j]]$gauges, 'try-error')) ) next

        # If gauges exist, they may or may not be in this priority domain
        lons = md[[j]]$gauges$lon
        lats = md[[j]]$gauges$lat

        if(md[[j]]$gauges$priority_gauges_only){
            keep = 1:length(lons)
        }else{
            # Search for the priority gauges
            # This can be slow, but for older files it was necessary
            nearest_d = sapply(1:length(lons), f<-function(x){
                nearest_point_in_multidomain(lons[x], lats[x], md)$closest_domain})
            keep = which(nearest_d == j)
        }
    
        if(length(keep) == 0) next
        # At this point, definitely we have some gauges of interest

        counter = counter + 1

        clean_gauges[[counter]] = md[[j]]$gauges

        # Get the "keep" subset of variables, corresponding to gauges in priority domain
        # Loop over 1-d variables lon/lat/ etc, and subset in space
        for(i in 1:length(clean_gauges[[counter]])){
            # Only work on arrays
            if(!is(clean_gauges[[counter]][[i]], 'array')) next
            # Do not subset time! 
            if(names(clean_gauges[[counter]])[i] == 'time') next
            clean_gauges[[counter]][[i]] = as.numeric(clean_gauges[[counter]][[i]][keep])
        }
        # Subset the time variables
        for(i in 1:length(clean_gauges[[counter]]$time_var)){

            # Non arrays are e.g. missing data -- ignore.
            if(!any(class(clean_gauges[[counter]]$time_var[[i]]) %in% c('matrix', 'array'))) next

            dd = dim(clean_gauges[[counter]]$time_var[[i]])
            if(is.null(dd) | length(dd) == 1){
                # There is only one gauge in the whole thing. Since length(keep) > 0, this must be it
                dim(clean_gauges[[counter]]$time_var[[i]]) = c(1, length(clean_gauges[[counter]]$time_var[[i]])) 
            }else{
                # There is more than one gauge. If we have
                # priority_gauges_only, then keep includes all rows and we can
                # do nothing
                if(!clean_gauges[[counter]]$priority_gauges_only){
                    clean_gauges[[counter]]$time_var[[i]] = clean_gauges[[counter]]$time_var[[i]][keep,,drop=FALSE] 
                }
            }
        }
        # Static variables
        for(i in 1:length(clean_gauges[[counter]]$static_var)){

            # Non arrays are e.g. missing data -- ignore?
            if(!any(class(clean_gauges[[counter]]$static_var[[i]]) %in% c('matrix', 'array'))) next
            # There is more than one gauge
            clean_gauges[[counter]]$static_var[[i]] = as.numeric(clean_gauges[[counter]]$static_var[[i]][keep])

        }

    }

    # Now, merge the "clean gauges" into a single data structure
    merged_gauges = clean_gauges[[1]]

    if(length(clean_gauges) > 1){

        for(counter in 2:length(clean_gauges)){

            # Combine single variables like lon/lat,gaugeID... but not time!
            for(var in 1:length(merged_gauges)){

                if(is(merged_gauges[[var]], 'list')) next

                if(names(clean_gauges[[counter]])[var] == 'time'){
                    # Treat time separately. But test that all times are the same
                    if(!all(clean_gauges[[counter]]$time == merged_gauges$time)){
                        stop('GAUGE TIMES DIFFER IN DIFFERENT DOMAINS')
                    }
                }else{
                    # Merge the variables
                    merged_gauges[[var]] = c(merged_gauges[[var]], clean_gauges[[counter]][[var]])
                }
            }
        }

        # Time variables are 2D -- so rbind should merge properly
        # Should be faster to use do.call than repeated rbind's, due to memory allocation
        for(var in 1:length(merged_gauges$time_var)){
            tmp_list = vector(mode='list', length=length(clean_gauges))
            for(counter in 1:length(clean_gauges)){
                stopifnot(ncol(merged_gauges$time_var[[var]]) == ncol( clean_gauges[[counter]]$time_var[[var]] ))
                tmp_list[[counter]] = clean_gauges[[counter]]$time_var[[var]]
            }
            merged_gauges$time_var[[var]] = do.call(rbind, tmp_list)
        }
        rm(tmp_list)

        for(var in 1:length(merged_gauges$static_var)){
            tmp_list = vector(mode='list', length=length(clean_gauges))
            for(counter in 1:length(clean_gauges)){
                tmp_list[[counter]] =  clean_gauges[[counter]]$static_var[[var]]
            }
            merged_gauges$static_var[[var]] = do.call(c, tmp_list)
        }
        rm(tmp_list)
    }

    # Bit of cleaning up. For missing variables, replace possibly multiple NAs with a single NA
    for(i in 1:length(merged_gauges$static_var)){
        if(all(is.na(merged_gauges$static_var[[i]]))) merged_gauges$static_var[[i]] = NA
    }
    for(i in 1:length(merged_gauges$time_var)){
        if(all(is.na(merged_gauges$time_var[[i]]))) merged_gauges$time_var[[i]] = NA
    }

    merged_gauges$multidomain_dir = multidomain_dir

    return(merged_gauges)
}


# Partition a set into 'k' groups with roughly equal sum
# This uses the naive algorithm.
partition_into_k<-function(vals, k){

    sums = rep(0, k)
    inds = vector(mode='list', length=k)

    sorted_vals = sort(vals, index.return=TRUE, decreasing=TRUE)
    for(i in 1:length(vals)){
        vi = sorted_vals$x[i]
        add_to = which.min(sums)
        inds[[add_to]] = c(inds[[add_to]], sorted_vals$ix[i])
        sums[add_to] = sums[add_to] + vi
    }
    return(list(inds, sums))
}

# A workhorse routine for partitioning with groups. See documentation below on
# the interface routine.
partition_into_k_with_grouping_WORK<-function(vals, k, vals_groups, random_ties=FALSE){

    cumulative_sum_vals_in_group = rep(0, k)
    inds_in_group = vector(mode='list', length=k)
    sorted_vals = sort(vals, index.return=TRUE, decreasing=TRUE)

    # Loop over each of the separate vals_groups, and partition that group
    # alone
    unique_vals_groups = unique(vals_groups)
    tmp = vector(mode='list', length = length(unique_vals_groups))
    counter = 0
    for(ug in unique_vals_groups){
        # Get the ones in this group
        p = which(vals_groups[sorted_vals$ix] == ug)

        # Group indices
        inds_p = sorted_vals$ix[p]
        # Group values
        vals_p = sorted_vals$x[p]
        # Split group values into k
        local_split = partition_into_k(vals_p, k)
        local_split[[1]] = lapply(local_split[[1]], f<-function(x) inds_p[x])
        counter = counter+1
        tmp[[counter]] = local_split
    }

    # Now loop over each group, ordered according to which has the largest value
    if(random_ties){
        between_group_order = rank(unlist(lapply(tmp, f<-function(x) max(x[[2]]))), ties='random')
    }else{
        between_group_order = rank(unlist(lapply(tmp, f<-function(x) max(x[[2]]))), ties='first')
    }
    # Flip
    between_group_order = max(between_group_order) + 1 - between_group_order

    for(g in between_group_order){
        # Within the groups, the 'vals' are already ordered into k groups, longest-time first
        within_group_inds = tmp[[g]][[1]]
        within_group_vals = tmp[[g]][[2]]
        # Here "ties='first'" would cause the test to fail
        # Experiments suggest it isn't useful to add extra randomness
        assign_to_order = rank(cumulative_sum_vals_in_group, ties='last')  
        for(i in 1:length(assign_to_order)){
            partition_ind = assign_to_order[i]
            inds_in_group[[partition_ind]] = c(inds_in_group[[partition_ind]], within_group_inds[[i]])
            cumulative_sum_vals_in_group[partition_ind] = cumulative_sum_vals_in_group[partition_ind] + 
                within_group_vals[i]
        }
    }
    
    return(list(inds_in_group, cumulative_sum_vals_in_group))
}



# Partition a set of values 'vals' into 'k' separate groups with roughly equal
# sum. Furthermore we assume the 'vals' are each associated with some
# 'vals_group', and the partitioning should make the sums rougly equal also
# within the vals_groups. 
#
# To make it concrete: suppose our 'vals' give model run-times for a
# partitioned multi-domain.  The 'vals_groups' might be chosen to give the
# "domain_index" of each, in which case all of the 'k' groups would have
# approximately equal values of 'vals' WITHIN each domain_index. 
#
# A practical situation where we might want to do this is for a model
# where some domains are only evolved for a fraction of the whole runtime,
# and we want to achieve good load-balancing for the entire model run.
#
# This basic method uses a naive partitioning algorithm. In some cases better
# performance can be obtained using random tie-breaking. The general interface
# (below) does this.
#
partition_into_k_with_grouping<-function(vals, k, vals_groups, ntries=1){

    # Deterministic case first -- so we never do worse than deterministic   
    best_try = partition_into_k_with_grouping_WORK(vals, k, vals_groups)
    if(ntries > 1){
        # Try to do better than deterministic case -- in some cases it works well
        for(i in 2:ntries){
            test = partition_into_k_with_grouping_WORK(vals, k, vals_groups, random_ties=TRUE)
            if(max(test[[2]]) < max(best_try[[2]])) best_try = test
        }
    }
    return(best_try)
}

# Basic test of partitioning
test_partition_into_k_with_grouping<-function(){

    # Example -- "nd" domains, each split into 32, each taking "time" ranging from 1:32.
    # 
    # By inspection, for this particular case we should be able to split the domains equally,
    # except for one domain which will have an extra unit of work

    nd = 8
    nsplit = 32
    vals = c()
    for(i in 1:nd) vals = c(vals, sample(1:nsplit, size=nsplit, replace=FALSE))
    group_inds = rep(1:nd, each=nsplit)

    # Test 1 -- this should partition perfectly
    test = partition_into_k_with_grouping(vals, k=nsplit, vals_groups=group_inds)
    max_val = max(test[[2]])
    if(sum(test[[2]] == max_val) == nsplit){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Test 2 -- break the 'equal partition' balance
    # By inspection, this should lead to one group that has one more element
    # than the others (when nd is even)
    vals[45] = vals[45] + 1 
    test = partition_into_k_with_grouping(vals, k=nsplit, vals_groups=group_inds)
    max_val = max(test[[2]])
    # There should be one value with max_val, and all others being max_val - 1
    if(sum(test[[2]] == (max_val-1)) == (nsplit-1)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Test 3 -- break it even more
    # Ideally we would have 2 values that differ from 1 by their optimal value.
    # But this naive algorithm does not achieve that.
    vals[145] = vals[145] + 1 
    testA = partition_into_k_with_grouping(vals, k=nsplit, vals_groups=group_inds)
    test = partition_into_k_with_grouping(vals, k=nsplit, vals_groups=group_inds, ntries=1000)
    # At least we check that the random tries led to improvement, which it seems to
    if(max(test[[2]]) < max(testA[[2]])){
        print('PASS')
    }else{
        print('FAIL')
    }
}


#' Make a load_balance_partition.txt file
#'
#' Suppose we have results from an existing multidomain model run, compiled with -DTIMER
#' so that the output files contain timing information for each domain. This function
#  can use those results to figure out a more balanced partitioning of the domains among images.
#' It makes a load balance partition file, using a 'greedy' approach to distribute
#' domains to images. 
#'
#' @param multidomain_dir directory containing the multidomain outputs
#' @param verbose if FALSE, suppress printing
#' @param domain_index_groups either an empty list, or a list with 2 or more vectors, each defining a group of domain indices.
#' In the latter case, the partition tries to be approximately equal WITHIN each group first, and then to combine the results
#' in a good way. This will generally be less efficient than not using groups, unless you have some other information
#' that tells you that the partition should be done in this way.
#' @return the function environment invisibly (so you have to use assignment to capture it)
#'
make_load_balance_partition<-function(multidomain_dir=NA, verbose=TRUE, domain_index_groups = list()){

    if(is.na(multidomain_dir)){
        multidomain_dir = rev(sort(Sys.glob('./OUTPUTS/RUN*')))[1]
    }

    if(verbose){
        cat('############ \n')
        cat(paste0('Generating load_balance_partition.txt file in ', multidomain_dir, '\n'))
        cat('    Future jobs which use this file should have the same number of images and threads \n')
        cat('\n')
    }

    md_files = Sys.glob(paste0(multidomain_dir, '/multi*.log'))

    # Get the runtime from the log file. 
    # This depends on the timer information having been printed
    md_runtime<-function(md_file){
        md_lines = readLines(md_file)

        domain_timer_starts = grep('Timer of md%domains(', md_lines, fixed=TRUE)
        total_wallclock_lines = grep('Total WALLCLOCK', md_lines)

        # Get the wallclock time as numeric
        nd = length(domain_timer_starts)
        total_wallclock = sapply(md_lines[total_wallclock_lines[1:nd]], 
            f<-function(x) as.numeric(strsplit(x, ':', fixed=TRUE)[[1]][2]))

        # Get the domain index and local index
        domain_index_and_local_index = sapply(md_lines[domain_timer_starts + 2], 
            f<-function(x) as.numeric(read.table(text=x)),
            simplify=FALSE)
        names(domain_index_and_local_index) = rep("", length(domain_index_and_local_index))
        domain_index_and_local_index = matrix(unlist(domain_index_and_local_index), ncol=2, byrow=TRUE)
    
        output = list(total_wallclock = as.numeric(total_wallclock), 
                      domain_index_and_local_index = domain_index_and_local_index) 
        #return(as.numeric(total_wallclock))
        return(output)
    }
    md_timer_data = sapply(md_files, md_runtime, simplify=FALSE)

    md_times = lapply(md_timer_data, f<-function(x) x$total_wallclock)
    md_domain_indices = lapply(md_timer_data, f<-function(x) x$domain_index_and_local_index[,1])
    md_local_indices = lapply(md_timer_data, f<-function(x) x$domain_index_and_local_index[,2])

    num_images = length(md_times)

    md_times_vec = unlist(md_times)
    md_domain_indices_vec = unlist(md_domain_indices)
    md_local_indices_vec = unlist(md_local_indices)

    if(length(domain_index_groups) == 0){
        # Make a partition of md_times_vec into num_images groups, with roughly equal sums
        new_times_and_grouping = partition_into_k(md_times_vec, num_images)
    }else{
        # Make a partition of md_times_vec into num_images groups, with roughly equal sums.
        # Further, ensure the sums are also rougly equal for images within the
        # same domain_index_group. 
        stopifnot(is.list(domain_index_groups))
        # Ensure there are no repeated domain indices
        tmp = unlist(domain_index_groups)
        stopifnot(length(tmp) == length(unique(tmp)))

        groups = md_domain_indices_vec * NA
        for(j in 1:length(domain_index_groups)){
            k = which(md_domain_indices_vec %in% domain_index_groups[[j]])
            groups[k] = j
        }

        new_times_and_grouping = partition_into_k_with_grouping(md_times_vec, num_images, groups)
    }

    range_new = range(new_times_and_grouping[[2]])
    dsplit = diff(range_new)
    if(verbose){
        cat(paste0('Range of partition total times: ', signif(range_new[1], 4), '-to-', signif(range_new[2], 4),  's \n'))
        cat(paste0('                  (difference): ', signif(dsplit, 4), 's \n'))
        cat(paste0('             (as a percentage): ', signif(dsplit/mean(new_times_and_grouping[[2]])*100, 4), '%\n'))
    }
    range_old = range(unlist(lapply(md_times, sum)))
    old_range = diff(range_old)
    if(verbose){
        cat(paste0('Previous time range: ', signif(range_old[1], 4), '-to-', signif(range_old[2], 4),  's \n'))
        cat(paste0('      (difference) : ', signif(old_range, 4), '\n'))
        cat(paste0('Potential run-time reduction (difference in max times):', signif(range_old[2] - range_new[2], 4), 's\n'))
        cat(paste0('      Potential run-time reduction (% of old max time):', 
                     signif((range_old[2] - range_new[2])/range_old[2] * 100, 4), '%\n'))
    }
    # 
    unique_domains = sort(unique(md_domain_indices_vec))
    if(! all(unique_domains == seq(1, max(unique_domains)))){
        stop('unique_domains is not a sequence of the form from 1, 2, ... max')
    }
    # Make the data for the load balance file
    load_balance_data = vector(mode='list', length=max(unique_domains))
    for(i in unique_domains){
        # Find local indices on this domain
        k = which(md_domain_indices_vec == i)
        local_inds = md_local_indices_vec[k]
        # Find which image gets this local index
        for(j in 1:max(local_inds)){
            n = which(local_inds == j)
            kn = k[n]
            # kn will only be in one entry of new_times_and_grouping[[1]] (the list of
            # indices on each images). Find it 
            target_domain = which(unlist(lapply(new_times_and_grouping[[1]], f<-function(x) kn%in%x)))
            load_balance_data[[i]] = c(load_balance_data[[i]], target_domain)
        }
    }

    output_file = paste0(multidomain_dir, '/load_balance_partition.txt')
    for(i in 1:length(load_balance_data)){
        if(i == 1){
            cat(load_balance_data[[i]], '\n', file=output_file)
        }else{
            cat(load_balance_data[[i]], '\n', file=output_file, append=TRUE)
        }
    }

    old_time_range = unlist(lapply(md_times, sum))
    new_time_range = new_times_and_grouping[[2]]
    if(verbose){
        cat('\n')
        cat('OLD MODEL DOMAIN TIMES (observed): Stem and leaf plot \n')
        stem(old_time_range)
        cat('\n')
        cat('LOAD BALANCED DOMAIN TIMES (prediction): Stem and leaf plot\n')
        stem(new_time_range)
    }
    return(invisible(environment()))
}


# Extract the domain_index from its folder name
#
# Domain folders from multidomains have names like here:
#    RUN_ID00000000590000000002_00070_20190529_094228.898
# Where the ID00... has first 10 digits being the image index (59 here),
# and the next 10 digits being the domain_index (2 here).
#
#
#stopifnot(domain_index_from_folder('RUN_ID00000000590000000002_00070_20190529_094228.898') == 2)
# 
domain_index_from_folder<-function(folder_name){
    
    ID_term = gsub("ID", "", as.character(strsplit(folder_name, '_')[[1]][2]))

    nch = nchar(ID_term)

    index_ID = substring(ID_term, nch-9, nch)

    return(as.numeric(index_ID))
}

# Extract the domain_image from its folder name
#
# Domain folders from multidomains have names like here:
#    RUN_ID00000000590000000002_00070_20190529_094228.898
# Where the ID00... has first 10 digits being the image index (59 here),
# and the next 10 digits being the domain_index (2 here).
#
#
#stopifnot(domain_image_from_folder('RUN_ID00000000590000000002_00070_20190529_094228.898') == 59)
# 
domain_image_from_folder<-function(folder_name){

    ID_term = gsub("ID", "", as.character(strsplit(folder_name, '_')[[1]][2]))

    nch = nchar(ID_term)

    image_ID = substring(ID_term, nch-19, nch-10)

    return(as.numeric(image_ID))
}

# Get the indices of all domains in a multidomain, using the folder name
#
# This can be used to avoid hard-coding the number of domains in scripts.
#
get_domain_indices_in_multidomain<-function(multidomain_dir){
    all_domain_run_folders = Sys.glob(paste0(multidomain_dir, '/RUN_ID*'))
    all_domain_indices = sapply(basename(all_domain_run_folders), domain_index_from_folder)
    return(unique(all_domain_indices))
}

# Get the interior bounding box (i.e. the domain box prior to halo padding) for
# each domain in a multidomain. 
#
# It is reported in both a 'merged form' (i.e. ignoring any parallel
# partitioning), and in a parallel partition
#
# @param multidomain_dir the directory with the multidomain
get_domain_interior_bbox_in_multidomain<-function(multidomain_dir){
    library(ncdf4)

    all_domain_run_folders = Sys.glob(paste0(multidomain_dir, '/RUN_ID*')) 
    if(length(all_domain_run_folders) == 0) stop('No domains found')

    all_domain_dx = lapply(all_domain_run_folders, get_dx)

    # Read the interior bounding box and dx/dy
    all_domain_ibb = vector(mode='list', length=length(all_domain_run_folders))
    for(i in 1:length(all_domain_run_folders)){
        nc_file = Sys.glob(paste0(all_domain_run_folders[i], '/Grid*.nc'))
        fid = nc_open(nc_file, readunlim=FALSE)
        all_att = ncatt_get(fid, varid=0)
        ibb = as.numeric(scan(text=(all_att$interior_bounding_box)[[1]], quiet=TRUE))
        nc_close(fid)
        dim(ibb) = c(4,2) 
        all_domain_ibb[[i]] = ibb
    }

    # Unsplit domain indices
    all_domain_indices = sapply(basename(all_domain_run_folders), domain_index_from_folder)
    all_domain_images = sapply(basename(all_domain_run_folders), domain_image_from_folder)

    # Unsplit domain bounding box
    split_domain_ibb = do.call(rbind, lapply(all_domain_ibb, c))
    lower_x = aggregate(split_domain_ibb[,1], by=list(all_domain_indices), min)
    upper_x = aggregate(split_domain_ibb[,2], by=list(all_domain_indices), max)
    lower_y = aggregate(split_domain_ibb[,5], by=list(all_domain_indices), min)
    upper_y = aggregate(split_domain_ibb[,7], by=list(all_domain_indices), max)
    unsplit_domain_ibb = cbind(lower_x$x, upper_x$x, upper_x$x, lower_x$x, 
                               lower_y$x, lower_y$x, upper_y$x, upper_y$x)
    merged_domain_indices = lower_x[,1]

    # Domain dx
    split_domain_dx = do.call(rbind, all_domain_dx)
    merged_domain_dx_1 = aggregate(split_domain_dx[,1], by=list(all_domain_indices), mean)
    merged_domain_dx_2 = aggregate(split_domain_dx[,1], by=list(all_domain_indices), mean)
    merged_domain_dx = cbind(merged_domain_dx_1$x, merged_domain_dx_2$x)

    # Useful to know the extent of the multidomain
    multidomain_xlim = c(min(lower_x$x), max(upper_x$x))
    multidomain_ylim = c(min(lower_y$x), max(upper_y$x))

    merged_domain_ibb = vector(mode='list', length=nrow(unsplit_domain_ibb))
    for(i in 1:length(merged_domain_ibb)){
        merged_domain_ibb[[i]] = matrix(unsplit_domain_ibb[i,], ncol=2)
    }
    names(merged_domain_ibb) = as.character(merged_domain_indices)

    outputs = list(domain_folders = all_domain_run_folders, 
                   domain_interior_bbox = all_domain_ibb,
                   domain_indices = all_domain_indices,
                   domain_images = all_domain_images,
                   domain_dx = all_domain_dx,
                   multidomain_xlim = multidomain_xlim,
                   multidomain_ylim = multidomain_ylim,
                   merged_domain_indices = merged_domain_indices,
                   merged_domain_interior_bbox = merged_domain_ibb,
                   merged_domain_dx = merged_domain_dx)
}


#'
#' Convenience function to get the timer-total-wallclock times for each
#' individual domain, as recorded in the md log file
#'
#' @param md_log_file filename
#' @param wallclock_time_line_spacing Number of lines after the RUN_ID00 filename that the wallclock-time line appears
#' @param a data.frame containing the domain ID info, the wallclock time, and the domain index according to the model setup
get_domain_wallclock_times_in_log<-function(md_log_file, wallclock_time_line_spacing = 15){
    x = readLines(md_log_file)
    # The timer information is preceeded by statement of the output directory, which will match this
    inds = grep('RUN_ID00', x)
    # Get the domain information
    domains = basename(x[inds])
    domains = unlist(lapply(strsplit(domains, '_'), f<-function(x) x[2]))
    # Get the time
    times = as.numeric(gsub('Total WALLCLOCK time: ', '', x[inds+wallclock_time_line_spacing]))
    # Get the domain number corresponding to the original
    domain_number = as.numeric(gsub('ID', '', domains))%%100000

    output = data.frame(domains=domains, times=times, index=domain_number)
    return(output)
}

#' Find the domain thas has a given point in its priority domain, 
#' while minimising the amount/time reading data
#'
#' @param xy_mat a 2-colum matrix with point coordinates
#' @param md result of get_multidomain or similar (i.e. list of domain objects)
#' @param multidomain_dir name of the directory containing the multidomain
#' 
find_domain_containing_point<-function(xy_mat, md=NULL, multidomain_dir=NULL){

    if(!is.null(multidomain_dir)){
        if(length(md) != 0) stop("cannot provide both multidomain_dir and md")

        # Read the multidomain, getting only the x/y coordinates
        md = get_multidomain(multidomain_dir, read_grids=FALSE, read_gauges=FALSE, 
            read_time_and_geometry=FALSE, always_read_xy=TRUE, quiet=TRUE)
    }
    
    # "xy_mat" should be a matrix -- but if it is a vector of length(2), then assume
    # this gives the xy coordinates of a single point, and coerce to a matrix
    if(length(dim(xy_mat)) == 0){
        if(length(xy_mat) == 2){
            dim(xy_mat) = c(1, 2)
        }
    }

    # Find dx/dy for each domain, and take the product as the 'cell area', with finer
    # cell-areas having higher priority domains
    md_cell_dx = unlist(lapply(md, f<-function(x) x$xs[2] - x$xs[1])) 
    md_cell_dy = unlist(lapply(md, f<-function(x) x$ys[2] - x$ys[1]))
    md_cell_areas = md_cell_dx * md_cell_dy 
    
    # Get the "interior bounding box" of each domain    
    md_bboxes = vector(mode='list', length=length(md))
    for(i in 1:length(md)){
        nc_file = Sys.glob(paste0(md[[i]]$output_folder, '/Grid_output_*.nc'))
        tmp_text = get_nc_global_attribute(nc_file, 'interior_bounding_box')
        md_bboxes[[i]] = matrix(as.numeric(scan(text=tmp_text$value, quiet=TRUE)), ncol=2)
    }


    md_cell_xmin = unlist(lapply(md_bboxes, f<-function(x) min(x[,1])))
    md_cell_xmax = unlist(lapply(md_bboxes, f<-function(x) max(x[,1])))
    md_cell_ymin = unlist(lapply(md_bboxes, f<-function(x) min(x[,2])))
    md_cell_ymax = unlist(lapply(md_bboxes, f<-function(x) max(x[,2])))

    # For each point, find which domain(s) it is in
    point_domain = rep(NA, nrow(xy_mat))
    point_dir = rep(NA, nrow(xy_mat))
    for(i in 1:nrow(xy_mat)){
        xy = xy_mat[i,1:2]
        point_in_domain = (xy[1] >= md_cell_xmin) & (xy[1] <= md_cell_xmax) & 
                          (xy[2] >= md_cell_ymin) & (xy[2] <= md_cell_ymax)

        if(!any(point_in_domain)){
            next
        }else{
            # Select the domain with the smallest cell area [which in SWALS is
            # always the priority domain].
            local_areas = md_cell_areas
            local_areas[!point_in_domain] = Inf
            point_domain[i] = which.min(local_areas)
            point_dir[i] = md[[point_domain[i]]]$output_folder
        }
    }

    output = list(points = xy_mat, domain_index=point_domain, domain_dir=point_dir)

    return(output)

}

#' Energy-type computations for leapfrog-type solvers that have "domain%linear_solver_is_truely_linear=.TRUE."
#'
#' FIXME: Update this to treat cases other than the linear shallow water equations. However there might not
#' be any point now that the logfiles output this information.
#'
#' @param domain_dir the domain directory
#' @param spherical Is the domain spherical?
#' @param radius_earth radius of the earth in m
#' @param gravity gravitational acceleration (m/s^2)
#' @param msl the mean-sea-level for the truely linear solver
#' @param dry_depth the depth at which velocities are set to zero
#' @param water_density density of (sea) water in kg/m^3
#' @param time_indices optional vector of time indices at which energy should be computed. By default compute all times.
#' @return List with the 'available energy' - including:
#'             potential energy (1/2 integral (h-msl)**2)
#'             kinetic energy   (1/2 integral (u^2 + v^2)*depth) 
#'             total energy (potential + kinetic)
#'        We use the 'available energy' (i.e. so a still ocean with stage=msl has energy = 0)
#
get_energy_truely_linear_domain<-function(domain_dir, spherical=TRUE, radius_earth=6371e+03, 
    gravity=9.8, msl=0.0, dry_depth=1.0e-05, water_density = 1024, time_indices=NULL){

    deg2rad = pi/180

    library(ncdf4)
    target_file = Sys.glob(paste0(domain_dir, "/Grid_output_ID*.nc"))
    if(length(target_file) != 1) stop(paste0('Did not find exactly one Grid_output_ID*.nc file in ', domain_dir))
    fid = nc_open(target_file)
    elev = ncvar_get(fid, 'elevation0')
    ys = ncvar_get(fid, 'y')
    xs = ncvar_get(fid, 'x')
    times = ncvar_get(fid, 'time')
    is_priority_domain = ncvar_get(fid, 'is_priority_domain')

    # If the priority domain is empty it seems filled with negative 127
    if(any(is_priority_domain < 0)){
        is_priority_domain = 1
    }

    # Depth below MSL -- avoiding singularity
    EPS = dry_depth
    depth = pmax(msl-elev, EPS)
    nx = dim(depth)
    is_wet = (elev < (msl-EPS))

    # Remove boundaries -- 2 cells.
    is_wet[1:2,] = 0
    is_wet[(nx[1]-1):nx[1],] = 0
    is_wet[,1:2] = 0
    is_wet[,(nx[2]-1):nx[2]] = 0

    # Get the grid cell area, assuming spherical coordinates
    cell_area = elev*0
    dy = mean(diff(ys))
    dx = mean(diff(xs))
    for(i in 1:ncol(cell_area)){
        if(spherical){
            # Spherical area
            cell_area[,i] = radius_earth**2 * cos(ys[i]*deg2rad) * dy * deg2rad * dx * deg2rad
        }else{
            cell_area[,i] = dx*dy
        }
    }

    # Cell area at i+1/2
    cell_area_ip = cell_area
    # Cell area at j+1/2
    cell_area_jp = cell_area
    cell_area_jp[, 1:(nx[2]-1)] = 0.5 * (cell_area[,2:nx[2]] + cell_area[,1:(nx[2]-1)])

    #
    depth_ip = depth
    depth_ip[1:(nx[1]-1),] = 0.5 * (depth[2:nx[1],] + depth[1:(nx[1]-1),])
    depth_jp = depth
    depth_jp[, 1:(nx[2]-1)] = 0.5 * (depth[,2:nx[2]] + depth[,1:(nx[2]-1)])

    # Below we'll either get results at a specific time, or all times.
    if(is.null(time_indices)){
        time_inds = 1:length(times)
    }else{
        time_inds = time_indices
        times = times[time_indices]
    }

    # Potential and kinetic energy
    energy_pot_on_rho = rep(NA, length(time_inds))
    energy_kin_on_rho = rep(NA, length(time_inds))
    for(i in time_inds){
        stg = ncvar_get(fid, 'stage', start=c(1, 1, i), count=c(nx[1], nx[2], 1))
        uh = ncvar_get(fid, 'uh', start=c(1, 1, i), count=c(nx[1], nx[2], 1))
        vh = ncvar_get(fid, 'vh', start=c(1, 1, i), count=c(nx[1], nx[2], 1))

        energy_pot_on_rho[i] = 0.5 * gravity * sum( (stg[,]*is_wet)**2 * cell_area * is_priority_domain)
        energy_kin_on_rho[i] = 0.5 * 
            sum( ((uh[,]**2/depth_ip)*cell_area + (vh[,]**2/depth_jp)* cell_area_jp) * is_priority_domain)
    }

    nc_close(fid)

    energy_tot = energy_pot_on_rho + energy_kin_on_rho

    output = list(times=times, energy_total=energy_tot*water_density, 
        energy_potential=energy_pot_on_rho*water_density, 
        energy_kinetic=energy_kin_on_rho*water_density)
    return(output)
}

#
# Apply get_energy_truely_linear_domain to a multidomain
# FIXME: Update this to cases other than linear
#
get_energy_truely_linear_multidomain<-function(multidomain_dir, mc_cores=1){

    all_domain_dirs = Sys.glob(paste0(multidomain_dir, '/RUN*'))

    if(length(all_domain_dirs) < 1){
        stop(paste0('Could not find any domain directories beginning with RUN* in ', multidomain_dir))
    }

    library(parallel)

    all_energies = mclapply(all_domain_dirs, get_energy_truely_linear_domain, mc.preschedule=FALSE, mc.cores=24)
    all_energy = all_energies[[1]]$energy_total
    for(i in 2:length(all_energies)) all_energy = all_energy + all_energies[[i]]$energy_total

    output = list(time = all_energies[[1]]$time, energy_total = all_energy, all_domain_energies=all_energies)
    return(output)
}
