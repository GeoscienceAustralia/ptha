#
# Scripts to work with SWALS output files. 
#
# It is not essential to use R to run SWALS. However, some kind of scripting
# language is useful to work with the outputs, especially if using a
# multidomain and distributed-memory parallel jobs. 
#
# In general any single domain will contain "halo regions" (i.e. where the flow
# state is received from other domains), and other regions where it is the
# "priority domain". Therefore, some care is required to combine results from
# multiple domains while only using "priority domain" cells.
#


#' Get name of folder most recently written to OUTPUTS (according to timestamp in name)
get_recent_output_folder<-function(){
    outputs = dir('OUTPUTS', pattern='RUN_')
    output_folder = sort(outputs, decreasing=TRUE)[1]
    output_folder = paste0('OUTPUTS/', output_folder)
    return(output_folder)
}

#' Get times at which output is written
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

get_lower_left_corner<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])

    # Get the model dimensions
    model_dim_line = grep('lower_left_corner:', file_metadata)
    model_dim = as.numeric(scan(text=file_metadata[model_dim_line], what='character', quiet=TRUE)[2:3])
    return(model_dim)
}

get_model_output_precision<-function(output_folder){
    file_metadata = readLines(Sys.glob(paste0(output_folder, '/', 'Domain_info*'))[1])
    model_dim_line = grep('output_precision', file_metadata)
    dp = as.numeric(scan(text=file_metadata[model_dim_line], what='character', quiet=TRUE)[2])
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
        time_var = time_series_var_list)

    return(output)        
}

#' Read gauges from netcdf file.
get_gauges_netcdf_format<-function(output_folder){
    library('ncdf4')

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
        if(static_names[i] %in% names(gauge_fid$var)){
            tmp = try(ncvar_get(gauge_fid, static_names[i]), silent=TRUE)
            if(class(tmp) == 'try-error') tmp = NA
        }else{
            tmp = NA
        }
        static_var[[static_names[i]]] = tmp
    }

    time_series_var = list()
    for(i in 1:length(time_series_names)){
        if(time_series_names[i] %in% names(gauge_fid$var)){
            tmp = try(ncvar_get(gauge_fid, time_series_names[i]), silent=TRUE)
            if(class(tmp) == 'try-error'){
                tmp = NA
            }
        }else{
            tmp = NA
        }
        time_series_var[[time_series_names[i]]] = tmp
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

    nc_file = Sys.glob(paste0(output_folder, '/Grid_output_*.nc'))

    # If netcdf output does not exist, then use home-brew binary
    if(length(nc_file) != 1){
    
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
        # Netcdf output
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
   
    if(drop_walls){
        output_stage_scan[1,,] = output_stage_scan[2,,]
        output_stage_scan[,1,] = output_stage_scan[,2,]
        output_stage_scan[model_dim[1],,] = output_stage_scan[model_dim[1]-1,,]
        output_stage_scan[,model_dim[2],] = output_stage_scan[, model_dim[2]-1,]
    }

    return(output_stage_scan) 
}

# Convenient function to read in pretty much all of the results
# Because this can sometimes be prohibitive due to memory issues,
# provide options to skip grids and gauges
get_all_recent_results<-function(output_folder=NULL, read_grids=TRUE, read_gauges=TRUE, 
                                 always_read_priority_domain=FALSE, quiet=FALSE){

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
    }

    if(always_read_priority_domain & !read_grids){
        is_priority_domain = try(get_gridded_variable(var='is_priority_domain', output_folder=output_folder), silent=TRUE)
    }

    # time is often not written until the sim is finished
    time = try(get_time(output_folder), silent=TRUE) 
    nx = try(get_model_dim(output_folder), silent=TRUE)
    dx = try(get_dx(output_folder), silent=TRUE)
    lower_left_corner = try(get_lower_left_corner(output_folder), silent=TRUE)

    # x and y --  get from netcdf if we can
    xs = try(get_gridded_variable(var='x', output_folder=output_folder), silent=TRUE)
    if(class(xs) != 'try_error'){
        ys = try(get_gridded_variable(var='y', output_folder=output_folder), silent=TRUE)
    }else{
        xs = lower_left_corner[1] + ((1:nx[1]) - 0.5) * dx[1]
        ys = lower_left_corner[2] + ((1:nx[2]) - 0.5) * dx[2]
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
        lower_left_corner=lower_left_corner, dx=dx))
}

## Nice plot example
#persp(seq(0,200, len=700), seq(0,200, len=700), X[,,11], phi=40, border=NA,
#    col='green', shade=1, axes=FALSE, box=FALSE, asp=1, scale=FALSE)


# X = get_all_recent_results(c(1100, 800))
# library('raster')
# E = raster(t(X$elev[,,1]), xmn=495000, xmx=(495000+11000), ymn=1610000, ymx=1610000+8000, 
#     crs=CRS("+init=epsg:3123"))
# E = flip(E, direction='y')
# library('rasterVis')
# plot3D(E)

#' Convert peak stage output to raster
#'
#' @param swals_out result of get_all_recent_results
#' @param proj4string
#' @param na_above_zero logical. If TRUE, then when elevation is > 0, set all stages to NA
#' @param return_elevation logical. If TRUE, return a list with max_stage AND elevation as rasters
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

#'
#' Make an image plot for a multidomain
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
#'
multidomain_image<-function(multidomain_dir, variable, time_index, xlim, ylim, zlim, cols, add=FALSE,
    var_transform_function = NULL, NA_if_stage_not_above_elev = FALSE){

    library('ncdf4')
    # Find all netcdf
    all_nc = Sys.glob(paste0(multidomain_dir, '/*/Grid*.nc')) 

    # Start a new plot
    if(!add) image(matrix(0, ncol=2, nrow=2), asp=1, col='white', xlim=xlim, ylim=ylim)

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

        # Set areas outside of priority domain to NA
        var[is_priority != 1] = NA

        if(NA_if_stage_not_above_elev){
            # Read stage and elevation, and set var to NA there
            if(variable == 'max_stage'){
                stage = var
            }else{
                stage = ncvar_get(fid, 'stage', start=c(1,1, time_index), count=c(-1,-1,1))
            }
            elevation = ncvar_get(fid, 'elev', start=c(1,1, time_index), count=c(-1,-1,1))
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
        image(xs2, ys2, var, zlim=zlim, col=cols, add=TRUE, useRaster=TRUE)
    }

}

# Given a domain object (i.e. output from 'get_all_recent_results'), find the index and
# distance of the point nearest to 'x', 'y'
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
# merge them. This can provide a method for working with coarray multidomains,
# by first combining each domain into a single output.
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
    if(all(is.null(nc_grid_files)) & (is.na(domain_index) | is.na(domain_index))){
        stop('Must provide EITHER nc_grid_files OR domain_index and multidomain_dir')
    }

    # Find the matching domain files
    if(all(is.null(nc_grid_files))){
        # NOTE: The number of '0' ahead of "domain_index" below protects us
        # against accidently matching other domains, but will eventually fail
        # if we have enough domains. For instance, if we have 10001 domains,
        # then domains '1' and '10001' would both match together, which
        # would be wrong. Some way off however!
        nc_grid_files = Sys.glob(paste0(multidomain_dir, '/RUN_ID0*000', domain_index, '_*/Grid*000', domain_index, '.nc'))
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
merge_multidomain_gauges<-function(md = NA, multidomain_dir=NA){
   
    # For all gauges, figure out if they are in the priority domain
    clean_gauges = list()
    counter = 0

    if(!is.na(multidomain_dir)){
        # We should read the multidomain information
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
        md = lapply(md_domains, 
                    f<-function(x) get_all_recent_results(x, read_grids=FALSE, always_read_priority_domain=TRUE) )
    }else{
        if(class(md) != 'list'){
            stop('Must provide named arguments with either "multidomain_dir" giving the multidomain directory, or a list "md" containing the output from get_all_recent_results for each domain in the multidomain')
        }
        multidomain_dir = dirname(md[[1]]$output_folder)
    }

    for(j in 1:length(md)){

        # Some domains might have no gauges -- skip them
        if( (length(md[[j]]$gauges) == 0) ) next
        if( (class(md[[j]]$gauges) == 'try-error') ) next

        # If gauges exist, they may or may not be in this priority domain
        lons = md[[j]]$gauges$lon
        lats = md[[j]]$gauges$lat
        nearest_d = sapply(1:length(lons), f<-function(x) nearest_point_in_multidomain(lons[x], lats[x], md)$closest_domain)
        
        keep = which(nearest_d == j)

        if(length(keep) == 0) next
        # At this point, definitely we have some gauges of interest

        counter = counter + 1

        clean_gauges[[counter]] = md[[j]]$gauges

        # Get the "keep" subset of variables, corresponding to gauges in priority domain
        # Loop over 1-d variables lon/lat/ etc, and subset in space
        for(i in 1:length(clean_gauges[[counter]])){
            if(class(clean_gauges[[counter]][[i]]) != 'array') next
            # Do not subset time! 
            if(names(clean_gauges[[counter]])[i] == 'time') next
            clean_gauges[[counter]][[i]] = as.numeric(clean_gauges[[counter]][[i]][keep])
        }
        # Subset the time variables
        for(i in 1:length(clean_gauges[[counter]]$time_var)){

            # Non arrays are e.g. missing data -- ignore.
            if(!(class(clean_gauges[[counter]]$time_var[[i]]) %in% c('matrix', 'array'))) next

            dd = dim(clean_gauges[[counter]]$time_var[[i]])
            if(is.null(dd) | length(dd) == 1){
                # There is only one gauge in the whole thing. Since length(keep) > 0, this must be it
                dim(clean_gauges[[counter]]$time_var[[i]]) = c(1, length(clean_gauges[[counter]]$time_var[[i]])) 
            }else{
                # There is more than one gauge
                clean_gauges[[counter]]$time_var[[i]] = clean_gauges[[counter]]$time_var[[i]][keep,,drop=FALSE] 
            }
        }
        # Static variables
        for(i in 1:length(clean_gauges[[counter]]$static_var)){

            # Non arrays are e.g. missing data -- ignore?
            if(!(class(clean_gauges[[counter]]$static_var[[i]]) %in% c('matrix', 'array'))) next
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

                if(class(merged_gauges[[var]]) == 'list') next

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
            # Time variables are 2D -- so rbind should merge properly
            for(var in 1:length(merged_gauges$time_var)){
                stopifnot(ncol(merged_gauges$time_var[[var]]) == ncol( clean_gauges[[counter]]$time_var[[var]] ))
                merged_gauges$time_var[[var]] = rbind(merged_gauges$time_var[[var]], clean_gauges[[counter]]$time_var[[var]])
            }
            for(var in 1:length(merged_gauges$static_var)){
                merged_gauges$static_var[[var]] = c(merged_gauges$static_var[[var]], clean_gauges[[counter]]$static_var[[var]])
            }

        }
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

#
# Make a load balance partition file, using a 'greedy' approach to distribute
# domains to images. The calculations assume that multidomain_dir did not
# use load balancing.
#
# @param multidomain_dir directory containing the multidomain outputs
# @param verbose if FALSE, suppress printing
# @return the function environment invisibly (so you have to use assignment to capture it)
#
make_load_balance_partition_DEFUNCT<-function(multidomain_dir=NA, verbose=TRUE){

    if(is.na(multidomain_dir)){
        multidomain_dir = rev(sort(Sys.glob('./OUTPUTS/RUN*')))[1]
    }

    print(paste0('Generating load_balance_partition.txt file in ', multidomain_dir))
    print('To be valid, this assumes that job used coarrays without any load balancing.')
    print(' Future jobs which use this file should have the same number of images and threads')

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

    md_times = sapply(md_files, f<-function(x) md_runtime(x)$total_wallclock)

    # 
    image_i_domains = matrix(NA, nrow=nrow(md_times), ncol=ncol(md_times))

    mean_domain_times = apply(md_times, 1, mean)

    counter = 0
    for(i in rev(order(mean_domain_times))){
        counter = counter+1
        if(counter == 1){
            # Assign the largest domain in the same order
            image_i_domains[i,] = 1:ncol(image_i_domains)
            total_time = md_times[i, image_i_domains[i,]]
        }else{
            # Assign domains to images based on their "total time" so far
            rank_next_domains = rank(md_times[i,], ties='random')
            rank_total_times = rank(-total_time, ties='random')
            image_i_domains[i,] = match(rank_total_times, rank_next_domains)
            total_time = total_time + md_times[i, image_i_domains[i,]]
        }

        if(length(unique(image_i_domains[i,])) != length(image_i_domains[i,])){
            stop('Some domains are missing / double ups')
        }
    }

    previous_sum_times = colSums(md_times)
    new_sum_times = total_time

    print(paste0('Previous time range: ', diff(range(previous_sum_times))))
    print(paste0('  Previous ', rev(c('min', 'max')), ' times : ', rev(range(previous_sum_times))))
    print(paste0('New time range: ', diff(range(new_sum_times))))
    print(paste0('  New ', rev(c('min', 'max')), ' times : ', rev(range(new_sum_times))))
    print(paste0('Theoretical time reduction: ', max(previous_sum_times) - max(new_sum_times)))

    write.table(image_i_domains, paste0(multidomain_dir, '/load_balance_partition.txt'), 
                sep=" ", row.names=FALSE, col.names=FALSE)

    return(invisible(environment()))
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


#
# Make a load balance partition file, using a 'greedy' approach to distribute
# domains to images. Unlike the old approach, there is no restruction that an equal
# number of domains ends up on each image (so long as their run-time is closer).
# This can improve the speed.
#
# @param multidomain_dir directory containing the multidomain outputs
# @param verbose if FALSE, suppress printing
# @param domain_index_groups either an empty list, or a list with 2 or more vectors, each defining a group of domain indices.
# In the latter case, the partition tries to be approximately equal WITHIN each group first, and then to combine the results
# in a good way. This will generally be less efficient than not using groups, unless you have some other information
# that tells you that the partition should be done in this way.
# @return the function environment invisibly (so you have to use assignment to capture it)
#
make_load_balance_partition<-function(multidomain_dir=NA, verbose=TRUE, domain_index_groups = list()){

    if(is.na(multidomain_dir)){
        multidomain_dir = rev(sort(Sys.glob('./OUTPUTS/RUN*')))[1]
    }

    print(paste0('Generating load_balance_partition.txt file in ', multidomain_dir))
    print(' Future jobs which use this file should have the same number of images and threads')

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
        splitter = partition_into_k(md_times_vec, num_images)
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

        splitter = partition_into_k_with_grouping(md_times_vec, num_images, groups)
    }

    range_new = range(splitter[[2]])
    dsplit = diff(range_new)
    print(paste0('Range of partition total times: ', signif(range_new[1], 4), '-to-', signif(range_new[2], 4),  's '))
    print(paste0('                  (difference): ', signif(dsplit, 4), 's '))
    print(paste0('             (as a percentage): ', signif(dsplit/mean(splitter[[2]])*100, 4),'%'))
    range_old = range(unlist(lapply(md_times, sum)))
    old_range = diff(range_old)
    print(paste0('Previous time range: ', signif(range_old[1], 4), '-to-', signif(range_old[2], 4),  's '))
    print(paste0('      (difference) : ', signif(old_range, 4)))
    print(paste0('Potential run-time reduction (difference in max times):', signif(range_old[2] - range_new[2], 4), 's'))
    print(paste0('      Potential run-time reduction (% of old max time):', 
                 signif((range_old[2] - range_new[2])/range_old[2] * 100, 4), '%'))
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
            # kn will only be in one entry of splitter[[1]] (the list of
            # indices on each images). Find it 
            target_domain = which(unlist(lapply(splitter[[1]], f<-function(x) kn%in%x)))
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

