# Avoid startup messages when calling library

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

    closest_domain = which.min(md_nearest_distances)

    out = md_nearest[[closest_domain]]
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
# all_nc = Sys.glob('RUN_ID000000000*00001_*/Grid*.nc')
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
        nc_grid_files = Sys.glob(paste0(multidomain_dir, '/RUN_ID0*000', domain_index, '_*/Grid*.nc'))
    }
    
    # Open all the files 
    fids = sapply(nc_grid_files, f<-function(x){nc_open(x, readunlim=FALSE)}, simplify=FALSE)

    # Get the 'x' and 'y' and 'time' dimension variables
    xs = lapply(fids, f<-function(x) ncvar_get(x, 'x'))
    ys = lapply(fids, f<-function(x) ncvar_get(x, 'y'))
    ts = lapply(fids, f<-function(x) ncvar_get(x, 'time'))

    # Check that times are compatible in all files
    if(length(ts) > 1){
        for(i in 2:length(ts)){
            if(!all(ts[[i]] == ts[[1]])){
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
#' @param md list with the multidomain info
#'
merge_multidomain_gauges<-function(md){
   
    # For all gauges, figure out if they are in the priority domain
    clean_gauges = list()
    counter = 0
    for(j in 1:length(md)){

        # Some domains might have no gauges -- skip them
        if(class(md[[j]]$gauges) == 'try-error') next

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


    return(merged_gauges)
}

