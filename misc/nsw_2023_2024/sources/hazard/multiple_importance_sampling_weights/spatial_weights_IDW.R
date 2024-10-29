## Weight raster resolution and extent
#regular_xy = as.matrix(expand.grid(seq(140, 180, by=1), seq(-55, -20, by=1)))
#r1 = rast(as.data.frame(cbind(regular_xy, regular_xy[,1]*0))) # Make this coarse
#
#for(sz in 1:length(source_zones)){
#
#    source_zone = source_zones[sz]
#    k = which(!is.na(
#    
#
#}


interp_weights_on_raster<-function(lonlat, weights, sample_names, raster_xs, raster_ys, radius=1){
    # Use inverse distance interpolation to interpolate the weights onto rasters
    library(terra)
    stopifnot(ncol(weights) == length(sample_names))

    # Template for raster outputs
    raster_xy = expand.grid(raster_xs, raster_ys)
    r1 = rast(cbind(raster_xy, raster_xy[,1]*0), crs="EPSG:4326")
    dm = dim(r1)
    
    default_weight = 1/ncol(weights) # Weights for pixels not near PTHA18 hazard points
    k = which(!is.na(rowSums(weights))) # Need to avoid NA PTHA18 hazard point curves.

    outrasts = vector(mode='list', length=length(sample_names)) # Store the rasts
    names(outrasts) = sample_names
    for(i in 1:ncol(weights)){
        r0 = interpIDW(r1, cbind(lonlat[k,], weights[k,i]), fill=default_weight, radius=radius)
        # Enforce default weights at boundary cells (since we will extend those to beyond the raster extent)
        r0[1,,1] = default_weight 
        r0[,1,1] = default_weight
        r0[dm[1],,1] = default_weight
        r0[,dm[2],1] = default_weight

        names(r0) = sample_names[i]

        outrasts[[i]] = r0
    }
    

    # Check the weights sum to give 1.
    rast_sum = outrasts[[1]]
    for(i in 2:length(outrasts)) rast_sum = rast_sum + outrasts[[i]]

    rs = as.matrix(rast_sum)
    if(any(abs(rs - 1) > 1e-06)) stop('Weights do not sum to give 1')
    
    return(outrasts)
}
