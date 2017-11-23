# Plot peak stage caused by each unit source, for a given site
library(rptha)

# Unit source files
uss_files = Sys.glob(
    '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/*/TSUNAMI_EVENTS/unit_source*.nc')

#
# Get lon/lat and filename for all unit sources
#
uss_lon_lat_filename = lapply(as.list(uss_files),
    f<-function(uss_file){
        fid = nc_open(uss_file, readunlim=FALSE)

        lon = ncvar_get(fid, 'lon_c')
        lat = ncvar_get(fid, 'lat_c')
        tide_gauge_file = ncvar_get(fid, 'tide_gauge_file')
        rake = ncvar_get(fid, 'rake')
 
        nc_close(fid)

        output = data.frame(lon=lon, lat=lat, tide_gauge_file=tide_gauge_file, rake=rake)

        return(output)
    })
names(uss_lon_lat_filename) = basename(dirname(dirname(uss_files)))

#
# Get hazard points
#
fid = nc_open(uss_lon_lat_filename[[1]]$tide_gauge_file[1], readunlim=FALSE)
hp = data.frame(
    lon = c(ncvar_get(fid, 'lon')),
    lat = c(ncvar_get(fid, 'lat')),
    elev = c(ncvar_get(fid, 'elevation0')),
    gaugeID = c(ncvar_get(fid, 'gaugeID'))
    )
nc_close(fid)


get_unit_source_wave_heights<-function(lon, lat, verbose=FALSE){

    # Find index of nearest point
    ni = lonlat_nearest_neighbours(cbind(lon, lat), cbind(hp$lon, hp$lat))

    output = uss_lon_lat_filename

    for(i in 1:length(uss_lon_lat_filename)){
        if(verbose) print(names(uss_lon_lat_filename)[i])
        # Store outputs here
        output[[i]] = cbind(output[[i]], 
            data.frame(maxstage=rep(NA, length=nrow(output[[i]])),
                minstage=rep(NA, length=nrow(output[[i]])))
            )

        for(j in 1:nrow(output[[i]])){
            # Read peak stage from the unit source file directly
            fid = nc_open(output[[i]]$tide_gauge_file[j], readunlim=FALSE)
            stg = ncvar_get(fid, 'stage', start=c(1, ni), count=c(-1, 1))
            output[[i]]$maxstage[j] = max(stg)
            output[[i]]$minstage[j] = min(stg)
            nc_close(fid)
        }
    }

    return(output)

}

plot_unit_source_wave_heights<-function(lon, lat, verbose=FALSE, bar_scale=5, 
    xlim=NULL, ylim=NULL, rake_range = c(-Inf, Inf), 
    unit_source_wave_heights=NULL){

    if(is.null(unit_source_wave_heights)){
        unit_source_wave_heights = get_unit_source_wave_heights(lon, lat, verbose=verbose)
    }

    plot(hp[,1:2], asp=1, pch='.', xlim=xlim, ylim=ylim)
    points(lon, lat, col='red', pch=19)

    all_us_info = do.call(rbind, unit_source_wave_heights)

    # Zero stages if the unit source rake is wrong
    keep_stage = (all_us_info$rake >= rake_range[1])&(all_us_info$rake <= rake_range[2])
    stage = all_us_info$maxstage
    stage[-which(keep_stage)] = 0

    stg_max = max(stage)
    stg_scaled = stage/stg_max

    colRamp = colorRamp(c('purple', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red'))
    stg_cols = colRamp(stg_scaled)
    stg_cols = rgb(stg_cols[,1], stg_cols[,2], stg_cols[,3], maxColorValue=255)

    ord = order(stg_scaled)

    arrows(all_us_info$lon[ord], all_us_info$lat[ord],
        all_us_info$lon[ord], (all_us_info$lat + stg_scaled*bar_scale)[ord],
        col = stg_cols[ord], length=0)
        
}
