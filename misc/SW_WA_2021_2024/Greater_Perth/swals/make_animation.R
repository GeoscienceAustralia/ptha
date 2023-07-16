#
# Make animations for Newcastle, and Gold Coast
#
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')
library(OpenStreetMap)
library(rgdal)
library(sp)
library(ncdf4)

# Multidomain directory
md_dir = './OUTPUTS/Production_with_animation_43731-full-ambient_sea_level_1.0/RUN_20210526_110810905/'
OUTPUT_BASEDIR = './animation/KT43731/'

# Get domains with high-resolution (domain indices known by inspection)
newcastle_files = sapply(paste0(md_dir, '/RUN*0000', c(153, 154, 158, 159)          , '_*/Grid*.nc'), Sys.glob, USE.NAMES=FALSE)
goldcoast_files = sapply(paste0(md_dir, '/RUN*0000', c(258, 259, 262, 263, 266, 267), '_*/Grid*.nc'), Sys.glob, USE.NAMES=FALSE)

get_output_times<-function(){
    fid = nc_open(newcastle_files[1])
    OUTPUT_TIMES = ncvar_get(fid, 'time')
    nc_close(fid)
    return(OUTPUT_TIMES)
}
OUTPUT_TIMES = get_output_times()

# Get a background image for the site, and save to a file, so that parallel
# plots can read it (without needing internet).
prepare_background_image<-function(site='newcastle', output_basedir=OUTPUT_BASEDIR){

    if(site == 'newcastle'){
        upperLeft  = c(-32.7, 151.65)
        lowerRight = c(-33.0, 151.95)
    }else if(site == 'goldcoast'){
        upperLeft  = c(-27.75, 153.3)
        lowerRight = c(-28.20, 153.6)
    }else if(site == 'goldcoast_south'){
        upperLeft  = c(-28.05, 153.35)
        lowerRight = c(-28.20, 153.6)
    }else{
        stop('unknown_site')
    }

    # Get the OSM tile from the web
    osm_backdrop = openmap(upperLeft = upperLeft, lowerRight= lowerRight, type='bing', zoom=14)
    # Reproject it
    osm_backdrop_reproj = openproj(osm_backdrop, proj4string(CRS("+init=epsg:4326")))

    output_dir = paste0(output_basedir, site, '/')
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    output_file = paste0(output_dir, 'bing_map_', site, '.RDS')
    saveRDS(osm_backdrop_reproj, output_file)

    return(invisible(0))
}

# Plot the wet regions at a single time
plot_time_index<-function(i, 
    site='newcastle', 
    output_basedir=OUTPUT_BASEDIR, 
    ZLIM=c(-6,6), 
    COLS=rainbow(255, alpha=0.5)[1:200],
    dry_threshold = 0.05, # Treat cells with smaller depth as dry
    remove_cells_with_maxima_below=1.05 # These cells never get seriously wet (aside from by the initial condition)
    ){

    if(site == 'newcastle'){
        site_grids = newcastle_files  
        image_ylim = c(-32.98, -32.85)
        image_xlim = c(151.7, 151.85)
        image_width = 6
        image_height = 6
        image_res = 200
    }else if(site == 'goldcoast'){
        site_grids = goldcoast_files

        image_ylim = c(-28.2, -27.88)
        image_xlim = c(153.3, 153.6)
        image_width = 6
        image_height = 7.6
        image_res = 200

    }else if(site == 'goldcoast_south'){
        site_grids = goldcoast_files

        image_ylim = c(-28.2, -28.08)
        image_xlim = c(153.4, 153.6)
        image_width = 6
        image_height = 5.1
        image_res = 200

    }else{
        return(try(log('unknown site_grids value')))
    }

    output_dir = paste0(output_basedir, site, '/')
    osm_backdrop = readRDS(paste0(output_dir, 'bing_map_', site, '.RDS'))
 
    # Merge the high-res tiles for a single time-step 
    if(!is.na(i)){
        # Typical case
        grid_stage = merge_domains_nc_grids(site_grids, desired_var='stage', desired_time_index=i)

        # Compute the time in a 'nice way'
        time_hours = floor(OUTPUT_TIMES[i]/3600)
        time_minutes = floor(OUTPUT_TIMES[i]/60 - time_hours*60)
        time_seconds = floor(OUTPUT_TIMES[i] - time_hours*3600 - time_minutes*60)
        time_string = paste0(
            paste0(formatC(c(time_hours, time_minutes, time_seconds), width=2, format='d', flag=0), collapse=":"),
            ' after earthquake')

    }else{
        # Shortcut to plot the max-stage
        grid_stage = merge_domains_nc_grids(site_grids, desired_var='max_stage')
        names(grid_stage)[3] = 'stage'
        time_string = 'Maximum inundation'
    }
    grid_elev = merge_domains_nc_grids(site_grids, desired_var='elevation0')
    grid_max_stage = merge_domains_nc_grids(site_grids, desired_var='max_stage')

    # Remove dry areas
    grid_stage$stage[grid_stage$stage - grid_elev$elevation0 <= dry_threshold] = NA
    # Remove initially wet areas that are never really affected by the tsunami
    grid_stage$stage[grid_max_stage$max_stage < remove_cells_with_maxima_below] = NA

    # Clip to ZLIM
    grid_stage$stage = pmax(pmin(grid_stage$stage, ZLIM[2]), ZLIM[1])


    #image_xlim = range(grid_stage$xs)
    #image_ylim = range(grid_stage$ys)

    if(!is.na(i)){
        output_file = paste0(output_dir, 'animation_tile_', 1e+06+i, '.png') 
    }else{
        output_file = paste0(output_dir, 'max_stage_tile.png')
    }
    png(output_file, width=image_width, height=image_height, units='in', res=image_res)

    # Make an empty frame
    plot(c(0, 1), c(0, 1), xlim=image_xlim, ylim=image_ylim, asp=1/cos(mean(image_ylim)/180*pi), 
        xlab="", ylab="", main=time_string, cex.main=1.8, xaxs='i', yaxs='i')

    # Plot background image
    plot(osm_backdrop, add=TRUE)

    # Append the gridded data
    image(grid_stage$xs, grid_stage$ys, grid_stage$stage, zlim=ZLIM, col=COLS, add=TRUE, useRaster=TRUE)

    dev.off()

    rm(grid_stage, grid_elev, grid_max_stage, osm_backdrop)
    gc()
    return(invisible(0))
}

make_movie<-function(png_dir){
    system(paste0("ffmpeg -f image2 -framerate 10 -pattern_type glob -i '", png_dir, 
        "/animation*.png' ", png_dir, "/scenario_movie.mp4"))
}

      
# This code downloads the image tiles for newcastle/goldcoast, and saves to files.
# It should only be run once, in serial, on a computer with an internet connection.
all_sites = c('newcastle', 'goldcoast', 'goldcoast_south')
for(i in 1:length(all_sites)){
    if( !file.exists(paste0(OUTPUT_BASEDIR, all_sites[i], '/bing_map_', all_sites[i], '.RDS'))){
        prepare_background_image(site=all_sites[i])
    }
}
      
# Plot the max-stage
plot_time_index(NA, site='newcastle')
plot_time_index(NA, site='goldcoast')

# Make the images. 
library(parallel)
MC_CORES=48
FIRST_TIME = 211 # Nothing happens for the first few hours, so better to do this.

# Newcastle
mclapply(FIRST_TIME:length(OUTPUT_TIMES), plot_time_index, site='newcastle', output_basedir=OUTPUT_BASEDIR,
    mc.cores=MC_CORES, mc.preschedule=TRUE)
make_movie(paste0(OUTPUT_BASEDIR, 'newcastle'))

# Gold Coast
mclapply(FIRST_TIME:length(OUTPUT_TIMES), plot_time_index, site='goldcoast', output_basedir=OUTPUT_BASEDIR,
    mc.cores=MC_CORES, mc.preschedule=TRUE)
make_movie(paste0(OUTPUT_BASEDIR, 'goldcoast'))

# Gold Coast South
mclapply(FIRST_TIME:length(OUTPUT_TIMES), plot_time_index, site='goldcoast_south', output_basedir=OUTPUT_BASEDIR,
    mc.cores=MC_CORES, mc.preschedule=TRUE)
make_movie(paste0(OUTPUT_BASEDIR, 'goldcoast_south'))
