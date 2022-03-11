#
# Plot 1/475 and 1/2475 inundation depth with background sea-level = 0m and 0.8m
#
library(rptha)
library(fields)

# Line shapefile representing the coast of Tonga, for the plot.
tonga_coast = readOGR('../../elevation/Tonga_coast/Tonga_coast_nearlon180.shp', 
    'Tonga_coast_nearlon180')

# Colours and x/y/depth limits
#COLZ = c('white', rev(rainbow(255)[1:200]))
library(cptcity)
COLZ = c('white', cpt('jjg_cbac_seq_cbacYlGnBu09', n=250)[51:250])
depth_limit = 6
XLIM = c(184.625, 185.0)
YLIM = c(-21.29, -21)
HIGH_TIDE=0.8

plot_panel<-function(raster_group_matching_strings, raster_group_tiles, xlim=NULL, ylim=NULL, 
    add_transparent_white_seabed=FALSE, add_white_outside_openstreetmap=FALSE){

    panel_titles = paste0(LETTERS[1:length(raster_group_tiles)], ')')

    for(j in 1:length(raster_group_matching_strings)){
        # Read the rasters in this group
        raster_set = raster_group_matching_strings[j]
        depth_files = Sys.glob(raster_set)
        depth_rasts = lapply(depth_files, raster)

        # Plot panel
        for(i in 1:length(depth_files)){

            if(grepl('percentile', depth_files[i])){
                # When all the percentile tifs are overlayed, we get some edge
                # artefacts (minor gaps between rasters) caused by only using
                # one pixel in each 3x3 window. Fix that here -- idea is to
                # buffer NA edges by 1 cell. Note this is not an issue for 
                # the non-percentile files, because they didn't sub-sample the grids.
                dr = focal(depth_rasts[[i]], w=matrix(1, ncol=3, nrow=3), fun=max, NAonly=TRUE, na.rm=TRUE)
                dr[is.nan(dr)] = NA
                depth_rasts[[i]] = dr
            }
            
            r1 = depth_rasts[[i]] * (depth_rasts[[i]] <  depth_limit) + 
                 depth_limit      * (depth_rasts[[i]] >= depth_limit)

            mean_lat = 0.5*(extent(r1)@ymin + extent(r1)@ymax)

            if(i == 1){
                if(is.null(xlim)) xlim=XLIM
                if(is.null(ylim)) ylim=YLIM
                plot(r1, asp=1/cos(mean_lat/180*pi), col=COLZ, xlim=xlim, ylim=ylim,
                     zlim=c(0, depth_limit), maxpixels=Inf, xaxs='i', yaxs='i', cex.axis=1.3, 
                     cex.lab=1.3, axis.args=list(cex.axis=1.2))
            }else{
                image(r1, asp=1/cos(mean_lat/180*pi), zlim = c(0, depth_limit), 
                      col=COLZ, add=TRUE, maxpixels=Inf)
            }
        }
        if(add_transparent_white_seabed){
            elevation_rasts = lapply(Sys.glob('initial_elevation*/*.tif'), raster)
            for(k in 1:length(elevation_rasts)){
                r1 = (elevation_rasts[[k]] < 0)
                cols = rgb(c(1, 1), c(1, 1), c(1, 1), alpha=c(0, 1))
                image(r1, add=TRUE, col=cols, maxpixels=Inf)
            }
        }

        if(add_white_outside_openstreetmap){
            land_water_rasts = lapply(Sys.glob('openstreetmap_land_water_grids/*.tif'), raster)
            for(k in 1:length(land_water_rasts)){
                r1 = (land_water_rasts[[k]] < 0)
                cols = rgb(c(1, 1), c(1, 1), c(1, 1), alpha=c(0, 1))
                image(r1, add=TRUE, col=cols, maxpixels=Inf)
            }

        }

        plot(tonga_coast, add=TRUE, lwd=1, col='black')
        title(raster_group_titles[j], cex.main=1.4)

        text(xlim[1] + 0.1 * diff(xlim), ylim[1] + 0.1*diff(ylim), 
             panel_titles[j], cex=2)

    }

}


#
# Plot of 1/475 and 1/2475 @ logic-tree mean
#

raster_group_matching_strings = c(
    'ptha18_tonga_MSL0_meshrefine4/depth_at_exceedance_rate_1_in_475_domain_*.tif',
    'ptha18_tonga_MSL0_meshrefine4_msl80cm/depth_at_exceedance_rate_1_in_475_domain_*.tif',
    'ptha18_tonga_MSL0_meshrefine4/depth_at_exceedance_rate_1_in_2475_domain_*.tif',
    'ptha18_tonga_MSL0_meshrefine4_msl80cm/depth_at_exceedance_rate_1_in_2475_domain_*.tif')
raster_group_titles = c(
    '1/475 exceedance-rate depth (sea-level = 0.0 m)',
    '1/475 exceedance-rate depth (sea-level = 0.8 m)',
    '1/2475 exceedance-rate depth (sea-level = 0.0 m)',
    '1/2475 exceedance-rate depth (sea-level = 0.8 m)')

png('depth_475_and_2475_multiple_sea_levels.png', width=10, height=8.2, 
    units='in', res=300)
    par(mfrow=c(2,2))
    par(mar=c(2,2,2,1))
    plot_panel(raster_group_matching_strings, raster_group_titles)
dev.off()

##
## Plot of 1/2475 @
## -- Logic-tree mean
## -- 84th percentile
## -- 16th percentile
## -- Logic-tree mean with higher sea-level
##
#raster_group_matching_strings = c(
#    'ptha18_tonga_MSL0_meshrefine4/depth_at_exceedance_rate_1_in_2475_domain_*.tif',
#    'ptha18_tonga_MSL0_meshrefine4/depth_rast_invexrate_2475_percentile_84_subsam_3_maxdpth_10_mindpth_1e-04_Nrand_10000_seed_123_domain_index_*.tif',
#    'ptha18_tonga_MSL0_meshrefine4/depth_rast_invexrate_2475_percentile_16_subsam_3_maxdpth_10_mindpth_1e-04_Nrand_10000_seed_123_domain_index_*.tif',
#    'ptha18_tonga_MSL0_meshrefine4_msl80cm/depth_at_exceedance_rate_1_in_2475_domain_*.tif')
#raster_group_titles = c(
#    '2% in 50y, logic-tree mean (sea-level = 0.0 m)',
#    '2% in 50y, 84th percentile (sea-level = 0.0 m)',
#    '2% in 50y, 16th percentile (sea-level = 0.0 m)',
#    '2% in 50y, logic-tree mean (sea-level = 0.8 m)')
#png('depth_2475_various_uncertaintiess.png', width=10, height=8.2, units='in', res=300)
#    par(mfrow=c(2,2))
#    par(mar=c(2,2,2,1))
#    plot_panel(raster_group_matching_strings, raster_group_titles)
#dev.off()

#
# Similar to above with better 'zoom', and only focus on 1/2475
#
raster_group_matching_strings = c(
    'ptha18_tonga_MSL0_meshrefine4/depth_at_exceedance_rate_1_in_2475_domain_*.tif',
    'ptha18_tonga_MSL0_meshrefine4/depth_rast_invexrate_2475_percentile_84_subsam_3_maxdpth_10_mindpth_1e-04_Nrand_10000_seed_123_domain_index_*.tif',
    'ptha18_tonga_MSL0_meshrefine4/depth_rast_invexrate_2475_percentile_16_subsam_3_maxdpth_10_mindpth_1e-04_Nrand_10000_seed_123_domain_index_*.tif',
    'ptha18_tonga_MSL0_meshrefine4_msl80cm/depth_at_exceedance_rate_1_in_2475_domain_*.tif')
raster_group_titles = c(
    'logic-tree mean (sea-level = 0.0 m), 2% in 50y',
    '84th percentile (sea-level = 0.0 m), 2% in 50y',
    '16th percentile (sea-level = 0.0 m), 2% in 50y',
    'logic-tree mean (sea-level = 0.8 m), 2% in 50y')
png('depth_2475_various_uncertainties_B.png', width=10, height=6.4, units='in', res=300)
    par(mfrow=c(2,2))
    par(mar=c(2,2,2,1))
    plot_panel(raster_group_matching_strings, raster_group_titles, 
               xlim=c(184.63, 184.99), ylim=c(-21.28, -21.06), 
               #add_transparent_white_seabed=TRUE)
               add_white_outside_openstreetmap=TRUE)
dev.off()


