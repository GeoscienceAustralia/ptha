#
# Quick raster plots of the inundation depth at a given exceedance-rate
#

library(rptha)
# Get a 'zero-contour' [actually this is from the openstreetmap coastline dataset]
zc = readOGR('../../elevation/Tonga_coast/Tonga_coast_nearlon180.shp', 'Tonga_coast_nearlon180')

# Workhorse plotting function.
# Pass it a set of tiffs (one per domain) giving the field to plot, as well as arguments that
# control the filename, the plot region (nukualofa or tongatapu), and any arguments to 'title'
plot_max_depth<-function(all_depth_tifs, output_file_name_tag, plot_region = 'nukualofa', ...){

    # Depth colours
    my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'))(1000)
    depth_zlim = c(1.0e-02, 10)

    # Choose xlim, ylim, and file controls differently depending on the plot region
    if(plot_region == 'nukualofa'){
       
        XLIM = c(184.75, 184.9)
        YLIM = c(-21.2, -21.13) 

        output_file = paste0(dirname(all_depth_tifs[1]), '/depth_image_nukualofa_', output_file_name_tag, '.png')
        ASP = 1/cos(mean(YLIM)/180*pi)
        png(output_file, width=10, height=10*diff(YLIM)/diff(XLIM) * 1/ASP * 1.5, res=300, units='in')

    }else if(plot_region == 'tongatapu'){

        XLIM = c(184.6, 184.96)
        YLIM = c(-21.28, -21.05) 

        output_file = paste0(dirname(all_depth_tifs[1]), '/depth_image_tongatapu_', output_file_name_tag, '.png')
        ASP = 1/cos(mean(YLIM)/180*pi)
        png(output_file, width=10, height=10*diff(YLIM)/diff(XLIM) * 1/ASP * 1.23, res=300, units='in')

    }else{
        stop('unknown plot_region')
    }

    # Plot the rasters, controlling the plot extent on the first pass
    plot_ext = extent(c(XLIM, YLIM))
    for(i in 1:length(all_depth_tifs)){
        r1 = raster(all_depth_tifs[i])
        r1[r1 < depth_zlim[1]] = NA
        r1[r1 > depth_zlim[2]] = depth_zlim[2]
        if(i == 1){
            plot(r1, col=my_col, asp=ASP, xlim=XLIM, ylim=YLIM, zlim=depth_zlim, ext=plot_ext, maxpixels=Inf)
        }else{
            image(r1, add=TRUE, zlim=depth_zlim, col=my_col, maxpixels=Inf)
        }
    }
    # Add the openstreetmap coast contour
    plot(zc, add=TRUE, col='slategrey', lwd=2)

    title(...)

    dev.off()
}

# Loop over the cases of interest, and make the plots
for(run_series_name in c("ptha18_tonga_MSL0", "ptha18_tonga_MSL0_meshrefine2", "ptha18_tonga_MSL0.8")){
    for(IS_weights in c("", "alternate_")){
        for(plot_region in c('nukualofa', 'tongatapu')){

            sea_level = ifelse(grepl("0.8", run_series_name, fixed=TRUE), 0.8, 0)
            IS_type = ifelse(IS_weights=="", 'selfNormalised', 'regular')

            plot_max_depth(
                Sys.glob(paste0(IS_weights, run_series_name, '/depth_at_exceedance_rate_1_in_2475_domain_*.tif')),
                output_file_name_tag='1_in_2475',
                plot_region = plot_region, 
                main=paste0('Depth with 2% chance of exceedance in 50 yrs \n Sea-level = ', sea_level, '; ', 
                            IS_type, ' importance sampling weights'), cex.main=1.5)

            plot_max_depth(
                Sys.glob(paste0(IS_weights, run_series_name, '/depth_at_exceedance_rate_1_in_475_domain_*.tif')),
                output_file_name_tag='1_in_475',
                plot_region = plot_region, 
                main=paste0('Depth with 10% chance of exceedance in 50 yrs \n Sea-level = ', sea_level, '; ', 
                            IS_type, ' importance sampling weights'), cex.main=1.5)
        }
    }
}
