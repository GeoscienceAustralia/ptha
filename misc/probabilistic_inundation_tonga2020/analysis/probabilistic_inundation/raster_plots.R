library(rptha)
library(cptcity)

# Name of folder in '../../swals/OUTPUTS/' where the random models are stored.
run_series_name = commandArgs(trailingOnly=TRUE)[1] #'ptha18_tonga_MSL0'

# Exceedance-rate with 2\% chance in 50 years. This is approximately 1/2475, and
# is often used as a reference exceedance-rate for risk-management purposes.
target_exrate = uniroot(f<-function(x){(1 - exp(-50*x)) - 0.02}, 
    lower=0, upper=10, tol=1e-12)$root

# Line shapefile representing the coast of Tonga, for the plot.
tonga_coast = readOGR('../../elevation/Tonga_coast/Tonga_coast_nearlon180.shp', 
    'Tonga_coast_nearlon180')

# We may need to create the 'merged' rasters that add 
#    0.5*(segmented + unsegmented source representations)
# This function does that.
make_combined_segunseg_tifs<-function(){
    # All the rasters with exceedance-rate [given depth-threshold and domain]
    # for this run_series
    all_rast_files = Sys.glob(paste0(run_series_name, 
        '/*_exceedance_rate_with_threshold*.tif'))

    # Extract key information from the filename
    all_rasts_domain_index_and_threshold = unlist(lapply(
        strsplit(basename(all_rast_files), '_domain_'), 
        f<-function(x) x[2]))
    # We need to add unsegmented and segmented source representations where
    # their domain index and threshold match
    unique_all_rasts_domain_i_and_t = 
        unique(all_rasts_domain_index_and_threshold)

    # Get all the rasters with the same domain index and threshold, and combine
    # them.
    for(i in 1:length(unique_all_rasts_domain_i_and_t)){

        k = which(all_rasts_domain_index_and_threshold == 
                  unique_all_rasts_domain_i_and_t[i])
        local_rasts = lapply(all_rast_files[k], f<-function(x){
            out = raster(x)
            out[is.na(out)] = 0
            return(out)
        })

        # To combine segmented and unsegmented we do 
        #    ( 0.5*unseg + 0.5*(sum_of_segments) )
        # which is the same as sum(0.5 * each_raster)
        summed_rast = local_rasts[[1]] * 0
        for(j in 1:length(local_rasts)){
            summed_rast = summed_rast + 0.5*local_rasts[[j]]
        }

        # Easiest to change 0 back to NA so we don't have zeros in halo regions.
        summed_rast[summed_rast == 0] = NA

        output_file = paste0(run_series_name, '/combined_segunseg_DOMAIN_', 
            unique_all_rasts_domain_i_and_t[i])
        writeRaster(summed_rast, output_file, options=c("COMPRESS=DEFLATE"), 
            overwrite=TRUE)
    }

    return(0)
}

# All the rasters with exceedance-rate [given depth-threshold and domain] for
# this run_series
all_rast_files = 
    Sys.glob(paste0(run_series_name, '/combined_segunseg_DOMAIN_*.tif'))

if(length(all_rast_files) == 0){
    # We haven't yet made the combined segmented/unsegmented files. Do it now.
    tmp = make_combined_segunseg_tifs()
    # If that worked, we can now define all_rast_files
    all_rast_files = 
        Sys.glob(paste0(run_series_name, '/combined_segunseg_DOMAIN_*.tif'))
    if(length(all_rast_files) == 0) stop('error finding combined rasters')
}

# The source model [should be 'combined_segunseg_']
all_rasts_source_model = unlist(lapply(
    strsplit(basename(all_rast_files), '_DOMAIN_'), f<-function(x) x[1]))
stopifnot(all(all_rasts_source_model == 'combined_segunseg'))

# The domain index [1-7]
all_rasts_domain_index = as.numeric(unlist(lapply(
    strsplit(basename(all_rast_files), '_DOMAIN_'), 
    f<-function(x) strsplit(x[2], '_')[[1]][1])))

# The depth threshold [0.001, 0.1, .... ]
all_rasts_threshold = as.numeric(unlist(lapply(
    strsplit(basename(all_rast_files), '_threshold_'), 
    f<-function(x) gsub('.tif', '', x[2], fixed=TRUE))))

# Read all the rasters
all_rasts = lapply(all_rast_files, raster)

# Make colours for the plot based on the depth threshold
unique_depth_thresholds = sort(unique(all_rasts_threshold))
depth_id = match(all_rasts_threshold, unique_depth_thresholds)
n = length(unique_depth_thresholds)
#colz_unique = colorRampPalette(
#    c('lightblue', 'darkblue', 'orange', 'red', 'purple', 'black'))(n)
#colz_unique = colorRampPalette(paste0('steelblue', 1:4))(n)
#colz_unique = c('white', cpt("jjg_serrate_seq_srtYlGnBu08", n=n-2), 'black')
colz_unique = c('white', cpt("cb_seq_YlGnBu_09", n=n-2), 'black')

colz_unique_labels = c(
    paste0(unique_depth_thresholds[1:(n-1)], ' - ', unique_depth_thresholds[2:n]),
    paste0('> ', unique_depth_thresholds[n]) )

colz = colz_unique[depth_id]

# Workhorse plotting function
plot_all_rasters<-function(
    threshold_exrate=target_exrate, 
    xlim=c(184.62, 185), 
    ylim=c(-21.25, -21.05)){

    # Plot the shallowest first so we can overlay the deeper ones
    plot_order = order(all_rasts_threshold)

    plot(xlim, ylim, col='white', asp=1/cos(mean(ylim)/180*pi), ann=FALSE)

    for(i in plot_order){
        r1 = (all_rasts[[i]] > threshold_exrate)
        r1[r1==0] = NA
        image(r1, add=TRUE, col=colz[i], maxpixels=Inf)
    }
}

png(paste0('Depth_exceedance_plot_Tongatapu_', run_series_name, '.png'), 
    width=10, height=7.5, res=400, units='in')
plot_all_rasters()
plot(tonga_coast, add=TRUE, col='red')
legend('bottomleft', colz_unique_labels[-1], fill=colz_unique[-1], cex=1.3, 
       title='Depth (m)')
title('Tsunami inundation depth: 2% chance of exceeding in 50 years', 
      cex.main=1.6)
dev.off()


png(paste0('Depth_exceedance_plot_Tongatapu_Zoom_', run_series_name, '.png'), 
    width=10, height=7.5, res=400, units='in')
plot_all_rasters(xlim=c(184.78, 184.845), ylim=c(-21.16, -21.12))
plot(tonga_coast, add=TRUE, col='red')
legend('bottomleft', colz_unique_labels[-1], fill=colz_unique[-1], cex=1.3, 
       title='Depth (m)')
title('Tsunami inundation depth: 2% chance of exceeding in 50 years', 
      cex.main=1.6)
dev.off()

