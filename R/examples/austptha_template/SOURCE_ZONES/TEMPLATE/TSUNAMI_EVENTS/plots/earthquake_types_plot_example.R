#
# Code to plot earthquake slip and vertical deformation
#

library(rgdal)
library(raster)
library(rgeos)
library(rptha)

# Inputs
unit_source_polygon_shapefile = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/EQ_SOURCE/unit_source_grid/kurilsjapan.shp'

stochastic_slip_events = '../kurilsjapan_tohoku_2011_03_11_Mw9.1_stochastic_stochastic/event_metadata.RDS' 

uniform_slip_events = '../kurilsjapan_tohoku_2011_03_11_Mw9.1_uniform_uniform/event_metadata.RDS' 

variable_uniform_slip_events = '../kurilsjapan_tohoku_2011_03_11_Mw9.1_variable_uniform_variable_uniform/event_metadata.RDS'

dem = raster('../../../../DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif')

zero_contour = '../../../../DATA/ELEV/merged_dem/zero_contour/contour.shp'

out_dir = 'source_animation'

# End inputs


# Read shapefiles
usp = readOGR(unit_source_polygon_shapefile, 
    layer=gsub('.shp', '', basename(unit_source_polygon_shapefile)))

zc = readOGR(zero_contour, gsub('.shp', '', basename(zero_contour)))

# Read tables with earthquake metadata
stoch_meta = readRDS(stochastic_slip_events)
unif_meta = readRDS(uniform_slip_events)
vari_unif_meta = readRDS(variable_uniform_slip_events)

# Unit source statistics
uss = stoch_meta$unit_source_statistics

# Make a plot
plot_region = extent(usp)
plot_poly = as(plot_region, 'SpatialPolygons')
smalldem = crop(dem, plot_region)

dir.create(out_dir, showWarnings=FALSE)


#
# Stochastic slip plots
#
# @param sim_meta  -- one of stoch_meta, unif_meta, vari_unif_meta
# @param base start of filename for each figure
# @param dem_zlim zlim for DEM in plot
# @param max_slip_val before plotting, clip all slip values to below this
# @param surface_def_max before plotting, clip all surface deformation values to below this
# @param surface_def_min before plotting, clip all surface deformation values to above this
# @param uniform_slip flag denoting uniform slip values
many_plots<-function(sim_meta, base='stochastic_slip',
    dem_zlim = c(-6500, 0), max_slip_val=50, 
    surface_def_max=12, surface_def_min=-8,
    uniform_slip=FALSE){
    
    for(i in 1:length(sim_meta$events_with_Mw[,1])){
        
        png_name = paste0(out_dir, '/', base, '_', 100000 + i, '.png')
        print(png_name)
    
        # Get included events
        included_events = get_unit_source_indices_in_event(sim_meta$events_with_Mw[i,])
        included_events_dd = uss$downdip_number[included_events]
        included_events_as = uss$alongstrike_number[included_events]
        included_events_usp = (usp$dwndp_n %in% included_events_dd & 
            usp$alngst_ %in% included_events_as)
    
        if(!uniform_slip){ 
            # Get slip
            included_events_slip = as.numeric(strsplit(
                sim_meta$events_with_Mw$event_slip_string[i], '_')[[1]]) 
        }else{ 
            included_events_slip = included_events * 0 + sim_meta$events_with_Mw$slip[i]
        }
    
        # Make deformation raster
        local_raster_files = uss$initial_condition_file[included_events]
    
        r1 = raster(local_raster_files[1]) * included_events_slip[1]
        l = length(local_raster_files)
        if(l>1){
            for(j in 2:l){
                r1 = r1 + raster(local_raster_files[j])*included_events_slip[j]
            }
        }
    
        # Map slip onto polygons
        included_slip_vals = rep(0, length(usp))
        for(j in 1:length(included_events)){
            matching_ind = which((usp$dwndp_n %in% included_events_dd[j] & 
                usp$alngst_ %in% included_events_as[j]))
            included_slip_vals[matching_ind] = included_events_slip[j]
        }
    
        png(png_name, width=8, height=16, res=75, units='in')
    
        par(mfrow=c(2,1))
    
        plot(usp, asp=1, axes=TRUE)
        image(smalldem, add=TRUE, col=grey(seq(0,1,len=100), alpha=0.3), 
            zlim=dem_zlim, maxpixels=1e+08)
        plot(zc, add=TRUE)
    
        #max_slip_val = 50 # For colours only, limit all slip values to this
        ncol = 100
        colz = c('green', rev(heat.colors(ncol)))
        slip2col = ceiling(pmin(included_slip_vals/max_slip_val, 1)*(ncol)) + 1
        plot(usp, add=TRUE, col=colz[slip2col], 
            density=50)
        title(paste0('Kurils-Japan stochastic slip events'), cex.main=2)
        # Add a legend to the plot
        plot(smalldem, add=TRUE, legend.only=TRUE, zlim=c(0, max_slip_val),
            col = colz)
            
    
        plot(usp, asp=1, axes=TRUE)
        #surface_def_max = 12
        #surface_def_min = -8
        plot(max(min(r1, surface_def_max), surface_def_min), 
            zlim=c(surface_def_min, surface_def_max), 
            col=rev(rainbow(255))[50:255], add=TRUE)
        title(main='Ocean Surface Deformation (m)', cex.main=2)
        plot(usp, add=TRUE)
        plot(zc, add=TRUE)
    
        dev.off()
    }

}

#stop()

many_plots(stoch_meta, base='stochastic_slip')
many_plots(unif_meta, base='uniform_slip', uniform_slip=TRUE)
many_plots(vari_unif_meta, base='variable_uniform_slip')
