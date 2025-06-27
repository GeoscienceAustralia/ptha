#
# Plot explaining how we define the scenario region
#
source('data_for_plots.R')
library(RFOC)
library(cptcity)

ptha18 = new.env()
source('/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R', chdir=TRUE, local=ptha18)

# Find the rasters
input_rasters = c(
    Sys.glob('set_range_of_mw_and_centroi*/chile1960*/*.tif'),
    Sys.glob('set_range_of_mw_and_centroi*/chile1960*/*.tif'))

OUTPUT_DIR = './plotting/chile1960_for_paper/'

dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Good-ish events that were also sampled
# HS: 141314
# FAUS: 9417
# VAUS: 139139


makeplot<-function(event_index, slip_type){

    event_data = ptha18$get_source_zone_events_data('southamerica', slip_type=slip_type, desired_event_rows=event_index)

    print(c(slip_type, event_index, event_data$events$Mw))

    # Open a PNG file
    output_file = paste0(OUTPUT_DIR, 'slip_chile1960_', slip_type, '_', event_index, '.png')
    png(output_file, width=5, height=7, units='in', res=100)

    # Setup the base plot
    plot(c(275, 300), c(-50, -25), col='white', asp=1/cos(-40/180*pi), xlab="", ylab="", cex.axis=1.4)
    plot(tm_wb, col='grey', border='grey', add=TRUE)
    plot(uss$southamerica, col='white', border='black', add=TRUE)

    event_info = get_unit_source_indices_in_event(event_data$events, also_return_slip=TRUE)
    event_slip = event_info$slip
    event_inds = event_info$inds

    event_uss_inds = sapply(event_inds, function(i){
        which(uss$southamerica$dwndp_n == event_data$unit_source_statistics$downdip_number[i] & 
              uss$southamerica$alngst_ == event_data$unit_source_statistics$alongstrike_number[i])})

    #col_val = colfun(pmin(sqrt(event_slip/50), 1))
    col_from_slip<-function(event_slip){
        # Map slip to 0-1
        col_val = approx(c(0, 5, 10, 25, 50, 9999), c(0, 1, 2, 3, 4, 5)/5, xout=event_slip)$y
        # These colours are meant to be robust for colour vision deficiencies,
        # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
        col_R4 = palette('R4') 
        colfun = colorRamp(col_R4[c(3, 4, 7, 2, 1)])
        tmp = colfun(col_val)
        rgb(tmp[,1], tmp[,2], tmp[,3], maxColorValue=255)
    }
    
    plot(uss$southamerica[event_uss_inds,], col=col_from_slip(event_slip), add=TRUE)
   
    # Add a legend 
    legend_slip = seq(0, 50)
    legend_y = rep(-30, length(legend_slip))
    x_scale = 10
    legend_x = 275 + x_scale*seq(0, 1, len=length(legend_slip))
    dx = legend_x[2]-legend_x[1]
    dy = 2
    rect(legend_x, legend_y, legend_x+dx, legend_y+dy, col=col_from_slip(legend_slip), border=NA)
    rect(min(legend_x), min(legend_y), max(legend_x)+dx, max(legend_y)+dy, border='black', col=NA)
    tick_vals = seq(0, 50, by=10)
    tick_x = approx(legend_slip, legend_x, xout=tick_vals)$y
    tick_dy = 0.5
    for(j in 1:length(tick_x)) points(c(tick_x[j], tick_x[j]), c(legend_y[1]-tick_dy, legend_y[1]), t='l')
    text(tick_x + dx/2, rep(legend_y[1], length=length(tick_x)) - tick_dy, tick_vals, pos=1, cex=1.2)
    text(mean(legend_x), max(legend_y)+dy, 'Slip (m)', pos=3, cex=1.6) 
    

    ## Index of Chile 1960
    #i = which(all_events$Mw_lower > 9 & grepl('1960', all_events$event_name))[1]

    ## The ISC-GEM catalogue gives it a rake of 140, like Kanamori's 2019 GJI paper. 
    #s = all_events_focal_mech$strk1[i]
    #d = all_events_focal_mech$dip1[i]
    #r = all_events_focal_mech$rake1[i]
    #lon = all_events_focal_mech$hypo_lon[i]
    #lon = lon + 360*(lon < -40)
    #lat = all_events_focal_mech$hypo_lat[i]
    #Mw = all_events_focal_mech$Mw[i]
    #size = (Mw - 6.5)*0.1
    #ball_col = c('orange', 'green', 'blue')[floor(Mw) - 6]
    #MECH = SDRfoc(s, d, r, PLOT=FALSE)
    #justfocXY(MECH, x=lon, y=lat, focsiz=size, fcol=ball_col)

    # Add band showing range of centroid, with an offset
    i = which(all_events$event_name == 'chile1960-batch2')
    lon_offset = -5
    p1 = c(all_events$event_point_1_lon[i] + lon_offset, all_events$event_point_1_lat[i])
    p2 = c(all_events$event_point_2_lon[i] + lon_offset, all_events$event_point_2_lat[i])
    #points(rbind(p1, p2), t='l', col='orange', lwd=2)
    arrows(p1[1], p1[2], p2[1], p2[2], code=3, angle=90, lwd=2, col='orange')
    text(0.5*(p1[1]+p2[1]) - 2.5, 0.5*(p1[2]+p2[2]) + 0.2, 'Centroid Range', col='orange', cex=1.4, srt=85)

    ## Add uplift/subsidence raster 
    ## raster_file = 'set_range_of_mw_and_centroid/chile1960/chile1960_heterogeneous_slip_136608_count_1.tif'
    #r1 = rast(raster_file)
    #ZLIM = c(-8,8)
    ## Clip raster for plotting
    #r1 = min(r1, ZLIM[2])
    #r1 = max(r1, ZLIM[1])
    ## Make color palette
    #colpal = "h5_dkbluered"
    #ncol = 101
    #COLS = cpt(pal=colpal, n=ncol)
    #COLS_rgb = t(col2rgb(COLS, alpha=TRUE))
    #COLS_rgb[ceiling(ncol/2), 4] = 0 # Middle colour is completely transparent
    #COLS_rgb[,4] = COLS_rgb[,4] * 0.75 # All colours are partly transparent
    #COLZ = rgb(COLS_rgb[,1], COLS_rgb[,2], COLS_rgb[,3], alpha=COLS_rgb[,4], maxColorValue=255)

    #image(r1, add=TRUE, col=COLZ, zlim=ZLIM)

    dev.off()
}

# Good-ish events
makeplot(141314, 'stochastic')
makeplot(9417, 'uniform')
makeplot(139139, 'variable_uniform') 

# Others
makeplot(136766, 'stochastic')
makeplot(9119, 'uniform')
makeplot(136656, 'variable_uniform')
