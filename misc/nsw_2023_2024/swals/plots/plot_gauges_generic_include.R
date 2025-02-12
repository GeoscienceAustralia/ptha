#
# Gauge plotting code that can be used at all sites after defining the following variables
#
#    target_sites_in_order = rev(c(
#        "Eden_1min_DPIE",
#        "PortKembla_1min_BOM",
#        "Sydney_BotanyBay_1min_PA",
#        "Sydney_FortDenison_1min_PA"
#        ))
#    png_name_start = 'Solomon2007_east_coast_'
#    png_width = 12
#    png_height = 6

inargs = commandArgs(trailingOnly=TRUE)

# First argument = name of RDS file containing the gauges
INPUT_RDS = inargs[1]

# Second argument (optional) gives the YLIM scale used for all plots. If not provided
# then each gauge uses its own y-scale. Both versions are useful.
if(length(inargs) > 1){
    YLIM_RANGE = as.numeric(inargs[2])
}else{
    YLIM_RANGE = NA
}

# Third and fourth arguments (optional) give the x-axis range in hours post earthquake
if(length(inargs) > 2){
    XLIM_RANGE_HRS = as.numeric(inargs[3:4])
}else{
    XLIM_RANGE_HRS = NA
}

x = readRDS(INPUT_RDS)

output_name_tag = paste0(basename(dirname(dirname(INPUT_RDS))), '_', basename(dirname(INPUT_RDS)))
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))
if(!is.na(XLIM_RANGE_HRS[1])) output_name_tag = paste0(output_name_tag, '_XLIM_', XLIM_RANGE_HRS[1], '_', XLIM_RANGE_HRS[2])

png(paste0(png_name_start, output_name_tag, '.png'), width=png_width, height=png_height, units='in', res=300)
par(mar=c(0, 3, 0, 0))
par(oma = c(4.5,1,0.5,0.2))
par(mfrow=c(length(target_sites_in_order), 1))
for(i in 1:length(target_sites_in_order)){

    nm = target_sites_in_order[i]
    obs = x[[nm]]$event_data$obs

    model_time = x[[nm]]$model_start_julian + x[[nm]]$model_time/(3600*24)
    model_stage = x[[nm]]$model_stage

    k0 = which( (obs$juliant > min(model_time)) & (obs$juliant < max(model_time)) )

    if(is.na(YLIM_RANGE)){
        YLIM = range(c(range(model_stage, na.rm=TRUE), range(obs$resid[k0], na.rm=TRUE)))
    }else{
        YLIM = c(-YLIM_RANGE, YLIM_RANGE)
    }

    if(is.na(XLIM_RANGE_HRS)[1]){
        XLIM = range(model_time)
    }else{
        XLIM = min(model_time) + XLIM_RANGE_HRS/24
    }

    plot(obs$juliant, obs$resid, t='l', ylim=YLIM, xlim=XLIM, las=1, axes=FALSE)
    points(model_time, model_stage, t='l', col='red')
    points(obs$juliant, obs$resid, t='l')
    #grid(col='orange')
    grid_dy = floor(max(YLIM)*10)/10 # Nearest 10cm
    #abline(h=seq(-grid_dy, grid_dy, length=3), col='orange', lty='dashed')
    abline(h=0, col='orange', lty='dashed')
    axes_locations = x[[nm]]$model_start_julian + seq(0, 60, by=1)/24
    abline(v=axes_locations, col='orange', lty='dotted')
    axis(side=2, at=c(0, grid_dy), line=-2, las=1, cex.axis=1.8)

    plot_title = strsplit(nm, '_')[[1]][1]
    plot_title = gsub('Alternate', '', plot_title)
    plot_title = gsub('InnerPumpoutJetty', '', plot_title)
    plot_title = gsub('FishermansWharf', '', plot_title)

    # Append the time spacing
    obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' min')
    plot_title = paste0(plot_title, ' (', obs_spacing, ')')
    # Remove time spacings matching (1 min) since most gauges have that, and we are short on space.
    plot_title = gsub('(1 min)', '', plot_title, fixed=TRUE) 

   
    text(XLIM[1], YLIM[1]*(0.4/0.7), plot_title, cex=2, pos=4)
}
axis(side=1, at=axes_locations[seq(1,61,5)], labels=seq(0, 60)[seq(1,61,5)], cex.axis=1.8)
mtext('Hours post earthquake', side=1, outer=TRUE, cex=2, line=3)
mtext('Stage (m)', side=2, outer=TRUE, cex=2, line=-1.2)
dev.off()


