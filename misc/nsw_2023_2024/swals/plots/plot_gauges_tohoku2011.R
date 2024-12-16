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

target_sites_in_order = rev(c(
    "TwofoldBay_1min_BOM",
    "Eden_1min_DPIE",
    "PortKembla_1min_BOM",
    "Sydney_BotanyBay_1min_PA",
    "Sydney_FortDenison_1min_PA",
    "Newcastle_east_1min_PA"
    ))

output_name_tag = paste0(basename(dirname(dirname(INPUT_RDS))), '_', basename(dirname(INPUT_RDS)))
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))
if(!is.na(XLIM_RANGE_HRS[1])) output_name_tag = paste0(output_name_tag, '_XLIM_', XLIM_RANGE_HRS[1], '_', XLIM_RANGE_HRS[2])

png(paste0('Tohoku2011_east_coast_', output_name_tag, '.png'), width=12, height=9, units='in', res=300)
par(mar=c(0, 3, 0, 0))
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

    plot(obs$juliant, obs$resid, t='l', ylim=YLIM, xlim=XLIM, las=1)
    points(model_time, model_stage, t='l', col='red')
    points(obs$juliant, obs$resid, t='l')
    grid(col='orange')

    plot_title = strsplit(nm, '_')[[1]][1]
    plot_title = gsub('Alternate', '', plot_title)

    obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' min')
   
    text(XLIM[1] + 0.5/24, YLIM[1]*(0.4/0.7), paste0(plot_title, ' (', obs_spacing, ')'), cex=2, pos=4)
}
dev.off()


