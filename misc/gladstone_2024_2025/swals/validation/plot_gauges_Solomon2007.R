inargs = commandArgs(trailingOnly=TRUE)

# First argument = name of RDS file containing the gauges, as created by
# _gauges_generic_include.R from running compute_gauges_tohoku_2011.R. Possibly like:
#  ../OUTPUTS/tohoku2011-full-ambient_sea_level_2.42/RUN.../
# gauge_RDS_old = inargs[1]
gauge_RDS_new = inargs[1]
model_sea_level <- 0.0

# Second argument (optional) gives the YLIM scale used for all plots. If not provided
# then each gauge uses its own y-scale. Both versions are useful.
if(length(inargs) == 2){
    YLIM_RANGE = as.numeric(inargs[2])
}else{
    YLIM_RANGE = NA
}

x_new = readRDS(gauge_RDS_new)
# x_old = readRDS(gauge_RDS_old)

# name: as in tide code_interface.R and value: is displayed name in plot
# target_sites_in_order = list(
#     "ROSSLYN_BAY_MSQ" = "Rosslyn Bay MSQ",
#     "SouthTreesStormSurge_MSQ" = "South Trees Storm Surge MSQ",
#     "PORT_ALMA_STORM_SURGE_DESI" = "Port Alma Storm Surge DESI",
#     "BURNETT_HEADS_STORM_SURGE_1min_DESI"="Burnett Heads Storm Surge DESI",
#     "BOWEN_STORM_SURGE_1min_DESI"="Bowen Storm Surge DESI"
# )

# for Coasts and Ports use different names
target_sites_in_order = list(
    "ROSSLYN_BAY_MSQ" = "Rosslyn Bay",
    "SouthTreesStormSurge_MSQ" = "South Trees",
    "PORT_ALMA_STORM_SURGE_DESI" = "Port Alma",
    "BURNETT_HEADS_STORM_SURGE_1min_DESI"="Burnett Heads",
    "BOWEN_STORM_SURGE_1min_DESI"="Bowen"
)


# target_sites_in_order = list(
#     "ROSSLYN_BAY_MSQ" = "Rosslyn Bay MSQ",
#     "SouthTreesStormSurge_MSQ" = "South Trees Storm Surge MSQ",
#     "GladstoneAucklandPoint_10min_MSQ" = "Gladstone Auckland Point MSQ",
#     "PORT_ALMA_STORM_SURGE_DESI" = "Port Alma Storm Surge DESI",
#     "BURNETT_HEADS_STORM_SURGE_1min_DESI"="Burnett Heads Storm Surge DESI",
#     "BOWEN_STORM_SURGE_1min_DESI"="Bowen Storm Surge DESI",
#     "MACKAY_STORM_TIDE_1min_DESI"="Mackay Storm Tide DESI"
# )

# target_sites_in_order = list(
#  "SouthTreesStormSurge_MSQ"="South Trees Storm Surge MSQ",
#  "GladstoneAucklandPoint_10min_MSQ"="Gladstone Auckland Point MSQ",
#  "ROSSLYN_BAY_MSQ"="Rosslyn Bay MSQ",
#  "PORT_ALMA_STORM_SURGE_DESI"="Port Alma Storm Surge DESI",
#  "BURNETT_HEADS_STORM_SURGE_1min_DESI"="Burnett Heads Storm Surge DESI",
#  "BOWEN_STORM_SURGE_1min_DESI"="Bowen Storm Surge DESI",
#  "MACKAY_STORM_TIDE_1min_DESI"="Mackay Storm Tide DESI",
#  "DART_21413_highres"="DART 21413 highres",
#  "DART_21414_highres"="DART 21414 highres",
#  "DART_32411_highres"="DART 32411 highres",
#  "DART_46409_highres"="DART 46409 highres",
#  "DART_46410_highres"="DART 46410 highres",
#  "DART_46413_highres"="DART 46413 highres",
#  "DART_46419_highres"="DART 46419 highres",
#  "DART_51407_highres"="DART 51407 highres",
#  "DART_52402_highres"="DART 52402 highres",
#  "DART_52403_highres"="DART 52403 highres",
#  "DART_52404_highres"="DART 52404 highres",
#  "DART_52405_highres"="DART 52405 highres"
# )

# check that target sites exist in names of x
assert <- all(names(x_new) %in% target_sites_in_order)
# assert <- all(names(x_old) %in% target_sites_in_order)


output_name_tag = paste0(basename(dirname(dirname(gauge_RDS_new))), '_', basename(dirname(gauge_RDS_new)))
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))

# png(paste0('gauge_plots/', output_name_tag, '_v1_v2.png'), width=9, height=12, units='in', res=300)
# png(paste0('gauge_plots/', output_name_tag, '_ALL.png'), width=9, height=18, units='in', res=300)
# png(paste0('gauge_plots/', output_name_tag, '_coasts_and_ports.png'), width=12, height=12, units='in', res=600)
png(paste0('gauge_plots/', output_name_tag, '_wide.png'), width=18, height=12, units='in', res=300)


par(mar=c(0, 5, 0, 0))
par(mfrow=c(length(target_sites_in_order), 1))
for (nm in names(target_sites_in_order)) {
    print(nm)
    obs = x_new[[nm]]$event_data$obs
    
    # model_time = x_old[[nm]]$model_start_julian + x_old[[nm]]$model_time/(3600*24)
    # model_resid = x_old[[nm]]$model_stage - model_sea_level

    model_time_new = x_new[[nm]]$model_start_julian + x_new[[nm]]$model_time/(3600*24)
    model_resid_new = x_new[[nm]]$model_stage - model_sea_level

    k0 = which( (obs$juliant > min(model_time_new)) & (obs$juliant < max(model_time_new)) )

    if(is.na(YLIM_RANGE)){
        # YLIM = range(c(range(model_resid, na.rm=TRUE), range(obs$resid[k0], na.rm=TRUE)), range(model_resid_new, na.rm=TRUE))
        YLIM = range(c(range(obs$resid[k0], na.rm=TRUE)), range(model_resid_new, na.rm=TRUE))
    }else{
        YLIM = c(-YLIM_RANGE, YLIM_RANGE)
    }
    # xlim = range(c(model_time, model_time_new), na.rm=TRUE)
    xlim = range(c(model_time_new), na.rm=TRUE)
    # limit to one day
    # xlim[2] = xlim[1] + 1

    plot(obs$juliant, obs$resid, t='l', ylim=YLIM, xlim=xlim, las=1, cex.axis=2, ylab="")
    # points(model_time, model_resid, t='l', col='red')
    points(obs$juliant[k0], obs$resid[k0], col='black', t='l')
    points(model_time_new, model_resid_new, t='l', col='blue', lty=2)
    # legend("topleft", legend=c("Expensive Model Residual at 2.42m", "New Model Residual at 2.42m", "Observations Residual"), col=c("red", "blue", "black"), lty=1)
    legend("bottomright", legend=c("Model Residual at 0.0m", "Observations Residual"), col=c("blue", "black"), lty=1, cex=2)

    # space the xTick at 2 hour intervals
    start_grid <- floor(as.numeric(xlim[1]))
    end_grid <- ceiling(as.numeric(xlim[2]))
    grid_time = seq(start_grid, end_grid, by=1/12)
    abline(v=grid_time, col='orange', lty=2)
    # horizontal grid lines
    grid(lty=2, col="orange")

    plot_title = target_sites_in_order[[nm]]
    plot_title = gsub('Alternate', '', plot_title)

    # get the minumum time difference between observations
    obs_spacing = paste0( round(min(diff(obs$juliant))*24*60,1), ' min')

    text(model_time_new[1] + 0.5/24, YLIM[1]*(0.4/0.7), paste0(plot_title, ' (', obs_spacing, ')'), cex=2, pos=4)
}
dev.off()


