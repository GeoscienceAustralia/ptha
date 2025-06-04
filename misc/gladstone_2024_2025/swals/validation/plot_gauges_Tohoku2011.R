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

target_sites_in_order = list(
    "ROSSLYN_BAY_BOM" = "Rosslyn Bay",
    "GladstoneFishermansLandingWharf_10min_MSQ" = "Fisherman's Wharf",
    "GladstoneAucklandPoint_10min_MSQ" = "Auckland Point",
    "SouthTreesStormSurge_MSQ" = "South Trees Storm Surge",
    "DART_55023" = "DART 850km NE of Townsville",
    "DART_55012" = "DART 1300km ENE of Townsville"
    )

# check that target sites exist in names of x
assert <- all(names(x_new) %in% target_sites_in_order)
# assert <- all(names(x_old) %in% target_sites_in_order)


output_name_tag = paste0(basename(dirname(dirname(gauge_RDS_new))), '_', basename(dirname(gauge_RDS_new)))
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))

# png(paste0('gauge_plots/', output_name_tag, '_v1_v2.png'), width=9, height=12, units='in', res=300)
png(paste0('gauge_plots/', output_name_tag, '.png'), width=12, height=9, units='in', res=300)

par(mar=c(0, 3, 0, 0))
par(mfrow=c(length(target_sites_in_order), 1))
for (nm in names(target_sites_in_order)) {

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

    plot(obs$juliant, obs$resid, t='l', ylim=YLIM, xlim=xlim, las=1, cex=2)
    # points(model_time, model_resid, t='l', col='red')
    points(model_time_new, model_resid_new, t='l', col='blue', lty=2)
    points(obs$juliant[k0], obs$resid[k0], col='black', t='l')
    grid(col='orange')
    # legend("topleft", legend=c("Expensive Model Residual at 2.42m", "New Model Residual at 2.42m", "Observations Residual"), col=c("red", "blue", "black"), lty=1)
    legend("topleft", legend=c("Model Residual at 0.0m", "Observations Residual"), col=c("blue", "black"), lty=1, cex=1.5)


    plot_title = target_sites_in_order[[nm]]
    plot_title = gsub('Alternate', '', plot_title)

    # if DART in nm
    if(grepl('DART', nm)){
                obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' sec')  # unconvinced this is correct

    } else {
        obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' min')
    }
   
    text(model_time_new[1] + 0.5/24, YLIM[1]*(0.4/0.7), paste0(plot_title, ' (', obs_spacing, ')'), cex=2, pos=4)
}
dev.off()


