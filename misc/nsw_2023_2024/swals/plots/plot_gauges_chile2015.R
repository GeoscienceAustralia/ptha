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

x = readRDS(INPUT_RDS)

target_sites_in_order = rev(c(
    "Eden_1min_DPIE",
    "TwofoldBay_1min_BOM",
    "Bermagui_1min_DPIE",
    "BatemansBay_PrincessJetty_1min_DPIE",
    "Ulladulla_1min_DPIE",
    "JervisBay_1min_DPIE",
    "CrookhavenHeads_1min_DPIE",
    "PortKembla_1min_BOM",
    "Bundeena_1min_DPIE",
    "Sydney_BotanyBay_1min_PA",
    "Sydney_FortDenison_1min_PA",
    "Sydney_MiddleHarbour_1min_DPIE",
    "Hawkesbury_Patonga_1min_DPIE",
    "Newcastle_east_1min_PA",
    "ShoalBay_1min_DPIE",
    "Forster_1min_DPIE",
    "CrowdyHeadFishermansWharf_1min_DPIE",
    "PortMacquarie_1min_DPIE",
    "CoffsHarbourInnerPumpoutJetty_1min_DPIE",
    "Yamba_1min_DPIE",
    "BallinaBreakwall_1min_DPIE",
    "BrunswickHeads_1min_DPIE",
    "TweedEntranceSouth_1min_DPIE",
    #"GoldCoastSandBypass",
    "LordHoweIsland_1min_DPIE"
    ))

output_name_tag = paste0(basename(dirname(dirname(INPUT_RDS))), '_', basename(dirname(INPUT_RDS)))
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))

png(paste0('Chile2015_east_coast_', output_name_tag, '.png'), width=12, height=15, units='in', res=300)
#png(paste0('Chile2015_east_coast_WIDE_', output_name_tag, '.png'), width=24, height=15, units='in', res=300)
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

    plot(obs$juliant, obs$resid, t='l', ylim=YLIM, xlim=range(model_time), las=1)
    points(model_time, model_stage, t='l', col='red')
    points(obs$juliant, obs$resid, t='l')
    grid(col='orange')

    plot_title = strsplit(nm, '_')[[1]][1]
    plot_title = gsub('Alternate', '', plot_title)

    obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' min')
   
    text(model_time[1] + 0.5/24, YLIM[1]*(0.4/0.7), paste0(plot_title, ' (', obs_spacing, ')'), cex=2, pos=4)
}
dev.off()


