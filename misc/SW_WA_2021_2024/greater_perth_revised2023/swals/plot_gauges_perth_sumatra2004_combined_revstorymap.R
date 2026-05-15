inargs = commandArgs(trailingOnly=TRUE)

# First argument = name of RDS file containing the gauges, as created by
# plots/plot_sumatra2004.R.
INPUT_RDS = inargs[1]
x = readRDS(INPUT_RDS)

# Second argument (optional) gives the YLIM scale used for all plots. If not provided
# then each gauge uses its own y-scale. Both versions are useful.
if(length(inargs) > 1){
    YLIM_RANGE = as.numeric(inargs[2])
}else{
    YLIM_RANGE = NA
}


target_sites_in_order = c(
    #"CocosIsland_BOM_1min_2004",
    "GeraldtonAlternate_WADoT_5min_2004",
    "JurianBay_WADoT_5min_2004",
    "Lancelin_WADoT_5min_2004",
    "Hillarys_BOM_1min_2004",
    "BarrackStreet_WADoT_5min_2004",
    "FreemantleAlternate_WADoT_5min_2004",
    "ManglesBay_WADoT_5min_2004",
    "MandurahFishermansJetty_WADoT_5min_2004",
    "PeelInlet_WADoT_5min_2004",
    "CapeBouvard_WADoT_5min_2004",
    "Caddup_WADoT_5min_2004",
    "Harvey_WADoT_5min_2004",
    "Bunbury_WADoT_5min_2004",
    "BunburyInner_WADoT_5min_2004") #,
    #"BusseltonPortGeographe_WADoT_5min_2004")

output_name_tag = paste0(basename(dirname(dirname(INPUT_RDS))), '_', basename(dirname(INPUT_RDS)))
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))

output_dir = 'plots/sumatra2004_gauges_and_separate_figures/'
dir.create(output_dir, showWarnings=FALSE)
png(paste0(output_dir, 'Perth_region_Sumatra2004_', output_name_tag, '.png'), width=9, height=12, units='in', res=300)
par(mar=c(0, 3, 0, 0))
par(oma=c(4,2,0,0))
par(mfrow=c(length(target_sites_in_order), 1))
for(i in 1:length(target_sites_in_order)){

    nm = target_sites_in_order[i]
    obs = x[[nm]]$event_data$obs

    model_time = x[[nm]]$model_start_julian + x[[nm]]$model_time/(3600*24)
    model_time_hours = x[[nm]]$model_time/(3600)
    model_stage = x[[nm]]$model_stage

    k0 = which( (obs$juliant > min(model_time)) & (obs$juliant < max(model_time)) )

    if(is.na(YLIM_RANGE)){
        YLIM = range(c(range(model_stage, na.rm=TRUE), range(obs$resid[k0], na.rm=TRUE)))
    }else{
        YLIM = c(-YLIM_RANGE, YLIM_RANGE)
    }

    obs_time_hours = (as.numeric(obs$juliant) - as.numeric(x[[nm]]$model_start_julian))*24
    plot(obs_time_hours, obs$resid, t='l', ylim=YLIM, xlim=range(model_time_hours), las=1,
        xlab='', ylab='', cex.axis=1.3, cex.lab=1.3)
    points(model_time_hours, model_stage, t='l', col='red')
    points(obs_time_hours, obs$resid, t='l')

    grid(col='orange')
    plot_title = strsplit(nm, '_')[[1]][1]
    plot_title = gsub('Alternate', '', plot_title)
    # Fix a few typos
    plot_title = gsub('Jurian', 'Jurien', plot_title)
    plot_title = gsub('Freemantle', 'Fremantle', plot_title)

    # Insert a space between repeated capitals
    tmp = gsub("([A-Z])", " \\1", plot_title)
    plot_lab_title = substring(tmp, 2, nchar(tmp)) 

    obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' min')
   
    #text(model_time[1] - 0.5/24, YLIM[1]*(0.4/0.7), paste0(plot_title, ' (', obs_spacing, ')'), cex=2, pos=4)
    text(0.5/24, YLIM[1]*(0.4/0.7), paste0(plot_lab_title), cex=2, pos=4)
    if(i == 1) legend('topleft', c('Observed', 'Modelled'), lty=c(1,1), col=c('black', 'red'), cex=1.1)

}
mtext('Hours post earthquake', side=1, line=1.9, cex=1.3, outer=TRUE)
mtext('Detided water level (m)', side=2, line=0.4, cex=1.3, outer=TRUE)

dev.off()


