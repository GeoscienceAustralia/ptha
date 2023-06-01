#INPUT_RDS = 'OUTPUTS/Fuji_andaman2004_24hrs-full-ambient_sea_level_0.0/RUN_20210812_183314483/gauges_plot_2004-12-26_Fuji_andaman2004_24hrs-full-ambient_sea_level_0.0.RDS'

#INPUT_RDS = "OUTPUTS/Fuji_andaman2004_24hrs_Mandurah2Geraldton-full-ambient_sea_level_0.0/RUN_20211014_143257535/gauges_plot_2004-12-26_Fuji_andaman2004_24hrs_Mandurah2Geraldton-full-ambient_sea_level_0.0.RDS"

INPUT_RDS = commandArgs(trailingOnly=TRUE)[1]


x = readRDS(INPUT_RDS)

target_sites_in_order = c(
    "CocosIsland_6min_BOM",
    "GeraldtonAlternate_WADoT_5min_2004",
    "JurianBay_WADoT_5min_2004",
    "Lancelin_WADoT_5min_2004",
    "BarrackStreet_WADoT_5min_2004",
    "FreemantleAlternate_WADoT_5min_2004",
    "ManglesBay_WADoT_5min_2004",
    "MandurahFishermansJetty_WADoT_5min_2004",
    "PeelInlet_WADoT_5min_2004",
    "CapeBouvard_WADoT_5min_2004",
    "Caddup_WADoT_5min_2004",
    "Harvey_WADoT_5min_2004",
    "Bunbury_WADoT_5min_2004",
    "BunburyInner_WADoT_5min_2004")

output_name_tag = basename(dirname(dirname(INPUT_RDS)))

png(paste0('plots/Perth_region_Sumatra2005_', output_name_tag, '.png'), width=9, height=12, units='in', res=300)
par(mar=c(0, 3, 0, 0))
par(mfrow=c(length(target_sites_in_order), 1))
for(i in 1:length(target_sites_in_order)){

    nm = target_sites_in_order[i]
    obs = x[[nm]]$event_data$obs

    model_time = x[[nm]]$model_start_julian + x[[nm]]$model_time/(3600*24)
    model_stage = x[[nm]]$model_stage

    k0 = which( (obs$juliant > min(model_time)) & (obs$juliant < max(model_time)) )

    #YLIM = range(c(range(model_stage, na.rm=TRUE), range(obs$resid[k0], na.rm=TRUE)))
    YLIM = c(-0.3, 0.3)

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


png('plots/Sumatra2004_Cocos_Island_6min_1min.png', width=12, height=5, units='in', res=300)
plot(x$CocosIsland_BOM_1min_2004$event_data$obs$juliant, x$CocosIsland_BOM_1min_2004$event_data$obs$resid, 
     t='l', col='red', xlim=c(12778, 12779), ylim=c(-0.4,0.4), xlab='Julian Day', ylab='Residual stage (m)',
     main='Cocos Island, Sumatra 2004: Comparison 6 min and 1 min data', cex.main=1.6, cex.lab=1.5, cex.axis=1.5)
points(x$CocosIsland_6min_BOM$event_data$obs$juliant, x$CocosIsland_6min_BOM$event_data$obs$resid, 
       t='l', col='black', lwd=2)
legend('topright', c('6 minute', '1 minute'), col=c('black', 'red'), pch=c(NA, NA), lty=c('solid', 'solid'), lwd=c(2,1), cex=1.5)
grid(col='orange')
dev.off()
