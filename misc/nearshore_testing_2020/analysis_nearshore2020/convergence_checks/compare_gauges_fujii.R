x_h = readRDS(Sys.glob('../gauge_RDS_files/OUTPUTS/Sumatra2004*/RUN*/gauges_plot_2004-12-26_Sumatra2004_FujiSatake2007_HIGHRES-risetime_0-full-linear_with_manning-0.035-highres_australia.RDS'))
x_l = readRDS(Sys.glob('../gauge_RDS_files/OUTPUTS/Sumatra2004*/RUN*/gauges_plot_2004-12-26_Sumatra2004_FujiSatake2007-risetime_0-full-linear_with_manning-0.035-highres_australia.RDS'))

low_maxima = rep(NA, length(x_h))
high_maxima = rep(NA, length(x_h))

pdf('compare_stage_ts_Fujii.pdf', width=12, height=7)
par(mfrow=c(2,1))
for(i in 1:length(x_h)){
    plot(x_h[[i]]$model_t, x_h[[i]]$model_stage, t='l', col='black')
    points(x_l[[i]]$model_t, x_l[[i]]$model_stage, t='l', col='red')
    title(names(x_h)[i], cex.main=1.5)
    grid()
    abline(h=0)
    legend('bottomleft', c('Highres', 'Lowres'), lty=c(1,1), col=c('black', 'red'), bty='n')

    low_maxima[i] = max(x_l[[i]]$model_stage)
    high_maxima[i] = max(x_h[[i]]$model_stage)
}
dev.off()


#
# Statistic for paper -- look at the max-stage at the nearshore gauges with high-resolution data 
convergence_stat = (high_maxima - low_maxima)/low_maxima
ni = 1:8
cbind(ni, names(x_h)[ni], convergence_stat[ni])
#     ni                                                       
#[1,] "1" "Sydney_FortDenison_1min_PA"    "-0.118968575334406" 
#[2,] "2" "Sydney_BotanyBay_1min_PA"      "0.00273697553161294"
#[3,] "3" "PortKembla_BOM_1min_2004"      "0.00581228779805328"
#[4,] "4" "Portland_BOM_1min_2004"        "0.210024911151178"  
#[5,] "5" "Hillarys_BOM_1min_2004"        "0.0540815996063954" 
#[6,] "6" "Freemantle_WADoT_5min_2004"    "0.020685658226721"  
#[7,] "7" "ManglesBay_WADoT_5min_2004"    "0.0680131954199823" 
#[8,] "8" "BarrackStreet_WADoT_5min_2004" "0.186619280188074"  
summary(convergence_stat[ni])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.118969  0.005043  0.037384  0.053626  0.097665  0.210025 

## HARD-CODED CHECKS OF OUR TWO CASES
convergence_stat_yamakazi = c(-0.0191915510709215, 0.0189574048947332, -0.116351717657874, 
0.0604646444782154)
convergence_stat_fujii = c(-0.118968575334406, 0.00273697553161294, 0.00581228779805328, 
    0.210024911151178, 0.0540815996063954, 0.020685658226721, 0.0680131954199823, 
    0.186619280188074)
summary(c(convergence_stat_yamakazi, convergence_stat_fujii))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.118969 -0.002745  0.019822  0.031074  0.062352  0.210025 


#
# Plot for paper. Port Kembla and Botany Bay give a reasonable sense of the convergence. 
#
png('Convergence_plot_PortKembla_Botany_Fujii.png', width=10, height=4.5, units='in', res=300)
ADD_OBS = FALSE
par(mfrow=c(2,1))
par(mar = c(1, 6, 3, 1))

i = which(grepl('PortKembla', names(x_h))) # Port Kembla
YLIM = c(-0.5, 0.5)
XLIM = c(10, 60)

plot(x_h[[i]]$model_t/3600, x_h[[i]]$model_stage, t='l', col='black', xlab="", ylab="",
     cex.axis=1.2, cex.lab=1.5, ylim=YLIM, las=1, xlim=XLIM)
#mtext('Time (hours post earthquake)', side=1, cex=1.5, line=2.2)
mtext('Stage (m)', side=2, cex=1.3, line=3.2)

points(x_l[[i]]$model_t/3600, x_l[[i]]$model_stage, t='l', col='red')
abline(h=0, col='orange')
#grid(col='orange')
abline(h=c(-1, 1)*0.3, col='orange', lty='dotted')

if(ADD_OBS){
    obs_time = as.numeric(x_l[[i]]$event_data$obs$juliant - x_l[[i]]$model_start_julian) * 24
    points(obs_time, x_l[[i]]$event_data$obs$resid, t='l', col='blue', lty='solid', lwd=0.4)
    legend('bottomleft', c('Observed'), lty=c(1), horiz=TRUE, 
           col=c('blue'), bty='n', cex=1, lwd=c(0.4))
}
legend('topleft', c('High resolution', 'Default resolution'), lty=c(1,1), horiz=TRUE, 
       col=c('black', 'red'), bty='n', cex=1, lwd=c(1, 1))

par(mar = c(4, 6, 1, 1))
i = which(grepl('BotanyBay', names(x_h))) # Botany Bay
YLIM = c(-0.15, 0.15)
plot(x_h[[i]]$model_t/3600, x_h[[i]]$model_stage, t='l', col='black', xlab="", ylab="",
     cex.axis=1.2, cex.lab=1.5, ylim=YLIM, las=1, xlim=XLIM)
mtext('Time (hours post earthquake)', side=1, cex=1.5, line=2.2)
mtext('Stage (m)', side=2, cex=1.3, line=3.2)

points(x_l[[i]]$model_t/3600, x_l[[i]]$model_stage, t='l', col='red')
abline(h=0, col='orange')
#grid(col='orange')

if(ADD_OBS){
    obs_time = as.numeric(x_l[[i]]$event_data$obs$juliant - x_l[[i]]$model_start_julian) * 24
    points(obs_time, x_l[[i]]$event_data$obs$resid, t='l', col='blue', lty='solid', lwd=0.4)
    #legend('bottomleft', c('Observed'), lty=c(1), 
    #       col=c('blue'), bty='n', cex=1, lwd=c(0.4))
}
#legend('topleft', c('High resolution', 'Default resolution'), lty=c(1,1), 
#       col=c('black', 'red'), bty='n', cex=1, lwd=c(1, 1))
abline(h=c(-1, 1)*0.1, col='orange', lty='dotted')

dev.off()

