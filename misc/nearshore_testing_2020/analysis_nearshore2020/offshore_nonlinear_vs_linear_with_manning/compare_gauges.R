x_h = readRDS(Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011*/RUN*/gauges_plot_2011-03-11_Tohoku2011_YamakaziEtAl2018-risetime_0-full-leapfrog_nonlinear-0.035-highres_NSW.RDS'))
x_l = readRDS(Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011*/RUN*/gauges_plot_2011-03-11_Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_manning-0.035-highres_NSW.RDS'))

simple_maxima = rep(NA, length(x_h))
nonlinear_maxima = rep(NA, length(x_h))

pdf('compare_stage_ts.pdf', width=12, height=7)
par(mfrow=c(2,1))
for(i in 1:length(x_h)){
    plot(x_h[[i]]$model_t, x_h[[i]]$model_stage, t='l', col='black')
    points(x_l[[i]]$model_t, x_l[[i]]$model_stage, t='l', col='red')
    title(names(x_h)[i], cex.main=1.5)
    grid()
    abline(h=0)
    legend('bottomleft', c('Fully Nonlinear', 'Linear+Manning'), lty=c(1,1), col=c('black', 'red'), bty='n')

    simple_maxima[i] = max(x_l[[i]]$model_stage)
    nonlinear_maxima[i] = max(x_h[[i]]$model_stage)
}
dev.off()

#
# Statistic for paper -- look at the max-stage at the nearshore gauges with high-resolution data 
convergence_stat = (nonlinear_maxima - simple_maxima)/simple_maxima
ni = 1:4
cbind(ni, names(x_h)[ni], convergence_stat[ni])
## RESULTS FOR EARLY VERSION -- BEFORE REVISIONS TO LEAPFROG_NONLINEAR NESTING
#     ni                                                      
#[1,] "1" "TwofoldBay_1min_BOM"        "0.0559981347769198"   
#[2,] "2" "Sydney_FortDenison_1min_PA" "-0.0964633169641134"  
#[3,] "3" "Sydney_BotanyBay_1min_PA"   "-0.0236560609074703"  
#[4,] "4" "PortKembla_1min_BOM"        "-0.000812829831840777"
## RESULTS FOR UPDATES WITH REVISIONS TO LEAPFROG_NONLINEAR NESTING
## VERY SMALL CHANGE WITH PREVIOUS VERSION
#     ni                                                     
#[1,] "1" "TwofoldBay_1min_BOM"        "0.0561076565256234"  
#[2,] "2" "Sydney_FortDenison_1min_PA" "-0.0979128502335137" 
#[3,] "3" "Sydney_BotanyBay_1min_PA"   "-0.0235289530659425" 
#[4,] "4" "PortKembla_1min_BOM"        "5.71139309226991e-05"



png('Linear_with_Manning_vs_Nonlinear_PortKembla_Botany.png', width=10, height=4.5, units='in', res=300)
ADD_OBS = FALSE
par(mfrow=c(2,1))
par(mar = c(1, 6, 3, 1))

i = which(grepl('PortKembla', names(x_h))) # PortKembla
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
legend('topleft', c('Nonlinear', 'Linear + Manning'), lty=c(1,1), horiz=TRUE, 
       col=c('black', 'red'), bty='n', cex=1.2, lwd=c(1, 1))

par(mar = c(4, 6, 1, 1))
i = which(grepl('BotanyBay', names(x_h))) # BotanyBay
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
