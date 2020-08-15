#
# Compare model types at a few gauges for illustration
#

#
# An example for Sumatra 2004 at Hillarys
#

Fujii_Sumatra2004_RDS = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/Sumatra2004_FujiSatake2007-risetime_0-full-linear_with_no_friction-0-highres_australia/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Sumatra2004_FujiSatake2007-risetime_0-full-linear_with_manning-0.035-highres_australia/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Sumatra2004_FujiSatake2007-risetime_0-full-linear_with_linear_friction-0-highres_australia/RUN*/gauge*.RDS')
    )

Fujii_sumatra = lapply(Fujii_Sumatra2004_RDS, f<-function(x) readRDS(x))
names(Fujii_sumatra) = Fujii_Sumatra2004_RDS
h = grep("Hillarys", names(Fujii_sumatra[[1]]))

Fujii_Hillarys = lapply(Fujii_sumatra, f<-function(x) x[[h]])

plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Fujii_Hillarys.png'

make_plot<-function(Fujii_Hillarys, plot_filename, plot_ylim, legend_xy=c(35, -0.25), extra_title = '', plot_xlim=c(0, 60)){

    obs_t = (Fujii_Hillarys[[1]]$event_data$obs$juliant - Fujii_Hillarys[[1]]$model_start_julian)*24
    obs_h = Fujii_Hillarys[[1]]$event_data$obs$resid

    png(plot_filename, , width=10, height=6, units='in', res=300)
    par(oma=c(4, 4, 0, 0))
    par(mfrow=c(3,1))
    par(mar=c(1,1,3,1))

    plot(Fujii_Hillarys[[1]]$model_t/3600, Fujii_Hillarys[[1]]$model_stage, t='l', col='red', lwd=0.7, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.axis=1.5, las=1)
    title(paste0('LSWE + frictionless on global grid', extra_title), cex.main=1.7, line=0.5)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    grid(col='orange')
    legend(x=legend_xy[1], y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    plot(Fujii_Hillarys[[2]]$model_t/3600, Fujii_Hillarys[[2]]$model_stage, t='l', col='red', lwd=1, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.main=1.7, cex.axis=1.5, las=1)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    title(paste0('LSWE + Manning-friction on global grid', extra_title), cex.main=1.7, line=0.5)
    grid(col='orange')
    legend(x=35, y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    plot(Fujii_Hillarys[[3]]$model_t/3600, Fujii_Hillarys[[3]]$model_stage, t='l', col='red', lwd=1.2, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.main=1.4, cex.axis=1.5, las=1)
    title(paste0('LSWE + Linear-friction on global grid', extra_title), cex.main=1.7, line=0.5)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    grid(col='orange')
    legend(x=35, y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    mtext('Time (hours after earthquake)', side=1, outer=TRUE, cex=1.4, line=2.5)
    mtext('Detided water-surface (m)', side=2, outer=TRUE, cex=1.4, line=2)
    dev.off()
}

make_plot(Fujii_Hillarys, plot_filename, plot_ylim, extra_title=' @ Hillarys, Sumatra 2004 tsunami')

h = grep("Portland", names(Fujii_sumatra[[1]]))
Fujii_Portland = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Fujii_Portland.png'
make_plot(Fujii_Portland, plot_filename, plot_ylim, extra_title=' @ Portland, Sumatra 2004 tsunami')

h = grep("Mangles", names(Fujii_sumatra[[1]]))
Fujii_ManglesBay = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Fujii_ManglesBay.png'
make_plot(Fujii_ManglesBay, plot_filename, plot_ylim, extra_title=' @ Mangles Bay, Sumatra 2004 tsunami')

h = grep("Barrack", names(Fujii_sumatra[[1]]))
Fujii_BarrackSt = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.1, 0.1)
plot_filename= 'Models_vs_data_Fujii_BarrackSt.png'
make_plot(Fujii_BarrackSt, plot_filename, plot_ylim, extra_title=' @ Barrack St, Sumatra 2004 tsunami')

h = grep("Freemantle", names(Fujii_sumatra[[1]]))
Fujii_Freemantle = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.7, 0.7)
plot_filename= 'Models_vs_data_Fujii_Freemantle.png'
make_plot(Fujii_Freemantle, plot_filename, plot_ylim, extra_title=' @ Freemantle, Sumatra 2004 tsunami')

h = grep("FortDenison", names(Fujii_sumatra[[1]]))
Fujii_FortDenison = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.3, 0.3)
plot_filename= 'Models_vs_data_Fujii_FortDenison.png'
make_plot(Fujii_FortDenison, plot_filename, plot_ylim, extra_title=' @ Fort Denison, Sumatra 2004 tsunami')

h = grep("PortKembla_BOM_1min", names(Fujii_sumatra[[1]]))
Fujii_PortKembla = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.4, 0.4)
plot_filename= 'Models_vs_data_Fujii_PortKembla.png'
make_plot(Fujii_PortKembla, plot_filename, plot_ylim, extra_title=' @ Port Kembla, Sumatra 2004 tsunami')

h = grep("BotanyBay", names(Fujii_sumatra[[1]]))
Fujii_BotanyBay = lapply(Fujii_sumatra, f<-function(x) x[[h]])
plot_ylim = c(-0.2, 0.2)
plot_filename= 'Models_vs_data_Fujii_BotanyBay.png'
make_plot(Fujii_BotanyBay, plot_filename, plot_ylim, extra_title=' @ Botany Bay, Sumatra 2004 tsunami')

#
# An example for Chile 1960 at Fort Denison
#

Ho_Chile1960_RDS = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/Chile1960_HoEtAl2019-risetime_0-full-linear_with_no_friction-0-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Chile1960_HoEtAl2019-risetime_0-full-linear_with_manning-0.035-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Chile1960_HoEtAl2019-risetime_0-full-linear_with_linear_friction-0-highres_NSW/RUN*/gauge*.RDS')
    )

Ho_Chile = lapply(Ho_Chile1960_RDS, f<-function(x) readRDS(x))
names(Ho_Chile) = Ho_Chile1960_RDS
h = grep("FortDenison", names(Ho_Chile[[1]]))

Ho_FortDenison = lapply(Ho_Chile, f<-function(x) x[[h]])

plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Ho_FortDenison.png'
make_plot(Ho_FortDenison, plot_filename, plot_ylim, extra_title=' @ Fort Denison, Chile 1960 tsunami', plot_xlim=c(10, 60))


#
# An example for Tohoku 2011 at Port Kembla and Botany Bay
#
Romano_Tohoku2011_RDS = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_RomanoEtAl2015-risetime_0-full-linear_with_no_friction-0-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_RomanoEtAl2015-risetime_0-full-linear_with_manning-0.035-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_RomanoEtAl2015-risetime_0-full-linear_with_linear_friction-0-highres_NSW/RUN*/gauge*.RDS')
    )

Romano_Tohoku = lapply(Romano_Tohoku2011_RDS, f<-function(x) readRDS(x))
names(Romano_Tohoku) = Romano_Tohoku2011_RDS
# 
h = grep("PortKembla_1min_BOM", names(Romano_Tohoku[[1]]))
Romano_PortKembla_1min_BOM = lapply(Romano_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Romano_PortKembla_1min_BOM.png'
make_plot(Romano_PortKembla_1min_BOM, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Port Kembla, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# Botany
h = grep("BotanyBay", names(Romano_Tohoku[[1]]))
Romano_BotanyBay = lapply(Romano_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.2, 0.2)
plot_filename= 'Models_vs_data_Romano_BotanyBay.png'
make_plot(Romano_BotanyBay, plot_filename, plot_ylim, legend_xy = c(10, -0.075), 
          extra_title= ' @ Botany Bay, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# FortDenison
h = grep("FortDenison", names(Romano_Tohoku[[1]]))
Romano_FortDenison = lapply(Romano_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.3, 0.3)
plot_filename= 'Models_vs_data_Romano_FortDenison.png'
make_plot(Romano_FortDenison, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Fort Denison, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# Twofold
h = grep("Twofold", names(Romano_Tohoku[[1]]))
Romano_FortDenison = lapply(Romano_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.5, 0.5)
plot_filename= 'Models_vs_data_Romano_Twofold.png'
make_plot(Romano_FortDenison, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Twofold Bay, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

#
# Yamakazi
#

Yamakazi_Tohoku2011_RDS = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_no_friction-0-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_manning-0.035-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_linear_friction-0-highres_NSW/RUN*/gauge*.RDS')
    )

Yamakazi_Tohoku = lapply(Yamakazi_Tohoku2011_RDS, f<-function(x) readRDS(x))
names(Yamakazi_Tohoku) = Yamakazi_Tohoku2011_RDS
# 
h = grep("PortKembla_1min_BOM", names(Yamakazi_Tohoku[[1]]))
Yamakazi_PortKembla_1min_BOM = lapply(Yamakazi_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Yamakazi_PortKembla_1min_BOM.png'
make_plot(Yamakazi_PortKembla_1min_BOM, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Port Kembla, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# Botany
h = grep("BotanyBay", names(Yamakazi_Tohoku[[1]]))
Yamakazi_BotanyBay = lapply(Yamakazi_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.2, 0.2)
plot_filename= 'Models_vs_data_Yamakazi_BotanyBay.png'
make_plot(Yamakazi_BotanyBay, plot_filename, plot_ylim, legend_xy = c(10, -0.075), 
          extra_title= ' @ Botany Bay, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# FortDenison
h = grep("FortDenison", names(Yamakazi_Tohoku[[1]]))
Yamakazi_FortDenison = lapply(Yamakazi_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.3, 0.3)
plot_filename= 'Models_vs_data_Yamakazi_FortDenison.png'
make_plot(Yamakazi_FortDenison, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Fort Denison, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# Twofold
h = grep("Twofold", names(Yamakazi_Tohoku[[1]]))
Yamakazi_FortDenison = lapply(Yamakazi_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.5, 0.5)
plot_filename= 'Models_vs_data_Yamakazi_Twofold.png'
make_plot(Yamakazi_FortDenison, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Twofold Bay, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

#
#
#
Satake_Tohoku2011_RDS = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_SatakeEtAl2013-risetime_0-full-linear_with_no_friction-0-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_SatakeEtAl2013-risetime_0-full-linear_with_manning-0.035-highres_NSW/RUN*/gauge*.RDS'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_SatakeEtAl2013-risetime_0-full-linear_with_linear_friction-0-highres_NSW/RUN*/gauge*.RDS')
    )

Satake_Tohoku = lapply(Satake_Tohoku2011_RDS, f<-function(x) readRDS(x))
names(Satake_Tohoku) = Satake_Tohoku2011_RDS
# 
h = grep("PortKembla_1min_BOM", names(Satake_Tohoku[[1]]))
Satake_PortKembla_1min_BOM = lapply(Satake_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.6, 0.6)
plot_filename= 'Models_vs_data_Satake_PortKembla_1min_BOM.png'
make_plot(Satake_PortKembla_1min_BOM, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Port Kembla, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# Botany
h = grep("BotanyBay", names(Satake_Tohoku[[1]]))
Satake_BotanyBay = lapply(Satake_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.2, 0.2)
plot_filename= 'Models_vs_data_Satake_BotanyBay.png'
make_plot(Satake_BotanyBay, plot_filename, plot_ylim, legend_xy = c(10, -0.075), 
          extra_title= ' @ Botany Bay, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# FortDenison
h = grep("FortDenison", names(Satake_Tohoku[[1]]))
Satake_FortDenison = lapply(Satake_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.3, 0.3)
plot_filename= 'Models_vs_data_Satake_FortDenison.png'
make_plot(Satake_FortDenison, plot_filename, plot_ylim, legend_xy = c(10, -0.25), 
          extra_title= ' @ Fort Denison, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

# Twofold
h = grep("Twofold", names(Satake_Tohoku[[1]]))
Satake_FortDenison = lapply(Satake_Tohoku, f<-function(x) x[[h]])
plot_ylim = c(-0.5, 0.5)
plot_filename= 'Models_vs_data_Satake_Twofold.png'
make_plot(Satake_FortDenison, plot_filename, plot_ylim, legend_xy = c(10, -0.23), 
          extra_title= ' @ Twofold Bay, Tohoku 2011 tsunami', plot_xlim=c(10, 60))

