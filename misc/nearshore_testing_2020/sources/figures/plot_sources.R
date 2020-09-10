#
# Plot the source models
#

library(rptha)
library(cptcity)

all_sources = c(
    '../Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM_KAJIURA_SMOOTHED.tif',
    # The Chile 1960 Ho-et-al source does not require Kajiura (it was anyway
    # derived with water-surface unit-sources)
    '../Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif',
    '../Sumatra2004/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Sumatra2004/LoritoEtAl2010/Lorito_sumatra2004_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Sumatra2004/PiatanesiLorito2007/Piatanesi_sumatra2004_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Tohoku2011/SatakeEtAl2013/Satake_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif',
    '../Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM_KAJIURA_SMOOTHED.tif',
    '../Chile2015/RomanoEtAl2016/Illapel_2015_Romano_KAJIURA_SMOOTHED.tif')

all_sources_names = c(
    'Chile 1960 F13',
    'Chile 1960 H19',
    'Sumatra 2004 F07',
    'Sumatra 2004 L10',
    'Sumatra 2004 P07',
    'Chile 2010 F13', 
    'Chile 2010 L11',
    'Tohoku 2011 S13',
    'Tohoku 2011 Y18',
    'Tohoku 2011 R14',
    'Chile 2015 W17',
    'Chile 2015 R16')

coastline = readOGR('zero_contour', layer='zero_contour')

rasts = lapply(all_sources, raster)

if(FALSE){
    #
    # Find the initial potential energy from the source models, using the multidomain file
    # So long as we get one log file per source, all good
    #

    reference_log_files = paste0('../../analysis_nearshore2020/gauge_RDS_files/OUTPUTS/',
                                 c(
        'Chile1960_FujiSatake2013-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200618_140206501/multidomain_log_image_00000000000000000001.log',
        'Chile1960_HoEtAl2019-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200618_154852777/multidomain_log_image_00000000000000000001.log',
        'Sumatra2004_FujiSatake2007-risetime_0-full-linear_with_delayed_linear_friction-0-highres_australia/RUN_20200618_173556019/multidomain_log_image_00000000000000000001.log',
        'Sumatra2004_LoritoEtAl2010-risetime_0-full-linear_with_delayed_linear_friction-0-highres_australia/RUN_20200618_193826090/multidomain_log_image_00000000000000000001.log',
        'Sumatra2004_PiatanesiLorito2007-risetime_0-full-linear_with_delayed_linear_friction-0-highres_australia/RUN_20200618_214110205/multidomain_log_image_00000000000000000001.log',
        'Chile2010_FujiSatake2013-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200618_234357549/multidomain_log_image_00000000000000000001.log',
        'Chile2010_LoritoEtAl2011-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200619_013059101/multidomain_log_image_00000000000000000001.log',
        'Tohoku2011_SatakeEtAl2013-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200619_031737093/multidomain_log_image_00000000000000000001.log',
        'Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200619_050353064/multidomain_log_image_00000000000000000001.log',
        'Tohoku2011_RomanoEtAl2015-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200619_065012059/multidomain_log_image_00000000000000000001.log',
        'Chile2015_WilliamsonEtAl2017-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200619_083635687/multidomain_log_image_00000000000000000001.log',
        'Chile2015_RomanoEtAl2016-risetime_0-full-linear_with_delayed_linear_friction-0-highres_NSW/RUN_20200619_102302795/multidomain_log_image_00000000000000000001.log'))

    initial_energy_on_rho = sapply(reference_log_files, f<-function(x){
        logfile = readLines(x)
        k = grep('Global energy-total', logfile)[1]
        energy0 = as.numeric(logfile[k+1])
        return(energy0)
    })

    names(initial_energy_on_rho) = substring(basename(dirname(dirname(reference_log_files))), 1, 20)

}else{
    # # Here I've copied the result derived from the code above, to reduce
    # # fragility of the file search in case more simulations are conducted in
    # # future
    # initial_energy_on_rho = c(FujiSatake2013 = 7403960839017, HoEtAl2019 = 6612244454705, 
    #         FujiSatake2007 = 2251862472431, LoritoEtAl2010 = 6341031803671, 
    #         PiatanesiLorito2007 = 3394111187122, FujiSatake2013 = 7403960839017, 
    #         LoritoEtAl2011 = 1049038594552, SatakeEtAl2013 = 1533431939409, 
    #         YamakaziEtAl2018 = 2450568199069, RomanoEtAl2015 = 3207686671737, 
    #         WilliamsonEtAl2017 = 63771049934.81, RomanoEtAl2016 = 50556024746.96
    #         )

    initial_energy_on_rho = c(Chile1960_FujiSatake = 7403960839017, Chile1960_HoEtAl2019 = 6612244454705, 
        Sumatra2004_FujiSata = 2249397410783, Sumatra2004_LoritoEt = 6341031803671, 
        Sumatra2004_Piatanes = 3394111187122, Chile2010_FujiSatake = 545012408712.3, 
        Chile2010_LoritoEtAl = 1049038594552, Tohoku2011_SatakeEtA = 1533431939409, 
        Tohoku2011_YamakaziE = 2450129131445, Tohoku2011_RomanoEtA = 3207686671737, 
        Chile2015_Williamson = 63771049934.81, Chile2015_RomanoEtAl = 50556024746.96
        )
}

initial_energy = initial_energy_on_rho * 1024 # Multiply by typical seawater density

make_plot<-function(add_energy_label=FALSE){

    if(add_energy_label){
        png('SourceModels_with_energy.png', width=12, height=9, units='in', res=300)
    }else{
        png('SourceModels.png', width=12, height=9, units='in', res=300)
    }
    par(mfrow=c(3,4))
    par(mar=c(2,2,2.5,3))
    mycol = cpt(pal='h5_dkbluered', n=101)
    mycol[ceiling(length(mycol)/2)] = 'white'
    for(i in 1:length(all_sources)){

        site = basename(dirname(dirname(all_sources[i]))) 

        M = switch(site,
                   "Chile1960" = 14,
                   "Tohoku2011" = 14,
                   "Sumatra2004" = 14,
                   "Chile2010" = 6,
                   "Chile2015" = 3)

        plot_ext = switch(site,
                   "Tohoku2011" = c(138, 146, 35, 42),
                   "Chile1960" = c(280, 290, -47, -36),
                   "Sumatra2004" = c(87, 97, 0, 16),
                   "Chile2010" = c(280, 290, -40, -32),
                   "Chile2015" = c(285, 290, -34, -28))

        XLIM = plot_ext[1:2]
        YLIM = plot_ext[3:4]
        ASP = 1/cos(mean(YLIM)/180*pi)

        # Clip the raster to the color scale range
        r = rasts[[i]]

        plot(r, zlim=c(-1, 1)*M, col=mycol, asp=ASP, xlim=XLIM, ylim=YLIM,
             smallplot = c(0.85, 0.88, 0.1, 0.9), cex.axis=1.4,
             axis.args=list(cex.axis=1.4))
        plot(coastline, add=TRUE, col='darkgreen')
        title(main=all_sources_names[i], line=0.3, cex.main=2)

        if(add_energy_label){
            title(main= gsub('e+', 'x10^', formatC(initial_energy[i], digits=2), fixed=TRUE), 
                  cex.main=2, line=-2, adj=0, 
                  col.main='black')
        }
    }
    dev.off()
}

make_plot(add_energy_label=TRUE)
make_plot()

