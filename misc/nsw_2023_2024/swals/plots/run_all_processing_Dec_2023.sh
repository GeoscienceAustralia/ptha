# [1] "../OUTPUTS/run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534754"         
# [2] "../OUTPUTS/run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535052"                 
# [4] "../OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159"               
# [5] "../OUTPUTS/run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535168"             
# [7] "../OUTPUTS/run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535311"              
# [8] "../OUTPUTS/run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535607"
#[10] "../OUTPUTS/run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023207"        
#[12] "../OUTPUTS/run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534998"                
#[13] "../OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028"              

#[1] "../OUTPUTS/run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534754/gauges_plot_1960-05-22_run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                  
#[2] "../OUTPUTS/run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535052/gauges_plot_1960-05-22_run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                                  
#[3] "../OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/gauges_plot_2010-02-27_run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                              
#[4] "../OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/gauges_plot_2015-09-16_run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                              
#[5] "../OUTPUTS/run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535311/gauges_plot_2021-03-04_run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                            
#[6] "../OUTPUTS/run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535607/gauges_plot_2021-02-10_run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0.RDS"
#[7] "../OUTPUTS/run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023207/gauges_plot_2009-07-15_run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                
#[8] "../OUTPUTS/run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534998/gauges_plot_2007-04-01_run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0.RDS"                                
#[9] "../OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/gauges_plot_2011-03-11_run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0.RDS"             


#Rscript process_gauges_chile2010.R ../OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/
#Rscript process_gauges_chile2015.R ../OUTPUTS/run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535168/
#Rscript process_gauges_kermadec2021.R ../OUTPUTS/run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535311/
#Rscript process_gauges_newhebrides2021.R ../OUTPUTS/run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535607
#Rscript process_gauges_puysegur2009.R ../OUTPUTS/run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023207/
#Rscript process_gauges_chile1960.R ../OUTPUTS/run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534754/
#Rscript process_gauges_chile1960.R ../OUTPUTS/run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535052/
#Rscript process_gauges_solomon2007.R ../OUTPUTS/run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534998/
#Rscript process_gauges_tohoku2011.R ../OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/

# Plots with variable ylim
#Rscript plot_gauges_chile2010.R ../OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/gauges_plot_2010-02-27_run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0.RDS
Rscript plot_gauges_kermadec2021.R ../OUTPUTS/run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535311/gauges_plot_2021-03-04_run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0.RDS 
#Rscript plot_gauges_newhebrides2021.R ../OUTPUTS/run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535607/gauges_plot_2021-02-10_run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0.RDS
#Rscript plot_gauges_puysegur2009.R ../OUTPUTS/run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023207/gauges_plot_2009-07-15_run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0.RDS
#Rscript plot_gauges_chile1960.R ../OUTPUTS/run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534754/gauges_plot_1960-05-22_run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0.RDS
#Rscript plot_gauges_chile1960.R ../OUTPUTS/run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535052/gauges_plot_1960-05-22_run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0.RDS
#Rscript plot_gauges_chile2015.R ../OUTPUTS/run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535168/gauges_plot_2015-09-16_run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0.RDS
#Rscript plot_gauges_solomon2007.R ../OUTPUTS/run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534998/gauges_plot_2007-04-01_run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0.RDS
#Rscript plot_gauges_tohoku2011.R ../OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/gauges_plot_2011-03-11_run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0.RDS

# Plots with fixed ylim
#Rscript plot_gauges_chile2010.R ../OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/gauges_plot_2010-02-27_run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.4
Rscript plot_gauges_kermadec2021.R ../OUTPUTS/run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535311/gauges_plot_2021-03-04_run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.4
#Rscript plot_gauges_newhebrides2021.R ../OUTPUTS/run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535607/gauges_plot_2021-02-10_run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.4
#Rscript plot_gauges_puysegur2009.R ../OUTPUTS/run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023207/gauges_plot_2009-07-15_run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.4
#Rscript plot_gauges_chile1960.R ../OUTPUTS/run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534754/gauges_plot_1960-05-22_run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0.RDS 1.0
#Rscript plot_gauges_chile1960.R ../OUTPUTS/run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535052/gauges_plot_1960-05-22_run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0.RDS 1.0
#Rscript plot_gauges_chile2015.R ../OUTPUTS/run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535168/gauges_plot_2015-09-16_run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.4
#Rscript plot_gauges_solomon2007.R ../OUTPUTS/run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102534998/gauges_plot_2007-04-01_run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.4
#Rscript plot_gauges_tohoku2011.R ../OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/gauges_plot_2011-03-11_run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0.RDS 0.5

