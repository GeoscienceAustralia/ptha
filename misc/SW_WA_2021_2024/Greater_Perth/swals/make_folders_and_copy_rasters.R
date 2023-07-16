dirs_to_make = paste0('OUTPUTS/', 
    c(
    'DEBUG_ptha18_random_scenarios_sunda2_row_0109366_Mw_95_HS_newProcessDataToSend-full-ambient_sea_level_0.6/RUN_20220303_104755705/',
    'VAUS_10231_SCALEDBY5_domains010322_newProcessDataToSend-full-ambient_sea_level_0.6/RUN_20220303_105641887/'
    ))

for(i in 1:length(dirs_to_make)) dir.create(dirs_to_make[i], showWarnings=FALSE, recursive=TRUE)

mydir = getwd()
for(i in 1:length(dirs_to_make)){
    setwd(dirs_to_make[i])
    copy_command = paste0('sshpass -f "/home/gareth/NCI_pass.txt" scp gxd547@gadi.nci.org.au:/g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/Greater_Perth/swals/', dirs_to_make[i], '/*.tif .')
    system(copy_command)
    # system("gdalbuildvrt -tr 6.858711e-05 6.858711e-05 all_max_stage.vrt max_stage*.tif")
    setwd(mydir)
}
