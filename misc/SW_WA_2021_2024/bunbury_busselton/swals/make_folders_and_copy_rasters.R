dirs_to_make = c(
    "OUTPUTS/extreme_source_24hrs_domain301122-full-ambient_sea_level_0.6/RUN_20230309_175044121/",
    "OUTPUTS/Fuji_andaman2004_24hrs_domain301122-full-ambient_sea_level_0.0/RUN_20230309_174201139/"
    )

for(i in 1:length(dirs_to_make)) dir.create(dirs_to_make[i], showWarnings=FALSE, recursive=TRUE)

mydir = getwd()
for(i in 1:length(dirs_to_make)){
    setwd(dirs_to_make[i])
    copy_command = paste0('sshpass -f "PATH_TO_LOCATION_WITH_MY_NCI_LOGIN_DETAILS.txt" scp gxd547@gadi.nci.org.au:/g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/bunbury_busselton/swals/', dirs_to_make[i], '/*.tif .')
    system(copy_command)
    # system("gdalbuildvrt -tr 6.858711e-05 6.858711e-05 all_max_stage.vrt max_stage*.tif")
    setwd(mydir)
}
