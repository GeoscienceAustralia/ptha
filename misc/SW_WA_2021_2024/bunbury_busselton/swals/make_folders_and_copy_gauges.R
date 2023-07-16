dirs_to_make = c(
    'OUTPUTS/Fuji_andaman2004_24hrs_domain301122_high_jetty_friction_time_varying-full-ambient_sea_level_0.0/RUN_20230310_231936232/',
    'OUTPUTS/Fuji_sumatra2005_24hrs_domain301122_high_jetty_friction-full-ambient_sea_level_0.0/RUN_20230310_201942011/'
    )

for(i in 1:length(dirs_to_make)) dir.create(dirs_to_make[i], showWarnings=FALSE, recursive=TRUE)

mydir = getwd()
for(i in 1:length(dirs_to_make)){
    setwd(dirs_to_make[i])
    copy_command = paste0('sshpass -f "PATH_TO_FILE_WITH_MY_NCI_LOGIN_INFO.txt" scp gxd547@gadi.nci.org.au:/g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/bunbury_busselton/swals/', dirs_to_make[i], '/gauge* .')
    system(copy_command)
    setwd(mydir)
}
