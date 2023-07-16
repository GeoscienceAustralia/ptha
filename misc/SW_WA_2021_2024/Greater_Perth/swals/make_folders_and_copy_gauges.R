dirs_to_make = paste0('OUTPUTS/', c(
    'Fuji_andaman2004_24hrs_domain010322_lowres_timevaryingRealistic-full-ambient_sea_level_0.0/RUN_20220312_212645770/'
    ))

for(i in 1:length(dirs_to_make)) dir.create(dirs_to_make[i], showWarnings=FALSE, recursive=TRUE)

mydir = getwd()
for(i in 1:length(dirs_to_make)){
    setwd(dirs_to_make[i])
    copy_command = paste0('sshpass -f "/home/gareth/NCI_pass.txt" scp gxd547@gadi.nci.org.au:/g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/Greater_Perth/swals/', dirs_to_make[i], '/gauge* .')
    system(copy_command)
    setwd(mydir)
}
