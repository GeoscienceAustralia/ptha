#
# Convenience script to copy rasters to my home machine
#
dirs_to_make = paste0(commandArgs(trailingOnly=TRUE), '/')

# Set this to the 'swals' folder containing the simulations of interest
gadi_swals_dir = '/g/data/w85/tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals'
# Point this to a file containing your NCI password
NCI_pass_info = '/home/gareth/NCI_pass.txt'

for(i in 1:length(dirs_to_make)) dir.create(dirs_to_make[i], showWarnings=FALSE, recursive=TRUE)

mydir = getwd()
for(i in 1:length(dirs_to_make)){
    setwd(dirs_to_make[i])
    copy_command = paste0('sshpass -f "', NCI_pass_info,'" scp gxd547@gadi.nci.org.au:', 
        gadi_swals_dir, '/', dirs_to_make[i], '/*.tif .')
    system(copy_command)
    # system("gdalbuildvrt -resolution highest all_max_stage.vrt max_stage*.tif")
    setwd(mydir)
}
