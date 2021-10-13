all_md_tar_files = Sys.glob('OUTPUTS/ptha18_tonga_MSL0*/ptha18*/RUN*.tar')

MYDIR = getwd()

convert_tar_to_folder<-function(tar_file){
    setwd(MYDIR)
    setwd(dirname(tar_file))

    unpacked_successfully = untar(basename(tar_file))
    # If the directory was unpacked properly, then delete the tar file
    if(unpacked_successfully == 0) unlink(basename(tar_file))

    setwd(MYDIR)
}

MC_CORES = 48
library(parallel)

mclapply(all_md_tar_files, convert_tar_to_folder, mc.cores=MC_CORES, mc.preschedule=FALSE)
