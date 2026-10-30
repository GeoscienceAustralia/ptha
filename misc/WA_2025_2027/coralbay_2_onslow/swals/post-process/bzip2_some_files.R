
all_files = Sys.glob('../OUTPUTS/ptha18-coralbay2onslow-hazard/random_*/*/RUN*.tar')

bzip2_file<-function(filename){
    cmd = paste0('bzip2 ', filename)  # compress
    #cmd = paste0('bzip2 -d ', filename)  # decompress
    result = system(cmd)
    return(result)
}

library(parallel)

mclapply(all_files, bzip2_file, mc.cores=48, mc.preschedule=FALSE)
