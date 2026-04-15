
all_files = Sys.glob('../OUTPUTS/ptha18-kalbarri2coralbay-hazar*/random_*/*/RUN*.tar')

bzip2_file<-function(filename){
    cmd = paste0('bzip2 ', filename)
    #cmd = paste0('bzip2 -d ', filename)
    result = system(cmd)
    return(result)
}

library(parallel)

mclapply(all_files, bzip2_file, mc.cores=48, mc.preschedule=FALSE)
