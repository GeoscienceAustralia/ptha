
#all_files = Sys.glob('OUTPUTS/ptha18-GreaterPerth-sealevel60cm/*/*/*.tar')
#all_files = Sys.glob('OUTPUTS/ptha18-GreaterPerth-sealevel60cm/*/*/raster_output*.tar.bz2')
#all_files = Sys.glob('OUTPUTS/ptha18-BunburyBusseltonRevised-sealevel60cm/*/*/RUN*.tar')
all_files = Sys.glob('OUTPUTS/ptha18-BunburyBusseltonShutFloodgateRevised-sealevel60cm/*/*/RUN*.tar')


bzip2_file<-function(filename){
    cmd = paste0('bzip2 ', filename)
    #cmd = paste0('bzip2 -d ', filename)
    result = system(cmd)
    return(result)
}

library(parallel)

mclapply(all_files, bzip2_file, mc.cores=48, mc.preschedule=FALSE)
