#all_md_dir = Sys.glob('OUTPUTS/*/RUN*')
#all_md_dir = Sys.glob('OUTPUTS/*load_balance*/RUN*')
#all_md_dir = Sys.glob('OUTPUTS/*HIGHRES*/RUN*')
all_md_dir = Sys.glob('OUTPUTS/*/RUN*')
#k = which(!grepl('ptha', all_md_dir))
#all_md_dir = all_md_dir[k]
MC_CORES=24

# Function to replace the directory 'md_dir' with the uncompressed tar archive md_dir.tar,
# and delete it if the operation is successful.
tar_and_remove_dir<-function(md_dir){
    # Only do something if we are working on a directory
    if(!file.info(md_dir)$isdir){
        print(paste0('skipping the following which is not a directory: ', md_dir))
        out = 'skip'
        return(out)
    }

    # Easiest to move to the directory containing md_dir, do the operation, and move back
    start_dir = getwd()
    setwd(dirname(md_dir))
    tar_worked = tar(paste0(basename(md_dir), '.tar'))
    # If the command worked ok we should get an exit code 0. In that case,
    # assume it's ok to delete the archive.
    if(tar_worked == 0) unlink(basename(md_dir), recursive=TRUE)
    # Back to where we were
    setwd(start_dir)
    return(tar_worked)
}

library(parallel)
mclapply(all_md_dir, tar_and_remove_dir, mc.cores=MC_CORES, mc.preschedule=FALSE)
