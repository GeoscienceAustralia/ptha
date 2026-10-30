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

# Take an incomplete model output folder name on the command line, and return the actual model output dir
model_dir_output_folder = Sys.glob(paste0(commandArgs(trailingOnly=TRUE)[1], '*/RUN_*'))
# Ensure we are not matching a tar file
keep = which( !endsWith(model_dir_output_folder, 'tar') )[1]
model_dir_output_folder = model_dir_output_folder[keep]

# Copy a logfile to the parent directory, so we can quickly check model statistics without untarring
md_logfile = Sys.glob(paste0(model_dir_output_folder, '/multidomain_*000001.log'))[1]
file.copy(md_logfile, dirname(model_dir_output_folder))

# Tar the output folder, and if successful, remove the folder
tar_and_remove_dir(model_dir_output_folder)
