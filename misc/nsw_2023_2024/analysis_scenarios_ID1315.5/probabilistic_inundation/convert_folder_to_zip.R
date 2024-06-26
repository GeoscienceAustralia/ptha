# Function to replace the directory 'my_dir' with the zip archive my_dir.zip,
# and delete 'my_dir' if the operation is successful.
zip_and_remove_dir<-function(my_dir){
    # Only do something if we are working on a directory
    if(!file.info(my_dir)$isdir){
        print(paste0('skipping the following which is not a directory: ', my_dir))
        out = 'skip'
        return(out)
    }

    # Easiest to move to the directory containing my_dir, do the operation, and move back
    start_dir = getwd()
    setwd(dirname(my_dir))
    zip_worked = zip(paste0(basename(my_dir), '.zip'), paste0(basename(my_dir)))
    # If the command worked ok we should get an exit code 0. In that case,
    # assume it's ok to delete the archive.
    if(zip_worked == 0) unlink(basename(my_dir), recursive=TRUE)
    # Back to where we were
    setwd(start_dir)
    return(zip_worked)
}
 
