#
# Make rasters in the multidomain directory. Run as:
#    Rscript create_tarred_rasters_from_tarred_multidomain.R string_matching_tarred_multidomain_dirs MYCHUNK NCHUNKS
#

input_args = commandArgs(trailingOnly=TRUE)
if(length(input_args) != 3){
    msg = paste0(
        'The code is called with 3 input arguments, like: \n ',
        '    Rscript create_tarred_rasters_from_tarred_multidomain.R string_matching_tarred_multidomain_dirs MYCHUNK NCHUNKS \n',
        'For example:\n',
        '    Rscript create_tarred_rasters_from_tarred_multidomain.R "OUTPUTS/swals_test_run_*/RUN_*.tar" 3 10 \n',
        'will find all tar files that match "OUTPUTS/swals_test_run_*/RUN_*.tar", then split them into 10 groups,\n',
        'and process the 3rd group of files.\n',
        'By splitting into groups it is easier to run the code in parallel.')
    stop(msg)
}

# Define the multidomain tar files we process 
MULTIDOMAIN_TAR_FILES = Sys.glob(input_args[1])

# We will replace 'MULTIDOMAIN_TAR_FILES' with a subset that is defined by the command-line arguments
NCHUNKS = as.numeric(input_args[3]) # Split into this many groups (e.g. 10)
MYCHUNK = as.numeric(input_args[2]) # Run only this group (ranging from 1, ... NCHUNKS)

# Parallel config
MC_CORES_REDUCED = 48 # When making rasters, reduce cores to avoid running out of memory
MC_CORES = 48 # Available cores when memory is not an issue

# The rasters will be packed in a tar archive to help manage the file counts.
RASTER_TAR_FILENAME_RELATIVE_TO_MULTIDOMAIN_TARFILE_DIR = 'raster_output_files_with_speed.tar'

STARTING_DIR = getwd() # Useful to prevent unexpected directory changes in case parts of the code fail.

# Get the SWALS post-processing scripts
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(file_nci)

#
# Check the input arguments
#

if(length(MULTIDOMAIN_TAR_FILES) == 0) stop('No matching MULTIDOMAIN_TAR_FILES were found')
if(!(MYCHUNK >= 1 & MYCHUNK <= NCHUNKS)) stop('Must have 1 <= MYCHUNK <= NCHUNKS')


# End checks

#' Untar a tarred multidomain directory 
#'
#' This produces a folder in the same subdirectory as the multidomain tar file.
#'
#' @param tarred_multidomain_dir The path to a .tar file created by tarring a multidomain directory.
#' @return TRUE if everything worked, or FALSE if it did not.
#'
untar_tarred_multidomain_dir<-function(tarred_multidomain_dir){

    setwd(STARTING_DIR) # Ensure we start here (will be true unless previous failures occurred to prevent some setwd commands)

    working_dir = STARTING_DIR

    if(!file.exists(tarred_multidomain_dir) | 
       !(endsWith(tarred_multidomain_dir, '.tar') | 
         endsWith(tarred_multidomain_dir, '.tar.bz2'))){
        print(paste0('Error: Could not find tar file ', tarred_multidomain_dir))
        return(FALSE)
    }

    # Go to the directory with the tar file, and untar it
    setwd(dirname(tarred_multidomain_dir))
    untarred_successfully = untar(basename(tarred_multidomain_dir))

    # Message if it didn't work
    untar_worked = (untarred_successfully == 0)
    if(!untar_worked) print(paste0('Could not untar the file ', tarred_multidomain_dir))

    setwd(working_dir)

    return(untar_worked)

}

#' Make raster for a single in a multidomain_directory.
#'
#' In practice we call this in parallel.
#'
#' @param i domain index for which we create the raster
#' @param multidomain_dir the multidomain directory
single_multidomain_raster_creator<-function(i, multidomain_dir){

    num_rasts = 0

    print(c('Making raster for domain ', i))
    # Max-stage
    ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='max_stage', 
        domain_index=i, return_raster=TRUE)
    # Elevation
    elev = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0', 
        domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, '/elevation0_domain_', i, '.tif')
    writeRaster(elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    num_rasts = num_rasts + 1

    # Difference between max-stage and elevation
    max_stg_less_elev = ms - elev
    max_stg_less_elev[max_stg_less_elev < 1.0e-03] = NA
    output_file = paste0(multidomain_dir, '/depth_as_max_stage_minus_elevation0_domain_', i, '.tif')
    writeRaster(max_stg_less_elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    num_rasts = num_rasts + 1
    rm(max_stg_less_elev)
    gc()

    ## Max stage, masked in dry areas
    ms[ms < elev + 1.0e-03] = NA
    output_file = paste0(multidomain_dir, '/max_stage_domain_', i, '.tif')
    writeRaster(ms, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    num_rasts = num_rasts + 1

    rm(ms, elev)
    gc()

    #for(desired_var in c('arrival_time', 'max_speed', 'max_flux')){
    for(desired_var in c('arrival_time', 'max_speed')){
        ## Arrival time
        gridded_var = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var=desired_var, 
            domain_index=i, return_raster=TRUE)
        output_file = paste0(multidomain_dir, '/', desired_var, '_domain_', i, '.tif')
        writeRaster(gridded_var, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        num_rasts = num_rasts + 1
        rm(gridded_var); gc()
    }

    ## Export UH at the last time.
    ## Idea is that we might notice nesting artefacts
    #all_times = get_multidomain_output_times(multidomain_dir)
    #N = length(all_times)
    #tm = round(all_times[N])
    #output_file = paste0(multidomain_dir, '/UH_time_', round(tm), 's_domain_', i, '.tif')
    #ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='uh',
    #    desired_time_index=N, domain_index=i, return_raster=TRUE)
    #writeRaster(ms, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    #rm(ms)
    #gc()

    return(num_rasts)
}

try_single_multidomain_raster_creator<-function(i, multidomain_dir){
    try(single_multidomain_raster_creator(i, multidomain_dir))
}


# Split into NCHUNKS -- and we only process the 'MYCHUNK' chunk
library(parallel)
all_chunks = splitIndices(length(MULTIDOMAIN_TAR_FILES), NCHUNKS)
MULTIDOMAIN_TAR_FILES = MULTIDOMAIN_TAR_FILES[all_chunks[[MYCHUNK]]]

# Ensure that the output raster tar files do not already exist. If they do, remove the associated multidomains
# from the vector naming those to be processed.
output_tarfiles = paste0(dirname(MULTIDOMAIN_TAR_FILES), '/', RASTER_TAR_FILENAME_RELATIVE_TO_MULTIDOMAIN_TARFILE_DIR)
k = which(file.exists(output_tarfiles))
if(length(k) > 0){
    print('Some output tarfiles (with rasters) already exist. Skipping the associated MULTIDOMAIN_TAR_FILES')
    MULTIDOMAIN_TAR_FILES = MULTIDOMAIN_TAR_FILES[-k]
}
rm(output_tarfiles, k)
if(length(MULTIDOMAIN_TAR_FILES) == 0) stop('All the multidomain tar files are completed')


#
# Untar the files
#
untar_worked = mclapply(MULTIDOMAIN_TAR_FILES, untar_tarred_multidomain_dir, mc.cores=MC_CORES)
if(any(! unlist(untar_worked))){
    print(paste0('WARNING: Some files did not untar in job ', MYCHUNK, '/', NCHUNKS))
    stop('Deliberate halt to investigate untar failure')
}

# Strip the '.tar' from the tar filename
#multidomain_dirs = paste0(substring(MULTIDOMAIN_TAR_FILES, 1, nchar(MULTIDOMAIN_TAR_FILES) - 4), '/')
remove_tarbz2 = gsub('.tar.bz2', '', MULTIDOMAIN_TAR_FILES, fixed=TRUE)
remove_tar = gsub('.tar', '', remove_tarbz2, fixed=TRUE)
multidomain_dirs = paste0(remove_tar, '/')
if(!all(file.exists(multidomain_dirs))){
    stop('Some multidomain directories do not exist -- failure with untar?')
}

#
# Make rasters 
#
make_the_rasters_parallel<-function(multidomain_dirs){
    all_domain_inds = get_domain_indices_in_multidomain(multidomain_dirs[1])
    all_inputs = expand.grid(all_domain_inds, multidomain_dirs, stringsAsFactors=FALSE)
    raster_jobs = mcmapply(try_single_multidomain_raster_creator, 
        i=all_inputs[,1], 
        multidomain_dir=all_inputs[,2],
        SIMPLIFY=FALSE,
        mc.cores=MC_CORES_REDUCED, # Choose this to avoid running out of memory
        mc.preschedule=TRUE) # Preschedule to reduce memory usage
    return(raster_jobs)
}
rasters_made = make_the_rasters_parallel(multidomain_dirs)
number_rasters_expected = sum(unlist(rasters_made))/length(multidomain_dirs)

#
# Now rasters should have been created in multidomain_dirs.
# Check they have
#
check_rasters_exist<-function(multidomain_dir, number_rasters_for_success){
    rasters_in_multidomain_dir = Sys.glob(paste0(multidomain_dir, '/*.tif'))
    did_it_work = (length(rasters_in_multidomain_dir) == number_rasters_for_success)
    return(did_it_work)
}
check = unlist(lapply(multidomain_dirs, check_rasters_exist, number_rasters_for_success=number_rasters_expected))
if(!all(check)){
    print('A surprising number of rasters were created')
    print(cbind(check, multidomain_dirs))
    stop()
}

#
# Put all the rasters in a tar file
#
tar_rasters_in_dir<-function(multidomain_dir){
    mydir = STARTING_DIR
    setwd(mydir) # In case of some failure earlier
    # Put the rasters in a tar file
    output_tar_filename = paste0(normalizePath(dirname(multidomain_dir)), '/', 
        RASTER_TAR_FILENAME_RELATIVE_TO_MULTIDOMAIN_TARFILE_DIR)

    if(file.exists(output_tar_filename)){
        print(paste0(output_tar_filename, ' already exists'))
    }
    
    # Move to the directory where we will make the tar file
    setwd(multidomain_dir)

    tif_files = Sys.glob('*.tif')
    ## For some reason this doesn't work
    #tarred_command = tar(output_tar_filename, tif_files, compression='none')
    ## But this system call works
    tarred_command = system(paste0('tar -cf ', output_tar_filename, ' *.tif'))

    # Should be 0 on success
    tarred_successfully = (tarred_command == 0)

    # If it worked, delete the rasters
    if(tarred_successfully) unlink(tif_files)

    # Move back
    setwd(mydir)

    return(tarred_successfully)
}

tarred_tiffs = unlist( mclapply(multidomain_dirs, tar_rasters_in_dir, mc.cores=MC_CORES) )
if(!all(tarred_tiffs)){
    print('Some tarring of tifs failed')
    print(cbind(tarred_tiffs, multidomain_dirs))
    stop()
}

#
# Delete the multidomain_dirs (noting we still have the tar files)
#
unlink(multidomain_dirs, recursive=TRUE)

