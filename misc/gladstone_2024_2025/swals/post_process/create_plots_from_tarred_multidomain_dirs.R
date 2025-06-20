#
# Make plots of max_flux in the multidomain directory. Run as:
#    Rscript create_plots_from_tarred_multidomain.R string_matching_tarred_multidomain_dirs
#

input_args = commandArgs(trailingOnly=TRUE)
if(length(input_args) != 1){
    msg = paste0('The code is called with 1 input arguments, like: \n ',
    '    Rscript create_plots_from_tarred_multidomain.R string_matching_tarred_multidomain_dirs')
    stop(msg)
}

# Define the multidomain tar files we process 
MULTIDOMAIN_TAR_FILES = Sys.glob(input_args[1])

# Parallel config
MC_CORES = 48 # Available cores when memory is not an issue

STARTING_DIR = getwd() # Useful to prevent unexpected directory changes in case parts of the code fail.

# Get the SWALS post-processing scripts
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(file_nci)

#
# Check the input arguments
#

if(length(MULTIDOMAIN_TAR_FILES) == 0) stop('No matching MULTIDOMAIN_TAR_FILES were found')

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
         endsWith(tarred_multidomain_dir, '.tar.bz2') | 
         endsWith(tarred_multidomain_dir, '.tar.gz')) ){
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

#' Make plot for a single in a multidomain_directory.
#'
#' In practice we call this in parallel.
#'
#' @param multidomain_dir the multidomain directory
single_plot_creator<-function(multidomain_dir){

    setwd(STARTING_DIR) # Ensure we always start here (protective if errors have occurred)

    setwd(multidomain_dir)
 
    bbox = get_domain_interior_bbox_in_multidomain('.', include_SpatialPolygonsDataFrame=TRUE)

    # Plot max-flux
    png('../max_flux_plot.png', 
        width=8, height=16, units='in', res=300)
    multidomain_image('.', 
        variable='max_flux', 
        time_index=NA, 
        xlim=c(115.55, 115.82), 
        ylim=c(-32.67, -31.4), 
        zlim=c(0.001, 10), 
        cols=rainbow(255)[1:200], 
        var_transform_function=sqrt)
    plot(bbox$all_domain_spdf, add=TRUE, col=NA, border='grey', lwd=0.5)
    title(main=multidomain_dir)
    dev.off()
    gc()

    # Plot max-speed
    png('../max_speed_plot.png', 
        width=8, height=16, units='in', res=300)
    multidomain_image('.', 
        variable='max_speed', 
        time_index=NA, 
        xlim=c(115.55, 115.82), 
        ylim=c(-32.67, -31.4), 
        zlim=c(sqrt(0.001), sqrt(20)), 
        cols=rainbow(255)[1:200], 
        var_transform_function=sqrt)
    plot(bbox$all_domain_spdf, add=TRUE, col=NA, border='grey', lwd=0.5)
    title(main=multidomain_dir)
    dev.off()
    gc()


    # Plot max-stage
    png('../max_stage_plot.png', 
        width=8, height=16, units='in', res=300)
    multidomain_image('.', 
        variable='max_stage', 
        time_index=NA, 
        xlim=c(115.55, 115.82), 
        ylim=c(-32.67, -31.4), 
        zlim=c(0.001, sqrt(10)), 
        cols=rainbow(255)[1:200], 
        #var_transform_function=sqrt)
        var_transform_function=function(x) sqrt(x-0.6),
        NA_if_stage_not_above_elev=TRUE)
    plot(bbox$all_domain_spdf, add=TRUE, col=NA, border='grey', lwd=0.5)
    title(main=multidomain_dir)
    dev.off()
    gc()

    setwd(STARTING_DIR)
    return(TRUE)
}

process_single_md<-function(multidomain_tar_file){

    # Extract the tar file
    untar_worked = untar_tarred_multidomain_dir(multidomain_tar_file)
    if(!untar_worked) stop()

    # Make the plot
    if(endsWith(multidomain_tar_file, '.tar')){
        multidomain_dir = paste0(substring(multidomain_tar_file, 1, nchar(multidomain_tar_file) - 4), '/')
    }else if(endsWith(multidomain_tar_file, '.tar.bz2')){
        multidomain_dir = paste0(substring(multidomain_tar_file, 1, nchar(multidomain_tar_file) - 8), '/')
    }else if(endsWith(multidomain_tar_file, '.tar.gz')){
        multidomain_dir = paste0(substring(multidomain_tar_file, 1, nchar(multidomain_tar_file) - 7), '/')
    }

    plot_worked = try(single_plot_creator(multidomain_dir))

    setwd(STARTING_DIR)

    # Whether or not the plot worked, if we got this far, we want to delete the directory.
    unlink(multidomain_dir, recursive=TRUE)

    return(plot_worked)
}

try_process_single_md<-function(multidomain_tar_file){
    try(process_single_md(multidomain_tar_file))
}

library(parallel)

mclapply(MULTIDOMAIN_TAR_FILES, try_process_single_md, mc.cores=MC_CORES)

