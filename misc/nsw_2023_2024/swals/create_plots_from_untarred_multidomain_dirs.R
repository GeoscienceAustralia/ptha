#
# Make plots of max_flux in the multidomain directory. Run as:
#    Rscript create_plots_from_tarred_multidomain.R string_matching_multidomain_dirs
#

input_args = commandArgs(trailingOnly=TRUE)
if(length(input_args) != 1){
    msg = paste0('The code is called with 1 input arguments, like: \n ',
    '    Rscript create_plots_from_untarred_multidomain.R string_matching_untarred_multidomain_dirs')
    stop(msg)
}

# Define the multidomain tar files we process 
MULTIDOMAIN_DIRS = Sys.glob(input_args[1])

# Parallel config
MC_CORES = 48 # Available cores when memory is not an issue

STARTING_DIR = getwd() # Useful to prevent unexpected directory changes in case parts of the code fail.

# Get the SWALS post-processing scripts
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(file_nci)

#
# Check the input arguments
#

if(length(MULTIDOMAIN_DIRS) == 0) stop('No matching MULTIDOMAIN_DIRS were found')

# End checks


#' Make plot for a single in a multidomain_directory.
#'
#' In practice we call this in parallel.
#'
#' @param multidomain_dir the multidomain directory
single_plot_creator<-function(multidomain_dir){

    setwd(STARTING_DIR) # Ensure we always start here (protective if errors have occurred)

    setwd(multidomain_dir)
 
    #bbox = get_domain_interior_bbox_in_multidomain('.', include_SpatialPolygonsDataFrame=TRUE)

    ## Plot max-flux
    #png('../max_flux_plot.png', 
    #    width=8, height=16, units='in', res=300)
    #multidomain_image('.', 
    #    variable='max_flux', 
    #    time_index=NA, 
    #    xlim=c(115.55, 115.82), 
    #    ylim=c(-32.67, -31.4), 
    #    zlim=c(0.001, 10), 
    #    cols=rainbow(255)[1:200], 
    #    var_transform_function=sqrt)
    #plot(bbox$all_domain_spdf, add=TRUE, col=NA, border='grey', lwd=0.5)
    #title(main=multidomain_dir)
    #dev.off()
    #gc()

    ## Plot max-speed
    #png('../max_speed_plot.png', 
    #    width=8, height=16, units='in', res=300)
    #multidomain_image('.', 
    #    variable='max_speed', 
    #    time_index=NA, 
    #    xlim=c(115.55, 115.82), 
    #    ylim=c(-32.67, -31.4), 
    #    zlim=c(sqrt(0.001), sqrt(20)), 
    #    cols=rainbow(255)[1:200], 
    #    var_transform_function=sqrt)
    #plot(bbox$all_domain_spdf, add=TRUE, col=NA, border='grey', lwd=0.5)
    #title(main=multidomain_dir)
    #dev.off()
    #gc()

    my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'))(1000)
    # Use a log10 transformation to stretch the max-stage colours
    stage_zlim = c(1.0e-03, 1)
    var_transform_fun<-function(x) log10(pmin(pmax(x, stage_zlim[1]), stage_zlim[2]))
    fields_axis_args = list(at=seq(log10(stage_zlim[1]),log10(stage_zlim[2])),
                            labels=round(10**(seq(log10(stage_zlim[1]),log10(stage_zlim[2]))), 3))


    # Plot max-stage
    png('../max_stage_plot.png', 
        width=17.15, height=8, units='in', res=200)
    multidomain_image('.', 
        variable='max_stage', 
        time_index=NA, 
        xlim=c(-40, 320), 
        ylim=c(-79, 68), 
        zlim=log10(stage_zlim), 
        cols=my_col, 
        clip_to_zlim=TRUE,
        use_fields=TRUE,
        buffer_is_priority_domain=TRUE,
        var_transform_function=var_transform_fun,
        fields_axis_args=fields_axis_args,
        NA_if_stage_not_above_elev=TRUE)
    #plot(bbox$all_domain_spdf, add=TRUE, col=NA, border='grey', lwd=0.5)
    title(main=multidomain_dir)
    dev.off()
    gc()

    setwd(STARTING_DIR)
    return(TRUE)
}

process_single_md<-function(multidomain_dir){

    plot_worked = try(single_plot_creator(multidomain_dir))

    setwd(STARTING_DIR)

    return(plot_worked)
}

try_process_single_md<-function(multidomain_dir){
    try(process_single_md(multidomain_dir))
}

library(parallel)

mclapply(MULTIDOMAIN_DIRS, try_process_single_md, mc.cores=MC_CORES)

