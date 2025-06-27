#
# Following the "modelling of source-inversions" work, I realised I have to tar the multidomain
# directories to manage file-inode-count issues. 
#
# This breaks the workflows that I used previously.
#
# We need a workflow where we can take a tarred multidomain directory as input, and then do:
#     untar it to a temporary location
#     create the gauge_RDS file and the pdf plot file
#     create the image plots of max-stage and elevation0
#     copy the gauge data, the plots, and one of the log-files to another location
#     delete the directory that I untarred.
#
# Beware this script makes various assumptions about the relative placement of directories,
# and it also creates and deletes files! So be careful
#

# Beginning of the working directory where we untar outputs, do processing, and
# then delete everything. Should be inside the current folder.
TEMP_DIR_BASE = './UNTAR_PROCESSING' 
# Location of the plotting scripts for each historic event -- should be a directory inside the current folder.
PLOT_SCRIPT_BASE = './plots/' 
# Location where we put the outputs of interest. Should be one-level-up from the current folder.
COPYOUT_DIR_BASE = '../analysis_ptha18_scenarios_2025/' 
# If there are more than this many folders being extracted in TEMP_DIR_BASE, throw an error. 
# This is to stop errors causing the runaway creation of very large file counts.
MAX_UNTARRED_FOLDERS = 200 

# Split a folder path into its directory names
split_folder_path <- function(x){
    if(dirname(x) == x){
        return(x)
    }else{
        out = c(basename(x),split_folder_path(dirname(x)))
        return(return(out))
    }
}

# Given a tarred multidomain directory, untar it in a new location derived by appending
# temp_dir_base to the full path name. The idea is that we work on the untarred directory,
# and then delete it later, without clobbering the original tar archive.
untar_multidomain_dir_and_return_dir_location<-function(md_tar_dir, temp_dir_base = TEMP_DIR_BASE){

    nc = nchar(md_tar_dir)

    # Remove './' from the start of the directory if it exists
    if( startsWith(md_tar_dir, './') ) md_tar_dir = substring(md_tar_dir, 3, nc)

    # The file could end in .tar or .tar.bz2
    is_tar = (endsWith(md_tar_dir, '.tar'))
    is_tar_bz2 = endsWith(md_tar_dir, '.tar.bz2')
    if(! (is_tar | is_tar_bz2) ) stop('Not a tar or tar.bz2 archive')

    nc = nchar(md_tar_dir)
    if(is_tar){
        # Put the tar directory in this location -- note we remove the '.tar' from the end
        tar_extraction_dir = paste0(temp_dir_base,  substring(md_tar_dir, 1, nc-4))
    }else if(is_tar_bz2){
        # Put the tar directory in this location -- note we remove the '.tar.bz2' from the end
        tar_extraction_dir = paste0(temp_dir_base,  substring(md_tar_dir, 1, nc-8))
    }else{
        stop('this should never be reached!!')
    }
    dir.create(tar_extraction_dir, showWarnings=FALSE, recursive=TRUE)

    # Check that we don't have too many untarred files on the go -- it could be a sign
    # that the code is failing, which could lead to runaway creation of untarred files (and 
    # fill out file system)
    existing_untar_folders = Sys.glob(paste0(dirname(dirname(tar_extraction_dir)), '/*/RUN*'))
    if(length(existing_untar_folders) > MAX_UNTARRED_FOLDERS){
        stop('Too many untarred folders, will not continue')
    }
   
    # Do the actual untarring (automatically deals with .tar.bz2)
    untar_exit_code = untar(md_tar_dir, exdir=dirname(tar_extraction_dir))
    if(untar_exit_code != 0) stop(paste0('Non-zero exit code when untarring ', md_dir_tar))

    return(tar_extraction_dir)

}


# Trim the tide-gauge observations in a gauge_data object to the time of the model.
# This is useful to prevent storage of too many large files.
trim_observations_to_model_time<-function(gauge_data){

    # Find the start/end times where we want data, in 'numeric julian days'
    start_time = as.numeric(gauge_data$model_start_julian)
    model_duration_days = diff(range(gauge_data$model_time))/(3600 * 24) # model_time is in seconds
    end_time = start_time + model_duration_days

    # Find the data subset matching the model time, plus a buffer
    data_keep = which( (gauge_data$event_data$obs$juliant > (start_time - 2/24)) &
                       (gauge_data$event_data$obs$juliant < (end_time   + 2/24)) )
    if(length(data_keep) > 0){
        # Reduce the data size
        gauge_data$event_data$obs = gauge_data$event_data$obs[data_keep,]
    }else{
        # Just return the head. At least we know this will not be large, and it will 
        # not cause big problems
        gauge_data$event_data$obs = head(gauge_data$event_data$obs)
    }
    return(gauge_data)
}

trim_RDS_file_and_zip_log_file<-function(md_dir){

    RDS_file = Sys.glob(paste0(md_dir, '/gauges_*.RDS'))
    if(length(RDS_file) != 1){
        stop(paste0('Non-unique RDS files: ', RDS_file, collapse=" "))
    }

    # Read the gauge data and trim the observations to the model time.    
    old_data = readRDS(RDS_file)
    new_data = lapply(old_data, trim_observations_to_model_time)
    # Overwrite the old RDS file
    saveRDS(new_data, RDS_file)

    # Make a zipped version of the 'mpi-rank-0' multidomain_log to the new location
    md_file = Sys.glob(paste0(dirname(RDS_file), '/multidomain*00001.log'))
    new_md_zip_file = paste0(md_file, '.zip')
    zip(new_md_zip_file, md_file)

    return(0)
}


# Given a multidomain directory, call the appropriate gauge plotting script which will create
# a gauge RDS file (which is crucial for later analysis) and maybe a pdf (depending on a hard-coded option
# in plots/plot_gauges_generic_include.R).
create_gauges_and_possibly_pdf_plot<-function(md_dir){

    nc = nchar(md_dir)
    if( endsWith(md_dir, '.tar') | endsWith(md_dir, '.tar.bz2') ){
        stop('Cannot create gauges on a tar archive (it might need to be unpacked)')
    }

    # We have a distinct script associated with each historic event, which knows about
    # the relevant data. The names should match part of the model output directory after
    # we convert to lower-case and remove underscores from the latter.
    event2script = list(
        'Chile1960'         = 'plot_southamerica1960.R',
        'Chile2010'         = 'plot_southamerica2010.R',
        'SouthAmerica2015'  = 'plot_southamerica2015.R',
        'Chile2015'         = 'plot_southamerica2015.R', # Different naming used for historical inversions
        'SouthAmerica2014'  = 'plot_southamerica2014.R',
        'Chile2014'         = 'plot_southamerica2014.R', # Different naming used for historical inversions
        'Sumatra2004'       = 'plot_sumatra2004.R',
        'Sumatra2005'       = 'plot_sumatra2005.R',
        'Java2006'          = 'plot_java2006.R',
        'Sumatra2007'       = 'plot_sumatra2007.R',
        'Tohoku2011'        = 'plot_tohoku2011.R',
        'Puysegur2009'      = 'plot_puysegur2009.R',
        'Solomon2007'       = 'plot_solomon2007.R',
        'Kermadec2021'      = 'plot_kermadectonga2021.R', 
        'NewHebrides2021'   = 'plot_newhebrides2021.R', 
        'Sandwich2021'      = 'plot_sandwich2021.R'
        )

    # Find the plotting script of interest, by matching the name in 'event2script'
    # with the model output directory
    plot_script = NA
    for(i in 1:length(event2script)){
        if(grepl(names(event2script)[i], gsub('_', '', md_dir, fixed=TRUE), ignore.case=TRUE)){
            plot_script = event2script[[i]]
            break
        }
    }
    if(is.na(plot_script)) stop(paste0('Could not find a plotting script matching ', md_dir))
    if(!file.exists(paste0(PLOT_SCRIPT_BASE, plot_script))){
        stop(paste0('Could not find the following script ', plot_script))
    }

    # Prepare to run the plotting script
    # First we have to make sure neither PLOT_SCRIPT_BASE nor md_dir start with
    # './', so I can cleanly compute the file paths
    #
    start_dir = getwd()
    if(startsWith(PLOT_SCRIPT_BASE, './')){
        PLOT_SCRIPT_BASE = substring(PLOT_SCRIPT_BASE, 3, nchar(PLOT_SCRIPT_BASE))
    }
    extra_indirection = paste0(rep('../', length(split_folder_path(PLOT_SCRIPT_BASE))-1), collapse="")
    if(startsWith(md_dir, './')) md_dir = substring(md_dir, 3, nchar(md_dir))

    # Change into the PLOT_SCRIPT_BASE directory, run the script, and change back.
    setwd(PLOT_SCRIPT_BASE)
    job_command = paste0('Rscript ', plot_script, ' ', extra_indirection, md_dir)
    print(job_command)
    system(job_command)
    setwd(start_dir)

    return(0)
}

# Make png plots of the max-stage and elevation for the model with multidomain_directory md_dir
plot_max_stage_and_maybe_elevation<-function(md_dir, plot_elevation=FALSE){

    source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')
    library(cptcity)
    library(rptha)

    # Add in domain boundaries
    md_bbox = get_domain_interior_bbox_in_multidomain(md_dir)
    add_bboxes<-function(LWD=1){

        for(i in 1:length(md_bbox$merged_domain_interior_bbox)){

            dx = md_bbox$merged_domain_dx[[i]][1]
            # Avoid boxes for the 1 arcmin domains
            if(dx > 0.9*1/60) next

            LTY=c('solid', 'solid') #c('dashed', 'solid')[(dx < 0.9*1/7*1/60) + 1]
            polygon(md_bbox$merged_domain_interior_bbox[[i]], border='grey', lty=LTY, lwd=LWD)
        }
    }

    #
    # Plot maximum stage
    #
    my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'))(1000)
    # Use a log10 transformation to stretch the max-stage colours
    stage_zlim = c(1.0e-03, 1)
    var_transform_fun<-function(x) log10(pmin(pmax(x, stage_zlim[1]), stage_zlim[2]))
    fields_axis_args = list(at=seq(log10(stage_zlim[1]),log10(stage_zlim[2])),
                            labels=round(10**(seq(log10(stage_zlim[1]),log10(stage_zlim[2]))), 3))
    out_file = paste0(md_dir, '/max_stage_full_model.png')
    png(out_file, width=17.1, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='max_stage', time_index=NULL,
        xlim=c(-40, 320), ylim=c(-79, 68), zlim=log10(stage_zlim),
        cols=my_col, clip_to_zlim=FALSE, use_fields=TRUE,
        NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE,
        var_transform_function=var_transform_fun,
        fields_axis_args=fields_axis_args)
    add_bboxes()
    dev.off()

    if(plot_elevation){
        #
        # Plot elevation
        #

        # Colours for elevation
        my_col = cpt("gmt_GMT_globe", n=1000)
        ZLIM = c(-11e+03, 6500) # This ensures the shoreline is close to a rapid transition in the colours

        out_file = paste0(md_dir, '/elevation0_full_model.png')
        png(out_file, width=17.1, height=8, units='in', res=212.504)
        multidomain_image(md_dir, variable='elevation0', time_index=NULL,
            xlim=c(-40, 320), ylim=c(-79, 68), zlim=ZLIM,
            cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
            NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
        add_bboxes()
        dev.off()
    }
    rm(md_bbox, my_col, add_bboxes, var_transform_fun, field_axis_args)
    gc()
    return(0)

}

# Once the gauges/plotting code have been run on the md directory,
# we should copy the required outputs to a new directory
copy_gauges_and_plots_to_analysis_directory<-function(md_dir){

    # Make sure the output location exists
    dir.create(COPYOUT_DIR_BASE, showWarnings=FALSE, recursive=TRUE)

    # Remove './' from the start of the directory if it exists
    nc = nchar(md_dir)
    if( startsWith(md_dir, './') ) md_dir = substring(md_dir, 3, nc)

    # Make a directory for these outputs in particular
    output_dir = paste0(COPYOUT_DIR_BASE, md_dir)
    dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

    # Find the outputs of interest
    files_of_interest = c(Sys.glob(paste0(md_dir, '/gauge*')),
                          #Sys.glob(paste0(md_dir, '/extra_gauges*')),
                          Sys.glob(paste0(md_dir, '/*.png')), 
                          Sys.glob(paste0(md_dir, '/multidomain*00001.log.zip')))
    nf = length(files_of_interest)
    if(nf == 0) stop(paste0('No files found to copy in ', md_dir))

    # Copy them to the output folder
    for(i in 1:nf) file.copy(files_of_interest[i], output_dir, overwrite=TRUE)

    return(0)
}

# Main driver routine -- this calls the above routines to do the job
extract_key_outputs_from_tarred_multidomain<-function(md_tar_dir){

    if(  length(md_tar_dir) != 1 ) stop('length(md_tar_dir) is not 1')
    if( !file.exists(md_tar_dir) ) stop(paste0('the input tarred directory does not exist: ', md_tar_dir))

    md_dir = try(untar_multidomain_dir_and_return_dir_location(md_tar_dir))
    if(is(md_dir, 'try-error')) return('FAILED TO UNTAR')
    gc()

    create_files_exit_code = try(create_gauges_and_possibly_pdf_plot(md_dir))
    if(is(create_files_exit_code, 'try-error')) return('FAILED TO CREATE GAUGES')
    gc()

    trim_RDS_files_exit_code = try(trim_RDS_file_and_zip_log_file(md_dir))
    if(is(trim_RDS_files_exit_code, 'try-error')) return('FAILED TO TRIM_RDS_FILES')
    gc()

    plot_exit_code = try(plot_max_stage_and_maybe_elevation(md_dir))
    if(is(plot_exit_code, 'try-error')) return('FAILED TO CREATE MAX_STAGE/ELEV PLOTS')
    gc()

    copy_files_exit_code = try(copy_gauges_and_plots_to_analysis_directory(md_dir))
    if(is(copy_files_exit_code, 'try-error')) return('FAILED TO COPY FILES')
    gc()

    # Finally remove the md_dir -- making sure it starts with RUN* to reduce the risk of accidental deletion
    if(startsWith(basename(md_dir), 'RUN')){
        unlink(md_dir, recursive=TRUE, force=TRUE)
    }
    gc()

    return(' ;) ')

}

#
# Driver program here
#

# Set up the parallel cluster
library(parallel)
MC_CORES = 48 # Beware on Gadi, detectCores() will find 96 cores on a 48 core machine
cl = makeCluster(MC_CORES)
clusterExport(cl, varlist=ls(all=TRUE))

PROCESS_RANDOM_SCENARIOS = TRUE # TRUE for random scenarios, FALSE for validation events 

if(PROCESS_RANDOM_SCENARIOS){
    # Find all the tarred multidomain files with random scenarios to operate on
    # I moved old runs [that covered smaller regions than newer runs] so the Sumatra events all use detailed SWWA info
    md_tar_dirs = unlist(lapply(
        c('OUTPUTS', 'OUTPUTS_SCRATCH', 'OUTPUTS_2022_new_events', 'OUTPUTS_2025_extend_WA', 'OUTPUTS_2025_NWWA'),
        #c('OUTPUTS_2025_NWWA'),
        function(x){
            matching_files = Sys.glob(paste0(x, '/ptha18_random_like_historic*/RUN*'))
            ## [The code below is now defunct, since we just move unwanted files to an OBSOLETE subdirectory]
            ## For events from Sumatra, we only want to keep models with 'SWWA' in the title, since
            ## previously the Sumatra models were run with a less extensive high-res zone in WA
            #k = which(grepl('umatra', matching_files) & !grepl('SWWA', matching_files))
            #if(length(k) > 0) matching_files = matching_files[-k]
            return(matching_files)
        })
    )
}else{
    # Validation events only
    md_tar_dirs =  Sys.glob('OUTPUTS_new_validation_events/*/RUN*')
}

# Run the calculations
all_runs = parLapplyLB(cl, md_tar_dirs, extract_key_outputs_from_tarred_multidomain)
print(all_runs)

stopCluster(cl)
