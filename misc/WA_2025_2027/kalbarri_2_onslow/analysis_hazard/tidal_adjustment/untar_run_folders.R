#
# Untar the tidal calculations that will be used
#

tars_to_extract = normalizePath(Sys.glob('../../swals/OUTPUTS/ptha18-kalbarri2onslow-tidal_testing/*/*/RUN*.tar'))
basedir = normalizePath(getwd())

# Some scenarios went unstable at one or more sea levels (because I wasn't
# careful enough to give them a short timestep, as was done in the production
# runs to work around these issues). Remove them -- we'll still have enough
# scenarios for this test.
remove_scenarios_that_failed_in_some_cases<-function(scenarios_sea_level_vary){

    scenarios_to_drop = c(
        'sunda2_row_0107985_Mw_94_HS',
        'sunda2_row_0108100_Mw_94_HS',
        'sunda2_row_0110660_Mw_96_HS',
        'sunda2_row_0110734_Mw_96_HS',
        'sunda2_row_0110907_Mw_96_HS',
        'sunda2_row_0109369_Mw_95_HS')

    scenario_keep_flag = rep(TRUE, length(scenarios_sea_level_vary))
    for(i in 1:length(scenarios_to_drop)){
        k = which(grepl(scenarios_to_drop[i], scenarios_sea_level_vary, fixed=TRUE))
        if(length(k) > 0) scenario_keep_flag[k] = FALSE
    }

    scenarios_sea_level_vary[scenario_keep_flag]
}
tars_to_extract = remove_scenarios_that_failed_in_some_cases(tars_to_extract)

# Make function to untar the folders
untar_unless_folder_already_extracted<-function(tarfile){
    on.exit(setwd(basedir)) # Ensure we always go back to the base directory
    setwd(dirname(tarfile)) # Move to the tar file's directory

    # If the expected output folder already exists, then quick exit
    if(file.exists(gsub(".tar", "", basename(tarfile)))){
        return('Already exists (quick exit)') # Error code
    }

    # Do the extraction
    untar_command = paste0("tar -xf ", basename(tarfile))
    untar_worked = system(untar_command)
    if(untar_worked != 0){
        return(paste0('Failed extracting ', tarfile))
    }else{
        return('Success')
    }
}
# Ensure failure in one case doesn't halt parallel execution of all
try_parallel_fun<-function(tarfile) try(untar_unless_folder_already_extracted(tarfile))

# Make a cluster
library(parallel)
MC_CORES=48
cl = makeCluster(MC_CORES)

# Export the workspace to the cluster
nothing = clusterExport(cl, ls(all=TRUE))

# Do the calculations
nothing = parLapplyLB(cl, tars_to_extract, try_parallel_fun, chunk.size=1)

print(unlist(nothing))

stopCluster(cl)
