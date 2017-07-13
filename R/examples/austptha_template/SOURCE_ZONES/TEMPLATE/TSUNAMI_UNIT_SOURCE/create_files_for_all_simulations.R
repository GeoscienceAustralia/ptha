# Batch produce files for all runs.

#
# INPUT PARAMETERS
#
config_env = new.env()
source('config.R', local=config_env)
#
# END INPUT
#


base_dir = getwd()
dir.create(config_env$all_runs_dir, recursive=TRUE)

for(file_ind in 1:length(config_env$initial_condition_files)){

    initial_condition_file = config_env$initial_condition_files[file_ind]
    initial_condition_file_basename_noext = 
        strsplit(basename(initial_condition_file), split='[.]')[[1]][1]

    # Step 1: Make a directory in 'all_runs_dir' with a name that matches the
    # initial condition
    run_dir = paste0(config_env$all_runs_dir, '/RUN_', 
        format(Sys.time(), format='%Y%m%d%H%M%S'), '_', 
        initial_condition_file_basename_noext
        )
    dir.create(run_dir, recursive=TRUE, showWarnings=FALSE)

    # Step 2: Copy initial condition file to the run directory
    file.copy(initial_condition_file, run_dir)

    # Step 3: Copy the template files to the output directory
    # Do this in a way that preserves symbolic links
    copy_command = paste0('cp -a ./template/* ', run_dir)
    system(copy_command)

    # Step 4: Change the names in jagurs input and output files to match the
    # input/output
    run_output_dir = paste0(config_env$all_runs_output_base_dir, '/', run_dir)
    dir.create(run_output_dir, recursive=TRUE, showWarnings=FALSE)

    setwd(run_dir)
    sed_output_replace_command = paste0('sed -i s:OUTPUTDIRREPLACEME:', 
        run_output_dir, ':g *')
    system(sed_output_replace_command)

    sed_input_deformation_replace_command = paste0(
        'sed -i s:STAGEDEFORMATIONFILEREPLACEME:', 
        initial_condition_file, ':g *')
    system(sed_input_deformation_replace_command)
  
    # Start fresh
    setwd(base_dir)
}


