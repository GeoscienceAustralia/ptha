#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
#
# The job commands look something like this:
## OMP_NUM_THREADS=13 mpiexec -np 16 --map-by ppr:4:socket:PE=13 ./model_sapphirerapids \
##    '../sources/like_historic/KermadecTonga2021/Romano2021/Kermadec2021_Romano.tif' \
##    'multidomain_design_control_NNL4_1arcminoffshore.nml' \
##    run_Kermadec2021_Romano_1arcminoffshore full 0.0 > outfile.log
##  
##

# Raster files containing Okada deformation
raster_files = Sys.glob('../sources/hazard/scenarios_ID4186.3/random_*/scenario_initial_conditions/*.tif')
sourcezone_dirs = basename(dirname(dirname(raster_files)))

model_names = paste0('ptha18-NSW2023-ID4186.3-sealevel110cm/', sourcezone_dirs, '/',
    'ptha18_random_scenarios_', gsub('.tif', '', basename(raster_files), fixed=TRUE))

# Folders for scenarios we have already run. If any of these scenarios are the same
# as the new scenarios then, to save computational effort, we can make a symbolic
# link to the old runs (rather than re-running the inundation calculation).
previous_scenarios_already_run_folders = c(
    Sys.glob('./OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_*/ptha18_*'),
    Sys.glob('./OUTPUTS/ptha18-NSW2023-ID1315.5-sealevel110cm/random_*/ptha18_*'))

match_to_previous_runs = sapply(
    basename(model_names), # Names of new model with prefix removed
    function(x){
        mtch = grep(x, previous_scenarios_already_run_folders)
        if(length(mtch) == 0){
            return(NA)
        }else{
            return(mtch[1])
        }}
    )

jobs_per_qsub_file = 3 # Run this many models in a single qsub file (or less), to stay within NCI walltime limits
number_of_qsub_files = ceiling( length(raster_files)/jobs_per_qsub_file )

# Mapping between the jobs and the qsub file. This cyclic distribution should make the work fairly even.
job_to_qsub_file_index = rep(1:number_of_qsub_files, times=jobs_per_qsub_file)[1:length(raster_files)]

# Each qsub script will begin like this. I request moderately generous run-times.
qsub_command_base = 
'#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=gdata/w85+scratch/w85

source SWALS_ifort_modules_2023_B.sh
# Load R as well (just for tarring directories)
module load R/4.3.1

'

#
# Make an evironment that will create our model run command. This requires many machine/job variables,
# so I put them in an environment to avoid polluting the namespace
#
model_run_creation_env = new.env()
with(model_run_creation_env, {

    # How many openmp threads
    omp_num_threads = 13

    # How many mpi processes
    mpi_np = 16

    # How many nodes we use on Gadi
    number_nodes = 2

    # Gadi machine data
    cores_per_node = 104
    sockets_per_node = 2
    memory_GB_per_node = 500

    # MSL
    ambient_sea_level = 1.1

    # Ensure we do the long runs at full resolution
    run_type = 'full'

    # Multidomain configuration file
    multidomain_design_control_namelist = 'multidomain_design_control_NNL4_1arcminoffshore.nml'

    number_core = number_nodes * cores_per_node
    ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

    # Logical checks
    if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
    if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')


    # Make the commands for each job
    make_job_command<-function(raster_file, model_name){

        slash_newline = " \\\n    "

        job_command = paste0(
            'mpiexec -np ', mpi_np, ' --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=', omp_num_threads, 
            ' -x OMP_PROC_BIND=true ',
            './model_sapphirerapids ', slash_newline, 
            raster_file, slash_newline,
            multidomain_design_control_namelist, slash_newline,
            model_name, slash_newline,
            run_type, slash_newline,
            ambient_sea_level,
            ' > outfile.log')
        #job_command_out = cat(job_command, sep="\n")

        return(job_command)
    }
})

# Make the files
for(i in 1:number_of_qsub_files){
    output_qsub_filename = paste0('run_ptha18_NSW2023_ID4186.3_sealevel110cm_', 1000 + i, '.sh')

    # Get all the individual model run commands
    jobs_to_run = which(job_to_qsub_file_index == i)
    job_run_commands = rep(NA, length(jobs_to_run))
    for(j in 1:length(jobs_to_run)){
        ind = jobs_to_run[j]

        if(is.na(match_to_previous_runs[ind])){
            # Typical case for which initial condition has not been previously simulated to inundation.
            # Create the commands that will do that.

            job_run_commands[j] = model_run_creation_env$make_job_command(raster_files[ind], model_names[ind])
            # Make a newline
            job_run_commands[j] = paste0(job_run_commands[j], '\n')
            # TAR THE FOLDER
            run_output_folder_match = paste0('./OUTPUTS/', model_names[ind])
            job_run_commands[j] = paste0(job_run_commands[j], '# Tar the results \n', 
                'Rscript tar_and_remove_matching_dir.R ', 
                run_output_folder_match, ' \n')

        }else if(is.finite(match_to_previous_runs[ind])){
            # Case where the initial condition has previously been simulated to inundation.
            # In this case we should just make a symbolic link to the previous run
            job_run_commands[j] = "" # No actual command

            new_run_output_folder = paste0('./OUTPUTS/', model_names[ind])
            previous_run_output_folder = previous_scenarios_already_run_folders[match_to_previous_runs[ind]]

            # Make the parent folder needed to store the simulation
            parent_folder_of_new_run = dirname(new_run_output_folder)
            dir.create(parent_folder_of_new_run, recursive=TRUE, showWarnings=FALSE)

            # Go into the simulation parent folder, make the symbolic link,
            # then move back to the current folder
            starting_dir = getwd()
            setwd(parent_folder_of_new_run)
            # The symbolic link needs to be relative to the
            # "parent_folder_of_new_run"
            folder_to_link_to = paste0('../../', gsub('./OUTPUTS/', '', previous_run_output_folder, fixed=TRUE))
            folder_to_link_to_exists = file.exists(folder_to_link_to) 
            if(!folder_to_link_to_exists){
                # Error handling
                setwd(starting_dir)
                stop(paste0('could not link to folder ', folder_to_link_to, ' from inside ', parent_folder_of_new_run))
            }
            # Make the link [this will work even if the link already exists]
            system(paste0('ln -sf ', folder_to_link_to))
            setwd(starting_dir)

        }else{
            stop('Problem with matching previous runs')
        }
    }

    all_run_commands = c(qsub_command_base, job_run_commands)

    cat(all_run_commands, file=output_qsub_filename, sep="\n")
}
