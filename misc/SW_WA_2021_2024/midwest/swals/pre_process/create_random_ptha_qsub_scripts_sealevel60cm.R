#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
#
# The job commands look something like this:
## mpirun -n 16 -map-by numa:SPAN -bind-to numa -x OMP_NUM_THREADS=13 -x OMP_PROC_BIND=TRUE \
##    ./model_sapphirerapids \
##    ../sources/like_historic/sumatra2005/Fuji_nias2005_unit_sources_SUM_KAJIURA_SMOOTHED.tif \
##    multidomain_design_control_NNL4_defaultres.nml \
##    Sumatra2005_FujiiEtAl2020 \
##    full \
##    0.0 > outfile.log
##  
##

# Paths for this script
path_swals = "../"
path_source = "../../sources/"

# path to swals for pbs script
path_swals_pbs = "./"

# Raster files containing Okada deformation
# This is relative to the current script
raster_files = Sys.glob(paste0(path_source, 'hazard/random_*/scenario_initial_conditions/*.tif'))
sourcezone_dirs = basename(dirname(dirname(raster_files)))

model_names = paste0('ptha18-midwest-sealevel60cm/', sourcezone_dirs, '/',
    'ptha18_random_scenarios_', gsub('.tif', '', basename(raster_files), fixed=TRUE))

jobs_per_qsub_file = 7 # Run this many models in a single qsub file (or less), to stay within NCI walltime limits
number_of_qsub_files = ceiling( length(raster_files)/jobs_per_qsub_file )

# Mapping between the jobs and the qsub file. This cyclic distribution should make the work fairly even.
job_to_qsub_file_index = rep(1:number_of_qsub_files, times=jobs_per_qsub_file)[1:length(raster_files)]

# Each qsub script will begin like this. I request moderately generous run-times.
qsub_command_pbs_vars = 
'#!/bin/bash
#PBS -P w85
#PBS -q normalsr
#PBS -l walltime=28:00:00
#PBS -lmem=1000GB
#PBS -lncpus=208
#PBS -l wd
#PBS -l storage=gdata/w85+scratch/w85
'

qsub_command_modules = 
'source SWALS_ifort_modules_2023_B.sh
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
    ambient_sea_level = 0.6

    # Ensure we do the long runs at full resolution
    run_type = 'full'

    # Multidomain configuration file
    multidomain_design_control_namelist = 'multidomain_design_control_NNL4_defaultres.nml'

    number_core = number_nodes * cores_per_node
    ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

    # Logical checks
    if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
    if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')


    # Make the commands for each job
    make_job_command<-function(raster_file, model_name){

        # mpi command
        newline_tab = " \\\n    "
        mpi_command = paste0(
            'mpirun -n ', mpi_np,
            ' -map-by numa:SPAN -bind-to numa',
            ' -x OMP_NUM_THREADS=', omp_num_threads,
            ' -x OMP_PROC_BIND=TRUE', newline_tab)

        # model command
        # NB: this path is relative to the swals directory
        log_file = paste0(path_swals_pbs, "run/log/$PBS_JOBID.log")
        job_command = paste0(
            mpi_command,
            paste0(path_swals_pbs, 'model_sapphirerapids '), newline_tab, 
            normalizePath(raster_file), newline_tab,
            paste0(path_swals_pbs, multidomain_design_control_namelist), newline_tab,
            model_name, newline_tab,
            run_type, newline_tab,
            ambient_sea_level,
            ' > ', log_file, '\n')

        return(job_command)
    }
})

# Make the files
for(i in 1:number_of_qsub_files){
    output_qsub_filename = paste0(path_swals, 'run/midwest_sealevel60cm_', 1000 + i, '.pbs')

    # Get all the individual model run commands
    jobs_to_run = which(job_to_qsub_file_index == i)
    job_run_commands = rep(NA, length(jobs_to_run))

    for(j in 1:length(jobs_to_run)){
        ind = jobs_to_run[j]
        job_run_commands[j] = model_run_creation_env$make_job_command(
            raster_files[ind],
            model_names[ind])

        # Make a newline
        job_run_commands[j] = paste0(job_run_commands[j], '\n')
        # FIXME: TAR THE FOLDER
        #     - SWALS appends '-full-ambient_sea_level_0.6' to the model_name to make the output folder
        run_output_folder_match = paste0(path_swals_pbs, 'OUTPUTS/', model_names[ind])
        tar_command = paste0(
            '# Tar the results \n',
            'Rscript ', path_swals_pbs, 'post_process/tar_and_remove_matching_dir.R ', run_output_folder_match, ' \n')
        job_run_commands[j] = paste0(job_run_commands[j], tar_command)
    }

    # Add in the pbs directives and source modules
    all_commands = c(qsub_command_pbs_vars, qsub_command_modules, job_run_commands)

    cat(all_commands, file=output_qsub_filename, sep="\n")
}
warnings()
