#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
#
# The job commands look something like this:
#     OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '../sources/like_historic/sumatra2005/Fuji_nias2005_unit_sources_SUM_KAJIURA_SMOOTHED.tif' Fuji_sumatra2005_domains181221 full '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 'load_balance_files/load_balance_181221_0.5_0.166666666666667_0.0333333333333333_24MPI.txt' 0.0 > outfile.log
#

# Raster files containing Okada deformation
raster_files = Sys.glob('../sources/hazard/random_*/scenario_initial_conditions/*.tif')
sourcezone_dirs = basename(dirname(dirname(raster_files)))

model_names = paste0('ptha18-GreaterPerth-sealevel60cm-lowres/', sourcezone_dirs, '/',
    'ptha18_random_scenarios_', gsub('.tif', '', basename(raster_files), fixed=TRUE))

jobs_per_qsub_file = 24 # Run this many models in a single qsub file (or less), to stay within NCI walltime limits
number_of_qsub_files = ceiling( length(raster_files)/jobs_per_qsub_file )

# Mapping between the jobs and the qsub file. This cyclic distribution should make the work fairly even.
job_to_qsub_file_index = rep(1:number_of_qsub_files, times=jobs_per_qsub_file)[1:length(raster_files)]

# Each qsub script will begin like this. I request moderately generous run-times.
qsub_command_base = 
'#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -lmem=570GB
#PBS -lncpus=144
#PBS -l wd
#PBS -l storage=gdata/w85+scratch/w85

source SWALS_ifort_modules_2022.sh
# Load R as well (just for tarring directories)
module load R/4.0.0 

'

#
# Make an evironment that will create our model run command. This requires many machine/job variables,
# so I put them in an environment to avoid polluting the namespace
#
model_run_creation_env = new.env()
with(model_run_creation_env, {

    # How many openmp threads
    omp_num_threads = 6

    # How many mpi processes
    mpi_np = 24

    # How many nodes we use on Gadi
    number_nodes = 3

    # Gadi machine data
    cores_per_node = 48
    sockets_per_node = 2
    memory_GB_per_node = 192

    # MSL
    ambient_sea_level = 0.6

    # Ensure we do the long runs at full resolution
    run_type = 'full'


    # Load balance file
    #load_balance_file = 'load_balance_files/load_balance_181221_0.5_0.166666666666667_0.0333333333333333_24MPI.txt'
    load_balance_file = 'load_balance_files/load_balance_120322_0.166666666666667_0.0333333333333333_24MPI_lowres.txt'

    # Three levels of nested domains
    #first_level_nesting_file = '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' 
    #second_level_nesting_file = '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' 
    #third_level_nesting_file = '../multidomain_design/domains_181221_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 
    # Updated domains
    first_level_nesting_file = '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/first_level_nesting_edited.csv' 
    second_level_nesting_file = '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/second_level_nesting_edited.csv' 
    third_level_nesting_file = '../multidomain_design/domains_010322_0.5_0.166666666666667_0.0333333333333333/third_level_nesting_edited.csv' 

    number_core = number_nodes * cores_per_node
    ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

    # Logical checks
    if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
    if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')
    if(number_nodes != 3 | run_type != 'full' ){
        stop("Need to make load-balance files for this case")
    }


    # Make the commands for each job
    make_job_command<-function(raster_file, model_name){


        job_command = paste0(
            'OMP_NUM_THREADS=', omp_num_threads, ' OMP_PROC_BIND=true ', 
            'mpiexec -np ', mpi_np, ' --map-by ppr:', ppr_socket, ':socket:PE=', omp_num_threads, ' ',
            './model ', raster_file, 
            ' ', model_name,
            ' ', run_type,
            ' ', first_level_nesting_file,
            ' ', second_level_nesting_file,
            ' ', third_level_nesting_file,
            ' ', load_balance_file,
            ' ', ambient_sea_level,
            ' > outfile.log')
        #job_command_out = cat(job_command, sep="\n")

        return(job_command)
    }
})

# Make the files
for(i in 1:number_of_qsub_files){
    output_qsub_filename = paste0('run_greaterperth_sealevel60cm_lowres_', 1000 + i, '.sh')

    # Get all the individual model run commands
    jobs_to_run = which(job_to_qsub_file_index == i)
    job_run_commands = rep(NA, length(jobs_to_run))
    for(j in 1:length(jobs_to_run)){
        ind = jobs_to_run[j]
        job_run_commands[j] = model_run_creation_env$make_job_command(raster_files[ind], model_names[ind])
        # Make a newline
        job_run_commands[j] = paste0(job_run_commands[j], '\n')
        # TAR THE FOLDER
        #     - SWALS appends '-full-ambient_sea_level_0.6' to the model_name to make the output folder
        run_output_folder_match = paste0('./OUTPUTS/', model_names[ind])
        job_run_commands[j] = paste0(job_run_commands[j], '# Tar the results \n', 'Rscript tar_and_remove_matching_dir.R ', 
            run_output_folder_match, ' \n')
    }

    all_run_commands = c(qsub_command_base, job_run_commands)

    cat(all_run_commands, file=output_qsub_filename, sep="\n")
}
