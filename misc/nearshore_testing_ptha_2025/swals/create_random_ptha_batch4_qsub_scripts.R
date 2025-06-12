#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
#
# The job commands look something like this:
# OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 'test_load_balance' 'load_balance_files/load_balance_default_australia_8nodes_32mpi.txt' linear_with_linear_friction 0.0 australia > outfile.log
#

# Raster files containing Okada deformation
raster_files = Sys.glob('../ptha18_scenarios_random/set_range_of_mw_and_centroid_batch4/*/*.tif')
model_names = paste0('ptha18_random_like_historic_', gsub('.tif', '', basename(raster_files), fixed=TRUE))

jobs_per_qsub_file = 12
number_of_qsub_files = ceiling( length(raster_files)/jobs_per_qsub_file )

# Mapping between the jobs and the qsub file. This cyclic distribution should make the work fairly even.
job_to_qsub_file_index = rep(1:number_of_qsub_files, times=jobs_per_qsub_file)[1:length(raster_files)]

# Each qsub script will begin like this. I request moderately generous run-times.
qsub_command_base = 
'#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=28:00:00
#PBS -lmem=1536GB
#PBS -lncpus=384
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2022.sh
ulimit -s unlimited

'

#
# Make an evironment that will create our model run command. This requires many machine/job variables,
# so I put them in an environment to avoid polluting the namespace
#
model_run_creation_env = new.env()
with(model_run_creation_env, {

    # How many openmp threads
    omp_num_threads = 12
    # How many mpi processes
    mpi_np = 32

    # How many nodes on Gadi
    number_nodes = 8

    # Gadi machine data
    cores_per_node = 48
    sockets_per_node = 2
    memory_GB_per_node = 192

    # Earthquake rise-time
    rise_time = 0.0

    # Ensure we do the long runs at full resolution
    run_type = 'full'

    # Model type
    offshore_model_type = 'linear_with_manning'

    number_core = number_nodes * cores_per_node
    ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

    # Logical checks
    if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
    if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')
    if(! any(offshore_model_type %in% c('linear_with_manning', 'linear_with_no_friction', 'linear_with_linear_friction'))){
        stop('unknown offshore model type')
    }
    if(number_nodes != 8 | run_type != 'full') stop("Need to make load-balance files for this case")

    if(grepl('manning', offshore_model_type)){
        offshore_manning = 0.035
    }else{
        offshore_manning = 0.0
    }


    # Make the commands for each job
    make_job_command<-function(raster_file, model_name){

        if(grepl('sandwich2021', model_name)){
            # For the sandwich simulations, include only perth at highres
            highres_region = 'perth'
        }else{
            # For other simulations, include only NSW at highres
            highres_region = 'NSW'
        }
        if(! any(highres_region %in% c('NSW', 'perth'))) stop('unknown highres region')


        # Choose load balance file based on highres_region and offshore model type
        if(highres_region == 'NSW'){
            load_balance_file = 'load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks_2022.txt'
        }else if(highres_region == 'perth'){
            load_balance_file = 'load_balance_files/load_balance_manning_offshore_perth_8nodes_32ranks_2022.txt'
        }else{
            stop('unknown highres_region')
        }

        job_command = paste0(
            'OMP_NUM_THREADS=', omp_num_threads, ' OMP_PROC_BIND=true ', 
            'mpiexec -np ', mpi_np, ' --map-by ppr:', ppr_socket, ':socket:PE=', omp_num_threads, ' ',
            './model ', raster_file, ' ', rise_time, ' ', model_name, ' ',
            run_type, ' ',  load_balance_file, ' ', offshore_model_type, ' ', offshore_manning, ' ',
            highres_region, ' > outfile.log' )
        #job_command_out = cat(job_command, sep="\n")

        return(job_command)
    }
})

# Make the files
for(i in 1:number_of_qsub_files){
    output_qsub_filename = paste0('run_ptha18_batch4_random_like_historic_batch_', 100 + i, '.sh')

    # Get all the individual model run commands
    jobs_to_run = which(job_to_qsub_file_index == i)
    job_run_commands = rep(NA, length(jobs_to_run))
    for(j in 1:length(jobs_to_run)){
        ind = jobs_to_run[j]
        job_run_commands[j] = model_run_creation_env$make_job_command(raster_files[ind], model_names[ind])
        # Make a newline
        job_run_commands[j] = paste0(job_run_commands[j], '\n')
    }

    all_run_commands = c(qsub_command_base, job_run_commands)

    cat(all_run_commands, file=output_qsub_filename, sep="\n")
}
