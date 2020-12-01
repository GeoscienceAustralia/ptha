#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
#
# The job commands look something like this:
#    OMP_NUM_THREADS=24 OMP_PROC_BIND=true mpiexec -np 4 --map-by ppr:1:socket:PE=24 ./model ../sources/like_historic/kermadectonga2_26849_variable_uniform_gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0.Rdata.tif 0 0.0 Tonga2006_PTHA18_VAUS_26849_animation 'full' 'load_balance_partition.txt' linear_with_manning 0.035 tonga regional animation > outfile.log
#

# Raster files containing Okada deformation
raster_files = Sys.glob('../sources/random/scenario_initial_conditions/*.tif')
model_names = paste0('ptha18_tonga_MSL0_meshrefine2/ptha18_random_scenarios_', gsub('.tif', '', basename(raster_files), fixed=TRUE))

jobs_per_qsub_file = 15 # Run this many models in a single qsub file (or less), to stay within NCI walltime limits
number_of_qsub_files = ceiling( length(raster_files)/jobs_per_qsub_file )

# Mapping between the jobs and the qsub file. This cyclic distribution should make the work fairly even.
job_to_qsub_file_index = rep(1:number_of_qsub_files, times=jobs_per_qsub_file)[1:length(raster_files)]

# Each qsub script will begin like this. I request moderately generous run-times.
# Note the NCI gadi page suggests we should only ask for 190GB/Node max memory [even though I understand they have 192GB]
# https://opus.nci.org.au/display/Help/Queue+Limits
running_on_scratch_only=TRUE # FALSE # Workaround for NCI gdata issues that are happening at the moment.
if(running_on_scratch_only){
qsub_command_base = 
'#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=380GB
#PBS -lncpus=96
#PBS -l wd
#PBS -l storage=scratch/n74

source SWALS_ifort_modules.sh
ulimit -s unlimited

'

}else{
qsub_command_base = 
'#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -lmem=380GB
#PBS -lncpus=96
#PBS -l wd
#PBS -l storage=scratch/n74+gdata/n74+gdata/w85

source SWALS_ifort_modules.sh
ulimit -s unlimited

'
}
#
# Make an evironment that will create our model run command. This requires many machine/job variables,
# so I put them in an environment to avoid polluting the namespace
#
model_run_creation_env = new.env()
with(model_run_creation_env, {

    # How many openmp threads
    omp_num_threads = 24

    # How many mpi processes
    mpi_np = 4 

    # How many nodes we use on Gadi
    number_nodes = 2

    # Gadi machine data
    cores_per_node = 48
    sockets_per_node = 2
    memory_GB_per_node = 192

    # Earthquake rise-time
    rise_time = 0.0

    # MSL
    ambient_sea_level = 0.0

    # Ensure we do the long runs at full resolution
    run_type = 'full'

    # Model type
    offshore_model_type = 'linear_with_manning'

    # Name of the highres region
    highres_region = 'tonga'

    # Control on the outer domain extent
    outer_domain_extent = 'regional'

    # Control how often gridded outputs are written
    output_style = 'few_grids'

    # Load balance file
    load_balance_file = 'load_balance_partition.txt'

    # Global model manning
    offshore_manning = 0.035

    number_core = number_nodes * cores_per_node
    ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

    # Logical checks
    if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
    if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')
    if(! any(offshore_model_type %in% c('linear_with_manning', 'linear_with_no_friction', 'linear_with_linear_friction'))){
        stop('unknown offshore model type')
    }
    if(number_nodes != 2 | run_type != 'full' | highres_region != 'tonga' | outer_domain_extent != 'regional'){
        stop("Need to make load-balance files for this case")
    }



    # Make the commands for each job
    make_job_command<-function(raster_file, model_name){


        job_command = paste0(
            'OMP_NUM_THREADS=', omp_num_threads, ' OMP_PROC_BIND=true ', 
            'mpiexec -np ', mpi_np, ' --map-by ppr:', ppr_socket, ':socket:PE=', omp_num_threads, ' ',
            './model ', raster_file, ' ', rise_time, ' ', ambient_sea_level, ' ', model_name, ' ',
            run_type, ' ',  load_balance_file, ' ', offshore_model_type, ' ', offshore_manning, ' ',
            highres_region, ' ', outer_domain_extent, ' ', output_style, ' > outfile.log' )
        #job_command_out = cat(job_command, sep="\n")

        return(job_command)
    }
})

# Make the files
for(i in 1:number_of_qsub_files){
    output_qsub_filename = paste0('run_ptha18_random_tonga_msl0_meshrefine2_scratch_only', 1000 + i, '.sh')

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
