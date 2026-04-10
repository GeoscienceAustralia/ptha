#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
##

# Raster files containing Okada deformation that we will use for tidal testing.
raster_files = c(
    "../sources/hazard/random_outerrisesunda/scenario_initial_conditions/outerrisesunda_row_0016470_Mw_81_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0098440_Mw_89_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0101901_Mw_90_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0105527_Mw_92_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0106688_Mw_93_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0107985_Mw_94_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0108069_Mw_94_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0108100_Mw_94_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0108285_Mw_94_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0108362_Mw_94_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0109368_Mw_95_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0109369_Mw_95_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0109374_Mw_95_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0109379_Mw_95_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0109700_Mw_95_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0110660_Mw_96_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0110693_Mw_96_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0110734_Mw_96_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0110907_Mw_96_HS.tif",
    "../sources/hazard/random_sunda2/scenario_initial_conditions/sunda2_row_0111004_Mw_96_HS.tif"
    )
sourcezone_dirs = basename(dirname(dirname(raster_files)))
model_group_name = 'ptha18-kalbarri2onslow-tidal_testing'
model_names = paste0(model_group_name, '/', sourcezone_dirs, '/',
    'ptha18_tidal_testing_', gsub('.tif', '', basename(raster_files), fixed=TRUE))

# Run at a set of fixed tidal levels as well as with the varying tides.
# All tidal levels (plus the varying tide model) are included in the same file.
fixed_tidal_levels = c(0, 0.65, 0.95, 1.17, 1.42) # NOTE: Tidal level 0 will automatically use the elevation adjustment

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

source SWALS_ifx_modules_2025_llvm.sh
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
    ambient_sea_level = 0.0

    # Ensure we do the long runs at full resolution
    run_type = 'full'

    # Multidomain configuration file
    multidomain_design_control_namelist_elevadjust = 'multidomain_design_control/multidomain_kalbarri2onslow_20251218_hazard.nml'
    multidomain_design_control_namelist_noelevadjust = 'multidomain_design_control/multidomain_kalbarri2onslow_20251218_hazard_noelevadjust.nml'

    number_core = number_nodes * cores_per_node
    ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

    # Logical checks
    if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
    if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')


    # Make the commands for each job
    make_job_command<-function(raster_file, model_name, background_sea_level){

        slash_newline = " \\\n    "

        if(background_sea_level == 0){
            multidomain_design_control_namelist = multidomain_design_control_namelist_elevadjust
        }else{
            multidomain_design_control_namelist = multidomain_design_control_namelist_noelevadjust
        }

        job_command = paste0(
            'mpiexec -np ', mpi_np, ' --map-by numa:SPAN --bind-to numa -x OMP_NUM_THREADS=', omp_num_threads, 
            ' -x OMP_PROC_BIND=true ',
            './model_build ', slash_newline, 
            raster_file, slash_newline,
            multidomain_design_control_namelist, slash_newline,
            model_name, slash_newline,
            run_type, slash_newline,
            background_sea_level,
            ' > outfile.log')
        #job_command_out = cat(job_command, sep="\n")

        return(job_command)
    }
})

# Make the files
for(i in 1:length(raster_files)){
    output_qsub_filename = paste0(model_group_name, '_', 1000 + i, '.sh')

    # This job will run raster_file[i] with a few different sea levels
    jobs_to_run =  i
    job_run_commands = rep(NA, length(fixed_tidal_levels))
    for(j in 1:length(job_run_commands)){
        job_run_commands[j] = model_run_creation_env$make_job_command(raster_files[i], model_names[i], fixed_tidal_levels[j])
        # Make a newline
        job_run_commands[j] = paste0(job_run_commands[j], '\n')
        # TAR THE FOLDER
        run_output_folder_match = paste0('./OUTPUTS/', model_names[i], '*ambient_sea_level_', fixed_tidal_levels[j])
        job_run_commands[j] = paste0(job_run_commands[j], '# Tar the results \n', 'Rscript tar_and_remove_matching_dir.R ', 
            '"', run_output_folder_match, '"', ' \n')
    }

    all_run_commands = c(qsub_command_base, job_run_commands)

    cat(all_run_commands, file=output_qsub_filename, sep="\n")
}
