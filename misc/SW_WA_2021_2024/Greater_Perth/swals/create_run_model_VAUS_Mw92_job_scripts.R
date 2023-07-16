#
# Make PBS submission scripts for a set of scenarios
#


JOB_SCRIPT_TEMPLATE = "
#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=5:00:00
#PBS -lmem=1140GB
#PBS -lncpus=288
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85

source SWALS_ifort_modules_2021.sh

OMP_NUM_THREADS=12 mpiexec -np 24 --map-by ppr:2:socket:PE=12 ./model '__SOURCE_FILE__' __SOURCE_NAME__ full '../multidomain_design/domains_0.5_0.125/first_level_nesting.csv' '../multidomain_design/domains_0.5_0.125/second_level_nesting.csv' 'load_balance_files/load_balance_120821_9x_0.5_0.125_24MPI.txt' 0.0 > outfile.log
"

SCENARIO_TIFS = Sys.glob('../sources/testing/VAUS_Mw92_sunda_arc/Mw92_scenarios_similar_size/*.tif')

for(i in 1:length(SCENARIO_TIFS)){
    scenario_tif_file = SCENARIO_TIFS[i]
    scenario_tif_name = paste0('VAUS_Mw92_sunda_arc_', gsub('.tif', '', basename(scenario_tif_file), fixed=TRUE))

    scenario_job_text = gsub('__SOURCE_FILE__', scenario_tif_file, JOB_SCRIPT_TEMPLATE, fixed=TRUE)
    scenario_job_text = gsub('__SOURCE_NAME__', scenario_tif_name, scenario_job_text, fixed=TRUE)

    scenario_job_filename = paste0('run_', scenario_tif_name, '.sh')
    writeLines(scenario_job_text, scenario_job_filename)
}
