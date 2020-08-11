#
# Automate (semi) the creation of commands to run the jobs.
# Less error prone.
#
# Beware this script doesn't create the full PBS script. It only makes the 'run commands', which look like this:
# OMP_NUM_THREADS=12 OMP_PROC_BIND=true mpiexec -np 32 --map-by ppr:2:socket:PE=12 ./model ../sources/Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif 0.0 Tohoku2011_Romano2015 'test_load_balance' 'load_balance_files/load_balance_default_australia_8nodes_32mpi.txt' linear_with_linear_friction 0.0 australia > outfile.log
#

# Raster files containing Okada deformation, which should have kajiura smoothing applied.
raster_base_name = '../sources/'
raster_files = c(
    'Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Chile1960/HoEtAl2019/Ho_Chile1960_initial_displacement.tif',
    'Sumatra2004/FujiSatake2007/Fuji_andaman2004_unit_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Sumatra2004/LoritoEtAl2010/Lorito_sumatra2004_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Sumatra2004/PiatanesiLorito2007/Piatanesi_sumatra2004_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Chile2010/LoritoEtAl2011/Lorito_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Tohoku2011/SatakeEtAl2013/Satake_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Tohoku2011/YamakaziEtAl2018/yamakazi18_Tohoku11_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Tohoku2011/RomanoEtAl2015/Tohoku2011_Romano_source_KAJIURA_SMOOTHED.tif',
    'Chile2015/WilliamsonEtAl2017/Williamson_chile_sources_SUM_KAJIURA_SMOOTHED.tif',
    'Chile2015/RomanoEtAl2016/Illapel_2015_Romano_KAJIURA_SMOOTHED.tif')
model_names = paste0(dirname(dirname(raster_files)), '_', basename(dirname(raster_files)))

# How many openmp threads
omp_num_threads = 12
# How many mpi processes
mpi_np = 32

# How many nodes on Gadi
number_nodes = 8

# Gadi machine data
cores_per_node = 48
sockets_per_node = 2
memory_GB_per_node = 96

# Earthquake rise-time
rise_time = 0.0

# Ensure we do the long runs at full resolution
run_type = 'full'

# Model type
offshore_model_type = 'linear_with_no_friction'

number_core = number_nodes * cores_per_node
ppr_socket = (cores_per_node/sockets_per_node)/omp_num_threads

# Logical checks
if( (omp_num_threads * mpi_np) != (cores_per_node * number_nodes)) stop('Not fully utilising cores')
if( round(ppr_socket) != ppr_socket ) stop('Need integer number of processes per socket')
if(! any(offshore_model_type %in% c('linear_with_manning', 'linear_with_no_friction', 'linear_with_linear_friction'))) stop('unknown offshore model type')
if(number_nodes != 8 | run_type != 'full') stop("Need to make load-balance files for this case")

if(grepl('manning', offshore_model_type)){
    offshore_manning = 0.035
}else{
    offshore_manning = 0.0
}


# Make the commands for each job
for(i in 1:length(model_names)){

    model_name = model_names[i]

    if(grepl('Sumatra2004', model_name)){
        # For the Sumatra simulations include NSW,Vic,WA at highres, because we have data.
        highres_region = 'australia'
    }else{
        # For simulations other than Sumatra, I only have data in NSW, so only refine grid there.
        highres_region = 'NSW'
    }
    if(! any(highres_region %in% c('NSW', 'australia', 'none'))) stop('unknown highres region')


    # Choose load balance file based on highres_region and offshore model type
    if(highres_region == 'NSW'){
        if(offshore_model_type == 'linear_with_manning'){
            load_balance_file = 'load_balance_files/load_balance_manning_offshore_NSW_8nodes_32ranks.txt'
        }else{
            load_balance_file = 'load_balance_files/load_balance_linear_offshore_NSW_8nodes_32ranks.txt'
        }
    }else if(highres_region == 'australia'){
        if(offshore_model_type == 'linear_with_manning'){
            load_balance_file = 'load_balance_files/load_balance_manning_offshore_australia_8nodes_32ranks.txt'
        }else{
            load_balance_file = 'load_balance_files/load_balance_linear_offshore_australia_8nodes_32ranks.txt'
        }
    }else if(highres_region == 'none'){
        load_balance_file = ''
    }else{
        stop('unknown highres_region')
    }

    job_command = paste0(
        'OMP_NUM_THREADS=', omp_num_threads, ' OMP_PROC_BIND=true ', 
        'mpiexec -np ', mpi_np, ' --map-by ppr:', ppr_socket, ':socket:PE=', omp_num_threads, ' ',
        './model ', raster_base_name, raster_files[i], ' ', rise_time, ' ', model_names[i], ' ',
        run_type, ' ',  load_balance_file, ' ', offshore_model_type, ' ', offshore_manning, ' ',
        highres_region, ' > outfile.log' )
    cat(job_command, sep="\n")
}
