regular_files = Sys.glob('OUTPUTS/run_kt43731_12h_final_NNL4_1arcminoffshore-full-ambient_sea_level_1.1/RUN_20241112_155633335/max_stage_domain_*.tif')
convergence_files = Sys.glob('OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/max_stage_domain_*.tif')

# Get total number of cells and non-NA cells
count_cells=function(x){
    library(terra)
    r1 = rast(x)
    total= prod(dim(r1)) # Beware this misses most halos from domain 1 (since it was merged to make the tif)
    total_nonNA = sum(!is.na(as.matrix(r1)))
    c(total, total_nonNA)
    }
regular_cell_counts = lapply(regular_files, count_cells)
regular_total_cells = sum(unlist(lapply(regular_cell_counts, function(x) x[1]))) # This doesn't capture halos on the global domains
regular_total_cells_noNA = sum(unlist(lapply(regular_cell_counts, function(x) x[2])))
print(paste0('Regular model uses ', regular_total_cells_noNA, ' cells in priority domain areas'))
print(paste0('Regular model uses more than ', regular_total_cells, ' cells in total (this calculation misses halos on domain 1)'))

convergence_cell_counts = lapply(convergence_files, count_cells)
convergence_total_cells = sum(unlist(lapply(convergence_cell_counts, function(x) x[1]))) # This doesn't capture halos on the global domains
convergence_total_cells_noNA = sum(unlist(lapply(convergence_cell_counts, function(x) x[2])))
print(paste0('convergence test model uses ', convergence_total_cells_noNA, ' cells in priority domain areas'))
print(paste0('convergence test model uses more than ', convergence_total_cells, ' cells in total (this calculation misses halos on domain 1)'))

