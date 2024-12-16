# Get the SWALS plot scripts
file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

#
# Cell counts from tif files on home machine
#
if(FALSE){
    # Cell counts for all tifs
    print('# Regular model')
    regular_files = Sys.glob('OUTPUTS/run_kt43731_12h_final_NNL4_1arcminoffshore-full-ambient_sea_level_1.1/RUN_20241112_155633335/elevation0_domain_*.tif')
    count_cells_all_files(regular_files)
# [1] "# Regular model"
# [1] "Total of 267433452 cells priority domain areas that did not have NA raster values"
# [1] "Total of 283608244 cells, except this calculation misses halos on merged domains"
# [1] "  Ratio: 0.942967835589434"
    
    print('')
    print('# CONVERGENCE MODEL')
    convergence_files = Sys.glob('OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/elevation0_domain_*.tif')
    count_cells_all_files(convergence_files)
# [1] ""
# [1] "# CONVERGENCE MODEL"
# [1] "Total of 1457837568 cells priority domain areas that did not have NA raster values"
# [1] "Total of 1527715964 cells, except this calculation misses halos on merged domains"
# [1] "  Ratio: 0.954259562872513"
}

#
# Cell counts from netcdf files on NCI -- these are better.
#
if(FALSE){
    #
    # Cell counts for netcdf, more accurate than above
    #

    all_netcdf_files = Sys.glob('OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/*/Grid*.nc')
    tmp = count_cells_all_files(all_netcdf_files)
# [1] "Total of 267433452 cells in priority domain areas"
# [1] "Total of 284580460 cells in total"
# [1] "  Ratio: 0.939746362065758"

    all_nc_files = Sys.glob('OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/*/Grid*.nc')
    tmp2 = count_cells_all_files(all_nc_files)
# [1] "Total of 1457837568 cells in priority domain areas"
# [1] "Total of 1533668044 cells in total"
# [1] "  Ratio: 0.950556134818964"

}
